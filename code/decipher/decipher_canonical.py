"""
Locally trained published-Decipher-marker gene-set surrogate. This is NOT the
licensed commercial Decipher Genomic Classifier and is not a reproduction of its
random forest; it is a separate L2-logistic model fit on a training-only
intersection of the published Erho et al. 2013 marker panel (Table 2, DOI
10.1371/journal.pone.0066855), read from config/decipher_marker_surrogate_v1.csv.
It does not read, refit, or modify the frozen 41-feature Met-Score.

The validation transform is transductive cohort-batch rank mapping to the
training reference, so the serialized artifact is suitable only for retrospective
cohort-batch scoring, not individual deployment.

Modes:
  --mode fit    Fit once on all training patients over the training-frozen panel,
                serialize the artifact, reload it, and generate JHU and Durham
                predictions from the stored coefficients + preprocessing state
                (not the live sklearn object). Writes a manifest last.
  --mode verify Clean process, no fit: load the artifact and authorized raw
                inputs, reconstruct both cohort scores, and compare to the saved
                predictions (max abs diff and class mismatches).
"""
import os
import sys
import csv
import json
import time
import argparse
import hashlib
import warnings
import numpy as np
import pandas as pd
import rdata


def _find_met_pca_root():
    env = os.environ.get("MET_PCA_ROOT")
    if env and os.path.isdir(env):
        return env
    d = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "code")) and (
            os.path.isdir(os.path.join(d, "outs")) or os.path.isdir(os.path.join(d, "output"))):
            return d
        d = os.path.dirname(d)
    raise FileNotFoundError("Could not locate Met_PCa project root; set MET_PCA_ROOT.")


ROOT = _find_met_pca_root()
OUT_DIR = os.path.join(ROOT, "outs", "Decipher")
os.makedirs(OUT_DIR, exist_ok=True)
CFG = os.path.join(ROOT, "config", "decipher_marker_surrogate_v1.csv")
JHU_RDA = os.path.join(ROOT, "outs", "MetastasisData_JHUOut.rda")
DUR_NPZ = os.path.join(OUT_DIR, "_durham_full_expr.npz")
DUR_CLIN = os.path.join(ROOT, "output", "Durham", "durham_metscore_batchcorrected.rda")

ART_NPZ = os.path.join(OUT_DIR, "decipher_surrogate_artifact.npz")
ART_JSON = os.path.join(OUT_DIR, "decipher_surrogate_artifact.json")
THR_TXT = os.path.join(OUT_DIR, "decipher_surrogate_threshold.txt")
MANIFEST = os.path.join(OUT_DIR, "decipher_surrogate_manifest.json")
JHU_PRED = os.path.join(OUT_DIR, "jhu_decipher_pred.csv")
DUR_PRED = os.path.join(OUT_DIR, "durham_decipher_pred.csv")

CV_SEED = 20260427
CV_FOLDS = 5
MAX_ITER = 20_000
CS_GRID = np.logspace(-7, 2, 10)   # documented decade grid 1e-7 .. 1e2
PREPROC_VERSION = "durham_555_before_rank_v1; fold-specific between-array QN in tuning"
N_TRAIN, N_JHU, N_DURHAM = 1000, 239, 555
CLASS_NO_METS, CLASS_METS = 694, 306


# ---- helpers -------------------------------------------------------------
def _sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def _sigmoid(z):
    # numerically stable logistic: no overflow for large |z|
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    e = np.exp(z[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def _short_jhu(s):
    import re
    m = re.search(r"JHU\d+", s)
    return m.group(0) if m else s


def _atomic_write_bytes(path, write_fn):
    tmp = path + ".tmp"
    write_fn(tmp)
    os.replace(tmp, path)


def _savez_atomic(path, **arrays):
    # np.savez appends .npz to a path string; write to a file handle instead so
    # the staged temp keeps its exact name, then replace atomically.
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez(fh, **arrays)
    os.replace(tmp, path)


def load_marker_config(path):
    mappable, nonmappable, seen = [], [], []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            row = {k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
            seen.append(str(row.get("marker_number", "")).strip())
            if str(row.get("gene_level_mappable", "")).upper() == "TRUE":
                aliases = [a.strip() for a in (row.get("gene_aliases") or "").split(";") if a.strip()]
                mappable.append({"canonical": row["canonical_gene"], "aliases": aliases,
                                 "reported": row["reported_name"]})
            else:
                nonmappable.append({"reported": row.get("reported_name", ""),
                                    "status": row.get("mapping_status", "")})
    # fail closed if the ledger drifts from the 22 Erho Table 2 rows plus the legacy exclusion
    expected = [str(i) for i in range(1, 23)] + ["legacy_excluded"]
    if seen != expected:
        raise ValueError(f"marker config rows drifted from locked ledger; got {seen}")
    return mappable, nonmappable


def resolve_panel(mappable, train_genes):
    """Fitted panel = unambiguous canonical-gene intersection with training only.
    Deterministic: for each mappable marker try canonical then aliases; first hit
    in training becomes the panel symbol; repeated canonical mappings kept once."""
    tg = set(train_genes)
    panel, excluded = [], []
    for m in mappable:
        hit = next((c for c in [m["canonical"]] + m["aliases"] if c in tg), None)
        if hit is None:
            excluded.append(m["canonical"])
        elif hit not in panel:
            panel.append(hit)
    return panel, excluded


def quantile_normalize_between_arrays(mat):
    """limma::normalizeBetweenArrays(method='quantile') equivalent."""
    mat = np.asarray(mat, dtype=float)
    n_genes, n_samples = mat.shape
    rank_mean = np.sort(mat, axis=0).mean(axis=1)
    out = np.empty_like(mat)
    for j in range(n_samples):
        ranks = pd.Series(mat[:, j]).rank(method="average").values
        out[:, j] = np.interp(ranks - 1, np.arange(n_genes), rank_mean)
    return out


def _fit_lr(X, y, C):
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning
    lr = LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=MAX_ITER,
                            tol=1e-10, fit_intercept=True)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)   # convergence failure is fatal
        lr.fit(X, y)
    if int(np.max(lr.n_iter_)) >= MAX_ITER:
        raise ValueError(f"solver hit max_iter at C={C}")
    return lr


def _bridge_rows(cohort_rows, ref_rows):
    """Rank-map each panel row of a cohort batch to the corresponding reference
    row (fold-training or full-training QN), mirroring cohort-batch deployment."""
    n = cohort_rows.shape[1]
    out = np.empty_like(cohort_rows, dtype=float)
    for i in range(cohort_rows.shape[0]):
        train_sorted = np.sort(ref_rows[i, :])
        ranks = pd.Series(cohort_rows[i, :]).rank(method="average").values
        out[i, :] = np.quantile(train_sorted, ranks / (n + 1), method="linear")
    return out


def tune_C(train_mat, y, tidx, panel):
    """Deterministic 5-fold tuning with fold-specific preprocessing (no leakage):
    within each fold, between-array QN is fit on that fold's training patients
    only and the held-out fold is rank-mapped as a cohort batch to that
    reference. Returns the grid, per-fold AUCs, mean AUC, and the selected
    interior C (fails if the optimum is on a grid boundary)."""
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    pidx = [tidx[g] for g in panel]
    skf = StratifiedKFold(CV_FOLDS, shuffle=True, random_state=CV_SEED)
    fold_auc = np.full((len(CS_GRID), CV_FOLDS), np.nan)
    for fi, (tr, te) in enumerate(skf.split(train_mat.T, y)):
        ref_panel = quantile_normalize_between_arrays(train_mat[:, tr])[pidx, :]  # fold-train only
        Xtr = ref_panel.T
        Xte = _bridge_rows(train_mat[pidx][:, te], ref_panel).T
        for ci, C in enumerate(CS_GRID):
            lr = _fit_lr(Xtr, y[tr], C)
            p = _sigmoid(Xte @ lr.coef_.ravel() + float(lr.intercept_[0]))
            fold_auc[ci, fi] = roc_auc_score(y[te], p)
    mean_auc = fold_auc.mean(axis=1)
    sel = int(np.argmax(mean_auc))                 # argmax mean AUC; lowest C on ties
    if sel in (0, len(CS_GRID) - 1):
        raise ValueError(f"tuning optimum on grid boundary (idx {sel}, C={CS_GRID[sel]:.3g})")
    return dict(grid=CS_GRID, fold_auc=fold_auc, mean_auc=mean_auc, sel_idx=sel,
                sel_C=float(CS_GRID[sel]))


def artifact_score(art, cohort_mat, cohort_genes):
    """Reconstruct cohort probabilities from the stored artifact only: per-gene
    rank map to the stored training reference, impute absent panel genes to the
    stored per-gene reference mean, then apply the stored logistic model."""
    panel = [str(g) for g in art["panel"]]
    coef = np.asarray(art["coef"], dtype=float).ravel()
    b = float(np.asarray(art["intercept"], dtype=float).ravel()[0])
    ref = np.asarray(art["train_qn_ref"], dtype=float)      # n_panel x n_train
    ref_mean = np.asarray(art["train_ref_mean"], dtype=float)
    lut = {g: i for i, g in enumerate(cohort_genes)}
    n = cohort_mat.shape[1]
    X = np.empty((n, len(panel)), dtype=float)
    imputed = []
    for j, g in enumerate(panel):
        ci = lut.get(g)
        if ci is None:
            X[:, j] = ref_mean[j]
            imputed.append(g)
        else:
            train_sorted = np.sort(ref[j, :])
            ranks = pd.Series(cohort_mat[ci, :]).rank(method="average").values
            X[:, j] = np.quantile(train_sorted, ranks / (n + 1), method="linear")
    return _sigmoid(X @ coef + b), imputed


def load_durham_555(clin_ids):
    # allow_pickle: the NPZ is produced by our Stage 1 and stores object-dtype
    # gene/sample arrays; it is a local, trusted artifact, not external input.
    d = np.load(DUR_NPZ, allow_pickle=True)
    dur_mat = d["mat"].astype(float)
    dur_genes = list(map(str, d["genes"]))
    dur_samples = list(map(str, d["samples"]))
    valid = set(clin_ids)
    keep_idx = [i for i, s in enumerate(dur_samples) if s in valid]
    kept = [dur_samples[i] for i in keep_idx]
    # exact clinical-cohort equality (555 before rank mapping)
    if set(kept) != valid:
        raise ValueError("Durham expression samples do not exactly cover the clin_valid IDs")
    if len(kept) != N_DURHAM or len(set(kept)) != N_DURHAM:
        raise ValueError(f"Durham analysis cohort must be {N_DURHAM} unique IDs, got {len(kept)}")
    return dur_mat[:, keep_idx], dur_genes, kept


def load_jhu():
    md = rdata.read_rda(JHU_RDA)
    test_genes = list(md["testMat"].coords[md["testMat"].dims[0]].values)
    test_mat = md["testMat"].values
    test_samples = list(md["testMat"].coords[md["testMat"].dims[1]].values)
    short = [_short_jhu(s) for s in test_samples]
    if len(short) != len(set(short)):
        raise ValueError("JHU short-ID collision")
    if len(short) != N_JHU:
        raise ValueError(f"JHU cohort must be {N_JHU}, got {len(short)}")
    return test_mat, test_genes, short


def load_clin_ids():
    clin = rdata.read_rda(DUR_CLIN)["clin_valid"]
    ids = [str(s) for s in clin["sample_id"].values]
    if len(ids) != N_DURHAM or len(set(ids)) != N_DURHAM:
        raise ValueError(f"clin_valid must be {N_DURHAM} unique IDs, got {len(ids)}/{len(set(ids))}")
    return ids


# ---- fit mode ------------------------------------------------------------
def run_fit():
    from sklearn.metrics import roc_curve
    import sklearn

    mappable, nonmappable = load_marker_config(CFG)
    print("Loading training + JHU + Durham …", flush=True)
    md = rdata.read_rda(JHU_RDA)
    train_genes = list(md["trainMat"].coords[md["trainMat"].dims[0]].values)
    train_mat = md["trainMat"].values
    y_train = np.asarray(md["trainGroup"].codes).astype(int)
    if train_mat.shape[1] != N_TRAIN:
        raise ValueError(f"training must be {N_TRAIN}, got {train_mat.shape[1]}")
    if int((y_train == 0).sum()) != CLASS_NO_METS or int((y_train == 1).sum()) != CLASS_METS:
        raise ValueError("training class counts != expected 694/306")
    if not np.all(np.isfinite(train_mat)):
        raise ValueError("nonfinite training expression")

    test_mat, test_genes, jhu_ids = load_jhu()
    clin_ids = load_clin_ids()
    dur_mat, dur_genes, dur_ids = load_durham_555(clin_ids)
    if not (np.all(np.isfinite(test_mat)) and np.all(np.isfinite(dur_mat))):
        raise ValueError("nonfinite validation expression")

    panel, excluded_absent = resolve_panel(mappable, train_genes)
    candidate_genes = [m["canonical"] for m in mappable]
    if len(set(panel)) != len(panel):
        raise ValueError("duplicated panel gene after resolution")
    print(f"Candidate gene-mappable markers: {len(candidate_genes)}; fitted panel "
          f"(training only): {len(panel)}")
    print(f"  Panel: {panel}")
    print(f"  Mappable but absent from training (excluded): {excluded_absent}")

    print("Quantile-normalising training (between-arrays) …", flush=True)
    t0 = time.time()
    train_qn = quantile_normalize_between_arrays(train_mat)
    print(f"  done in {time.time()-t0:.1f}s")
    tidx = {g: i for i, g in enumerate(train_genes)}
    X_train = train_qn[[tidx[g] for g in panel], :].T
    train_qn_ref = train_qn[[tidx[g] for g in panel], :]
    train_ref_mean = train_qn_ref.mean(axis=1)

    print("Tuning C: deterministic 5-fold, fold-specific QN (leakage-free) …", flush=True)
    tune = tune_C(train_mat, y_train, tidx, panel)
    best_C = tune["sel_C"]
    cv_auc_tuning = float(tune["mean_auc"][tune["sel_idx"]])
    print(f"  Selected C = {best_C:.10g} (grid idx {tune['sel_idx']}/{len(CS_GRID)-1}, interior)")
    print(f"  Mean fold AUC at selected C (tuning diagnostic only): {cv_auc_tuning:.7f}")
    print("Final L2 fit on all 1000 fully-normalized training patients …", flush=True)
    lr = _fit_lr(X_train, y_train, best_C)
    print(f"  Final fit n_iter_ (max) = {int(np.max(lr.n_iter_))} of max_iter = {MAX_ITER}")
    coef = lr.coef_.ravel()
    intercept = float(lr.intercept_[0])

    p_tr = lr.predict_proba(X_train)[:, 1]
    fpr, tpr, thrs = roc_curve(y_train, p_tr)
    threshold = float(thrs[np.argmax(tpr - fpr)])
    if not (0.0 <= threshold <= 1.0):
        raise ValueError("threshold out of [0,1]")
    print(f"  Secondary categorization threshold (Youden J): {threshold:.10g}")

    versions = {"python": sys.version.split()[0], "numpy": np.__version__,
                "pandas": pd.__version__, "sklearn": sklearn.__version__}
    input_hashes = {"MetastasisData_JHUOut.rda": _sha256(JHU_RDA),
                    "_durham_full_expr.npz": _sha256(DUR_NPZ),
                    "durham_metscore_batchcorrected.rda": _sha256(DUR_CLIN),
                    "config": _sha256(CFG)}

    # ---- serialize the artifact atomically ------------------------------
    _savez_atomic(
        ART_NPZ, panel=np.array(panel, dtype=object), coef=coef,
        intercept=np.array([intercept]), C=np.array([best_C]),
        cv_seed=np.array([CV_SEED]), n_folds=np.array([CV_FOLDS]),
        threshold=np.array([threshold]), train_qn_ref=train_qn_ref,
        train_ref_mean=train_ref_mean,
        cs_grid=tune["grid"], fold_auc=tune["fold_auc"], mean_fold_auc=tune["mean_auc"],
        sel_idx=np.array([tune["sel_idx"]]),
        candidate_markers=np.array(candidate_genes, dtype=object),
        excluded_absent_markers=np.array(excluded_absent, dtype=object))
    meta = {
        "model": "locally_trained_published_decipher_marker_gene_set_surrogate",
        "not_licensed_commercial_decipher_gc": True,
        "not_a_reproduction_of_erho_random_forest": True,
        "marker_source": "Erho et al. 2013 PLoS ONE Table 2 (DOI 10.1371/journal.pone.0066855)",
        "panel": panel, "n_panel": len(panel),
        "candidate_gene_mappable_markers": candidate_genes,
        "excluded_absent_from_training": excluded_absent,
        "nonmappable_markers": nonmappable,
        "coef": coef.tolist(), "intercept": intercept, "C": best_C,
        "cv_seed": CV_SEED, "n_folds": CV_FOLDS, "fold_seed": CV_SEED,
        "cs_grid": tune["grid"].tolist(), "fold_auc": tune["fold_auc"].tolist(),
        "mean_fold_auc": tune["mean_auc"].tolist(), "selected_index": tune["sel_idx"],
        "selection_rule": "argmax mean 5-fold AUC (lowest C on ties); fail if optimum on grid boundary",
        "preprocessing_description": "fold-specific between-array QN on fold-training only; held-out fold rank-mapped as a cohort batch; final fit on all-1000 full-transcriptome QN",
        "cv_auc_tuning_diagnostic_only": cv_auc_tuning,
        "threshold": threshold, "threshold_role": "secondary fixed categorization rule",
        "preprocessing_version": PREPROC_VERSION,
        "validation_transform": "transductive cohort-batch rank mapping to training reference",
        "cohort_batch_limitation": "retrospective cohort-batch scoring only; not individual deployment",
        "sample_counts": {"training": N_TRAIN, "JHU": N_JHU, "Durham": N_DURHAM},
        "class_counts_training": {"No_Mets": CLASS_NO_METS, "Mets": CLASS_METS},
        "input_hashes": input_hashes, "software_versions": versions,
    }
    _atomic_write_bytes(ART_JSON, lambda p: open(p, "w").write(json.dumps(meta, indent=2)))
    # full-precision repr so NPZ/JSON/text thresholds agree exactly
    _atomic_write_bytes(THR_TXT, lambda p: open(p, "w").write(f"{threshold!r}\n"))
    print("  Artifact + metadata + threshold written.")

    # stored-coefficient round trip: reloaded coef reproduce live model on X_train.
    # allow_pickle: artifact NPZ is written by this script (local, trusted).
    art = dict(np.load(ART_NPZ, allow_pickle=True))
    p_stored = _sigmoid(X_train @ np.asarray(art["coef"], float).ravel()
                        + float(np.asarray(art["intercept"], float).ravel()[0]))
    rt = float(np.max(np.abs(p_stored - p_tr)))
    print(f"  Stored-coefficient round trip max abs error (training): {rt:.3e}")
    if rt > 1e-12:
        raise ValueError("stored-coefficient round trip exceeds 1e-12")

    # ---- predictions from the RELOADED artifact (not the live lr) -------
    p_jhu, imp_jhu = artifact_score(art, test_mat, test_genes)
    p_dur, imp_dur = artifact_score(art, dur_mat, dur_genes)
    for name, p in (("JHU", p_jhu), ("Durham", p_dur)):
        if not (np.all(np.isfinite(p)) and p.min() >= 0.0 and p.max() <= 1.0):
            raise ValueError(f"{name} probabilities not finite in [0,1]")
    _atomic_write_bytes(JHU_PRED, lambda p: pd.DataFrame(
        {"sample_id": jhu_ids, "decipher_surrogate_prob": p_jhu}).to_csv(p, index=False))
    _atomic_write_bytes(DUR_PRED, lambda p: pd.DataFrame(
        {"sample_id": dur_ids, "decipher_surrogate_prob": p_dur}).to_csv(p, index=False))
    print(f"  JHU n={len(jhu_ids)} (imputed {len(imp_jhu)}); Durham n={len(dur_ids)} (imputed {len(imp_dur)})")

    # ---- manifest last ---------------------------------------------------
    manifest = {"artifact_npz": _sha256(ART_NPZ), "artifact_json": _sha256(ART_JSON),
                "threshold_txt": _sha256(THR_TXT), "jhu_pred": _sha256(JHU_PRED),
                "durham_pred": _sha256(DUR_PRED)}
    _atomic_write_bytes(MANIFEST, lambda p: open(p, "w").write(json.dumps(manifest, indent=2)))
    print("  Manifest written.")
    print(f"FIT DONE: panel={len(panel)} C={best_C:.6g} intercept={intercept:.6g} "
          f"threshold={threshold:.6g}")
    print(f"  JHU High/Low={int((p_jhu>threshold).sum())}/{int((p_jhu<=threshold).sum())}; "
          f"Durham High/Low={int((p_dur>threshold).sum())}/{int((p_dur<=threshold).sum())}")


# ---- verify mode (clean process, no fit) ---------------------------------
def run_verify():
    # clean-process, no-fit, fail-closed contract check.
    # 1) manifest hashes match current files
    manifest = json.load(open(MANIFEST))
    for key, path in {"artifact_npz": ART_NPZ, "artifact_json": ART_JSON, "threshold_txt": THR_TXT,
                      "jhu_pred": JHU_PRED, "durham_pred": DUR_PRED}.items():
        if _sha256(path) != manifest.get(key):
            raise ValueError(f"manifest hash mismatch: {key}")
    # 2) current raw-input and config hashes match the stored metadata
    meta = json.load(open(ART_JSON))
    cur = {"MetastasisData_JHUOut.rda": _sha256(JHU_RDA), "_durham_full_expr.npz": _sha256(DUR_NPZ),
           "durham_metscore_batchcorrected.rda": _sha256(DUR_CLIN), "config": _sha256(CFG)}
    for k, v in cur.items():
        if meta["input_hashes"].get(k) != v:
            raise ValueError(f"current input/config hash mismatch vs metadata: {k}")
    # 3) NPZ/JSON agreement + exact threshold-text agreement
    art = dict(np.load(ART_NPZ, allow_pickle=True))   # local, trusted artifact
    npz_panel = [str(g) for g in art["panel"]]
    if npz_panel != list(meta["panel"]):
        raise ValueError("panel order NPZ vs JSON mismatch")
    if not np.array_equal(np.asarray(art["coef"], float).ravel(), np.asarray(meta["coef"], float)):
        raise ValueError("coefficients NPZ vs JSON mismatch")
    thr_npz = float(np.asarray(art["threshold"], float).ravel()[0])
    thr_json = float(meta["threshold"])
    thr_txt = float(open(THR_TXT).read().strip())
    if not (thr_npz == thr_json == thr_txt):
        raise ValueError("threshold disagreement among NPZ / JSON / text")
    if float(np.asarray(art["intercept"], float).ravel()[0]) != float(meta["intercept"]):
        raise ValueError("intercept NPZ vs JSON mismatch")
    if float(np.asarray(art["C"], float).ravel()[0]) != float(meta["C"]):
        raise ValueError("C NPZ vs JSON mismatch")
    if int(np.asarray(art["cv_seed"]).ravel()[0]) != int(meta["cv_seed"]):
        raise ValueError("cv_seed NPZ vs JSON mismatch")
    if meta.get("preprocessing_version") != PREPROC_VERSION:
        raise ValueError("preprocessing metadata mismatch")
    ref = np.asarray(art["train_qn_ref"], float)
    if ref.shape != (len(npz_panel), N_TRAIN):
        raise ValueError(f"train_qn_ref shape {ref.shape} != ({len(npz_panel)}, {N_TRAIN})")
    if len(npz_panel) != len(set(npz_panel)):
        raise ValueError("duplicate panel gene in artifact")
    if not (0.0 <= thr_npz <= 1.0):
        raise ValueError("threshold out of [0,1]")
    # 4) cohorts: counts, unique IDs/features, finite inputs
    test_mat, test_genes, jhu_ids = load_jhu()
    clin_ids = load_clin_ids()
    dur_mat, dur_genes, dur_ids = load_durham_555(clin_ids)
    if len(jhu_ids) != N_JHU or len(dur_ids) != N_DURHAM:
        raise ValueError("unexpected cohort sample counts")
    if len(set(jhu_ids)) != N_JHU or len(set(dur_ids)) != N_DURHAM:
        raise ValueError("duplicate cohort sample IDs")
    if len(set(test_genes)) != len(test_genes) or len(set(dur_genes)) != len(dur_genes):
        raise ValueError("duplicate cohort features")
    if not (np.all(np.isfinite(test_mat)) and np.all(np.isfinite(dur_mat))):
        raise ValueError("nonfinite cohort expression")
    # 5) reconstruct scores (no fit) and compare to saved predictions
    p_jhu, _ = artifact_score(art, test_mat, test_genes)
    p_dur, _ = artifact_score(art, dur_mat, dur_genes)
    for nm, p in (("JHU", p_jhu), ("Durham", p_dur)):
        if not (np.all(np.isfinite(p)) and p.min() >= 0.0 and p.max() <= 1.0):
            raise ValueError(f"{nm} probabilities not finite in [0,1]")
    saved_j = pd.read_csv(JHU_PRED); saved_d = pd.read_csv(DUR_PRED)
    if list(saved_j["sample_id"].astype(str)) != list(map(str, jhu_ids)):
        raise ValueError("JHU sample-ID order mismatch vs saved predictions")
    if list(saved_d["sample_id"].astype(str)) != list(map(str, dur_ids)):
        raise ValueError("Durham sample-ID order mismatch vs saved predictions")
    dj = float(np.max(np.abs(p_jhu - saved_j["decipher_surrogate_prob"].values)))
    dd = float(np.max(np.abs(p_dur - saved_d["decipher_surrogate_prob"].values)))
    cj = int(np.sum((p_jhu > thr_npz) != (saved_j["decipher_surrogate_prob"].values > thr_npz)))
    cd = int(np.sum((p_dur > thr_npz) != (saved_d["decipher_surrogate_prob"].values > thr_npz)))
    print("VERIFY: manifest hashes OK; input/config hashes OK; NPZ/JSON/text agree; "
          f"counts {N_JHU}/{N_DURHAM}")
    print(f"VERIFY: JHU max|diff|={dj:.3e} class_mismatch={cj}; "
          f"Durham max|diff|={dd:.3e} class_mismatch={cd}")
    print(f"VERIFY summary hash JHU={hashlib.sha256(np.round(p_jhu,12).tobytes()).hexdigest()[:16]} "
          f"Durham={hashlib.sha256(np.round(p_dur,12).tobytes()).hexdigest()[:16]}")
    if max(dj, dd) > 1e-12 or (cj + cd) > 0:
        raise ValueError("artifact-only reconstruction disagrees with saved predictions")
    print("VERIFY PASS (fail-closed contract)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["fit", "verify"], default="fit")
    args = ap.parse_args()
    if args.mode == "fit":
        run_fit()
    else:
        run_verify()


if __name__ == "__main__":
    main()

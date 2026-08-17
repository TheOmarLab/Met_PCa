"""
Figure 3: Gleason score 7 (GS7) locked Met-Score validation, plotting consumer.

Panels:
  (a) locked-probability distribution by ever-metastasis status (descriptive,
      censoring-naive), with Benjamini-Hochberg q from the producer.
  (b) locked-probability ROC (descriptive, censoring-naive), curve and AUC/95% CI
      from the producer. Untitled.
  (c) adjusted Cox summary associations over observed follow-up: All-GS7 primary
      continuous-score HR and exploratory GG2/GG3 HRs (per full-cohort SD of the
      locked Met-Score), with Lin-Wei robust 95% CIs. No subgroup significance
      stars. Proportional-hazards diagnostics and the robust-Wald Met-Score x
      Grade-Group interaction are in GS7_model_diagnostics.csv / GS7_cohort_summary.csv
      and belong in the figure caption, not on the plot.
  (d) GS7 incremental cause-specific concordance: Grade-Group-only vs
      Grade-Group + continuous Met-Score C-index with 95% CIs (one dumbbell per
      cohort) and the paired delta-C, read verbatim from
      GS7_incremental_concordance.csv. Binary-class rows are not plotted.

load_and_verify still asserts the complete six-fit diagnostics (present,
converged, estimable, finite PH/EPV/dfbeta) and that each robust CI brackets its
point HR, even though those values are not drawn.

This script fits no model and computes no inferential p-value or confidence
interval. It reconstructs the GS7 ledger only to assert the producer's ledger
hash, n/events, subgroup counts, and point AUCs, then reads all reported values
from outs/Figure3/GS7_*.csv (produced by
code/survival_analysis/Met_PCa_Survival_Multivariate.R). The locked threshold is
read from config. Values describe association and descriptive discrimination in a
GS7 subgroup; they are not a subgroup-specific validation, equivalence, clinical-
utility, or recalibration claim.
"""

import os
import sys
import hashlib


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
_CODE_DIR = os.path.join(ROOT, "code")
if os.path.isdir(_CODE_DIR) and _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.lines import Line2D
import rdata

OUT_BASE   = os.path.join(ROOT, "figures")
FIG3_DIR   = os.path.join(ROOT, "outs/Figure3")
SUMMARY    = os.path.join(FIG3_DIR, "GS7_cohort_summary.csv")
ASSOC      = os.path.join(FIG3_DIR, "GS7_adjusted_associations.csv")
DIAG_CSV   = os.path.join(FIG3_DIR, "GS7_model_diagnostics.csv")
ROC_CSV    = os.path.join(FIG3_DIR, "GS7_ROC_curves.csv")
CONC_CSV   = os.path.join(FIG3_DIR, "GS7_incremental_concordance.csv")
LOCKED_META = os.path.join(ROOT, "config/metscore_locked_v1_metadata.csv")

# CSV cohort label -> figure key; expected GS7 counts (structural asserts).
CSV_COHORT = {"JHU": "JHU", "Durham VA": "Durham"}
EXPECT = {  # n, events, GG2 n/ev, GG3 n/ev
    "JHU":    dict(n=132, events=23, gg2_n=86, gg2_e=12, gg3_n=46, gg3_e=11),
    "Durham": dict(n=422, events=25, gg2_n=336, gg2_e=11, gg3_n=86, gg3_e=14),
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.labelsize": 9.5, "axes.titlesize": 10.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
    "axes.edgecolor": "#222222", "axes.labelcolor": "#000000",
    "xtick.color": "#222222", "ytick.color": "#222222", "text.color": "#000000",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "savefig.edgecolor": "white",
    "axes.unicode_minus": False,
})
COHORT_COLORS = {"JHU": "#0072B2", "Durham": "#009E73"}


def load_rda(p):
    return rdata.conversion.convert(rdata.parser.parse_file(p),
                                    default_encoding="utf-8", force_default_encoding=True)


def load_locked_threshold():
    meta = pd.read_csv(LOCKED_META, dtype=str)
    thr = float(meta.loc[meta["field"] == "threshold", "value"].iloc[0])
    if not (0.0 <= thr <= 1.0):
        raise ValueError(f"locked threshold out of range: {thr}")
    return thr


THRESHOLD = load_locked_threshold()


def style_axis(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#222222"); ax.spines["bottom"].set_color("#222222")
    ax.grid(False)


# ---- ledger reconstruction (mirrors the R build_gs7; for asserts only) ----
def _gg_from(total, primary):
    total = pd.to_numeric(total, errors="coerce")
    primary = pd.to_numeric(primary, errors="coerce")
    gg = np.full(len(total), None, dtype=object)
    gg[total <= 6] = "GG1"
    gg[(total == 7) & (primary == 3)] = "GG2"
    gg[(total == 7) & (primary == 4)] = "GG3"
    gg[total == 8] = "GG4"
    gg[total >= 9] = "GG5"
    return gg


def _pt_collapse(x):
    s = pd.Series(x).astype(str).str.strip().str.upper()
    out = np.full(len(s), None, dtype=object)
    out[s.str.startswith("T2")] = "T2"
    out[s.str.startswith("T3")] = "T3"
    out[s.str.startswith("T4")] = "T4"
    return out


def reconstruct_ledger():
    """Rebuild the common complete-case GS7 ledger per cohort from the
    scored-object primary variables (locked probability, GG, log2PSA, pT, time, event)."""
    jhu = load_rda(os.path.join(ROOT, "outs/coxdata.rda"))["CoxData_jhu"]
    dj = pd.DataFrame({
        "GG": _gg_from(jhu["Pathological GS"], jhu["pathgs_p"]),
        "log2PSA": np.log2(pd.to_numeric(jhu["preop_psa"], errors="coerce") + 1),
        "pT": _pt_collapse(jhu["pstage"]),
        "event": pd.to_numeric(jhu["met"], errors="coerce"),
        "time": pd.to_numeric(jhu["met_time"], errors="coerce"),
        "ms_prob": pd.to_numeric(jhu["Met-Score prob"], errors="coerce"),
    })
    dur = load_rda(os.path.join(ROOT, "output/Durham/durham_metscore_batchcorrected.rda"))["clin_valid"]
    dd = pd.DataFrame({
        "GG": _gg_from(dur["PathGleason"], dur["pogl1"]),
        "log2PSA": np.log2(pd.to_numeric(dur["psapresurg"], errors="coerce") + 1),
        "pT": _pt_collapse(dur["stg"]),
        "event": pd.to_numeric(dur["mets"], errors="coerce"),
        "time": pd.to_numeric(dur["surgmets"], errors="coerce"),
        "ms_prob": pd.to_numeric(dur["MetScore_prob"], errors="coerce"),
    })

    def gs7(df):
        d = df[df["GG"].isin(["GG2", "GG3"])].copy()
        d = d.dropna(subset=["GG", "log2PSA", "pT", "ms_prob"])
        d["event"] = d["event"].astype(int)
        return d.reset_index(drop=True)

    return {"JHU": gs7(dj), "Durham": gs7(dd)}


def ledger_hash(d):
    key = [f"{g}|{p}|{int(e)}|{lp:.6f}|{ms:.6f}|{t:.4f}"
           for g, p, e, lp, ms, t in zip(d["GG"], d["pT"], d["event"],
                                         d["log2PSA"], d["ms_prob"], d["time"])]
    return hashlib.sha256("\n".join(sorted(key)).encode()).hexdigest()


def emp_auc(score, y):
    r = pd.Series(score).rank().values
    n1 = int((y == 1).sum()); n0 = int((y == 0).sum())
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


# displayed continuous fits that must carry complete diagnostics
DISPLAYED = [("JHU", "All-GS7"), ("JHU", "GG2"), ("JHU", "GG3"),
             ("Durham VA", "All-GS7"), ("Durham VA", "GG2"), ("Durham VA", "GG3")]
DIAG_COLS = ["n", "events", "n_coef", "events_per_coef", "converged", "ms_z_estimable",
             "ph_global_p", "ph_ms_z_p", "dfbeta_ms_z_max_abs", "n_warnings", "warning_disposition"]


# ---- producer outputs + fail-closed asserts ----
def load_and_verify():
    for f in (SUMMARY, ASSOC, DIAG_CSV, ROC_CSV):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"missing producer output: {f}")
    summ = pd.read_csv(SUMMARY)
    assoc = pd.read_csv(ASSOC)
    diag = pd.read_csv(DIAG_CSV)
    roc = pd.read_csv(ROC_CSV)
    if summ.duplicated(subset=["cohort"]).any():
        raise ValueError("duplicate cohort rows in GS7_cohort_summary.csv")
    # complete diagnostics for every displayed continuous fit (fail-closed)
    for col in DIAG_COLS:
        if col not in diag.columns:
            raise ValueError(f"GS7_model_diagnostics.csv missing column: {col}")
    for coh, strat in DISPLAYED:
        dr = diag[(diag["Cohort"] == coh) & (diag["Stratum"] == strat)]
        if len(dr) != 1:
            raise ValueError(f"diagnostics: expected one row for {coh}/{strat}, got {len(dr)}")
        dr = dr.iloc[0]
        if not (bool(dr["converged"]) and bool(dr["ms_z_estimable"])):
            raise ValueError(f"diagnostics: {coh}/{strat} not converged/estimable")
        if not np.isfinite(float(dr["events_per_coef"])):
            raise ValueError(f"diagnostics: {coh}/{strat} nonfinite events_per_coef")
        # PH / dfbeta may be NA when unavailable (case-cohort risk sets); require finite only when present
        for col in ("ph_global_p", "ph_ms_z_p", "dfbeta_ms_z_max_abs"):
            v = dr[col]
            if pd.notna(v) and not np.isfinite(float(v)):
                raise ValueError(f"diagnostics: {coh}/{strat} nonfinite {col}")
    # robust CIs present and bracket the point HR for every displayed fit
    for coh, strat in DISPLAYED:
        ar = assoc[(assoc["Cohort"] == coh) & (assoc["Stratum"] == strat)]
        if len(ar) != 1:
            raise ValueError(f"associations: expected one row for {coh}/{strat}, got {len(ar)}")
        ar = ar.iloc[0]
        hr, lo, hi = float(ar["HR"]), float(ar["CI_lo_robust"]), float(ar["CI_hi_robust"])
        if not (np.isfinite(hr) and np.isfinite(lo) and np.isfinite(hi) and lo <= hr <= hi):
            raise ValueError(f"associations: {coh}/{strat} robust CI does not bracket HR")
    led = reconstruct_ledger()
    out = {}
    for csv_name, key in CSV_COHORT.items():
        s = summ[summ["cohort"] == csv_name]
        if len(s) != 1:
            raise ValueError(f"{csv_name}: expected one summary row, got {len(s)}")
        s = s.iloc[0]
        d = led[key]; ex = EXPECT[key]
        n, ev = len(d), int(d["event"].sum())
        g2, g3 = d[d["GG"] == "GG2"], d[d["GG"] == "GG3"]
        got = dict(n=n, events=ev, gg2_n=len(g2), gg2_e=int(g2["event"].sum()),
                   gg3_n=len(g3), gg3_e=int(g3["event"].sum()))
        if got != ex:
            raise ValueError(f"{key}: reconstructed GS7 counts {got} != expected {ex}")
        if any(got[k] != int(s[c]) for k, c in
               [("n", "n"), ("events", "events"), ("gg2_n", "GG2_n"), ("gg2_e", "GG2_events"),
                ("gg3_n", "GG3_n"), ("gg3_e", "GG3_events")]):
            raise ValueError(f"{key}: reconstructed counts disagree with producer summary")
        h = ledger_hash(d)
        if h != str(s["ledger_sha256"]):
            raise ValueError(f"{key}: ledger hash mismatch (renderer {h[:12]} vs producer {str(s['ledger_sha256'])[:12]})")
        auc_chk = emp_auc(d["ms_prob"].values, d["event"].values)
        if not np.isfinite(float(s["desc_AUC"])) or abs(auc_chk - float(s["desc_AUC"])) > 5e-4:
            raise ValueError(f"{key}: point AUC {auc_chk:.4f} != producer {s['desc_AUC']}")
        rc = roc[roc["Cohort"] == csv_name][["fpr", "tpr"]].to_numpy()
        if rc.shape[0] < 3 or not np.isfinite(rc).all():
            raise ValueError(f"{key}: ROC coordinates missing/nonfinite")
        out[key] = dict(summary=s, ledger=d, roc=rc)
    return out, assoc, diag


# ---- Panel a: locked-probability distribution by ever-metastasis status ----
def panel_a(ax, V):
    cohorts = ["JHU", "Durham"]
    positions_no, positions_yes = [], []
    width, pad, inner = 0.35, 1.0, 0.10
    for i in range(len(cohorts)):
        positions_no.append(i * pad - inner); positions_yes.append(i * pad + inner)
    allp = np.concatenate([V[c]["ledger"]["ms_prob"].values for c in cohorts])
    ymin, ymax = float(allp.min()), float(allp.max()); span = ymax - ymin
    for i, c in enumerate(cohorts):
        d = V[c]["ledger"]; z = d["ms_prob"].values; y = d["event"].values
        z0, z1 = z[y == 0], z[y == 1]
        v0 = ax.violinplot(z0, positions=[positions_no[i]], widths=width * 1.4,
                           showmeans=False, showmedians=False, showextrema=False)
        for pc in v0["bodies"]:
            vv = pc.get_paths()[0].vertices; vv[:, 0] = np.minimum(vv[:, 0], positions_no[i])
            pc.set_facecolor("#B6BDC6"); pc.set_edgecolor("none"); pc.set_alpha(0.55)
        v1 = ax.violinplot(z1, positions=[positions_yes[i]], widths=width * 1.4,
                           showmeans=False, showmedians=False, showextrema=False)
        for pc in v1["bodies"]:
            vv = pc.get_paths()[0].vertices; vv[:, 0] = np.maximum(vv[:, 0], positions_yes[i])
            pc.set_facecolor(COHORT_COLORS[c]); pc.set_edgecolor("none"); pc.set_alpha(0.65)
        bp = ax.boxplot([z0, z1], positions=[positions_no[i], positions_yes[i]], widths=0.085,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", linewidth=1.1),
                        boxprops=dict(linewidth=0.6), whiskerprops=dict(color="black", linewidth=0.6),
                        capprops=dict(color="black", linewidth=0.6))
        for b in bp["boxes"]:
            b.set_facecolor("#FFFFFF"); b.set_edgecolor("black")
        q = float(V[c]["summary"]["MW_BH_q"])
        lab = f"q = {q:.1e}".replace("e-0", "e-") if q < 1e-3 else f"q = {q:.3g}"
        ax.text((positions_no[i] + positions_yes[i]) / 2, ymax + 0.05 * span, lab,
                ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks([i * pad for i in range(len(cohorts))])
    ax.set_xticklabels(["JHU Nat. History", "Durham VA"], fontsize=9.0)
    ax.set_ylabel("Met-Score probability", fontsize=9.5)
    ax.axhline(THRESHOLD, color="#BBBBBB", linewidth=0.5, linestyle=(0, (1, 2)), zorder=0)
    ax.set_ylim(ymin - 0.05 * span, ymax + 0.18 * span)
    style_axis(ax)
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#B6BDC6", alpha=0.7, ec="none"),
               plt.Rectangle((0, 0), 1, 1, fc="#444444", alpha=0.7, ec="none")]
    ax.legend(handles, ["No metastasis", "Metastasis"], loc="upper center",
              bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False, fontsize=8.5,
              handlelength=1.0, handletextpad=0.45, columnspacing=1.6, borderaxespad=0.0)


# ---- Panel b: descriptive ROC (producer coordinates + AUC/CI) ----
def panel_b(ax, V):
    for c in ["JHU", "Durham"]:
        rc = V[c]["roc"]; s = V[c]["summary"]
        order = np.argsort(rc[:, 0], kind="stable")
        ax.step(rc[order, 0], rc[order, 1], where="post", color=COHORT_COLORS[c],
                linewidth=2.0, solid_capstyle="round",
                label=f"{'JHU' if c=='JHU' else 'Durham'}  {float(s['desc_AUC']):.2f} "
                      f"({float(s['AUC_lo']):.2f}–{float(s['AUC_hi']):.2f})")
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=0.8, color="#BBBBBB")
    ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks(np.arange(0, 1.001, 0.2)); ax.set_yticks(np.arange(0, 1.001, 0.2))
    ax.set_xlabel("False positive rate (1 – specificity)", fontsize=9.5)
    ax.set_ylabel("True positive rate (sensitivity)", fontsize=9.5)
    style_axis(ax)
    leg = ax.legend(title="AUC (95% CI)", loc="lower right", bbox_to_anchor=(0.99, 0.05),
                    frameon=False, fontsize=8.5, title_fontsize=8.5,
                    handlelength=1.2, handletextpad=0.45, labelspacing=0.3, borderaxespad=0)
    leg.get_title().set_fontweight("bold")


# ---- Panel c: adjusted continuous-score associations (robust CI) ----
def panel_c(ax, assoc):
    order = [("JHU", "All-GS7"), ("JHU", "GG2"), ("JHU", "GG3"),
             ("Durham VA", "All-GS7"), ("Durham VA", "GG2"), ("Durham VA", "GG3")]
    rowlab = {"All-GS7": "All GS7", "GG2": "GG2 (3+4)", "GG3": "GG3 (4+3)"}
    yvals = [5, 4, 3, 2, 1, 0]
    yticklabels, xmax = [], 1.0
    for (coh, strat), y in zip(order, yvals):
        r = assoc[(assoc["Cohort"] == coh) & (assoc["Stratum"] == strat)].iloc[0]
        # Lin-Wei robust (sandwich) CI is used uniformly in the forest
        hr, lo, hi = float(r["HR"]), float(r["CI_lo_robust"]), float(r["CI_hi_robust"])
        key = "JHU" if coh == "JHU" else "Durham"; col = COHORT_COLORS[key]
        primary = (strat == "All-GS7")
        ax.errorbar(hr, y, xerr=[[hr - lo], [hi - hr]], fmt="o", color=col, ecolor=col,
                    markersize=7 if primary else 5.5, markerfacecolor=col,
                    elinewidth=1.5 if primary else 1.1, capsize=3.5, zorder=3)
        yticklabels.append(f"{rowlab[strat]}\n(n={int(r['n'])}, {int(r['events'])} ev)")
        # HR (robust 95% CI) just above the right edge (upper CI) of each bar
        ax.text(hi, y + 0.18, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", ha="right", va="bottom",
                fontsize=9.0, clip_on=False)
        xmax = max(xmax, hi)
    ax.axvline(1.0, color="#888888", linewidth=0.9, linestyle="--", zorder=1)
    ax.axhline(2.5, color="#BBBBBB", linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log"); ax.set_xlim(0.5, max(8.0, xmax * 1.15))
    ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4, 8]))
    ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
    ax.set_yticks(yvals); ax.set_yticklabels(yticklabels, fontsize=8.0)
    ax.set_ylim(-0.6, 5.6)
    ax.set_xlabel("Cox summary association over observed follow-up "
                  "(per full-cohort SD of locked Met-Score)", fontsize=8.5)
    ax.spines["left"].set_visible(False)
    style_axis(ax); ax.tick_params(axis="y", length=0)
    # cohort named once (vertical), left of the y tick labels; fixed points offset
    # keeps the label clear of the tick text regardless of panel width
    for key, yc in [("JHU", 4.0), ("Durham", 1.0)]:
        ax.annotate("JHU" if key == "JHU" else "Durham",
                    xy=(0.0, yc), xycoords=ax.get_yaxis_transform(),
                    xytext=(-80, 0), textcoords="offset points",
                    rotation=90, va="center", ha="center", fontsize=10, fontweight="bold",
                    color=COHORT_COLORS[key], annotation_clip=False)


# ---- Panel d input: accepted GS7 incremental concordance (fail-closed) ----
CONC_NEED = ["cohort", "model", "score_form", "corrected_or_frozen_C", "C_lo", "C_hi",
             "dC_vs_GG_only", "dC_lo", "dC_hi", "boot_attempted", "boot_success", "boot_failed"]


def load_concordance():
    """Read the accepted concordance aggregate. Plots only the primary continuous
    raw-probability comparison; performs no computation of its own."""
    if not os.path.isfile(CONC_CSV):
        raise FileNotFoundError(f"missing concordance aggregate: {CONC_CSV}")
    df = pd.read_csv(CONC_CSV)
    miss = [c for c in CONC_NEED if c not in df.columns]
    if miss:
        raise ValueError(f"GS7_incremental_concordance.csv missing columns: {miss}")
    if df.duplicated(subset=["cohort", "model", "score_form"]).any():
        raise ValueError("GS7_incremental_concordance.csv has duplicate cohort/model rows")
    if (pd.to_numeric(df["boot_failed"], errors="coerce").fillna(-1) != 0).any():
        raise ValueError("GS7_incremental_concordance.csv reports bootstrap failures")
    if (pd.to_numeric(df["boot_success"], errors="coerce")
            != pd.to_numeric(df["boot_attempted"], errors="coerce")).any():
        raise ValueError("GS7_incremental_concordance.csv: bootstrap success != attempted")
    out = {}
    for csv_coh, key in CSV_COHORT.items():             # JHU, Durham VA -> JHU, Durham
        sub = df[df["cohort"] == csv_coh]
        gg = sub[sub["model"] == "Grade Group only"]
        ms = sub[(sub["model"] == "Grade Group + Met-Score") &
                 (sub["score_form"] == "continuous raw probability")]
        if len(gg) != 1 or len(ms) != 1:
            raise ValueError(f"concordance: expected one GG-only and one continuous row for {csv_coh}")
        gg, ms = gg.iloc[0], ms.iloc[0]
        vals = {
            "gg": (float(gg["corrected_or_frozen_C"]), float(gg["C_lo"]), float(gg["C_hi"])),
            "ms": (float(ms["corrected_or_frozen_C"]), float(ms["C_lo"]), float(ms["C_hi"])),
            "dC": (float(ms["dC_vs_GG_only"]), float(ms["dC_lo"]), float(ms["dC_hi"])),
        }
        for grp in ("gg", "ms"):
            c, lo, hi = vals[grp]
            if not (np.isfinite(c) and np.isfinite(lo) and np.isfinite(hi) and lo <= c <= hi):
                raise ValueError(f"concordance: {csv_coh} {grp} CI does not bracket C")
        d, dlo, dhi = vals["dC"]
        if not (np.isfinite(d) and np.isfinite(dlo) and np.isfinite(dhi) and dlo <= d <= dhi):
            raise ValueError(f"concordance: {csv_coh} delta-C CI malformed")
        out[key] = vals
    return out


def _fmt_signed(x):
    s = f"{abs(x):.3f}"
    return ("−" + s) if x < 0 else s


def fmt_dc(d, lo, hi):
    dc = ("+" + f"{d:.3f}") if d >= 0 else ("−" + f"{abs(d):.3f}")
    return f"ΔC {dc} ({_fmt_signed(lo)}–{_fmt_signed(hi)})"


# ---- Panel d: GS7 incremental cause-specific concordance (dumbbell) ----
def panel_d(ax, C):
    rows = [("JHU", 1.0), ("Durham", 0.0)]              # JHU top, Durham bottom
    for key, y in rows:
        col = COHORT_COLORS[key]
        gc, glo, ghi = C[key]["gg"]; mc, mlo, mhi = C[key]["ms"]; d, dlo, dhi = C[key]["dC"]
        ax.plot([gc, mc], [y, y], color=col, linewidth=2.2, solid_capstyle="round", zorder=2)
        # GG-only: open marker; GG + Met-Score: filled marker (fill distinguishes model)
        ax.errorbar(gc, y, xerr=[[gc - glo], [ghi - gc]], fmt="o", color=col, ecolor=col,
                    markersize=8, markerfacecolor="white", markeredgecolor=col,
                    markeredgewidth=1.7, elinewidth=1.3, capsize=3.2, zorder=3)
        ax.errorbar(mc, y, xerr=[[mc - mlo], [mhi - mc]], fmt="o", color=col, ecolor=col,
                    markersize=8, markerfacecolor=col, markeredgecolor=col,
                    markeredgewidth=1.7, elinewidth=1.3, capsize=3.2, zorder=4)
        ax.text((gc + mc) / 2.0, y + 0.28, fmt_dc(d, dlo, dhi), ha="center", va="bottom",
                fontsize=8.5, clip_on=False)
    ax.set_yticks([1.0, 0.0]); ax.set_yticklabels(["JHU", "Durham VA"], fontsize=9.5)
    for tl, key in zip(ax.get_yticklabels(), ["JHU", "Durham"]):
        tl.set_color(COHORT_COLORS[key]); tl.set_fontweight("bold")
    ax.set_ylim(-0.6, 1.78)
    ax.set_xlim(0.45, 0.98)
    ax.xaxis.set_major_locator(FixedLocator([0.5, 0.6, 0.7, 0.8, 0.9]))
    ax.set_xlabel("Cause-specific concordance", fontsize=9.5)
    style_axis(ax); ax.tick_params(axis="y", length=0)
    h_gg = Line2D([0], [0], marker="o", color="#555555", markerfacecolor="white",
                  markeredgecolor="#555555", markeredgewidth=1.7, markersize=8, linestyle="none")
    h_ms = Line2D([0], [0], marker="o", color="#555555", markerfacecolor="#555555",
                  markeredgecolor="#555555", markersize=8, linestyle="none")
    ax.legend([h_gg, h_ms], ["Grade Group", "Grade Group + Met-Score"],
              loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False,
              fontsize=8.5, handletextpad=0.4, columnspacing=1.1, borderaxespad=0.0)


def build():
    V, assoc, _diag = load_and_verify()   # _diag asserted inside load_and_verify
    C = load_concordance()
    # canvas widened modestly so the split bottom row (c + new panel d) is not cramped;
    # height and bottom margin trimmed together so the panels keep their size while the
    # excess whitespace below the x-axis labels is removed
    fig = plt.figure(figsize=(12.8, 7.62), dpi=600, facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[1.28, 1.0], hspace=0.50,
                             left=0.115, right=0.965, top=0.947, bottom=0.072)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.24)
    ax_a = fig.add_subplot(top[0, 0]); ax_b = fig.add_subplot(top[0, 1])
    panel_a(ax_a, V); panel_b(ax_b, V); ax_b.set_anchor("N")
    # bottom row: panel c ~62% width, panel d ~38%
    bottom = outer[1].subgridspec(1, 2, width_ratios=[0.62, 0.38], wspace=0.30)
    ax_c = fig.add_subplot(bottom[0, 0]); panel_c(ax_c, assoc)
    ax_d = fig.add_subplot(bottom[0, 1]); panel_d(ax_d, C)

    fig.canvas.draw()
    _lab = dict(fontsize=15, fontweight="black", family="DejaVu Sans",
                va="bottom", ha="left", color="#000000")
    # a/b share the top-row baseline; c/d share the bottom-row baseline
    y_top = ax_a.get_position().y1 + 0.012
    y_bot = max(ax_c.get_position().y1, ax_d.get_position().y1) + 0.012
    placements = ((ax_a, "a", y_top), (ax_b, "b", y_top),
                  (ax_c, "c", y_bot), (ax_d, "d", y_bot))
    for ax_, letter, ylab in placements:
        bb = ax_.get_position()
        fig.text(max(0.01, bb.x0 - 0.045), ylab, letter, **_lab)

    os.makedirs(OUT_BASE, exist_ok=True)
    out_pdf = os.path.join(OUT_BASE, "Figure3_GS7.pdf")
    out_tif = os.path.join(OUT_BASE, "Figure3_GS7.tiff")
    out_png = os.path.join(OUT_BASE, "Figure3_GS7.png")
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_tif, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_png, format="png", dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")
    print(f"Saved {out_tif}")


if __name__ == "__main__":
    build()

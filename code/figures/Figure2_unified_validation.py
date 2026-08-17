"""
Unified Figure 2 — locked Met-Score validation across two cohorts (JHU Nat.
History and Durham VA).

Layout (three rows, five panels):
  (a) locked-probability distribution by ever-metastasis status, both cohorts.
  (b) locked-probability ROC (descriptive, ever-metastasis label), both cohorts.
  (c) JHU  Aalen-Johansen metastasis cumulative incidence by the locked HIGH/LOW
      threshold, with pointwise CI and at-risk table.
  (d) Durham Aalen-Johansen metastasis cumulative incidence by the locked HIGH/LOW
      threshold, with pointwise CI and at-risk table.
  (e) fixed-score time-dependent AUC at 5 and 10 years (point + 95% CI), both cohorts.

The figure renders descriptive distributions/ROC and plots pre-computed
Aalen-Johansen curves; it fits no model. The threshold is read from the locked
config; the Aalen-Johansen curves/CIs/at-risk counts and the inferential values
(high-vs-low HR/CI/Cox p, 10-year HIGH-LOW CIF gap, time-dependent AUCs) are read
from outs/RiskStratification_LockedMetScore_AllCohorts.csv and
outs/RiskStratification_AJ_curves_AllCohorts.csv, produced by
code/data_preparation/Calibration_LockedLR_AllCohorts.R.

Design: Okabe-Ito palette (JHU #0072B2, Durham #009E73); Helvetica/DejaVu Sans;
only left/bottom spines; saved as vector PDF and 600-dpi LZW TIFF.
"""

import os
import sys

# Resolve project root: MET_PCA_ROOT env var, else walk up to the dir holding
# code/ and outs|output.
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

# Make _stats_utils.py importable regardless of cwd.
_CODE_DIR = os.path.join(ROOT, "code")
if os.path.isdir(_CODE_DIR) and _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import rdata
from sklearn.metrics import roc_curve
from scipy.stats import mannwhitneyu

OUT_BASE = os.path.join(ROOT, "figures")
SUMMARY_CSV = os.path.join(ROOT, "outs/RiskStratification_LockedMetScore_AllCohorts.csv")
AJ_CURVES_CSV = os.path.join(ROOT, "outs/RiskStratification_AJ_curves_AllCohorts.csv")
LOCKED_META = os.path.join(ROOT, "config/metscore_locked_v1_metadata.csv")

# Expected full-cohort sizes; the renderer fails closed if the loaded data or
# the producer CSV disagree with these.
EXPECT_N = {"JHU": 239, "Durham": 555}
# CSV cohort label -> figure cohort key.
CSV_COHORT = {"JHU": "JHU", "Durham VA": "Durham"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "axes.labelsize": 9.5,
    "axes.titlesize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "axes.edgecolor": "#000000",
    "axes.labelcolor": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
    "text.color": "#000000",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.unicode_minus": False,
})

COHORT_COLORS = {"JHU": "#0072B2", "Durham": "#009E73"}
RISK_COLOR = {"Low": "#4C78A8", "High": "#F58518"}


def load_rda(path):
    return rdata.conversion.convert(rdata.parser.parse_file(path),
                                    default_encoding="utf-8", force_default_encoding=True)


def load_locked_threshold():
    """Read the fixed decision threshold from the locked config (no hard-coding)."""
    meta = pd.read_csv(LOCKED_META, dtype=str)
    row = meta.loc[meta["field"] == "threshold", "value"]
    if row.empty:
        raise ValueError(f"threshold field missing from {LOCKED_META}")
    thr = float(row.iloc[0])
    if not np.isfinite(thr) or not (0.0 <= thr <= 1.0):
        raise ValueError(f"locked threshold out of range: {thr}")
    return thr


THRESHOLD = load_locked_threshold()


# ---------------- data ----------------
def get_cohort_data():
    """Per-cohort locked probability, ever-metastasis label, and MFS KM inputs.
    The locked LR probability is used, never a per-cohort refit."""
    out = {}

    jhu_pheno = load_rda(os.path.join(ROOT, "outs/jhu_pheno_filter_MetaScorer_Zscore.rda"))["pheno_jhu"]
    jhu_cox   = load_rda(os.path.join(ROOT, "outs/coxdata.rda"))["CoxData_jhu"]
    jhu_df = jhu_pheno[["Metastasis", "score", "sample_id"]].dropna(subset=["score"]).copy()
    jhu_df["y"] = (jhu_df["Metastasis"].astype(str).str.strip() == "Mets").astype(int)
    jhu_locked = jhu_cox[["sample_id", "Met-Score prob"]].rename(columns={"Met-Score prob": "pred_locked"})
    jhu_df = jhu_df.merge(jhu_locked, on="sample_id", how="left")
    jhu_surv = jhu_cox[["met", "met_time", "MetScoreClass"]].copy()
    jhu_surv = jhu_surv[jhu_surv["MetScoreClass"].notna() & jhu_surv["met"].notna() & jhu_surv["met_time"].notna()]
    jhu_surv["met"] = jhu_surv["met"].astype(float).astype(int)
    jhu_surv["met_time"] = jhu_surv["met_time"].astype(float)
    jhu_surv["risk"] = (jhu_surv["MetScoreClass"].astype(str)
                        .str.replace(" risk", "", regex=False).str.strip())
    out["JHU"] = {
        "y": jhu_df["y"].values,
        "pred_locked": pd.to_numeric(jhu_df["pred_locked"], errors="coerce").values,
        "label": "JHU Nat. History",
        "surv_t": jhu_surv["met_time"].values,
        "surv_e": jhu_surv["met"].values,
        "surv_g": jhu_surv["risk"].values,
    }

    dur_obj  = load_rda(os.path.join(ROOT, "output/Durham/durham_metscore_batchcorrected.rda"))
    dur_full = dur_obj["clin_valid"]
    # Drop pogl==0 sentinel rows (placeholder records) if the rda predates the
    # canonical fix in Durham_MetScore_Validation_BatchCorrected.R.
    if "pogl" in dur_full.columns:
        _pogl = pd.to_numeric(dur_full["pogl"], errors="coerce")
        _n_before = len(dur_full)
        dur_full = dur_full[~(_pogl.isna() | (_pogl == 0))].reset_index(drop=True)
        _n_dropped = _n_before - len(dur_full)
        if _n_dropped > 0:
            print(f"[Durham] dropped {_n_dropped} pogl==0 sentinel rows (of {_n_before}).")
    dur_df = dur_full[["mets", "MetScore_prob"]].dropna(subset=["mets", "MetScore_prob"]).copy()
    dur_df["mets"] = dur_df["mets"].astype(float).astype(int)
    dur_surv = dur_full[["mets", "surgmets", "MetScore_class"]].copy()
    dur_surv = dur_surv[dur_surv["MetScore_class"].notna() & dur_surv["mets"].notna() & dur_surv["surgmets"].notna()]
    dur_surv["mets"] = dur_surv["mets"].astype(float).astype(int)
    dur_surv["surgmets"] = dur_surv["surgmets"].astype(float)
    dur_surv["risk"] = dur_surv["MetScore_class"].astype(str)
    out["Durham"] = {
        "y": dur_df["mets"].values,
        "pred_locked": pd.to_numeric(dur_df["MetScore_prob"], errors="coerce").values,
        "label": "Durham VA",
        "surv_t": dur_surv["surgmets"].values,
        "surv_e": dur_surv["mets"].values,
        "surv_g": dur_surv["risk"].values,
    }
    return out


def read_summary(data):
    """Read producer-owned inferential values (fail closed). Returns
    {cohort: {hr, lo, hi, p, auc5(+ci), auc10(+ci), n_total, group counts}}."""
    if not os.path.isfile(SUMMARY_CSV):
        raise FileNotFoundError(f"missing risk-stratification CSV: {SUMMARY_CSV}")
    df = pd.read_csv(SUMMARY_CSV)
    need = ["Cohort", "RiskGroup", "N", "Events", "n_total", "HR_high_vs_low", "HR_lo",
            "HR_hi", "Cox_p", "AJ_diff_10y", "AJ_diff_lo", "AJ_diff_hi",
            "AUC_5y", "AUC_5y_lo", "AUC_5y_hi",
            "AUC_10y", "AUC_10y_lo", "AUC_10y_hi"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"risk-stratification CSV missing columns: {miss}")
    if df.duplicated(subset=["Cohort", "RiskGroup"]).any():
        raise ValueError("duplicate (Cohort, RiskGroup) rows in risk-stratification CSV")
    summ = {}
    for csv_name, key in CSV_COHORT.items():
        sub = df[df["Cohort"] == csv_name]
        if set(sub["RiskGroup"]) != {"Low", "High"}:
            raise ValueError(f"{csv_name}: expected exactly Low and High rows, got {sorted(sub['RiskGroup'])}")
        r = sub.iloc[0]
        vals = dict(
            hr=float(r["HR_high_vs_low"]), lo=float(r["HR_lo"]), hi=float(r["HR_hi"]), p=float(r["Cox_p"]),
            aj_diff=float(r["AJ_diff_10y"]), aj_diff_lo=float(r["AJ_diff_lo"]), aj_diff_hi=float(r["AJ_diff_hi"]),
            auc5=float(r["AUC_5y"]), auc5lo=float(r["AUC_5y_lo"]), auc5hi=float(r["AUC_5y_hi"]),
            auc10=float(r["AUC_10y"]), auc10lo=float(r["AUC_10y_lo"]), auc10hi=float(r["AUC_10y_hi"]),
            n_total=int(r["n_total"]),
            n_low=int(sub[sub["RiskGroup"] == "Low"]["N"].iloc[0]),
            n_high=int(sub[sub["RiskGroup"] == "High"]["N"].iloc[0]),
        )
        if not all(np.isfinite([vals["hr"], vals["lo"], vals["hi"], vals["p"],
                                vals["aj_diff"], vals["aj_diff_lo"], vals["aj_diff_hi"],
                                vals["auc5"], vals["auc5lo"], vals["auc5hi"],
                                vals["auc10"], vals["auc10lo"], vals["auc10hi"]])):
            raise ValueError(f"{csv_name}: non-finite inferential value in CSV")
        if vals["n_total"] != EXPECT_N[key]:
            raise ValueError(f"{csv_name}: CSV n_total {vals['n_total']} != expected {EXPECT_N[key]}")
        summ[key] = vals

    # Cross-check the loaded per-patient data against the producer counts.
    for key in ("JHU", "Durham"):
        g = np.asarray(data[key]["surv_g"])
        nlow = int(np.sum(g == "Low")); nhigh = int(np.sum(g == "High"))
        if nlow + nhigh != EXPECT_N[key]:
            raise ValueError(f"{key}: loaded KM n {nlow + nhigh} != expected {EXPECT_N[key]}")
        if (nlow, nhigh) != (summ[key]["n_low"], summ[key]["n_high"]):
            raise ValueError(f"{key}: loaded group sizes ({nlow},{nhigh}) != CSV "
                             f"({summ[key]['n_low']},{summ[key]['n_high']})")
        pl = np.asarray(data[key]["pred_locked"], dtype=float)
        yy = np.asarray(data[key]["y"])
        if not np.isfinite(pl[np.isfinite(pl)]).all() or np.sum(np.isfinite(pl)) == 0:
            raise ValueError(f"{key}: locked probability all non-finite")
        if not set(np.unique(yy[np.isin(yy, [0, 1])])).issubset({0, 1}):
            raise ValueError(f"{key}: ever-metastasis label not 0/1")
    return summ


def read_aj_curves():
    """Load the producer's Aalen-Johansen CIF curves (fail closed). Returns
    {cohort: {group: DataFrame(time, cif, cif_lo, cif_hi, n_at_risk)}}."""
    if not os.path.isfile(AJ_CURVES_CSV):
        raise FileNotFoundError(f"missing AJ curve CSV: {AJ_CURVES_CSV}")
    df = pd.read_csv(AJ_CURVES_CSV)
    need = ["Cohort", "RiskGroup", "time", "cif", "cif_lo", "cif_hi", "n_at_risk"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"AJ curve CSV missing columns: {miss}")
    curves = {}
    for csv_name, key in CSV_COHORT.items():
        curves[key] = {}
        for grp in ("Low", "High"):
            sub = df[(df["Cohort"] == csv_name) & (df["RiskGroup"] == grp)].sort_values("time")
            if sub.empty or not (sub[["time", "cif", "cif_lo", "cif_hi"]].to_numpy()
                                 == sub[["time", "cif", "cif_lo", "cif_hi"]].to_numpy()).all():
                raise ValueError(f"{csv_name}/{grp}: empty or non-finite AJ curve")
            if 120 not in set(sub["time"].astype(int)):
                raise ValueError(f"{csv_name}/{grp}: AJ curve missing the 120-month horizon")
            curves[key][grp] = sub.reset_index(drop=True)
    return curves


# ---------------- helpers ----------------
def style_axis(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#222222")
    ax.tick_params(colors="#222222")
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_color("#222222")
    ax.xaxis.label.set_color("#222222")
    ax.yaxis.label.set_color("#222222")
    ax.grid(False)


from utils._stats_utils import auc_with_ci as _auc_with_ci_shared


def auc_with_bootstrap_ci(z, y, n_boot=2000, seed=20260427, label=""):
    return _auc_with_ci_shared(z, y, n_boot=n_boot, seed=seed, label=label)


# ---------------- panel a: distribution by ever-metastasis status ----------------
def panel_a_distribution(ax, data):
    """Boxplot + half-violins of the locked probability per cohort, split by
    ever-metastasis status."""
    cohorts = ["JHU", "Durham"]
    positions_no, positions_yes = [], []
    width, pad, inner_off = 0.35, 1.0, 0.10
    for i, c in enumerate(cohorts):
        center = i * pad
        positions_no.append(center - inner_off)
        positions_yes.append(center + inner_off)

    def _pl(c):
        v = np.asarray(data[c]["pred_locked"], dtype=float)
        y = np.asarray(data[c]["y"], dtype=int)
        m = np.isfinite(v) & np.isin(y, [0, 1])
        return v[m], y[m]

    all_p = np.concatenate([_pl(c)[0] for c in cohorts])
    global_ymin, global_ymax = float(all_p.min()), float(all_p.max())
    global_span = global_ymax - global_ymin
    p_y_anchor = global_ymax + 0.05 * global_span
    ax_top = global_ymax + 0.18 * global_span

    p_strs = []
    for i, c in enumerate(cohorts):
        z, y = _pl(c)
        z0, z1 = z[y == 0], z[y == 1]

        v0 = ax.violinplot(z0, positions=[positions_no[i]], widths=width * 1.4,
                           showmeans=False, showmedians=False, showextrema=False)
        for pc in v0["bodies"]:
            verts = pc.get_paths()[0].vertices
            verts[:, 0] = np.minimum(verts[:, 0], positions_no[i])
            pc.set_facecolor("#B6BDC6"); pc.set_edgecolor("none"); pc.set_alpha(0.55)

        v1 = ax.violinplot(z1, positions=[positions_yes[i]], widths=width * 1.4,
                           showmeans=False, showmedians=False, showextrema=False)
        for pc in v1["bodies"]:
            verts = pc.get_paths()[0].vertices
            verts[:, 0] = np.maximum(verts[:, 0], positions_yes[i])
            pc.set_facecolor(COHORT_COLORS[c]); pc.set_edgecolor("none"); pc.set_alpha(0.65)

        bp = ax.boxplot(
            [z0, z1], positions=[positions_no[i], positions_yes[i]],
            widths=0.085, showfliers=False, patch_artist=True,
            medianprops=dict(color="black", linewidth=1.1),
            boxprops=dict(linewidth=0.6),
            whiskerprops=dict(color="black", linewidth=0.6),
            capprops=dict(color="black", linewidth=0.6),
        )
        for box in bp["boxes"]:
            box.set_facecolor("#FFFFFF"); box.set_edgecolor("black")

        p_strs.append((i, mannwhitneyu(z0, z1, alternative="two-sided").pvalue))

    # BH-correct (Benjamini-Hochberg) the cohort-level Wilcoxon p-values.
    raw_ps = np.array([p for _, p in p_strs], dtype=float)
    order, m_ = np.argsort(raw_ps), len(raw_ps)
    bh_ps = np.empty_like(raw_ps); cummin = 1.0
    for rank, idx in enumerate(order[::-1]):
        cummin = min(cummin, raw_ps[idx] * m_ / (m_ - rank))
        bh_ps[idx] = min(cummin, 1.0)
    for (i, _), p_adj in zip(p_strs, bh_ps):
        label = f"q = {p_adj:.1e}".replace("e-0", "e-") if p_adj < 1e-3 else f"q = {p_adj:.3g}"
        ax.text((positions_no[i] + positions_yes[i]) / 2, p_y_anchor, label,
                ha="center", va="bottom", fontsize=8.5, color="#000000")

    ax.set_xticks([i * pad for i in range(len(cohorts))])
    ax.set_xticklabels([data[c]["label"] for c in cohorts], fontsize=9.0)
    ax.set_xlabel("")
    ax.set_ylabel("Met-Score probability", fontsize=9.5)
    # Faint reference at the locked decision threshold, dichotomising HIGH/LOW.
    ax.axhline(THRESHOLD, color="#BBBBBB", linewidth=0.5, linestyle=(0, (1, 2)), zorder=0)
    ax.set_ylim(global_ymin - 0.05 * global_span, ax_top)
    style_axis(ax)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc="#B6BDC6", alpha=0.7, ec="none"),
        plt.Rectangle((0, 0), 1, 1, fc="#444444", alpha=0.7, ec="none"),
    ]
    ax.legend(handles, ["No metastasis", "Metastasis"],
              loc="upper center", bbox_to_anchor=(0.5, -0.10),
              ncol=2, frameon=False, fontsize=8.5,
              handlelength=1.0, handletextpad=0.45,
              columnspacing=1.6, borderaxespad=0.0)


# ---------------- panel b: descriptive ROC (ever-metastasis) ----------------
def panel_b_roc(ax, data):
    """ROC of the locked probability for ever-metastasis status (descriptive,
    not a time-to-event estimate), one curve per cohort."""
    for c in ["JHU", "Durham"]:
        z = np.asarray(data[c]["pred_locked"], dtype=float)
        y = np.asarray(data[c]["y"], dtype=int)
        m = np.isfinite(z) & np.isin(y, [0, 1])
        z, y = z[m], y[m]
        fpr, tpr, _ = roc_curve(y, z)
        auc, lo, hi = auc_with_bootstrap_ci(z, y)
        ax.plot(fpr, tpr, color=COHORT_COLORS[c], linewidth=2.0,
                label=f"{c}  {auc:.2f} ({lo:.2f}–{hi:.2f})", solid_capstyle="round")
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=0.8, color="#BBBBBB")
    ax.set_aspect("equal")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks(np.arange(0, 1.001, 0.2)); ax.set_yticks(np.arange(0, 1.001, 0.2))
    ax.set_xlabel("False positive rate (1 – specificity)", fontsize=9.5)
    ax.set_ylabel("True positive rate (sensitivity)", fontsize=9.5)
    style_axis(ax)
    leg = ax.legend(title="AUC (95% CI)", loc="lower right", bbox_to_anchor=(0.99, 0.05),
                    frameon=False, fontsize=8.5, title_fontsize=8.5,
                    handlelength=1.2, handletextpad=0.45, labelspacing=0.3, borderaxespad=0)
    leg.get_title().set_color("#000000"); leg.get_title().set_fontweight("bold")
    for txt in leg.get_texts():
        txt.set_color("#000000")


# ---------------- panels c/d: Aalen-Johansen CIF; curves/CIs/at-risk from the CSV ----------------
def panel_aj_cif(ax, curve, summ, cohort_label, cohort_color, *, ylabel=True, show_legend=True):
    """Aalen-Johansen metastasis cumulative incidence by locked HIGH/LOW, with a
    pointwise CI band and numbers at risk. Curves, CIs, at-risk counts, and the
    annotated statistics are all read from the producer CSV; nothing is estimated
    here and no 1-KM curve is drawn."""
    ymax = 0.0
    for grp, color in [("Low", RISK_COLOR["Low"]), ("High", RISK_COLOR["High"])]:
        sub = curve[grp]
        t   = sub["time"].to_numpy(dtype=float)
        cif = sub["cif"].to_numpy(dtype=float)
        clo = sub["cif_lo"].to_numpy(dtype=float)
        chi = sub["cif_hi"].to_numpy(dtype=float)
        ax.fill_between(t, clo, chi, step="post", color=color, alpha=0.15, linewidth=0)
        ax.step(t, cif, where="post", color=color, linewidth=1.7, solid_joinstyle="round")
        ymax = max(ymax, float(np.nanmax(chi)))

    ytop = min(1.0, np.ceil(ymax / 0.1) * 0.1 + 0.05)
    ax.set_title(cohort_label, fontsize=10, fontweight="bold", color=cohort_color, pad=6)
    ax.set_xlabel("Time (months)", fontsize=9.5)
    ax.set_ylabel("Metastasis cumulative incidence", fontsize=9.5) if ylabel else ax.set_ylabel("")
    ax.set_ylim(0, ytop)
    ax.set_xlim(0, 120)
    ax.set_xticks([0, 60, 120])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)

    # Producer-owned annotation: cause-specific HR and the 10-year HIGH-LOW CIF gap.
    p_str = "p<0.0001" if summ["p"] < 1e-4 else f"p={summ['p']:.3g}"
    txt = (f"Cause-specific HR {summ['hr']:.2f} ({summ['lo']:.2f}–{summ['hi']:.2f}), {p_str}\n"
           f"Δ10y {summ['aj_diff']*100:.1f}% "
           f"({summ['aj_diff_lo']*100:.1f}–{summ['aj_diff_hi']*100:.1f}%)")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="#000000")

    if show_legend:
        # same style as panel b: colour-line handles, no frame
        leg_handles = [plt.Line2D([0], [0], color=RISK_COLOR["Low"], lw=2.0),
                       plt.Line2D([0], [0], color=RISK_COLOR["High"], lw=2.0)]
        ax.legend(leg_handles, ["Low-risk", "High-risk"], loc="upper left",
                  bbox_to_anchor=(0.03, 0.82), frameon=False, fontsize=8.5,
                  handlelength=1.2, handletextpad=0.45, labelspacing=0.3, borderaxespad=0)

    # Numbers at risk: count columns sit beneath the 0 / 60 / 120-month ticks
    # (axis fractions 0.00 / 0.50 / 1.00), edge-aligned at the extremes so they
    # stay under the ticks without overflowing. The "No. at risk" header and the
    # group labels occupy a separate left position outside the count columns and
    # are drawn unclipped so they and the time-zero counts remain fully visible.
    def nar(grp):
        s = curve[grp]
        return [int(s.loc[s["time"] == tp, "n_at_risk"].iloc[0]) for tp in (0, 60, 120)]
    x_lab = -0.20
    cols = ((0.0, "left"), (0.5, "center"), (1.0, "right"))
    ax.text(x_lab, -0.20, "No. at risk", transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="#444444", style="italic", clip_on=False)
    for gi, (grp, name) in enumerate((("Low", "Low-risk"), ("High", "High-risk"))):
        yo = -0.30 - gi * 0.09
        color = RISK_COLOR[grp]
        ax.text(x_lab, yo, name, transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, color=color, fontweight="bold", clip_on=False)
        for (xf, cha), n in zip(cols, nar(grp)):
            ax.text(xf, yo, str(n), transform=ax.transAxes, ha=cha, va="top",
                    fontsize=7.5, color=color, clip_on=False)


# ---------------- panel e: fixed-score time-dependent AUC forest (5y, 10y) ----------------
def panel_e_td_auc(ax, summ):
    """Grouped horizontal forest of fixed-score time-dependent AUC (point + 95% CI):
    cohort named once (JHU top, Durham bottom), with 5- and 10-year rows and a
    dashed separator between cohorts. Values read from the CSV."""
    # y: JHU 5y=3, JHU 10y=2, Durham 5y=1, Durham 10y=0
    rows = [("JHU", "5 y", "auc5", 3), ("JHU", "10 y", "auc10", 2),
            ("Durham", "5 y", "auc5", 1), ("Durham", "10 y", "auc10", 0)]
    # faint vertical gridlines and the 0.50 reference
    for xt in (0.6, 0.7, 0.8, 0.9, 1.0):
        ax.axvline(xt, color="#F0F0F0", linewidth=0.6, zorder=0)
    ax.axvline(0.50, color="#888888", linewidth=0.9, linestyle="--", zorder=1)
    # dashed separator between the JHU and Durham blocks
    ax.axhline(1.5, color="#BBBBBB", linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    for coh, hz, key, y in rows:
        v, lo, hi = summ[coh][key], summ[coh][key + "lo"], summ[coh][key + "hi"]
        c = COHORT_COLORS[coh]
        ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", color=c, ecolor=c,
                    markersize=6.0, elinewidth=1.5, capsize=3.5, zorder=3)
        # numeric result just outside the right spine (axis-fraction x, data y)
        ax.text(1.03, y, f"{v:.2f} ({lo:.2f}–{hi:.2f})", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=8.5, color="#000000", clip_on=False)
    ax.set_yticks([r[3] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=9.0)
    # cohort name once, vertical, centred on its block, in cohort color
    ax.text(-0.085, 2.5, "JHU", transform=ax.get_yaxis_transform(), rotation=90,
            va="center", ha="center", fontsize=10.5, fontweight="bold",
            color=COHORT_COLORS["JHU"])
    ax.text(-0.085, 0.5, "Durham", transform=ax.get_yaxis_transform(), rotation=90,
            va="center", ha="center", fontsize=10.5, fontweight="bold",
            color=COHORT_COLORS["Durham"])
    ax.set_ylim(-0.5, 3.5)
    ax.set_xlim(0.50, 1.00)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Time-dependent AUC (95% CI)", fontsize=9.5)
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def build():
    data = get_cohort_data()
    summ = read_summary(data)
    curves = read_aj_curves()

    fig = plt.figure(figsize=(9.0, 8.4), dpi=600, facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.40,
                             left=0.12, right=0.985, top=0.955, bottom=0.44)
    bot_outer = fig.add_gridspec(1, 1, left=0.20, right=0.72, top=0.29, bottom=0.11)

    top = outer[0].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.18)
    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(top[0, 1])
    panel_a_distribution(ax_a, data)
    panel_b_roc(ax_b, data)
    ax_b.set_anchor("S")

    mid = outer[1].subgridspec(1, 2, wspace=0.46)
    ax_c = fig.add_subplot(mid[0, 0])
    ax_d = fig.add_subplot(mid[0, 1])
    panel_aj_cif(ax_c, curves["JHU"], summ["JHU"], "JHU Nat. History",
                 COHORT_COLORS["JHU"], ylabel=True, show_legend=True)
    panel_aj_cif(ax_d, curves["Durham"], summ["Durham"], "Durham VA",
                 COHORT_COLORS["Durham"], ylabel=False, show_legend=False)

    ax_e = fig.add_subplot(bot_outer[0, 0])
    panel_e_td_auc(ax_e, summ)

    # Panel letters anchored to each axes' bounding box (immediately above its
    # upper-left corner). Draw once first so the aspect-locked ROC panel reports
    # its final position.
    fig.canvas.draw()
    _lab = dict(fontsize=16, fontweight="bold", family="DejaVu Sans",
                va="bottom", ha="left", color="#000000")
    for _ax, _letter in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"), (ax_e, "e")):
        _bb = _ax.get_position()
        fig.text(_bb.x0 - 0.035, _bb.y1 + 0.010, _letter, **_lab)

    os.makedirs(OUT_BASE, exist_ok=True)
    out_pdf = os.path.join(OUT_BASE, "Figure2_unified.pdf")
    out_tif = os.path.join(OUT_BASE, "Figure2_unified.tiff")
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_tif, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Saved {out_pdf}")
    print(f"Saved {out_tif}")


if __name__ == "__main__":
    build()

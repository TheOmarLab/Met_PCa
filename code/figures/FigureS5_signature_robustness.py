"""
Figure S5 renderer: Met-Score signature robustness.

Calculation-free. Reads only the identifier-free aggregates written to
outs/FigureS5/ by code/signature_discovery/Met_Score_Gene_Consistency.R
(--figure-s5-only) and code/signature_discovery/Met_PCa_LASSO_vs_Ridge.R
(--figure-s5-only). It performs no fitting, AUC estimation, confidence-interval
calculation, resampling, or testing.

  a: 45-gene pooled Hedges g with 95% prediction intervals
  b: leave-one-cohort-out reselection membership
  c: nested internal-external discovery-cohort AUC
  d: 10-year external time-dependent AUC in JHU and Durham
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator


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
AGG = os.path.join(ROOT, "outs", "FigureS5")
OUT_BASE = os.path.join(ROOT, "figures")

INK = "#1a1a1a"
DIR_COLOR = {"POS": "#B0443E", "NEG": "#2E5A8C"}   # up in mets / down in mets
BAR_COLOR = "#41608C"
DOT_ON, DOT_OFF = "#1f3b63", "#C7CDD4"             # retained / dropped
# panel c models and panel d models share the four alternative colours; the
# M1 slot is navy (Full-panel ridge in c, Frozen Met-Score in d).
MODEL_COLOR = {
    "Full-panel ridge (IECV)": "#2E5A8C", "Frozen Met-Score": "#2E5A8C",
    "Top-10 ridge": "#B0443E", "Top-20 ridge": "#E08214",
    "LASSO (alpha=1)": "#4D9221", "Elastic net (alpha=0.5)": "#7B3294",
    "LOCO-reselected signature": "#111111"}
LOSO_LABEL = "LOCO-reselected signature"
C_ORDER = ["Full-panel ridge (IECV)", "Top-10 ridge", "Top-20 ridge", "LASSO (alpha=1)", "Elastic net (alpha=0.5)"]
C_ORDER_PLUS = C_ORDER + [LOSO_LABEL]
D_ORDER = ["Frozen Met-Score", "Top-10 ridge", "Top-20 ridge", "LASSO (alpha=1)", "Elastic net (alpha=0.5)"]
COH_SHORT = {"GSE116918": "116918", "GSE41408": "41408", "GSE46691": "46691",
             "GSE51066": "51066", "GSE55935": "55935", "GSE70769": "70769"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.9,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9, "xtick.major.size": 3.2, "ytick.major.size": 3.0,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.unicode_minus": False})


def panel_a(ax, a):
    a = a.sort_values("pooled_hedges_g", ascending=True).reset_index(drop=True)
    y = np.arange(len(a))
    ax.axvline(0, color="#BBBBBB", lw=0.8, ls="--", zorder=1)
    for i, r in a.iterrows():
        col = DIR_COLOR[r["direction"]]
        ax.plot([r["pi_lo"], r["pi_hi"]], [i, i], color=col, lw=1.0, alpha=0.85, zorder=2)
        ax.scatter([r["pooled_hedges_g"]], [i], s=22, marker="o",
                   facecolor=col, edgecolor=col, linewidth=1.0, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(a["gene"], fontsize=5.2)
    ax.set_ylim(-1, len(a))
    ax.set_xlabel("pooled Hedges g (95% prediction interval)", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=DIR_COLOR["POS"],
                  markeredgecolor=DIR_COLOR["POS"], markersize=6, label="up in mets (POS)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=DIR_COLOR["NEG"],
                  markeredgecolor=DIR_COLOR["NEG"], markersize=6, label="down in mets (NEG)")]
    ax.legend(handles=leg, loc="lower right", fontsize=7, frameon=False, handletextpad=0.3)


def panel_b(ax_bar, ax_mat, b):
    order = pd.read_csv(os.path.join(AGG, "FigureS5_panelA_heterogeneity.csv")).sort_values(
        "pooled_hedges_g", ascending=True)["gene"].tolist()
    b = b.set_index("gene").loc[order].reset_index()
    cols = [c for c in b.columns if c.startswith("omit_")]
    x = np.arange(len(cols))
    kept = np.array([int(b[c].sum()) for c in cols])
    # top: number of the 45 genes reselected when each cohort is omitted
    ax_bar.bar(x, kept, color=BAR_COLOR, width=0.68, zorder=3)
    for xi, k in zip(x, kept):
        ax_bar.text(xi, k + 1.2, str(k), ha="center", va="bottom", fontsize=7, color=INK)
    ax_bar.set_xlim(-0.6, len(cols) - 0.4); ax_bar.set_ylim(0, 45)
    ax_bar.set_xticks([]); ax_bar.set_yticks([0, 20, 40]); ax_bar.tick_params(labelsize=7)
    ax_bar.set_ylabel("genes kept (of 45)", fontsize=8)
    for s in ("top", "right"):
        ax_bar.spines[s].set_visible(False)
    # matrix: retained (dark) vs dropped (grey)
    y = np.arange(len(b))
    for j, c in enumerate(cols):
        on = b[c].values.astype(int) == 1
        ax_mat.scatter(np.full(on.sum(), j), y[on], s=15, marker="o", color=DOT_ON, edgecolor="none", zorder=3)
        ax_mat.scatter(np.full((~on).sum(), j), y[~on], s=11, marker="o", color=DOT_OFF, edgecolor="none", zorder=2)
    ax_mat.set_xticks(x)
    ax_mat.set_xticklabels([COH_SHORT[c.replace("omit_", "")] for c in cols], fontsize=7, rotation=45, ha="right")
    ax_mat.set_yticks(y); ax_mat.set_yticklabels(b["gene"], fontsize=5.2)
    ax_mat.set_ylim(-1, len(b)); ax_mat.set_xlim(-0.6, len(cols) - 0.4)
    ax_mat.set_xlabel("cohort omitted", fontsize=9)
    for s in ("top", "right"):
        ax_mat.spines[s].set_visible(False)
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=DOT_ON, markersize=6, label="reselected"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=DOT_OFF, markersize=6, label="dropped")]
    ax_bar.legend(handles=leg, loc="upper left", fontsize=6.5, frameon=False, handletextpad=0.3,
                  bbox_to_anchor=(0.0, 1.0), ncol=2, columnspacing=1.0)


def _forest(ax, df, groups, group_key, order, xlab, xlim, y_gap=1.0, markers=None):
    n = len(order)
    yc = []
    cur = (len(groups) * (n + y_gap))
    labels = []
    for g in groups:
        base = cur
        for k, m in enumerate(order):
            yy = base - k
            row = df[(df[group_key] == g) & (df["model"] == m)]
            if len(row):
                r = row.iloc[0]
                col = MODEL_COLOR[m]; mk = (markers or {}).get(m, "o")
                ax.plot([r["lo"], r["hi"]], [yy, yy], color=col, lw=1.1, zorder=2)
                ax.scatter([r["auc"]], [yy], s=24, marker=mk, color=col, edgecolor="white", linewidth=0.6, zorder=4)
            yc.append(yy)
        labels.append((base - (n - 1) / 2.0, g))
        cur = base - n - y_gap
    ax.axvline(0.5, color="#CCCCCC", lw=0.7, zorder=1)
    ax.set_yticks([p for p, _ in labels])
    ax.set_yticklabels([l for _, l in labels], fontsize=8)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlab, fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_c(ax, c, loso):
    folds = ["GSE116918", "GSE41408", "GSE46691", "GSE51066", "GSE55935", "GSE70769"]
    d = c.rename(columns={"fold_omitted": "grp", "auc_lo": "lo", "auc_hi": "hi"})[["grp", "model", "auc", "lo", "hi"]].copy()
    ls = loso.rename(columns={"cohort": "grp", "ci_low": "lo", "ci_high": "hi"}).copy()
    ls["model"] = LOSO_LABEL
    d = pd.concat([d, ls[["grp", "model", "auc", "lo", "hi"]]], ignore_index=True)
    markers = {LOSO_LABEL: "D"}
    # inclusive label: five nested model strategies + the LOCO-reselected
    # signature validation (which has no inner tuning loop)
    _forest(ax, d, folds, "grp", C_ORDER_PLUS, "held-out discovery-cohort AUC", (0.38, 1.0), markers=markers)
    ax.set_yticklabels([COH_SHORT[f] for f in folds], fontsize=8)
    # reserve a blank band above the top fold so the legend clears every point/CI
    lo_y, hi_y = ax.get_ylim(); ax.set_ylim(lo_y, hi_y + 0.26 * (hi_y - lo_y))
    leg = [Line2D([0], [0], marker=markers.get(m, "o"), color="w", markerfacecolor=MODEL_COLOR[m],
                  markersize=6, label=m) for m in C_ORDER_PLUS]
    ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=3,
              fontsize=5.8, frameon=False, handletextpad=0.2, labelspacing=0.2,
              columnspacing=0.8, borderpad=0.2)


def panel_d(ax, dpanel):
    d = dpanel.rename(columns={"ci_lo": "lo", "ci_hi": "hi"})
    cohorts = ["JHU", "Durham"]
    _forest(ax, d, cohorts, "cohort", D_ORDER, "10-year time-dependent AUC", (0.5, 0.95))
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=MODEL_COLOR[m],
                  markersize=6, label=m) for m in D_ORDER]
    ax.legend(handles=leg, loc="lower left", fontsize=6.3, frameon=False, handletextpad=0.2,
              labelspacing=0.25, borderpad=0.2)


def build():
    a = pd.read_csv(os.path.join(AGG, "FigureS5_panelA_heterogeneity.csv"))
    b = pd.read_csv(os.path.join(AGG, "FigureS5_panelB_membership.csv"))
    c = pd.read_csv(os.path.join(AGG, "FigureS5_panelC_outer_auc.csv"))
    loso = pd.read_csv(os.path.join(AGG, "FigureS5_panelC_LOSO_signature_auc.csv"))
    d = pd.read_csv(os.path.join(AGG, "FigureS5_panelD_timeAUC.csv"))

    fig = plt.figure(figsize=(8.9, 10.9), dpi=600, facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[2.7, 1.0], width_ratios=[1.0, 1.0],
                          hspace=0.13, wspace=0.20, top=0.965, bottom=0.055, left=0.09, right=0.985)
    ax_a = fig.add_subplot(gs[0, 0])
    gs_b = gs[0, 1].subgridspec(2, 1, height_ratios=[1.0, 8.0], hspace=0.04)
    ax_bbar = fig.add_subplot(gs_b[0, 0]); ax_bmat = fig.add_subplot(gs_b[1, 0])
    ax_c = fig.add_subplot(gs[1, 0]); ax_d = fig.add_subplot(gs[1, 1])
    panel_a(ax_a, a); panel_b(ax_bbar, ax_bmat, b); panel_c(ax_c, c, loso); panel_d(ax_d, d)

    lab = dict(fontsize=14, fontweight="semibold", family="DejaVu Sans", va="bottom", ha="left", color=INK)
    for ax, ch in ((ax_a, "a"), (ax_bbar, "b"), (ax_c, "c"), (ax_d, "d")):
        p = ax.get_position()
        fig.text(p.x0 - 0.028, p.y1 + 0.006, ch, **lab)

    os.makedirs(OUT_BASE, exist_ok=True)
    out_pdf = os.path.join(OUT_BASE, "FigureS5_signature_robustness.pdf")
    out_tif = os.path.join(OUT_BASE, "FigureS5_signature_robustness.tiff")
    out_png = os.path.join(OUT_BASE, "FigureS5_signature_robustness.png")
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_tif, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_png, format="png", dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")
    print(f"Saved {out_tif}")
    print(f"Saved {out_png}")


if __name__ == "__main__":
    build()

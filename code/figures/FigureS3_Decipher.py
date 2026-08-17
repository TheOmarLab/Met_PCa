"""
Figure S3 renderer: JHU / Durham Met-Score vs Decipher-marker surrogate.

Calculation-free. Reads only the identifier-free aggregates written by
code/ancillary/Met_PCa_Survival_DECIPHER.R to outs/FigureS3/. It performs no
survival fitting, AUC estimation, bootstrapping, or resampling.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
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
AGG_DIR = os.path.join(ROOT, "outs", "FigureS3")
OUT_BASE = os.path.join(ROOT, "figures")

MS_LABEL, DEC_LABEL = "Met-Score", "Decipher-marker surrogate"
SIG_COLORS = {MS_LABEL: "#2E5A8C", DEC_LABEL: "#B0443E"}
SIG_ORDER = [MS_LABEL, DEC_LABEL]
COHORT_ORDER = ["JHU", "Durham"]
COHORT_LABEL = {"JHU": "JHU Nat. History", "Durham": "Durham VA"}
SHAPE = "s"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0, "ytick.major.width": 0.0,
    "xtick.major.size": 3.4, "ytick.major.size": 0.0,
    "xtick.direction": "out",
    "axes.labelsize": 10.0, "xtick.labelsize": 9.0, "ytick.labelsize": 0.001,
    "axes.edgecolor": "#1a1a1a", "axes.labelcolor": "#1a1a1a",
    "xtick.color": "#1a1a1a", "text.color": "#1a1a1a",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "axes.unicode_minus": False,
})

# display rows top-to-bottom: All-patients block (JHU, Durham), then GS7 block
ROW_SPEC = [("All patients", "JHU"), ("All patients", "Durham"),
            ("GS7", "JHU"), ("GS7", "Durham")]


def _y_positions(gap=1.0):
    y, cur = [], (len(ROW_SPEC) + gap)
    for i in range(2):
        y.append(cur); cur -= 1
    cur -= gap
    for i in range(2):
        y.append(cur); cur -= 1
    return np.array(y, dtype=float)


def _lookup(df, subset, cohort, score, cols):
    r = df[(df["subset"] == subset) & (df["cohort"] == cohort) & (df["score"] == score)]
    if len(r) != 1:
        raise ValueError(f"expected one row for {subset}/{cohort}/{score}, got {len(r)}")
    return tuple(float(r.iloc[0][c]) for c in cols)


def _row_forest(ax, df, cols, y_pos, div_y):
    off = {MS_LABEL: +0.16, DEC_LABEL: -0.16}
    for (subset, cohort), y in zip(ROW_SPEC, y_pos):
        for sig in SIG_ORDER:
            v, lo, hi = _lookup(df, subset, cohort, sig, cols)
            yy = y + off[sig]
            ax.plot([lo, hi], [yy, yy], color=SIG_COLORS[sig], linewidth=1.4, zorder=2)
            for cap in (lo, hi):
                ax.plot([cap, cap], [yy - 0.07, yy + 0.07], color=SIG_COLORS[sig], linewidth=1.0)
            ax.scatter([v], [yy], marker=SHAPE, s=56, color=SIG_COLORS[sig],
                       edgecolor="white", linewidth=0.9, zorder=5)


def build():
    panelA = pd.read_csv(os.path.join(AGG_DIR, "FigureS3_panelA_timeAUC.csv"))
    panelB = pd.read_csv(os.path.join(AGG_DIR, "FigureS3_panelB_hr.csv"))
    y_pos = _y_positions()
    div_y = (y_pos[1] + y_pos[2]) / 2.0

    fig = plt.figure(figsize=(9.6, 4.9), dpi=600, facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.5,
                          top=0.90, bottom=0.20, left=0.20, right=0.88)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)

    ax_a.axhline(div_y, color="#BBBBBB", linewidth=0.6, linestyle="--", zorder=1)
    ax_b.axhline(div_y, color="#BBBBBB", linewidth=0.6, linestyle="--", zorder=1)
    ax_b.axvline(1.0, color="#BBBBBB", linewidth=1.0, zorder=1)

    _row_forest(ax_a, panelA, ("auc", "ci_lo", "ci_hi"), y_pos, div_y)
    _row_forest(ax_b, panelB, ("hr", "ci_lo_robust", "ci_hi_robust"), y_pos, div_y)

    # left-margin row identifiers (cohort names + block headers)
    ax_a.set_yticks([])
    for (subset, cohort), y in zip(ROW_SPEC, y_pos):
        ax_a.text(-0.30, y, COHORT_LABEL[cohort], transform=ax_a.get_yaxis_transform(),
                  ha="left", va="center", fontsize=9.0, color="#1a1a1a", clip_on=False)
    hdr = dict(transform=ax_a.get_yaxis_transform(), ha="left", fontsize=10.0,
               color="#1a1a1a", family="DejaVu Sans", fontweight="bold", clip_on=False)
    ax_a.text(-0.34, float(y_pos[0]) + 0.62, "All patients", va="bottom", **hdr)
    ax_a.text(-0.34, div_y, "Gleason 7 subset", va="center", **hdr)

    ax_a.set_xlim(0.55, 1.0)
    ax_a.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0])
    ax_a.set_xlabel("ten-year time-dependent AUC", fontsize=10.0)
    ax_b.set_xscale("log")
    ax_b.set_xlim(0.8, 8.5)
    ax_b.set_xticks([1, 2, 4, 8])
    ax_b.set_xticklabels(["1", "2", "4", "8"])
    ax_b.xaxis.set_minor_locator(NullLocator())
    ax_b.set_xlabel("cause-specific summary HR per cohort SD", fontsize=10.0)
    for ax in (ax_a, ax_b):
        ax.set_ylim(min(y_pos) - 0.7, max(y_pos) + 1.3)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color("#1a1a1a")
        ax.spines["bottom"].set_linewidth(1.0)
        ax.set_yticks([])

    # lowercase panel labels directly above each panel
    _lab = dict(fontsize=13, fontweight="semibold", family="DejaVu Sans",
                va="bottom", ha="left", color="#1a1a1a")
    fig.text(ax_a.get_position().x0, 0.925, "a", **_lab)
    fig.text(ax_b.get_position().x0, 0.925, "b", **_lab)

    # compact legend (marker identity only)
    ax_leg = fig.add_axes([0.0, 0.02, 1.0, 0.06]); ax_leg.set_axis_off()
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
    for label, x in ((MS_LABEL, 0.30), (DEC_LABEL, 0.50)):
        ax_leg.scatter([x], [0.5], marker=SHAPE, s=55, color=SIG_COLORS[label],
                       edgecolor="white", linewidth=0.7, transform=ax_leg.transAxes)
        ax_leg.text(x + 0.013, 0.5, label, ha="left", va="center", fontsize=9.0,
                    color="#1a1a1a", transform=ax_leg.transAxes)

    os.makedirs(OUT_BASE, exist_ok=True)
    out_pdf = os.path.join(OUT_BASE, "FigureS3_Decipher.pdf")
    out_tif = os.path.join(OUT_BASE, "FigureS3_Decipher.tiff")
    out_png = os.path.join(OUT_BASE, "FigureS3_Decipher.png")
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_tif, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_png, format="png", dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")
    print(f"Saved {out_tif}")
    print(f"Saved {out_png}")


if __name__ == "__main__":
    build()

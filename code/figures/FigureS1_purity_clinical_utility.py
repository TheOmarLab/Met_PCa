#!/usr/bin/env python3
"""Figure S1: tumour-purity robustness and clinical utility. Calculation-free.

Reads pre-computed aggregates and draws them; it computes nothing itself.

  a  before/after ESTIMATE-purity Met-Score HR (per cohort SD) forest, JHU + Durham
  b  10-year two-arm decision-curve analysis, JHU and Durham as internal facets
  c  JHU 10-year risk-band reclassification flow (clinical vs + Met-Score)
  d  Durham 10-year risk-band reclassification flow

Inputs:
  outs/purity/FigureS1a_purity_robustness.csv                          (panel a)
  outs/DCA/dca_curves.csv, outs/DCA/delta_net_benefit.csv              (panel b)
  outs/DCA/risk_band_transitions_10y.csv                               (panels c/d)
Outputs (figures/): FigureS1_purity_clinical_utility.{pdf,tiff,png}.
The CAPRA-S selective-reflex curve is a separate supporting candidate, not in Figure S1.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Patch


def _find_root():
    env = os.environ.get("MET_PCA_ROOT")
    if env and os.path.isdir(env):
        return env
    d = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "code")) and os.path.isdir(os.path.join(d, "outs")):
            return d
        d = os.path.dirname(d)
    raise FileNotFoundError("Could not locate project root; set MET_PCA_ROOT.")


ROOT = _find_root()
OUTD = os.path.join(ROOT, "outs", "DCA")
PURD = os.path.join(ROOT, "outs", "purity")
FIGD = os.path.join(ROOT, "figures")
FIGD_DCA = os.path.join(ROOT, "figures", "DCA")
os.makedirs(FIGD, exist_ok=True)
os.makedirs(FIGD_DCA, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.9,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "axes.unicode_minus": False,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "axes.labelsize": 9.5,
})

COHORT_COLORS = {"JHU": "#0072B2", "Durham": "#009E73"}
COHORT_LABELS = {"JHU": "JHU Nat. History", "Durham": "Durham VA"}
PANEL_KW = dict(fontsize=16, fontweight="bold", family="DejaVu Sans",
                va="bottom", ha="left", color="#000000")


def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, **PANEL_KW)


# ---------------------------------------------------------------------------
# Data.
# ---------------------------------------------------------------------------
purity = pd.read_csv(os.path.join(PURD, "FigureS1a_purity_robustness.csv"))
mecorr = pd.read_csv(os.path.join(PURD, "FigureS1a_microenvironment_correlations.csv"))
dca = pd.read_csv(os.path.join(OUTD, "dca_curves.csv"))
delta = pd.read_csv(os.path.join(OUTD, "delta_net_benefit.csv"))
trans = pd.read_csv(os.path.join(OUTD, "risk_band_transitions_10y.csv"))
capras = pd.read_csv(os.path.join(OUTD, "capras_rule_utility.csv"))

LABS = ["<10%", "10-20%", ">20%"]
IDX = {b: i for i, b in enumerate(LABS)}
EXPECT_N = {"JHU": 235, "Durham": 555}


def _matrix(coh):
    d = trans[trans.cohort == coh]
    M = np.zeros((3, 3), dtype=int)
    for _, r in d.iterrows():
        M[IDX[r.clinical_band], IDX[r.full_band]] = int(r.n)
    return M


for coh in ("JHU", "Durham"):
    if _matrix(coh).sum() != EXPECT_N[coh]:
        raise ValueError(f"{coh} transition total {_matrix(coh).sum()} != {EXPECT_N[coh]}")


def build_streams(coh):
    d = trans[trans.cohort == coh]
    rows = []
    for _, r in d.iterrows():
        n = int(r.n)
        if n == 0:
            continue
        di = ("stable" if IDX[r.clinical_band] == IDX[r.full_band]
              else ("down" if IDX[r.full_band] < IDX[r.clinical_band] else "up"))
        rows.append({"cohort": coh, "from": r.clinical_band, "to": r.full_band, "n": n, "direction": di})
    return pd.DataFrame(rows)


def build_totals(coh):
    d = trans[trans.cohort == coh]
    return pd.DataFrame([{"cohort": coh, "band": b,
                          "clin_n": int(d[d.clinical_band == b].n.sum()),
                          "comb_n": int(d[d.full_band == b].n.sum())} for b in LABS])


streams = pd.concat([build_streams("JHU"), build_streams("Durham")], ignore_index=True)
bt = pd.concat([build_totals("JHU"), build_totals("Durham")], ignore_index=True)

C_METS = "#0072B2"   # Clinical + Met-Score
C_CLIN = "#E69F00"   # Clinical
C_REF = "#9AA0A6"    # Treat-all / Treat-none references
C_REFLEX = "#6A3D9A"  # CAPRA-S selective reflex


# ---------------------------------------------------------------------------
# Panel a (composite): left = paired before/after-purity HR dumbbell (one row per
# cohort); right = Met-Score microenvironment correlation dot-whisker.
# ---------------------------------------------------------------------------
def dumbbell_panel(ax):
    yc = {"JHU": 1.0, "Durham": 0.0}; dy = 0.16
    ax.axvline(1.0, color="#888888", ls=(0, (4, 3)), lw=1.1, zorder=1)
    for coh in ("JHU", "Durham"):
        col = COHORT_COLORS[coh]; y = yc[coh]
        b = purity[(purity.cohort == coh) & (purity.model == "before")].iloc[0]
        a = purity[(purity.cohort == coh) & (purity.model == "after")].iloc[0]
        ax.plot([b.HR_perSD, a.HR_perSD], [y + dy, y - dy], color=col, lw=1.2, alpha=0.55, zorder=2)
        for r, yy, fill in ((b, y + dy, "white"), (a, y - dy, col)):
            ax.plot([r.CI_lo, r.CI_hi], [yy, yy], color=col, lw=1.6, zorder=3, solid_capstyle="round", alpha=0.9)
            ax.scatter([r.HR_perSD], [yy], s=58, facecolor=fill, edgecolor=col, linewidth=1.4, zorder=4)
        # HR text sits just above the right-most edge (max upper CI) of the two dumbbell bars
        xr = max(float(b.CI_hi), float(a.CI_hi))
        ax.text(xr, y + dy + 0.12, f"{b.HR_perSD:.2f} $\\rightarrow$ {a.HR_perSD:.2f}",
                ha="right", va="bottom", fontsize=9.2, color="#222222")
    ax.set_yticks([yc["JHU"], yc["Durham"]]); ax.set_yticklabels(["JHU", "Durham"], fontsize=10.5)
    for lab, coh in (("JHU", "JHU"), ("Durham", "Durham")):
        ax.get_yticklabels()[0 if coh == "JHU" else 1].set_color(COHORT_COLORS[coh])
    ax.set_xscale("log"); ax.set_xlim(0.9, 2.75)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5]); ax.set_xticklabels(["1.0", "1.5", "2.0", "2.5"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Met-Score HR per 1 SD (95% CI)", fontsize=10.5)
    ax.tick_params(axis="x", labelsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(1.2); ax.tick_params(axis="y", length=0)


def mecorr_panel(ax):
    rows = ["Tumour purity", "Stromal score", "Immune score"]; ypos = {s: (len(rows) - 1 - i) for i, s in enumerate(rows)}
    off = {"JHU": 0.16, "Durham": -0.16}
    ax.axvline(0.0, color="#888888", ls=(0, (4, 3)), lw=1.1, zorder=1)
    for coh in ("JHU", "Durham"):
        col = COHORT_COLORS[coh]
        for sc in rows:
            r = mecorr[(mecorr.cohort == coh) & (mecorr.score == sc)].iloc[0]
            y = ypos[sc] + off[coh]
            ax.plot([r.ci_lo, r.ci_hi], [y, y], color=col, lw=1.6, zorder=3, solid_capstyle="round")
            ax.scatter([r.rho], [y], s=48, color=col, edgecolor="white", linewidth=0.8, zorder=4)
    ax.set_yticks([ypos[s] for s in rows]); ax.set_yticklabels(rows, fontsize=10.0)
    ax.set_xlim(-0.34, 0.42); ax.set_xticks([-0.2, 0.0, 0.2, 0.4])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("Spearman ρ with Met-Score", fontsize=10.5)
    ax.tick_params(axis="x", labelsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(1.2); ax.tick_params(axis="y", length=0)


# ---------------------------------------------------------------------------
# Panel b: two-arm DCA facets (no Durham ribbon).
# ---------------------------------------------------------------------------
def dca_panel(ax, cohort, basis, title, ylim, yticks):
    d = dca[(dca.cohort == cohort) & (dca.horizon == 120)]
    clin = d[(d.basis == basis) & (d.model == "clin")].sort_values("threshold")
    full = d[(d.basis == basis) & (d.model == "full")].sort_values("threshold")
    ta = d[d.model == "treat_all"].sort_values("threshold")
    ax.axhline(0, color=C_REF, ls=(0, (1, 1.8)), lw=1.6, zorder=1, label="Treat none")
    ax.plot(ta.threshold, ta.nb, color=C_REF, ls=(0, (5, 2.5)), lw=1.8, zorder=2, label="Treat all")
    ax.plot(clin.threshold, clin.nb, color=C_CLIN, lw=3.2, zorder=3, solid_capstyle="round", label="Clinical model")
    ax.plot(full.threshold, full.nb, color=C_METS, lw=3.0, zorder=4, solid_capstyle="round", label="Clinical model + Met-Score")
    ax.set_xlim(0.01, 0.20); ax.set_ylim(*ylim)
    ax.set_xticks([0.01, 0.05, 0.10, 0.15, 0.20]); ax.set_xticklabels(["1%", "5%", "10%", "15%", "20%"])
    ax.set_yticks(yticks)
    ax.set_xlabel("Threshold probability", fontsize=12); ax.set_ylabel("Net benefit per 100", fontsize=12)
    ax.tick_params(axis="both", labelsize=10.5)
    ax.set_title(title, fontsize=13, color=COHORT_COLORS[cohort], fontweight="bold", pad=6)
    ax.spines["left"].set_linewidth(1.2); ax.spines["bottom"].set_linewidth(1.2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------------------
# Panels c/d: reclassification Sankey.
# ---------------------------------------------------------------------------
C_DOWN, C_UP, C_STABLE = "#3D8DBF", "#D55E00", "#BFC4CC"
NODE_FILL, NODE_EDGE, TXT = "#EAECEF", "#4A4F57", "#22262B"
X_L, X_R = 0.17, 0.83
NODE_W = 0.045
GAP = 0.055
TOP, BOT = 0.83, 0.13


def band_layout(counts):
    total = counts.sum()
    span = (TOP - BOT) - GAP * (len(LABS) - 1)
    y = TOP
    out = {}
    for b in LABS[::-1]:
        h = span * (counts[b] / total) if total > 0 else 0.0
        out[b] = (y - h, y)
        y = y - h - GAP
    return out


def ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha):
    cx = (x0 + x1) / 2
    verts = [(x0, y0a), (cx, y0a), (cx, y1a), (x1, y1a), (x1, y1b),
             (cx, y1b), (cx, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def sankey_panel(ax, cohort):
    s = streams[streams.cohort == cohort].copy()
    b = bt[bt.cohort == cohort].set_index("band")
    clin_counts, comb_counts = b["clin_n"], b["comb_n"]
    L = band_layout(clin_counts)
    R = band_layout(comb_counts)

    def unit(counts, layout, band):
        y0, y1 = layout[band]
        c = counts[band]
        return (y1 - y0) / c if c > 0 else 0.0

    s["order"] = s.direction.map({"stable": 0, "down": 1, "up": 1})
    s = s.sort_values(["order", "n"], ascending=[True, False]).reset_index(drop=True)
    left_run = {k: v[1] for k, v in L.items()}
    right_run = {k: v[1] for k, v in R.items()}
    mids = []
    for _, row in s.iterrows():
        fb, tb = row["from"], row["to"]
        hL = unit(clin_counts, L, fb) * row.n
        hR = unit(comb_counts, R, tb) * row.n
        y0a = left_run[fb]; y0b = y0a - hL
        y1a = right_run[tb]; y1b = y1a - hR
        left_run[fb] = y0b; right_run[tb] = y1b
        col = {"down": C_DOWN, "up": C_UP, "stable": C_STABLE}[row.direction]
        al = 0.30 if row.direction == "stable" else 0.72
        ribbon(ax, X_L + NODE_W / 2, X_R - NODE_W / 2, y0a, y0b, y1a, y1b, col, al)
        if row.direction != "stable":
            mids.append((((y0a + y0b) / 2 + (y1a + y1b) / 2) / 2, row))

    for col_x, layout, counts, side in [(X_L, L, clin_counts, "left"), (X_R, R, comb_counts, "right")]:
        for band in LABS:
            y0, y1 = layout[band]
            ax.add_patch(patches.Rectangle((col_x - NODE_W / 2, y0), NODE_W, y1 - y0,
                         facecolor=NODE_FILL, edgecolor=NODE_EDGE, lw=1.3, zorder=3))
            lab = f"{band}\n(n={int(counts[band])})"
            dx = -(NODE_W / 2 + 0.016) if side == "left" else (NODE_W / 2 + 0.016)
            ax.text(col_x + dx, (y0 + y1) / 2, lab, ha=("right" if side == "left" else "left"),
                    va="center", fontsize=9.6, color=TXT, zorder=4)

    mids.sort(key=lambda t: t[0])
    min_gap = 0.098
    ys = [m[0] for m in mids]
    for i in range(1, len(ys)):
        if ys[i] - ys[i - 1] < min_gap:
            ys[i] = ys[i - 1] + min_gap
    if ys and ys[-1] > TOP:
        ys = [y - (ys[-1] - TOP) for y in ys]
    for (_, row), my in zip(mids, ys):
        col = C_DOWN if row.direction == "down" else C_UP
        txt = f"{row['from']} to {row['to']}   n={int(row.n)}"
        ax.text(0.5, my, txt, ha="center", va="center", fontsize=8.9, color="white",
                zorder=6, bbox=dict(boxstyle="round,pad=0.30", fc=col, ec="none", alpha=0.96))

    ax.text(X_L, TOP + 0.055, "Clinical\nmodel", ha="center", va="bottom", fontsize=10.5, color=TXT, fontweight="bold")
    ax.text(X_R, TOP + 0.055, "+ Met-Score", ha="center", va="bottom", fontsize=10.5, color=TXT, fontweight="bold")
    ax.set_title(f"{COHORT_LABELS[cohort]}  (n={int(clin_counts.sum())})",
                 fontsize=13, color=COHORT_COLORS[cohort], fontweight="bold", pad=6)
    ax.set_xlim(0, 1); ax.set_ylim(0.09, 0.97); ax.axis("off")


# ---------------------------------------------------------------------------
# Figure S1 assembly: a forest (full width) / b DCA facets / c-d flows.
# ---------------------------------------------------------------------------
def render_figS1():
    # Explicit placement so the a->b and b->(c,d) vertical gaps are exactly equal.
    fig = plt.figure(figsize=(10.2, 11.4), dpi=600)

    LEFT, RIGHT = 0.115, 0.972
    WIDTH = RIGHT - LEFT
    G_H = 0.075                       # horizontal gutter between the two facets
    COLW = (WIDTH - G_H) / 2.0
    X2 = LEFT + COLW + G_H

    GAP = 0.050                       # equal vertical gap: a->b and b-block->(c,d)
    A_XLAB, B_XLAB = 0.036, 0.040     # room below a / b axes for x-tick + axis title
    LEG_H = 0.024                     # DCA legend band height

    ya1 = 0.955
    Ha = 0.140
    ya0 = ya1 - Ha

    yb1 = ya0 - A_XLAB - GAP
    Hb = 0.185
    yb0 = yb1 - Hb

    leg_top = yb0 - B_XLAB
    leg_y = leg_top - LEG_H / 2.0     # DCA legend centre, just below b's x labels

    ycd1 = (leg_top - LEG_H) - GAP    # equal gap below the full b block (incl. legend)
    Hcd = 0.315                       # taller reclassification facets
    ycd0 = ycd1 - Hcd

    # panel a = composite (dumbbell + microenvironment dot-whisker), one letter.
    # Row a uses a slightly wider internal gap than b / c-d while keeping the outer
    # edges aligned at LEFT and RIGHT (the two are different plot types, not facets).
    G_Ha = G_H + 0.032
    COLWa = (WIDTH - G_Ha) / 2.0
    X2a = LEFT + COLWa + G_Ha
    ax_al = fig.add_axes([LEFT, ya0, COLWa, Ha])
    ax_ar = fig.add_axes([X2a, ya0, COLWa, Ha])
    dumbbell_panel(ax_al)
    mecorr_panel(ax_ar)
    from matplotlib.lines import Line2D
    lh = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLORS["JHU"], markeredgecolor=COHORT_COLORS["JHU"], markersize=9, label="JHU"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLORS["Durham"], markeredgecolor=COHORT_COLORS["Durham"], markersize=9, label="Durham"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#555555", markersize=9, label="before purity"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", markeredgecolor="#555555", markersize=9, label="after purity")]
    fig.legend(handles=lh, loc="center", ncol=4, frameon=False, fontsize=9.3,
               bbox_to_anchor=(0.5, 0.980), columnspacing=1.6, handletextpad=0.3)

    ax_bj = fig.add_axes([LEFT, yb0, COLW, Hb])
    ax_bd = fig.add_axes([X2, yb0, COLW, Hb])
    dca_panel(ax_bj, "JHU", "optimism_corrected", COHORT_LABELS["JHU"] + " · 10-year", (0, 20.5), [0, 5, 10, 15, 20])
    dca_panel(ax_bd, "Durham", "frozen", COHORT_LABELS["Durham"] + " · 10-year", (0, 6.5), [0, 2, 4, 6])

    ax_c = fig.add_axes([LEFT, ycd0, COLW, Hcd])
    ax_d = fig.add_axes([X2, ycd0, COLW, Hcd])
    sankey_panel(ax_c, "JHU")
    sankey_panel(ax_d, "Durham")

    # DCA legend directly below the DCA facets; flow legend at the very bottom
    h, l = ax_bj.get_legend_handles_labels()
    order = ["Clinical model + Met-Score", "Clinical model", "Treat all", "Treat none"]
    hl = {lab: hh for hh, lab in zip(h, l)}
    fig.legend([hl[o] for o in order], order, loc="center", ncol=4, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, leg_y), handlelength=2.2, columnspacing=1.9)
    leg = [Patch(fc=C_DOWN, alpha=0.72, label="Down-classified to lower risk (fewer treated)"),
           Patch(fc=C_UP, alpha=0.72, label="Up-classified to higher risk (more treated)"),
           Patch(fc=C_STABLE, alpha=0.5, label="Unchanged risk band")]
    fig.legend(handles=leg, loc="center", ncol=3, frameon=False, fontsize=10.5,
               bbox_to_anchor=(0.5, ycd0 - 0.045))

    for ax, letter in [(ax_al, "a"), (ax_bj, "b"), (ax_c, "c"), (ax_d, "d")]:
        p = ax.get_position()
        panel_letter(fig, p.x0 - 0.058, p.y1 + 0.008, letter)

    base = os.path.join(FIGD, "FigureS1_purity_clinical_utility")
    fig.savefig(base + ".pdf", format="pdf", metadata={"CreationDate": None})
    fig.savefig(base + ".tiff", format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base + ".png", format="png", dpi=300)
    plt.close(fig)
    print("wrote", base + ".{pdf,tiff,png}")


# ---------------------------------------------------------------------------
# CAPRA-S supporting figure, not in Figure S1.
# ---------------------------------------------------------------------------
def render_capras():
    cr = capras[(capras.type == "rule_dca") & (capras.horizon == 120)].sort_values("threshold")
    fig, ax = plt.subplots(figsize=(5.0, 4.4), dpi=600)
    ax.axhline(0, color=C_REF, ls=(0, (1, 1.8)), lw=1.2, zorder=1, label="Treat none")
    ax.plot(cr.threshold, cr.nb_treat_all, color=C_REF, ls=(0, (5, 2.5)), lw=1.4, zorder=2, label="Treat all")
    ax.plot(cr.threshold, cr.nb_capras_high, color=C_CLIN, lw=2.4, zorder=3, solid_capstyle="round", label="CAPRA-S high")
    ax.plot(cr.threshold, cr.nb_reflex, color=C_REFLEX, lw=2.6, zorder=4, solid_capstyle="round", label="CAPRA-S selective reflex")
    ax.set_xlim(0.01, 0.20); ax.set_ylim(0, 20)
    ax.set_xticks([0.01, 0.05, 0.10, 0.15, 0.20]); ax.set_xticklabels(["1%", "5%", "10%", "15%", "20%"])
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit per 100")
    ax.set_title(f"{COHORT_LABELS['JHU']} CAPRA-S · 10-year", fontsize=11, fontweight="bold", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    order = ["CAPRA-S selective reflex", "CAPRA-S high", "Treat all", "Treat none"]
    hl = {lab: hh for hh, lab in zip(handles, labels)}
    ax.legend([hl[o] for o in order], order, loc="lower left", frameon=False, fontsize=8.5)
    fig.tight_layout()
    base = os.path.join(FIGD_DCA, "CAPRAS_candidate")
    fig.savefig(base + ".pdf", format="pdf", metadata={"CreationDate": None})
    fig.savefig(base + ".tiff", format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base + ".png", format="png", dpi=150)
    plt.close(fig)
    print("wrote", base + ".{pdf,tiff,png}")


if __name__ == "__main__":
    render_figS1()
    render_capras()

#!/usr/bin/env python3
"""Supplementary Figure S6: Met-Score robustness montage (six panels, 3x2).

Calculation-free renderer. Reads only the aggregate CSVs written by
code/ancillary/MetScore_Sensitivity.R (outs/FigureS6/); it fits, resamples, and
estimates nothing.

  a  follow-up robustness: high-vs-low HR, full / 5y / 10y admin truncation
  b  12-month salvage landmark: high-vs-low HR, landmark and + salvage-by-12mo
  c  Met-Score x race interaction (Durham): adjusted Met-Score HR per SD, White and Black
  d  multivariable covariate forest: GG, log2(PSA+1), pT, Met-Score per SD
  e  Met-Score HR under common-clinical vs + CCP gene-set adjustment
  f  alternative-cutoff sensitivity: locked Youden vs median / Sens90 / Spec90
     development-derived cutoffs; upper = cutoffs on the probability scale,
     lower = JHU/Durham 10-year sensitivity-specificity tradeoff

Outputs (figures/): FigureS6_sensitivity.{pdf,tiff,png}.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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
AGG = os.path.join(ROOT, "outs", "FigureS6")
FIGD = os.path.join(ROOT, "figures")
os.makedirs(FIGD, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.9,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "axes.unicode_minus": False,
})
COHORT_COLORS = {"JHU": "#0072B2", "Durham": "#009E73"}
COHORT_LABEL = {"JHU": "JHU", "Durham": "Durham"}
PANEL_KW = dict(fontsize=15, fontweight="bold", family="DejaVu Sans", va="bottom", ha="left", color="#000000")


def draw_forest(ax, rows, xlabel, xlim, xticks, right_x=0.82, header_gap=0.55, sep=False, text_beside=False):
    # rows: dicts with label, hr, lo, hi, color, header(optional), sep_before(optional)
    ypos = []; y = 0.0; sep_yy = None
    for i, r in enumerate(rows):
        if r.get("header") and i != 0:
            sep_yy = y + header_gap / 2.0; y += header_gap
        ypos.append(y); y += 1.0
    ymax = y
    ax.set_xscale("log")
    ax.axvline(1.0, color="#999999", lw=0.9, ls=(0, (4, 3)), zorder=1)
    if sep and sep_yy is not None:
        ax.axhline(ymax - sep_yy - 1.0, color="#CBCBCB", lw=0.8, ls=(0, (3, 3)), zorder=1, xmin=0.02, xmax=0.98)
    for r, yy in zip(rows, ypos):
        yp = ymax - yy - 1.0; col = r["color"]
        ax.plot([r["lo"], r["hi"]], [yp, yp], color=col, lw=2.1, zorder=3, solid_capstyle="round")
        ax.scatter([r["hr"]], [yp], s=44, color=col, zorder=4, edgecolor="white", linewidth=0.8)
        ax.text(-0.02, yp, r["label"], transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=9.6, color="#222222")
        if r.get("header"):
            ax.text(-0.02, yp + 0.42, r["header"], transform=ax.get_yaxis_transform(), ha="right",
                    va="bottom", fontsize=11.0, color=col, fontweight="bold")
        htxt = f"{r['hr']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})"
        if text_beside:
            # HR text immediately beside its bar (right of the upper CI, centred on the bar)
            ax.text(r["hi"] * 1.04, yp, htxt, ha="left", va="center", fontsize=8.2, color="#333333", clip_on=False)
        else:
            # HR text just above the right-most edge of its own bar, clamped to the axis
            ax.text(min(r["hi"], xlim[1]), yp + 0.26, htxt, ha="right", va="bottom",
                    fontsize=8.6, color="#333333", clip_on=False)
    ax.set_ylim(-0.6, ymax - 0.4); ax.set_yticks([])
    ax.set_xlim(*xlim); ax.set_xticks(xticks); ax.set_xticklabels([("%g" % t) for t in xticks])
    ax.set_xlabel(xlabel, fontsize=10.0); ax.tick_params(axis="x", labelsize=9.5)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)


# ---- load aggregates ------------------------------------------------------
pa = pd.read_csv(os.path.join(AGG, "panelA_followup_robustness.csv"))
pb = pd.read_csv(os.path.join(AGG, "panelB_salvage_landmark.csv"))
pc = pd.read_csv(os.path.join(AGG, "panelC_race_interaction.csv"))
pd_ = pd.read_csv(os.path.join(AGG, "panelD_multivariable.csv"))
pe = pd.read_csv(os.path.join(AGG, "panelE_ccp_adjustment.csv"))
pf_cuts = pd.read_csv(os.path.join(AGG, "panelF_cutoff_definitions.csv"))
pf_val = pd.read_csv(os.path.join(AGG, "panelF_threshold_strategies.csv"))

SCEN = ["Full follow-up", "Administrative truncation at 5 years", "Administrative truncation at 10 years"]
SCEN_LAB = {"Full follow-up": "Full follow-up", "Administrative truncation at 5 years": "Truncate 5 y",
            "Administrative truncation at 10 years": "Truncate 10 y"}


def rows_a():
    out = []
    for coh in ("JHU", "Durham"):
        first = True
        for sc in SCEN:
            r = pa[(pa.cohort == coh) & (pa.scenario == sc)].iloc[0]
            out.append(dict(label=SCEN_LAB[sc], hr=float(r.hr), lo=float(r.ci_lo), hi=float(r.ci_hi),
                            color=COHORT_COLORS[coh], header=COHORT_LABEL[coh] if first else None))
            first = False
    return out


def rows_b():
    out = []
    for coh in ("JHU", "Durham"):
        sub = pb[(pb.cohort == coh) & (pb.available == True)]
        first = True
        for _, r in sub.iterrows():
            lab = "Landmark 12 mo" if r.model == "Landmark 12 months" else "+ salvage by 12 mo"
            out.append(dict(label=lab, hr=float(r.hr), lo=float(r.ci_lo), hi=float(r.ci_hi),
                            color=COHORT_COLORS[coh], header=COHORT_LABEL[coh] if first else None))
            first = False
    return out


def rows_c():
    out = []
    for i, (_, r) in enumerate(pc.iterrows()):
        out.append(dict(label=str(r.label), hr=float(r.hr), lo=float(r.ci_lo), hi=float(r.ci_hi),
                        color=COHORT_COLORS["Durham"], header="Durham" if i == 0 else None))
    return out


def rows_d():
    out = []
    for coh in ("JHU", "Durham"):
        sub = pd_[(pd_.cohort == coh) & (pd_.estimable == True) & (pd_.display == True)]
        first = True
        for _, r in sub.iterrows():
            out.append(dict(label=str(r.label), hr=float(r.hr), lo=float(r.ci_lo), hi=float(r.ci_hi),
                            color=COHORT_COLORS[coh], header=COHORT_LABEL[coh] if first else None))
            first = False
    return out


def rows_e():
    out = []
    for coh in ("JHU", "Durham"):
        first = True
        for adj in ("Common clinical", "+ CCP gene-set score"):
            r = pe[(pe.cohort == coh) & (pe.adjustment == adj)].iloc[0]
            out.append(dict(label=adj, hr=float(r.hr), lo=float(r.ci_lo), hi=float(r.ci_hi),
                            color=COHORT_COLORS[coh], header=COHORT_LABEL[coh] if first else None))
            first = False
    return out


fig = plt.figure(figsize=(10.6, 12.0), dpi=600)
gs = fig.add_gridspec(3, 2, wspace=0.92, hspace=0.34, left=0.145, right=0.86, top=0.955, bottom=0.065)
ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0]); ax_d = fig.add_subplot(gs[1, 1])
ax_e = fig.add_subplot(gs[2, 0]); ax_f = fig.add_subplot(gs[2, 1])

draw_forest(ax_a, rows_a(), "High-vs-low metastasis HR (95% CI)", (0.7, 60), [1, 2, 5, 10, 20], right_x=0.80)
draw_forest(ax_b, rows_b(), "High-vs-low HR (95% CI), 12-mo landmark", (0.8, 90), [1, 2, 5, 10, 20, 40], right_x=0.80)
draw_forest(ax_c, rows_c(), "Adjusted Met-Score HR per SD (95% CI)", (0.8, 3.2), [1, 1.5, 2, 2.5, 3], right_x=0.82)
# formal Met-Score x race interaction (estimate/interval, from the aggregate)
_ic = pc.iloc[0]
ax_c.text(0.98, 0.05, f"Met-Score × race interaction  HR {float(_ic.interaction_ratio):.2f} "
          f"({float(_ic.interaction_lo):.2f}–{float(_ic.interaction_hi):.2f}), p = {float(_ic.interaction_p):.2f}",
          transform=ax_c.transAxes, ha="right", va="bottom", fontsize=8.2, color="#333333")
draw_forest(ax_d, rows_d(), "Adjusted HR (95% CI), multivariable model", (0.3, 44), [0.5, 1, 5, 20], right_x=0.88, sep=True, text_beside=True)
draw_forest(ax_e, rows_e(), "Met-Score HR per SD (95% CI)", (0.9, 2.6), [1.0, 1.5, 2.0, 2.5], right_x=0.86, sep=True)

# panel f — alternative-cutoff sensitivity: cutoff strip (upper) + 10y Se/Sp tradeoff (lower)
_pf = ax_f.get_position(); ax_f.remove()
FU_H = 0.052; LGAP = 0.052
ax_fu = fig.add_axes([_pf.x0, _pf.y1 - FU_H, _pf.width, FU_H])            # upper cutoff strip
ax_fl = fig.add_axes([_pf.x0, _pf.y0, _pf.width, (_pf.y1 - FU_H - LGAP) - _pf.y0])  # lower Se/Sp

# upper: four development-derived cutoffs on the probability scale, Youden distinguished
XPROB = (0.20, 0.42)
ax_fu.axhline(0.0, color="#666666", lw=1.0, zorder=1)
# (label_x, level) spread horizontally so labels clear in the 0.24-0.28 cluster
_LAB = {"Sens90": (0.226, 0), "Youden (locked)": (0.262, 1), "Median": (0.302, 0), "Spec90": (0.400, 0)}
for _, r in pf_cuts.iterrows():
    x = float(r.cutoff_value_out); prim = bool(r.is_primary); nm = str(r.cutoff_name)
    mk = "D" if prim else "o"; sz = 70 if prim else 40
    face = "#000000" if prim else "white"
    ax_fu.scatter([x], [0], marker=mk, s=sz, facecolor=face, edgecolor="#000000",
                  linewidth=1.1, zorder=4, clip_on=False)
    lx, lv = _LAB.get(nm, (x, 0)); ytxt = 1.0 + 1.15 * lv
    ax_fu.plot([x, lx], [0.07, ytxt - 0.08], color="#9A9A9A", lw=0.6, zorder=2)  # angled leader
    ax_fu.text(lx, ytxt, f"{nm}\n{x:.3f}", ha="center", va="bottom", fontsize=7.6,
               color=("#000000" if prim else "#333333"), fontweight=("bold" if prim else "normal"),
               linespacing=0.95)
ax_fu.set_xlim(*XPROB); ax_fu.set_ylim(-0.5, 3.6)
ax_fu.set_yticks([]); ax_fu.set_xticks([0.20, 0.25, 0.30, 0.35, 0.40])
ax_fu.tick_params(axis="x", labelsize=8.4, length=3)
ax_fu.set_xlabel("Development cutoff (predicted probability)", fontsize=8.8, labelpad=1.5)
for s in ("top", "right", "left"):
    ax_fu.spines[s].set_visible(False)
ax_fu.tick_params(axis="y", length=0)

# lower: 10-year (120-mo) sensitivity/specificity for the four strategies in JHU and Durham
val10 = pf_val[pf_val.horizon_months == 120]
STRAT_ORDER = ["Sens90", "Youden (locked)", "Median", "Spec90"]
band = {s: (len(STRAT_ORDER) - 1 - i) for i, s in enumerate(STRAT_ORDER)}
DY = {"JHU": 0.20, "Durham": -0.20}
ax_fl.axvline(0.5, color="#DADADA", lw=0.8, ls=(0, (2, 2)), zorder=1)
for s in STRAT_ORDER:
    yb = band[s]
    for coh in ("JHU", "Durham"):
        row = val10[(val10.cutoff_name == s) & (val10.cohort == coh)]
        if row.empty:
            continue
        row = row.iloc[0]; col = COHORT_COLORS[coh]; y = yb + DY[coh]
        # sensitivity (filled) and specificity (open) with 95% CI whiskers
        ax_fl.plot([float(row.Se_lo), float(row.Se_hi)], [y, y], color=col, lw=1.6, zorder=3, solid_capstyle="round")
        ax_fl.scatter([float(row.Se)], [y], marker="o", s=34, facecolor=col, edgecolor="white", linewidth=0.6, zorder=4)
        ax_fl.plot([float(row.Sp_lo), float(row.Sp_hi)], [y, y], color=col, lw=1.6, zorder=3, solid_capstyle="round")
        ax_fl.scatter([float(row.Sp)], [y], marker="o", s=34, facecolor="white", edgecolor=col, linewidth=1.3, zorder=4)
for s in STRAT_ORDER:
    ax_fl.text(-0.02, band[s], s.replace(" (locked)", "\n(locked)"), transform=ax_fl.get_yaxis_transform(),
               ha="right", va="center", fontsize=8.6, color="#222222", linespacing=0.95)
ax_fl.set_ylim(-0.6, len(STRAT_ORDER) - 0.4); ax_fl.set_yticks([])
ax_fl.set_xlim(0.0, 1.0); ax_fl.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0]); ax_fl.tick_params(axis="x", labelsize=9.0)
ax_fl.set_xlabel("Sensitivity / specificity at 10 years (95% CI)", fontsize=9.6)
for s in ("top", "right", "left"):
    ax_fl.spines[s].set_visible(False)
ax_fl.tick_params(axis="y", length=0)
# in-panel metric legend (Se filled / Sp open)
mleg = [Line2D([0], [0], marker="o", color="#555555", lw=0, markersize=7, markerfacecolor="#555555", label="Sensitivity"),
        Line2D([0], [0], marker="o", color="#555555", lw=0, markersize=7, markerfacecolor="white", label="Specificity")]
ax_fl.legend(handles=mleg, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=2, frameon=False,
             fontsize=8.2, handletextpad=0.3, columnspacing=1.2)

# shared cohort legend (bottom)
handles = [Line2D([0], [0], marker="o", color=COHORT_COLORS[c], lw=0, markersize=8, label=COHORT_LABEL[c])
           for c in ("JHU", "Durham")]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10.5,
           bbox_to_anchor=(0.5, 0.012), columnspacing=2.4)

fig.canvas.draw()
for ax, lab in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"), (ax_e, "e"), (ax_fu, "f")]:
    p = ax.get_position()
    fig.text(max(0.004, p.x0 - 0.10), min(0.996, p.y1 + 0.012), lab, **PANEL_KW)

base = os.path.join(FIGD, "FigureS6_sensitivity")
fig.savefig(base + ".pdf", format="pdf", metadata={"CreationDate": None})
fig.savefig(base + ".tiff", format="tiff", dpi=500, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(base + ".png", format="png", dpi=300)
plt.close(fig)
print("wrote", base + ".{pdf,tiff,png}")

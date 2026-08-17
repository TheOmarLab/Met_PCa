"""
Figure S2: Durham secondary endpoints. Plotting consumer only.

Panel a: PCSM Aalen-Johansen cumulative incidence by locked Met-Score class with
         95% bands and a 0/60/120-month at-risk table.
Panel b: BCR, OS, PCSM high-vs-low Cox robust hazard-ratio forest (log scale).
Panel c: fixed-marker MFS 5/10-year time-dependent AUC (Gleason vs Met-Score).
Panel d: fixed-marker 10-year IPCW concordance forest across MFS/BCR/OS/PCSM
         (Gleason vs Met-Score).

Reads only the five aggregate CSVs written by the Durham producer. It does not
fit, bootstrap, threshold, or calculate statistics, and reads no RDA or
patient-level file.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.lines import Line2D


def _find_root():
    env = os.environ.get("MET_PCA_ROOT")
    if env and os.path.isdir(env):
        return env
    d = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "code")) and os.path.isdir(os.path.join(d, "output")):
            return d
        d = os.path.dirname(d)
    raise FileNotFoundError("Could not locate project root; set MET_PCA_ROOT.")


ROOT = _find_root()
DUR = os.path.join(ROOT, "output", "Durham")
CURVE = os.path.join(DUR, "FigureS2_Durham_PCSM_CIF_curve.csv")
SUMMARY = os.path.join(DUR, "FigureS2_Durham_PCSM_summary.csv")
HR_CSV = os.path.join(DUR, "FigureS2_Durham_secondary_HR.csv")
AUC_CSV = os.path.join(DUR, "FigureS2_Durham_MFS_timeAUC.csv")
CIDX_CSV = os.path.join(DUR, "FigureS2_Durham_fixed_marker_Cindex.csv")
OUT_DIR = os.path.join(ROOT, "figures")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.labelsize": 9.5, "axes.titlesize": 10.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.0,
    "axes.edgecolor": "#222222", "text.color": "#000000",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "axes.unicode_minus": False,
})
CLS_COL = {"Low risk": "#2b2eb5", "High risk": "#d55e00"}
# three-series markers with Priyanka/Itzel colours (gray / blue / muted green)
MK_ORDER = ["Gleason", "Met-Score", "Combined"]
MK_COL = {"Gleason": "#6F6F6F", "Met-Score": "#2E5A8C", "Combined": "#6B9080"}


def style_axis(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(False)


def load_and_verify():
    for f in (CURVE, SUMMARY, HR_CSV, AUC_CSV, CIDX_CSV):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"missing aggregate CSV: {f}")
    curve = pd.read_csv(CURVE); summ = pd.read_csv(SUMMARY)
    hr = pd.read_csv(HR_CSV); auc = pd.read_csv(AUC_CSV); cidx = pd.read_csv(CIDX_CSV)
    if not {"class", "month", "cif", "cif_lo", "cif_hi", "n_risk"}.issubset(curve.columns):
        raise ValueError("CIF curve CSV missing columns")
    if set(curve["class"].unique()) != {"Low risk", "High risk"}:
        raise ValueError("CIF curve classes unexpected")
    for cl in ("Low risk", "High risk"):
        m = curve.loc[curve["class"] == cl, "month"].to_numpy()
        if not (len(m) == 121 and (m == np.arange(0, 121)).all()):
            raise ValueError(f"{cl}: months must be 0..120 ascending")
    if not np.isfinite(curve[["cif", "cif_lo", "cif_hi", "n_risk"]].to_numpy()).all():
        raise ValueError("nonfinite in curve")
    s = dict(zip(summ["metric"], summ["value"]))
    for k in ("n_total", "n_low", "n_high", "gray_p", "cif_10y_low", "cif_10y_high",
              "diff_10y_high_minus_low", "diff_10y_boot_lo", "diff_10y_boot_hi"):
        if k not in s or not np.isfinite(float(s[k])):
            raise ValueError(f"summary missing/nonfinite: {k}")
    if int(s["n_total"]) != 555 or int(s["n_low"]) + int(s["n_high"]) != 555:
        raise ValueError("summary tally mismatch")
    if list(hr["endpoint"]) != ["BCR", "OS", "PCSM"]:
        raise ValueError("HR endpoints must be BCR, OS, PCSM in order")
    if not np.isfinite(hr[["HR", "ci_lo_robust", "ci_hi_robust", "p_robust", "n", "events"]].to_numpy()).all():
        raise ValueError("nonfinite in HR")
    if set(auc["marker"].unique()) != set(MK_ORDER) or set(auc["horizon_months"]) != {60, 120}:
        raise ValueError("AUC table markers/horizons unexpected")
    if not np.isfinite(auc[["auc", "ci_lo", "ci_hi"]].to_numpy()).all():
        raise ValueError("nonfinite in AUC")
    if set(cidx["endpoint"]) != {"MFS", "BCR", "OS", "PCSM"} or set(cidx["marker"].unique()) != set(MK_ORDER):
        raise ValueError("C-index endpoints/markers unexpected")
    if not np.isfinite(cidx[["c_index", "ci_lo", "ci_hi"]].to_numpy()).all():
        raise ValueError("nonfinite in C-index")
    return curve, s, hr, auc, cidx


def panel_a(ax, curve, s):
    for cl in ("Low risk", "High risk"):
        d = curve[curve["class"] == cl].sort_values("month"); col = CLS_COL[cl]
        ax.fill_between(d["month"], d["cif_lo"], d["cif_hi"], color=col, alpha=0.15, linewidth=0, step="post")
        lab = f"{'Low-risk' if cl=='Low risk' else 'High-risk'} (n={int(s['n_low']) if cl=='Low risk' else int(s['n_high'])})"
        ax.step(d["month"], d["cif"], where="post", color=col, linewidth=2.0, label=lab)
    ax.set_xlim(0, 120); ax.set_ylim(0, 0.10)
    ax.xaxis.set_major_locator(FixedLocator([0, 30, 60, 90, 120]))
    ax.yaxis.set_major_locator(FixedLocator(np.arange(0, 0.1001, 0.02)))
    ax.set_xlabel("Months from surgery", fontsize=9.5)
    ax.set_ylabel("PCSM cumulative incidence", fontsize=9.5)
    style_axis(ax)
    ax.legend(loc="upper left", frameon=False, fontsize=8.2, handlelength=1.4,
              bbox_to_anchor=(0.0, 1.02), borderaxespad=0.0)
    lab = (f"Gray p = {float(s['gray_p']):.3f}\n"
           f"10-y CIF diff = {float(s['diff_10y_high_minus_low']):.3f} "
           f"({float(s['diff_10y_boot_lo']):.3f} - {float(s['diff_10y_boot_hi']):.3f})")
    ax.text(0.97, 0.60, lab, transform=ax.transAxes, ha="right", va="top", fontsize=8.0)
    ax.text(-0.03, -0.22, "No. at risk", transform=ax.transAxes, ha="right", va="top",
            fontsize=7.8, fontweight="bold", clip_on=False)
    for row, cl in enumerate(("Low risk", "High risk")):
        y = -0.29 - row * 0.075
        d = curve[curve["class"] == cl]
        ax.text(-0.03, y, "Low-risk" if cl == "Low risk" else "High-risk",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.8, color=CLS_COL[cl], clip_on=False)
        for m, ha in ((0, "left"), (60, "center"), (120, "right")):
            nr = int(d.loc[d["month"] == m, "n_risk"].iloc[0])
            ax.text(m / 120.0, y, str(nr), transform=ax.transAxes, ha=ha, va="top", fontsize=7.8, clip_on=False)


def panel_b(ax, hr):
    order = ["BCR", "OS", "PCSM"]; yvals = [2, 1, 0]; xmax = 1.0
    for ep, y in zip(order, yvals):
        r = hr[hr["endpoint"] == ep].iloc[0]
        h, lo, hi = float(r["HR"]), float(r["ci_lo_robust"]), float(r["ci_hi_robust"])
        ax.errorbar(h, y, xerr=[[h - lo], [hi - h]], fmt="o", color="#222222", ecolor="#222222",
                    markersize=6, elinewidth=1.4, capsize=3.5, zorder=3)
        ax.text(h, y + 0.30, f"{h:.2f} ({lo:.2f} - {hi:.2f}), p = {float(r['p_robust']):.3f}",
                ha="center", va="bottom", fontsize=8.0)
        xmax = max(xmax, hi)
    ax.axvline(1.0, color="#888888", linewidth=0.9, linestyle="--", zorder=1)
    ax.set_xscale("log"); ax.set_xlim(0.5, xmax * 1.5)
    ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20]))
    ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20"])
    ax.set_yticks(yvals)
    ax.set_yticklabels([f"{ep}\n({int(hr[hr['endpoint']==ep]['events'].iloc[0])} events)" for ep in order], fontsize=8.5)
    ax.set_ylim(-0.5, 2.8)
    ax.set_xlabel("Robust Cox HR summary (high vs low), log scale", fontsize=9.5)
    ax.spines["left"].set_visible(False)
    style_axis(ax); ax.tick_params(axis="y", length=0)


def shared_marker_legend(fig, y):
    h = [Line2D([0], [0], marker="o", color=MK_COL[m], linestyle="none", markersize=7, label=m)
         for m in MK_ORDER]
    fig.legend(handles=h, loc="lower center", frameon=False, fontsize=8.5, handletextpad=0.3,
               ncol=3, columnspacing=1.8, bbox_to_anchor=(0.5, y))


def panel_c(ax, auc):
    # grouped vertical dot-and-CI at 5 and 10 years; Gleason / Met-Score / Combined
    xc = {60: 0.0, 120: 1.0}; off = {"Gleason": -0.22, "Met-Score": 0.0, "Combined": 0.22}
    for yref in (0.6, 0.7, 0.8, 0.9):
        ax.axhline(yref, color="#EEEEEE", linewidth=0.7, zorder=0)
    ax.axhline(0.5, color="#999999", linewidth=0.7, linestyle=(0, (3, 3)), zorder=1)
    for mk in MK_ORDER:
        for h in (60, 120):
            r = auc[(auc["marker"] == mk) & (auc["horizon_months"] == h)].iloc[0]
            x = xc[h] + off[mk]; a, lo, hi = float(r["auc"]), float(r["ci_lo"]), float(r["ci_hi"])
            ax.plot([x, x], [lo, hi], color=MK_COL[mk], linewidth=1.6, alpha=0.9, zorder=2)
            ax.scatter([x], [a], marker="o", s=44, color=MK_COL[mk], edgecolor="white", linewidth=0.9, zorder=3)
    ax.set_xlim(-0.7, 1.7); ax.set_ylim(0.5, 1.0)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["5 year", "10 year"])
    ax.yaxis.set_major_locator(FixedLocator(np.arange(0.5, 1.001, 0.1)))
    ax.set_ylabel("Time-dependent AUC (MFS)", fontsize=9.5)
    style_axis(ax)
    ax.text(0.5, -0.12, "Durham VA", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=10.0, color="#1a1a1a")


def panel_d(ax, cidx):
    # horizontal concordance forest: endpoints as rows; Gleason / Met-Score / Combined
    eps = ["MFS", "BCR", "OS", "PCSM"]
    ypos = {e: (len(eps) - i) for i, e in enumerate(eps)}   # MFS top .. PCSM bottom
    off = {"Gleason": 0.24, "Met-Score": 0.0, "Combined": -0.24}
    for xref in (0.6, 0.7, 0.8, 0.9):
        ax.axvline(xref, color="#EEEEEE", linewidth=0.7, zorder=0)
    ax.axvline(0.5, color="#999999", linewidth=0.7, linestyle=(0, (3, 3)), zorder=1)
    for mk in MK_ORDER:
        for ep in eps:
            r = cidx[(cidx["endpoint"] == ep) & (cidx["marker"] == mk)].iloc[0]
            y = ypos[ep] + off[mk]; c, lo, hi = float(r["c_index"]), float(r["ci_lo"]), float(r["ci_hi"])
            ax.plot([lo, hi], [y, y], color=MK_COL[mk], linewidth=1.6, alpha=0.9, zorder=2)
            ax.scatter([c], [y], marker="o", s=38, color=MK_COL[mk], edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_xlim(0.45, 1.0); ax.set_ylim(0.4, len(eps) + 0.6)
    ax.xaxis.set_major_locator(FixedLocator(np.arange(0.5, 1.001, 0.1)))
    ax.set_yticks([ypos[e] for e in eps]); ax.set_yticklabels(eps, fontsize=9.0)
    ax.set_xlabel("10-year concordance index", fontsize=9.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.grid(False); ax.tick_params(axis="y", length=0)


def build():
    curve, s, hr, auc, cidx = load_and_verify()
    fig = plt.figure(figsize=(9.8, 8.4), dpi=600, facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.0], height_ratios=[1.0, 1.0],
                          wspace=0.30, hspace=0.68, left=0.10, right=0.975, top=0.95, bottom=0.12)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0]); ax_d = fig.add_subplot(gs[1, 1])
    panel_a(ax_a, curve, s); panel_b(ax_b, hr); panel_c(ax_c, auc); panel_d(ax_d, cidx)
    shared_marker_legend(fig, 0.02)

    fig.canvas.draw()
    _lab = dict(fontsize=15, fontweight="black", family="DejaVu Sans", va="bottom", ha="left")
    for ax_, letter in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d")):
        bb = ax_.get_position()
        fig.text(max(0.004, bb.x0 - 0.058), bb.y1 + 0.012, letter, **_lab)

    os.makedirs(OUT_DIR, exist_ok=True)
    out_pdf = os.path.join(OUT_DIR, "FigureS2_Durham.pdf")
    out_tif = os.path.join(OUT_DIR, "FigureS2_Durham.tiff")
    out_png = os.path.join(OUT_DIR, "FigureS2_Durham.png")
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_tif, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_png, format="png", dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")
    print(f"Saved {out_tif}")
    print(f"Saved {out_png}")


if __name__ == "__main__":
    build()

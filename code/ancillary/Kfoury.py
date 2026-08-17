
# solid metastatic tissue (Tumor)
# liquid BM at the vertebral level of spinal cord compression (Involved)
# liquid BM from a different vertebral body distant from the tumor site but within the surgical field (Distal)

import os
import re
import tarfile
import gzip
from pathlib import Path
import scanpy as sc
from matplotlib.pyplot import ion
import pandas as pd
import scipy.sparse as sp
import loompy as lp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
from scipy.io import mmread
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp


# ============================================================
# CONFIG
# ============================================================
# resolve project root portably: MET_PCA_ROOT env var, else walk up to the dir holding code/ + outs|output
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
ROOT = Path(_find_met_pca_root())

BASE_DIR = ROOT / "data" / "kfoury"
RAW_DIR  = BASE_DIR / "RAW"
META_PATH = BASE_DIR / "GSE143791_cell.annotation.human.csv"

OUTS_DIR = ROOT / "outs" / "kfoury"
FIG_DIR  = ROOT / "figures" / "scRNAseq" / "kfoury"

for d in [BASE_DIR, RAW_DIR, OUTS_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir = str(FIG_DIR)
sc.set_figure_params(dpi=150, dpi_save=400, format="png", transparent=False)

# Met-Score gene lists
POS_PATH = OUTS_DIR / "MetScore_pos_genes.txt"
NEG_PATH = OUTS_DIR / "MetScore_neg_genes.txt"

# ============================================================
# HELPERS
# ============================================================
def parse_sample_from_filename(fname: str):
    """
    Returns dict with:
      sample_id: e.g., GSM4274678_BMET1-Tumor
      gsm: GSM4274678
      patient: BMET1 / BMM4 / Naive_BM1 / etc
      fraction: Tumor / Involved / Distal / Benign / Naive / Other
      cohort: BMET / BMM / Naive / Other
    """
    base = fname.replace(".counts.csv.gz","").replace(".count.csv.gz","")
    # base examples:
    # GSM4274678_BMET1-Tumor
    # GSM4274705_BMM4-Benign
    # GSM5549078_Naive_BM1
    # GSM4490340_SCG_MBMMet1_KO1
    # GSM5551112_Lung.Met1

    gsm = base.split("_")[0]

    # BMET / BMM style
    m = re.match(r"^(GSM\d+)_(BMET\d+|BMM\d+)\-(Tumor|Involved|Distal|Benign)$", base)
    if m:
        gsm, patient, frac = m.group(1), m.group(2), m.group(3)
        cohort = "BMET" if patient.startswith("BMET") else "BMM"
        return {"sample_id": base, "gsm": gsm, "patient": patient, "fraction": frac, "cohort": cohort}

    # Naive_BM style
    m = re.match(r"^(GSM\d+)_Naive_(BM\d+)$", base)
    if m:
        gsm, bm = m.group(1), m.group(2)
        return {"sample_id": base, "gsm": gsm, "patient": f"Naive_{bm}", "fraction": "Naive", "cohort": "Naive"}

    # Everything else
    return {"sample_id": base, "gsm": gsm, "patient": "Other", "fraction": "Other", "cohort": "Other"}

def make_unique(names):
    seen = {}
    out = []
    for x in names:
        if x not in seen:
            seen[x] = 0
            out.append(x)
        else:
            seen[x] += 1
            out.append(f"{x}.{seen[x]}")
    return out

def _read_tsv(path: Path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return pd.read_csv(f, sep="\t", header=None, dtype=str)
    return pd.read_csv(path, sep="\t", header=None, dtype=str)

def _read_mtx(path: Path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return mmread(f).tocsr()
    return mmread(path).tocsr()

def present_genes(adata_obj, genes, use_raw=True):
    genes_up = set([str(g).upper() for g in genes])
    if use_raw and adata_obj.raw is not None:
        var_up = pd.Index([g.upper() for g in adata_obj.raw.var_names])
        var_orig = pd.Index(adata_obj.raw.var_names)
    else:
        var_up = pd.Index([g.upper() for g in adata_obj.var_names])
        var_orig = pd.Index(adata_obj.var_names)

    m = var_up.isin(list(genes_up))
    return var_orig[m].tolist()

def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() > 0:
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out

def fisher_enrichment(up_genes, background_genes, met_genes):
    bg = {g.upper() for g in background_genes}
    up = {g.upper() for g in up_genes} & bg
    met = {g.upper() for g in met_genes} & bg

    a = len(up & met)
    b = len(up - met)
    c = len(met - up)
    d = len(bg - up - met)

    table = np.array([[a, b],
                      [c, d]], dtype=int)
    OR, p = stats.fisher_exact(table, alternative="greater")
    return OR, p, table, dict(a=a, b=b, c=c, d=d, up=len(up), met=len(met), bg=len(bg))

def make_rank_series(de_df, padj_col="pvals_adj", lfc_col="logfoldchanges"):
    df = de_df.copy()
    if lfc_col not in df.columns and "logfoldchange" in df.columns:
        lfc_col = "logfoldchange"
    df = df.dropna(subset=["gene", padj_col, lfc_col]).copy()

    p = df[padj_col].astype(float).replace(0, np.nextafter(0, 1))
    lfc = df[lfc_col].astype(float)
    r = np.sign(lfc) * (-np.log10(p))

    s = pd.Series(r.values, index=df["gene"].astype(str).values)
    s = s.groupby(s.index).apply(lambda x: x.iloc[np.argmax(np.abs(x.values))])
    s = s.sort_values(ascending=False)
    return s

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p

# ============================================================
# 1) LOAD METADATA (and infer sample names, cohort, patient, fraction from sample IDs)
# ============================================================
metadata = pd.read_csv('./data/kfoury/GSE143791_cell.annotation.human.csv', sep=',', header=0)
metadata['cells'].value_counts()
metadata.columns


kfoury_dir = 'data/kfoury/GSE143791_RAW'
filenames = [sample for sample in os.listdir(kfoury_dir) if sample.endswith('.csv.gz')]

meta_rows = [parse_sample_from_filename(f) for f in filenames]
meta_files = pd.DataFrame(meta_rows)

# Keep prostate BM niche relevant cohorts
meta_files_keep = meta_files[meta_files["cohort"].isin(["BMET","BMM","Naive"])].copy()

print(meta_files_keep["cohort"].value_counts())
print(meta_files_keep.groupby(["cohort","fraction"])["sample_id"].nunique())
print(meta_files_keep.head())

# ============================================================
# 2) LOAD RAW MATRICES FOR ALL SAMPLES
# ============================================================
RAW_DIR = Path("data/kfoury/GSE143791_RAW")

adatas = []
for _, r in meta_files_keep.iterrows():
    fpath = RAW_DIR / (r["sample_id"] + ".count.csv.gz")
    if not fpath.exists():
        fpath = RAW_DIR / (r["sample_id"] + ".counts.csv.gz")
    if not fpath.exists():
        raise FileNotFoundError(f"Missing file for {r['sample_id']}")

    a = sc.read_csv(str(fpath)).transpose()   # cells x genes

    # --------- IMPORTANT (from old code + fixes)
    # 2) now make unique (after standardization)
    a.var_names_make_unique()

    # 3) store a stable gene column like your old code
    a.var["gene"] = a.var_names.astype(str)

    # 4) reset var index like your old code (helps avoid merge quirks)
    a.var.reset_index(drop=True, inplace=True)

    # --------- attach sample-level metadata
    a.obs["sample_id"] = r["sample_id"]
    a.obs["gsm"] = r["gsm"]
    a.obs["patient"] = r["patient"]
    a.obs["fraction"] = r["fraction"]
    a.obs["cohort"] = r["cohort"]

    # --------- unique cell IDs across samples
    a.obs_names = [f"{r['sample_id']}_{cid}" for cid in a.obs_names]

    adatas.append(a)

adata = ad.concat(adatas, join="outer",  label = 'sample', merge="first", fill_value=0)
adata.var_names = adata.var["gene"].astype(str).values
adata.var_names_make_unique()
print(adata)
print(adata.obs[["cohort","fraction"]].value_counts())

print("concat n_vars:", adata.n_vars)
print("concat n_obs:", adata.n_obs)
print("X min/max:", float(adata.X.min()), float(adata.X.max()))
print("Example genes present:",
      [g for g in ["SPP1","APOE","C1QA","C1QB","C1QC","FCN1","S100A8","S100A9"] if g in adata.var_names])

# ============================================================
# save the integrated raw adata
# ============================================================


# Write
adata.write('./outs/kfoury/kfoury_raw.h5ad')

## read the adata
adata = sc.read('./outs/kfoury/kfoury_raw.h5ad')
adata


# ============================================================
# Integrate with cell labels from the paper
# ============================================================

adata.obs["cell_bc"] = adata.obs_names.to_series().str.split("_").str[-1].values
adata.obs["sample_prefix_for_meta"] = np.where(
    adata.obs["cohort"].astype(str).eq("Naive"),
    adata.obs["patient"].astype(str),  # "Naive_BM1"
    adata.obs["patient"].astype(str) + "-" + adata.obs["fraction"].astype(str)  # "BMET6-Tumor"
)
adata.obs["barcode_key"] = adata.obs["sample_prefix_for_meta"] + "_" + adata.obs["cell_bc"]
print(adata.obs["barcode_key"].head())
print(metadata["barcode"].head())

metadata2 = metadata.copy()
metadata2["barcode"] = metadata2["barcode"].astype(str)
metadata2 = metadata2.drop_duplicates(subset=["barcode"]).set_index("barcode")

adata.obs["celltype_ref"] = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
hit = adata.obs["barcode_key"].isin(metadata2.index)

adata.obs.loc[hit, "celltype_ref"] = metadata2.loc[
    adata.obs.loc[hit, "barcode_key"], "cells"
].astype(str).values
print("Mapped celltype_ref fraction:", float(np.mean(pd.notna(adata.obs["celltype_ref"]))))
adata = adata[pd.notna(adata.obs["celltype_ref"])].copy()

print(adata.obs["celltype_ref"].value_counts(dropna=False))
print(metadata["cells"].value_counts(dropna=False))


# categorical for plotting
adata.obs["celltype_ref"] = adata.obs["celltype_ref"].astype("category")
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
adata.obs["patient"] = adata.obs["patient"].astype("category")
adata.obs["fraction"] = adata.obs["fraction"].astype("category")
adata.obs["cohort"] = adata.obs["cohort"].astype("category")

#####################################################
## Preprocessing
#####################################################

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var['mt'].value_counts()
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

min_genes = 200
max_genes = 6000
max_mt = 20.0

sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :].copy()
adata = adata[adata.obs["pct_counts_mt"] <= max_mt, :].copy()

sc.pp.filter_genes(adata, min_cells=3)
print("After QC:", adata)

# Store raw counts in a layer before normalization
adata.layers["counts"] = adata.X.copy()

# Normalize + log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()


# ComBat batch correction?
#sc.pp.combat(adata, key='sample')

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=4000,
    layer="counts",
    subset=False,
)

# Regress out effects of total counts per cell
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# scale the data to unit variance.
sc.pp.scale(adata, max_value=10)

# pca
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable = True)

# Computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)

# umap
sc.tl.umap(adata)
sc.pl.umap(adata, color = ['sample'])

# clustering
sc.tl.leiden(adata, resolution=0.3)
adata.obs['leiden'].value_counts()
sc.pl.umap(adata, color = ['leiden', 'celltype_ref'])

# DE genes
sc.tl.rank_genes_groups(adata, 'leiden', method="wilcoxon", use_raw=True)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# ============================================================
# save the annotated adata
# ============================================================

# Write
adata.write('./outs/kfoury/kfoury_annot.h5ad')

## read the adata
adata = sc.read('./outs/kfoury/kfoury_annot.h5ad')

pd.crosstab(adata.obs['patient'], adata.obs['cohort'])
pd.crosstab(adata.obs['cohort'], adata.obs['fraction'])

# ============================================================
# MET-SCORE (ALL CELLS; and later subset by cell types)
# met_score below is the directional 45-gene module score computed per cell as
# met_score_pos - met_score_neg (scanpy score_genes over the positive and negative
# signature gene lists). It is not the frozen 41-feature ridge-logistic Met-Score
# classifier: the bulk logistic probability is not defined for a single cell, so
# the directional module score is used as the single-cell surrogate.
# ============================================================
pos_path = ROOT / "outs" / "MetScore_pos_genes.txt"
neg_path = ROOT / "outs" / "MetScore_neg_genes.txt"

PositiveGenes = pd.read_csv(pos_path, header=None, sep="\t").iloc[:, 0].astype(str).str.upper().tolist()
NegativeGenes = pd.read_csv(neg_path, header=None, sep="\t").iloc[:, 0].astype(str).str.upper().tolist()

pos_in = present_genes(adata, PositiveGenes, use_raw=True)
neg_in = present_genes(adata, NegativeGenes, use_raw=True)
MetScoreGenes_in = pos_in + neg_in


print(f"Met-Score overlap in this dataset: pos={len(pos_in)}, neg={len(neg_in)}")
if len(pos_in) < 5 or len(neg_in) < 5:
    print("[WARN] Low overlap with Met-Score genes; results may be noisy.")

sc.tl.score_genes(adata, pos_in, score_name="met_score_pos", use_raw=True)
sc.tl.score_genes(adata, neg_in, score_name="met_score_neg", use_raw=True)
adata.obs["met_score"] = adata.obs["met_score_pos"] - adata.obs["met_score_neg"]

#################
# Summary of met-score per cell type and fraction
df = adata.obs[["celltype_ref", "fraction", "cohort", "patient",
                "met_score_pos", "met_score_neg"] + (["met_score"] if "met_score" in adata.obs.columns else [])].copy()

df["celltype_ref"] = df["celltype_ref"].astype(str)
df["fraction"] = df["fraction"].astype(str)
if "cohort" in df.columns:
    df["cohort"] = df["cohort"].astype(str)
if "patient" in df.columns:
    df["patient"] = df["patient"].astype(str)

for c in ["met_score_pos", "met_score_neg"] + (["met_score"] if "met_score" in df.columns else []):
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["celltype_ref", "fraction", "met_score_pos", "met_score_neg"])

# Celltype_ref × fraction summary
agg_dict = {
    "n_cells": ("met_score_pos", "size"),

    "met_pos_median": ("met_score_pos", "median"),
    "met_pos_mean":   ("met_score_pos", "mean"),

    "met_neg_median": ("met_score_neg", "median"),
    "met_neg_mean":   ("met_score_neg", "mean"),
}

if "met_score" in df.columns:
    agg_dict.update({
        "met_score_median": ("met_score", "median"),
        "met_score_mean":   ("met_score", "mean"),
    })

summary_ct_frac = (
    df.groupby(["celltype_ref", "fraction"], observed=True)
      .agg(**agg_dict)
      .reset_index()
)

# optional: add derived deltas so you can see which arm dominates
summary_ct_frac["pos_minus_neg_median"] = summary_ct_frac["met_pos_median"] - summary_ct_frac["met_neg_median"]
summary_ct_frac["pos_minus_neg_mean"]   = summary_ct_frac["met_pos_mean"]   - summary_ct_frac["met_neg_mean"]

# sort to make browsing easier
summary_ct_frac = summary_ct_frac.sort_values(["celltype_ref", "fraction"]).reset_index(drop=True)

out1 = OUTS_DIR / "GSE143791_met_score_pos_neg_by_celltype_by_fraction.csv"
summary_ct_frac.to_csv(out1, index=False)
print("Wrote:", out1)
print(summary_ct_frac.head(10))

# Celltype_ref-only summary (collapsed across fractions)
summary_ct = (
    df.groupby(["celltype_ref"], observed=True)
      .agg(
          n_cells=("met_score_pos", "size"),
          met_pos_median=("met_score_pos", "median"),
          met_pos_mean=("met_score_pos", "mean"),
          met_neg_median=("met_score_neg", "median"),
          met_neg_mean=("met_score_neg", "mean"),
          **({"met_score_median": ("met_score", "median"),
              "met_score_mean": ("met_score", "mean")} if "met_score" in df.columns else {})
      )
      .reset_index()
)

summary_ct["pos_minus_neg_median"] = summary_ct["met_pos_median"] - summary_ct["met_neg_median"]
summary_ct["pos_minus_neg_mean"]   = summary_ct["met_pos_mean"]   - summary_ct["met_neg_mean"]

summary_ct = summary_ct.sort_values("n_cells", ascending=False).reset_index(drop=True)

out2 = OUTS_DIR / "GSE143791_met_score_pos_neg_by_celltype.csv"
summary_ct.to_csv(out2, index=False)
print("Wrote:", out2)
print(summary_ct.head(10))

# Fraction-only summary (collapsed across cell types)
# -----------------------------
summary_frac = (
    df.groupby(["fraction"], observed=True)
      .agg(
          n_cells=("met_score_pos", "size"),
          met_pos_median=("met_score_pos", "median"),
          met_pos_mean=("met_score_pos", "mean"),
          met_neg_median=("met_score_neg", "median"),
          met_neg_mean=("met_score_neg", "mean"),
          **({"met_score_median": ("met_score", "median"),
              "met_score_mean": ("met_score", "mean")} if "met_score" in df.columns else {})
      )
      .reset_index()
)
summary_frac["pos_minus_neg_median"] = summary_frac["met_pos_median"] - summary_frac["met_neg_median"]
summary_frac["pos_minus_neg_mean"]   = summary_frac["met_pos_mean"]   - summary_frac["met_neg_mean"]

out3 = OUTS_DIR / "GSE143791_met_score_pos_neg_by_fraction.csv"
summary_frac.to_csv(out3, index=False)
print("Wrote:", out3)
print(summary_frac)

####
# save as multi-sheet excel file for supplementary material
import openpyxl
# Helper for Excel-safe sheet names
def _safe_sheet_name(s: str) -> str:
    # Excel limits: 31 chars, cannot contain: : \ / ? * [ ]
    bad = [":", "\\", "/", "?", "*", "[", "]"]
    for b in bad:
        s = s.replace(b, " ")
    return s[:31]

# -----------------------------
# Write workbook
# -----------------------------
out_xlsx = OUTS_DIR / "GSE143791_MetScore_pos_neg_summaries.xlsx"
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    summary_ct_frac.to_excel(writer, sheet_name=_safe_sheet_name("ByCelltype_ByFraction"), index=False)
    summary_ct.to_excel(writer, sheet_name=_safe_sheet_name("ByCelltype"), index=False)
    summary_frac.to_excel(writer, sheet_name=_safe_sheet_name("ByFraction"), index=False)

print("Wrote Excel workbook:", out_xlsx)

# ============================================================
# ALL-CELLS FIGURES + PAIRED FRACTION TESTS (Involved vs Distal)
# ============================================================
# UMAPs
plt.figure(figsize=(7, 6))
sc.pl.umap(adata, color=["fraction", "celltype_ref"], wspace=0.4, frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_fraction_cells.png", dpi=600, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8, 7))
sc.pl.umap(adata, color=["celltype_ref"],
           legend_loc="on data",
           legend_fontsize=6,
           size=5,
           frameon=False,
           show=False
           )
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_cellTypes.png", dpi=400, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8, 7))
sc.pl.umap(adata, color=["fraction"],
           legend_loc="upper left",
           legend_fontsize=6,
           size=5,
           frameon=False,
           show=False
           )
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_fraction.png", dpi=400, bbox_inches="tight")
plt.close()


plt.figure(figsize=(7, 6))
sc.pl.umap(adata, color=["met_score"], color_map="magma", frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_met_score_all_cells.png", dpi=600, bbox_inches="tight")
plt.close()

plt.figure(figsize=(7, 6))
sc.pl.umap(adata, color=["met_score_pos"], color_map="magma", frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_met_score_pos_all_cells.png", dpi=600, bbox_inches="tight")
plt.close()

################
# dotplot of met-score genes per celltype
dp = sc.pl.DotPlot(adata, var_names = MetScoreGenes_in, groupby='celltype_ref', cmap = 'Reds', use_raw=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig(FIG_DIR / "dotplot_met_score_all_cells.png", dpi=600, bbox_inches="tight")
################

# Per-patient, per-fraction summary across ALL cells
tmp = adata.obs[["cohort", "patient", "fraction", "met_score"]].copy()
tmp["cohort"] = tmp["cohort"].astype(str)
tmp["patient"]   = tmp["patient"].astype(str)
tmp["fraction"]  = tmp["fraction"].astype(str)
tmp["met_score"] = pd.to_numeric(tmp["met_score"], errors="coerce")
tmp = tmp.dropna(subset=["patient", "fraction", "met_score"])

# keep BMET only for these paired marrow/tumor comparisons
tmp_bmet = tmp[tmp["cohort"].eq("BMET")].copy()


all_summary = (
    tmp.groupby(["patient", "fraction"])
       .agg(
           n_cells=("met_score", "size"),
           met_median=("met_score", "median"),
           met_mean=("met_score", "mean"),
       )
       .reset_index()
)

all_summary.to_csv(OUTS_DIR / "GSE143791_allcells_met_score_per_patient_fraction.csv", index=False)

# ----------------------------
# Helper to run paired test + plot
# ----------------------------
def paired_test_and_plot(df_summary, frac_A, frac_B, label, out_png):
    # wide table with both fractions
    wide = df_summary.pivot(index="patient", columns="fraction", values="met_median")
    wide_n = df_summary.pivot(index="patient", columns="fraction", values="n_cells")

    # keep only patients with both fractions present
    wide = wide.dropna(subset=[frac_A, frac_B]).copy()
    wide_n = wide_n.loc[wide.index].copy()

    # paired deltas
    d = (wide[frac_A] - wide[frac_B]).astype(float).values
    d_nonzero = d[d != 0]

    if len(d_nonzero) >= 5:
        W, p = stats.wilcoxon(d_nonzero, alternative="two-sided")
    else:
        W, p = np.nan, np.nan

    print(f"{label}: n paired patients = {wide.shape[0]}")
    print(f"{label}: paired Wilcoxon p = {p}")

    # long form for plot
    plot_df = wide.reset_index().melt(
        id_vars="patient",
        value_vars=[frac_A, frac_B],
        var_name="fraction",
        value_name="met_median"
    )

    # box + points (no lines)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    sns.boxplot(
        data=plot_df,
        x="fraction", y="met_median",
        order=[frac_B, frac_A],   # show “control” first (e.g., Involved/Distal) then Tumor
        width=0.35,
        showcaps=True,
        boxprops={"facecolor":"white", "edgecolor":"black", "linewidth":1.5},
        medianprops={"color":"black", "linewidth":1.8},
        whiskerprops={"color":"black", "linewidth":1.2},
        capprops={"color":"black", "linewidth":1.2},
        ax=ax
    )
    sns.stripplot(
        data=plot_df,
        x="fraction", y="met_median",
        order=[frac_B, frac_A],
        color="black", alpha=0.6, size=3.5, jitter=0.12, ax=ax
    )

    ax.set_xlabel("")
    ax.set_ylabel("Per-patient median 45-gene signature module score (all cells)")
    ax.set_title(label)

    # p-value bracket
    y_max = plot_df["met_median"].max()
    y_min = plot_df["met_median"].min()
    yr = (y_max - y_min) if y_max > y_min else 1.0
    y_bar = y_max + 0.08 * yr
    h = 0.03 * yr

    ax.plot([0, 0, 1, 1], [y_bar, y_bar + h, y_bar + h, y_bar],
            lw=1.5, c="black", clip_on=False)
    ax.text(0.5, y_bar + h + 0.02 * yr,
            f"paired Wilcoxon p = {p:.2e}" if np.isfinite(p) else "paired Wilcoxon p = NA",
            ha="center", va="bottom", fontsize=11)

    ax.set_ylim(y_min - 0.05 * yr, y_bar + h + 0.12 * yr)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close()

    # also save the paired table for reproducibility
    wide_out = wide.copy()
    wide_out["delta"] = wide_out[frac_A] - wide_out[frac_B]
    wide_out = wide_out.reset_index()
    wide_out.to_csv(OUTS_DIR / f"paired_{frac_A}_vs_{frac_B}_allcells.csv", index=False)

    return wide_out, p

# ----------------------------
# Run paired tests: comparing met-score per fraction (all cells)
# ----------------------------
# Tumor vs Involved (paired within BMET patients that have both)
wide_TI, p_TI = paired_test_and_plot(
    all_summary,
    frac_A="Tumor",
    frac_B="Involved",
    label="Tumor vs Involved marrow (paired BMET patients)",
    out_png=FIG_DIR / "paired_Tumor_vs_Involved_allcells.png"
)

# Tumor vs Distal (paired within BMET patients that have both)
wide_TD, p_TD = paired_test_and_plot(
    all_summary,
    frac_A="Tumor",
    frac_B="Distal",
    label="Tumor vs Distal marrow (paired BMET patients)",
    out_png=FIG_DIR / "paired_Tumor_vs_Distal_allcells.png"
)

# Distal vs Involved (paired within BMET patients that have both)
wide_DI, p_DI = paired_test_and_plot(
    all_summary,
    frac_A="Distal",
    frac_B="Involved",
    label="Distal vs Involved marrow (paired BMET patients)",
    out_png=FIG_DIR / "paired_Distal_vs_Involved_allcells.png"
)

# =====================================================================================================================
# MYELOID / TAM-FOCUSED ANALYSIS
# =====================================================================================================================
cells_unique = set(adata.obs["celltype_ref"].astype(str).unique())
print("Unique metadata cell labels:", list(sorted(cells_unique)))

# myeloid/TAMs
MYELOID_CELLTYPES = ["Mono1", "Mono2", "Mono3", "Monocyte prog", "TAM", "TIM", "mDC", "pDC", "Osteoclasts"]

myeloid_mask = adata.obs["celltype_ref"].astype(str).isin(MYELOID_CELLTYPES)
adata_my = adata[myeloid_mask].copy()
print("Myeloid cells:", adata_my.n_obs)
print(adata_my.obs["celltype_ref"].value_counts())

# --- Re-embed myeloid subset using COUNTS -> log1p
adata_my.X = adata_my.layers["counts"].copy()
sc.pp.normalize_total(adata_my, target_sum=1e4)
sc.pp.log1p(adata_my)
adata_my.raw = adata_my

#sc.pp.highly_variable_genes(adata_my, n_top_genes=4000, flavor="seurat_v3", subset=True)
sc.pp.scale(adata_my, max_value=10)
sc.tl.pca(adata_my, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata_my, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata_my)
sc.tl.leiden(adata_my, resolution=0.8, key_added="leiden_my")

plt.figure(figsize=(7, 6))
sc.pl.umap(adata_my, color=["celltype_ref", "leiden_my"], legend_loc="on data", legend_fontsize=8, frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "myeloid_umap_cells_leiden.png", dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Define TAM programs (SPP1/lipid, osteoclast-like, inflammatory)
# ============================================================
SPP1_LIPID_TAM_POS = ["SPP1","TREM2","APOE","GPNMB","LPL","CTSB","CTSD","LGALS3","SLC40A1","MARCO","MSR1","CSTA"]
INFLAM_MONO_POS    = ["FCN1","S100A8","S100A9","VCAN","CCR2"]
OSTEOCLAST_POS     = ["ACP5","CTSK","CALCR","MMP9","SPI1","TYROBP","FCER1G","LST1"]

# Score modules (using raw)
for name, genes in {
    "spp1_lipid": SPP1_LIPID_TAM_POS,
    "inflam_mono": INFLAM_MONO_POS,
    "osteoclast": OSTEOCLAST_POS
}.items():
    gl = present_genes(adata_my, genes, use_raw=True)
    sc.tl.score_genes(adata_my, gl, score_name=f"ms_{name}", use_raw=True)

# Define a TAM mask using either labels or myeloid markers:
def is_tam(lbl: str) -> bool:
    s = str(lbl).lower()
    return ("tam" in s) or ("macro" in s)

adata_my.obs["is_tam"] = [is_tam(x) for x in adata_my.obs["celltype_ref"].astype(str)]
tam_mask = adata_my.obs["is_tam"].astype(bool)
print("TAM-labeled cells:", int(tam_mask.sum()))

# Define SPP1-high / lipid-associated TAMs within TAMs:
# Use a robust within-(patient,fraction) threshold (top quartile) to avoid global cutoffs.
adata_my.obs["SPP1hi_TAM"] = False

if tam_mask.sum() >= 20:
    tmp = adata_my.obs.loc[tam_mask, ["patient", "fraction", "ms_spp1_lipid"]].copy()
    tmp["ms_spp1_lipid"] = pd.to_numeric(tmp["ms_spp1_lipid"], errors="coerce")
    tmp = tmp.dropna()

    # compute q75 per patient×fraction
    q = tmp.groupby(["patient","fraction"])["ms_spp1_lipid"].quantile(0.75).rename("q75").reset_index()
    tmp = tmp.merge(q, on=["patient","fraction"], how="left")
    hit_idx = tmp.index[tmp["ms_spp1_lipid"] >= tmp["q75"]]
    # map back to adata_my obs indices
    # tmp index is aligned with adata_my.obs[tam_mask].index after .loc; reconstruct:
    tam_obs_idx = adata_my.obs.index[tam_mask]
    adata_my.obs.loc[tam_obs_idx[hit_idx], "SPP1hi_TAM"] = True

# myeloid_state label without changing "cells"
def my_state(row):
    if bool(row["is_tam"]):
        return r"SPP1$^{hi}$ TAM" if bool(row["SPP1hi_TAM"]) else "Other TAM"
    return "Other myeloid"

adata_my.obs["myeloid_state"] = adata_my.obs.apply(my_state, axis=1).astype("category")

# UMAP overlays
plt.figure(figsize=(7, 6))
sc.pl.umap(adata_my, color=["myeloid_state", "ms_spp1_lipid", "ms_inflam_mono", "ms_osteoclast"],
           wspace=0.4, frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "myeloid_umap_states_programs.png", dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Is Met-Score higher in SPP1hi TAMs? (patient×fraction paired)
# ============================================================
# Ensure met_score exists in adata_my (it should, inherited from parent)
if "met_score" not in adata_my.obs.columns:
    raise ValueError("met_score missing in myeloid subset; compute on full adata before subsetting.")

df_tam = adata_my.obs.loc[tam_mask, ["patient","fraction","myeloid_state","met_score"]].copy()
df_tam["met_score"] = pd.to_numeric(df_tam["met_score"], errors="coerce")
df_tam = df_tam.dropna()

# Per patient×fraction×state medians
tmp = adata_my.obs[["patient", "fraction", "met_score", "myeloid_state"]].copy()
tmp["patient"]   = tmp["patient"].astype(str)
tmp["fraction"]  = tmp["fraction"].astype(str)
tmp["myeloid_state"]  = tmp["myeloid_state"].astype(str)
tmp["met_score"] = pd.to_numeric(tmp["met_score"], errors="coerce")
tmp = tmp.dropna(subset=["patient", "fraction", "met_score", "myeloid_state"])

tam_summary = (
    tmp.groupby(["patient","fraction","myeloid_state"])
          .agg(n_cells=("met_score","size"),
               met_median=("met_score","median"),
               met_mean=("met_score","mean"),
          )
           .reset_index()
)

tam_summary.to_csv(OUTS_DIR / "GSE143791_TAM_met_score_per_patient_fraction_state.csv", index=False)

# Wide for paired within each patient×fraction
wide_tam = tam_summary.pivot_table(
    index=["patient","fraction"],
    columns="myeloid_state",
    values="met_median",
    aggfunc="first"
).reset_index()

if r"SPP1$^{hi}$ TAM" in wide_tam.columns and "Other TAM" in wide_tam.columns:
    wide_tam = wide_tam.dropna(subset=[r"SPP1$^{hi}$ TAM", "Other TAM"]).copy()
    wide_tam["delta_med"] = wide_tam[r"SPP1$^{hi}$ TAM"] - wide_tam["Other TAM"]
else:
    print("[WARN] Could not form paired SPP1hi vs Other TAM table (missing columns).")

wide_tam.to_csv(OUTS_DIR / "GSE143791_TAM_met_score_WIDE_patient_fraction.csv", index=False)

# Paired Wilcoxon within patient×fraction rows (this is the clean paired unit here)
if "delta_med" in wide_tam.columns and wide_tam.shape[0] >= 10:
    d = wide_tam["delta_med"].astype(float).values
    d_nonzero = d[d != 0]
    W, p_tam = stats.wilcoxon(d_nonzero) if len(d_nonzero) >= 5 else (np.nan, np.nan)
else:
    p_tam = np.nan

# Plot (no lines)
if "delta_med" in wide_tam.columns and wide_tam.shape[0] >= 5:
    plot_long = wide_tam.melt(
        id_vars=["patient","fraction"],
        value_vars=["Other TAM", r"SPP1$^{hi}$ TAM"],
        var_name="TAM_state",
        value_name="met_median"
    )
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    sns.boxplot(
        data=plot_long,
        x="TAM_state", y="met_median",
        order=["Other TAM", r"SPP1$^{hi}$ TAM"],
        color="white",
        boxprops={"edgecolor":"black", "linewidth":1.5},
        medianprops={"color":"black", "linewidth":1.8},
        whiskerprops={"color":"black", "linewidth":1.2},
        capprops={"color":"black", "linewidth":1.2},
        ax=ax
    )
    sns.stripplot(
        data=plot_long,
        x="TAM_state", y="met_median",
        order=["Other TAM", r"SPP1$^{hi}$ TAM"],
        color="black", alpha=0.6, size=3.5, jitter=0.12, ax=ax
    )
    ax.set_xlabel("")
    ax.set_ylabel("Per patient×fraction median 45-gene signature module score (TAMs)")
    ax.set_title("45-gene signature module score in SPP1-high vs other TAMs")

    y_max = plot_long["met_median"].max()
    y_min = plot_long["met_median"].min()
    yr = y_max - y_min if y_max > y_min else 1.0
    y_bar = y_max + 0.10*yr
    h = 0.03*yr
    ax.plot([0,0,1,1],[y_bar,y_bar+h,y_bar+h,y_bar], lw=1.5, c="black", clip_on=False)
    if np.isfinite(p_tam):
        ax.text(0.5, y_bar+h+0.02*yr, f"paired Wilcoxon p = {p_tam:.2e}", ha="center", va="bottom", fontsize=11)
    ax.set_ylim(y_min - 0.05*yr, y_bar+h+0.15*yr)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "TAM_met_score_paired_patient_fraction.png", dpi=600, bbox_inches="tight")
    plt.close()


# ============================================================
# Cell-type–stratified per-patient summaries + paired tests
#   For each patient × fraction, compute median Met-Score within:
#     - TAM/TIM
#     - SPP1-high / lipid-TAM program-high TAM/TIM
#     - Osteoclasts
#     - Monocytes (Mono1/2/3 / Monocyte prog)
#   Then paired tests within patient:
#     - Involved vs Distal
#     - Tumor vs Involved
#     - Tumor vs Distal
#
# ============================================================
# ----------------------------
# 0) Build a clean analysis table (BMET only)
# ----------------------------
req = ["cohort","patient","fraction","celltype_ref","met_score"]
missing = [c for c in req if c not in adata.obs.columns]
if missing:
    raise KeyError(f"Missing obs columns: {missing}")

tmp = adata.obs[req].copy()
tmp["cohort"] = tmp["cohort"].astype(str)
tmp["patient"] = tmp["patient"].astype(str)
tmp["fraction"] = tmp["fraction"].astype(str)
tmp["celltype_ref"] = tmp["celltype_ref"].astype(str)
tmp["met_score"] = pd.to_numeric(tmp["met_score"], errors="coerce")
tmp = tmp.dropna(subset=["patient","fraction","celltype_ref","met_score"])

# For fraction contrasts Tumor/Involved/Distal, only BMET actually has these fractions.
tmp_bmet = tmp[(tmp["cohort"] == "BMET") & (tmp["fraction"].isin(["Tumor","Involved","Distal"]))].copy()

print("BMET fraction counts:\n", tmp_bmet["fraction"].value_counts())
print("BMET patients:", tmp_bmet["patient"].nunique())
print("Celltypes:", tmp_bmet["celltype_ref"].nunique())


# ----------------------------
# 1) Per-patient × fraction × celltype_ref medians
# ----------------------------
ct_summary = (
    tmp_bmet
    .groupby(["patient","fraction","celltype_ref"], observed=True)
    .agg(
        n_cells=("met_score","size"),
        met_median=("met_score","median"),
        met_mean=("met_score","mean"),
    )
    .reset_index()
)
ct_summary.to_csv(OUTS_DIR / "GSE143791_BMET_met_score_per_patient_fraction_celltype.csv", index=False)

# ----------------------------
# 2) Paired Wilcoxon per cell type
#    (uses per-patient medians; avoids pseudo-replication)
# ----------------------------
CONTRASTS = [("Involved","Distal"), ("Tumor","Involved"), ("Tumor","Distal")]
MIN_CELLS_PER_PTFRAC = 10

def paired_wilcoxon_for_celltype(df, celltype, fracA, fracB, min_cells=10):
    d = df[df["celltype_ref"].eq(celltype)].copy()
    d = d[d["n_cells"] >= min_cells].copy()

    wide = d.pivot_table(index="patient", columns="fraction", values="met_median", aggfunc="first")
    if (fracA not in wide.columns) or (fracB not in wide.columns):
        return None

    wide = wide.dropna(subset=[fracA, fracB]).copy()
    if wide.shape[0] < 5:  # too few paired patients
        return {"celltype_ref": celltype, "contrast": f"{fracA} vs {fracB}",
                "n_paired": int(wide.shape[0]), "W": np.nan, "p": np.nan, "median_delta": np.nan}

    delta = (wide[fracA] - wide[fracB]).astype(float).values
    delta_nz = delta[delta != 0]

    if len(delta_nz) < 5:
        W, p = np.nan, np.nan
    else:
        W, p = stats.wilcoxon(delta_nz, alternative="two-sided")

    return {"celltype_ref": celltype, "contrast": f"{fracA} vs {fracB}",
            "n_paired": int(wide.shape[0]), "W": W, "p": p, "median_delta": float(np.median(delta))}

results = []
celltypes = sorted(ct_summary["celltype_ref"].unique().tolist())

for ct in celltypes:
    for a, b in CONTRASTS:
        out = paired_wilcoxon_for_celltype(ct_summary, ct, a, b, min_cells=MIN_CELLS_PER_PTFRAC)
        if out is not None:
            results.append(out)

res = pd.DataFrame(results)
# BH adjust within each contrast
res["p_adj"] = np.nan
for c in res["contrast"].unique():
    m = res["contrast"].eq(c) & res["p"].notna()
    if m.sum() > 0:
        res.loc[m, "p_adj"] = multipletests(res.loc[m, "p"].values, method="fdr_bh")[1]

res.to_csv(OUTS_DIR / "GSE143791_BMET_paired_tests_by_celltype.csv", index=False)

print("Saved:")
print(" -", OUTS_DIR / "GSE143791_BMET_met_score_per_patient_fraction_celltype.csv")
print(" -", OUTS_DIR / "GSE143791_BMET_paired_tests_by_celltype.csv")

# ----------------------------
# 3) Plot heatmap of median deltas per celltype (quick read)
# ----------------------------
# Make a delta table (Tumor-Involved, Tumor-Distal, Involved-Distal) by celltype (median of per-patient deltas)
delta_rows = []
for ct in celltypes:
    dct = ct_summary[(ct_summary["celltype_ref"]==ct) & (ct_summary["n_cells"]>=MIN_CELLS_PER_PTFRAC)].copy()
    wide = dct.pivot_table(index="patient", columns="fraction", values="met_median", aggfunc="first")

    def med_delta(a,b):
        if a in wide.columns and b in wide.columns:
            w = wide.dropna(subset=[a,b])
            if w.shape[0] >= 5:
                return float(np.median((w[a]-w[b]).values))
        return np.nan

    delta_rows.append({
        "celltype_ref": ct,
        "Tumor-Involved": med_delta("Tumor","Involved"),
        "Tumor-Distal": med_delta("Tumor","Distal"),
        "Involved-Distal": med_delta("Involved","Distal"),
    })

delta_df = pd.DataFrame(delta_rows).set_index("celltype_ref")
delta_df.to_csv(OUTS_DIR / "GSE143791_BMET_median_delta_by_celltype.csv")

# Plot heatmap
plt.figure(figsize=(6.5, 0.28*len(delta_df) + 2.0))
sns.heatmap(delta_df, center=0, cmap="RdBu_r", linewidths=0.0, cbar_kws={"label":"Median Δ 45-gene signature module score"})
plt.title("")
plt.ylabel("")
plt.xlabel("")
plt.xticks(rotation=45)
plt.grid(False)
plt.tight_layout()
plt.savefig(FIG_DIR / "BMET_celltype_delta_heatmap.png", dpi=600, bbox_inches="tight")
plt.close()

print("Heatmap saved:", FIG_DIR / "BMET_celltype_delta_heatmap.png")

# Plot top hits per contrast
def paired_boxplot_celltype(ct_summary, celltype, fracA="Tumor", fracB="Distal",
                            min_cells=20, title_prefix="", out_path=None, p_override=None):
    d = ct_summary[(ct_summary["celltype_ref"] == celltype) & (ct_summary["n_cells"] >= min_cells)].copy()
    if d.empty:
        return None

    wide = d.pivot_table(index="patient", columns="fraction", values="met_median", aggfunc="first")
    if (fracA not in wide.columns) or (fracB not in wide.columns):
        return None
    wide = wide.dropna(subset=[fracA, fracB]).copy()
    if wide.shape[0] < 5:
        return None

    delta = (wide[fracA] - wide[fracB]).astype(float).values
    delta_nz = delta[delta != 0]
    if p_override is not None and np.isfinite(p_override):
        p = float(p_override)
    else:
        p = stats.wilcoxon(delta_nz, alternative="two-sided")[1] if len(delta_nz) >= 5 else np.nan

    plot_long = wide.reset_index().melt(
        id_vars="patient",
        value_vars=[fracB, fracA],  # order control then test
        var_name="fraction",
        value_name="met_median"
    )

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    sns.boxplot(
        data=plot_long, x="fraction", y="met_median",
        order=[fracB, fracA],
        color="white",
        width=0.38,
        boxprops={"edgecolor":"black", "linewidth":1.5},
        medianprops={"color":"black", "linewidth":1.8},
        whiskerprops={"color":"black", "linewidth":1.2},
        capprops={"color":"black", "linewidth":1.2},
        ax=ax
    )
    sns.stripplot(
        data=plot_long, x="fraction", y="met_median",
        order=[fracB, fracA],
        color="black", alpha=0.65, size=3.5, jitter=0.12, ax=ax
    )

    ax.set_xlabel("")
    ax.set_ylabel("Per-patient median 45-gene signature module score")
    ax.set_title(f"{title_prefix}{celltype} ({fracA} vs {fracB})")

    # p-value bracket
    y_max = plot_long["met_median"].max()
    y_min = plot_long["met_median"].min()
    yr = (y_max - y_min) if y_max > y_min else 1.0
    y_bar = y_max + 0.10 * yr
    h = 0.03 * yr
    ax.plot([0, 0, 1, 1], [y_bar, y_bar + h, y_bar + h, y_bar], lw=1.3, c="black", clip_on=False)
    ax.text(
        0.5, y_bar + h + 0.02 * yr,
        f"paired Wilcoxon p = {p:.2e}" if np.isfinite(p) else "paired Wilcoxon p = NA",
        ha="center", va="bottom", fontsize=10
    )
    ax.set_ylim(y_min - 0.05 * yr, y_bar + h + 0.15 * yr)
    sns.despine(ax=ax)
    plt.tight_layout()

    if out_path is not None:
        plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()

    return {"celltype_ref": celltype, "fracA": fracA, "fracB": fracB, "n_paired": int(wide.shape[0]), "p": p,
            "median_delta": float(np.median(delta))}

# -----------------------
# Choose which contrast(s) to plot
# -----------------------
PLOT_CONTRASTS = [
    ("Tumor vs Distal", "Tumor", "Distal"),
    ("Tumor vs Involved", "Tumor", "Involved"),
]

# Select BH-significant cell types for each contrast and plot them
plot_summaries = []
for contrast_label, fracA, fracB in PLOT_CONTRASTS:
    sig = res[(res["contrast"] == contrast_label) & (res["p_adj"].notna()) & (res["p_adj"] < 0.05)].copy()
    sig = sig.sort_values("p_adj")
    print(f"{contrast_label}: {sig.shape[0]} cell types with BH p_adj < 0.05")

    for _, r in sig.iterrows():
        ct = r["celltype_ref"]
        out = FIG_DIR / f"paired_{fracA}_vs_{fracB}_{ct.replace('/','_').replace(' ','_')}.png"
        s = paired_boxplot_celltype(
            ct_summary,
            ct,
            fracA=fracA,
            fracB=fracB,
            min_cells=MIN_CELLS_PER_PTFRAC,
            title_prefix="GSE143791 BMET: ",
            out_path=out,
            p_override=float(r["p"]) if np.isfinite(r["p"]) else None  # use your computed p
        )
        if s is not None:
            s["contrast"] = contrast_label
            s["p_adj"] = float(r["p_adj"]) if np.isfinite(r["p_adj"]) else np.nan
            plot_summaries.append(s)

plot_summaries_df = pd.DataFrame(plot_summaries)
plot_summaries_df.to_csv(Path(OUTS_DIR) / "GSE143791_BMET_paired_plots_generated.csv", index=False)
print("Saved plot list:", Path(OUTS_DIR) / "GSE143791_BMET_paired_plots_generated.csv")
print("Plots saved in:", FIG_DIR)



###############################
#  Paired plots (no lines) for the significant cell types
MIN_CELLS_PER_PTFRAC = 20  # must match what you used when generating res/ct_summary
CONTRASTS = [("Tumor", "Distal"), ("Tumor", "Involved"), ("Involved", "Distal")]

def _median_delta_by_celltype(ct_summary, celltype, fracA, fracB, min_cells=20, min_pairs=5):
    d = ct_summary[(ct_summary["celltype_ref"] == celltype) & (ct_summary["n_cells"] >= min_cells)].copy()
    if d.empty:
        return np.nan
    wide = d.pivot_table(index="patient", columns="fraction", values="met_median", aggfunc="first")
    if (fracA not in wide.columns) or (fracB not in wide.columns):
        return np.nan
    wide = wide.dropna(subset=[fracA, fracB])
    if wide.shape[0] < min_pairs:
        return np.nan
    return float(np.median((wide[fracA] - wide[fracB]).values))

# Build delta matrix
celltypes = sorted(ct_summary["celltype_ref"].unique().tolist())
delta_rows = []
for ct in celltypes:
    delta_rows.append({
        "celltype_ref": ct,
        "Tumor–Distal": _median_delta_by_celltype(ct_summary, ct, "Tumor", "Distal",
                                                 min_cells=MIN_CELLS_PER_PTFRAC),
        "Tumor–Involved": _median_delta_by_celltype(ct_summary, ct, "Tumor", "Involved",
                                                    min_cells=MIN_CELLS_PER_PTFRAC),
        "Involved–Distal": _median_delta_by_celltype(ct_summary, ct, "Involved", "Distal",
                                                     min_cells=MIN_CELLS_PER_PTFRAC),
    })
delta_df = pd.DataFrame(delta_rows).set_index("celltype_ref")

# Attach p_adj into labels (for the 3 contrasts)
# res.contrast strings are like "Tumor vs Distal" etc
pmap = {}
for _, r in res.dropna(subset=["p_adj"]).iterrows():
    pmap[(r["celltype_ref"], r["contrast"])] = float(r["p_adj"])

def _fmt_p(x):
    if x is None or (not np.isfinite(x)):
        return "NA"
    # concise: 2 sig fig in exponent form when small
    return f"{x:.2g}" if x >= 0.001 else f"{x:.1e}"

# Create a row label that shows the best (smallest) adjusted p among the three contrasts
row_labels = []
for ct in delta_df.index:
    p_td = pmap.get((ct, "Tumor vs Distal"), np.nan)
    p_ti = pmap.get((ct, "Tumor vs Involved"), np.nan)
    p_id = pmap.get((ct, "Involved vs Distal"), np.nan)
    best = np.nanmin([p_td, p_ti, p_id])
    best_str = _fmt_p(best) if np.isfinite(best) else "NA"
    row_labels.append(f"{ct}  (min BH p={best_str})")

delta_plot = delta_df.copy()
delta_plot.index = row_labels

# Plot heatmap
plt.figure(figsize=(7.2, 0.30 * len(delta_plot) + 2.0))
sns.heatmap(
    delta_plot,
    center=0,
    cmap="RdBu_r",
    linewidths=0.0,
    cbar_kws={"label": "Median Δ 45-gene signature module score (paired, per-patient)"},
)
plt.title(f"GSE143791 BMET: per-patient median Δ 45-gene signature module score by cell type (min n_cells={MIN_CELLS_PER_PTFRAC})")
plt.ylabel("")
plt.xlabel("")
plt.xticks(rotation=35, ha="right")
plt.tight_layout()
out_heat = FIG_DIR / "BMET_celltype_delta_heatmap_with_padj.png"
plt.savefig(out_heat, dpi=600, bbox_inches="tight")
plt.close()

print("Saved:", out_heat)
delta_df.to_csv(Path(OUTS_DIR) / "GSE143791_BMET_median_delta_by_celltype.csv")

# ============================================================
# 9) PSEUDOBULK DE + Fisher + GSEA
#    DE: SPP1hi TAMs vs Other TAMs
# ============================================================
# Build pseudobulk counts per patient×fraction×TAM_state
# Use counts from parent object for stability; we stored counts in adata.layers["counts"].
# We'll subset to myeloid TAM cells and sum counts.

# Build pseudobulk counts per patient×fraction×TAM_state
# Use counts from parent object for stability; we stored counts in adata.layers["counts"].
# We'll subset to myeloid TAM cells and sum counts.

# Create a TAM-only AnnData with counts as X
adata_tam = adata_my[tam_mask].copy()
adata_tam.X = adata_tam.layers["counts"].copy()
adata_tam.X.min()
adata_tam.X.max()

grp_hi = r"SPP1$^{hi}$ TAM"
grp_ot = "Other TAM"

# Only keep the two groups
adata_tam = adata_tam[adata_tam.obs["myeloid_state"].isin([r"SPP1$^{hi}$ TAM", "Other TAM"])].copy()
print("TAM cells for pseudobulk:", adata_tam.n_obs)
print(adata_tam.obs["myeloid_state"].value_counts())
print(adata_tam.obs["fraction"].value_counts())

# Aggregate counts
def pseudobulk_sum(adata_in, group_cols):
    X = adata_in.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    gdf = adata_in.obs[group_cols].astype(str).copy()
    gdf["__grp__"] = gdf.apply(lambda r: "|".join(r.values.tolist()), axis=1)
    groups = gdf["__grp__"].values
    uniq = pd.Index(pd.unique(groups))

    codes = pd.Categorical(groups, categories=uniq).codes
    G = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                      shape=(len(codes), len(uniq)))

    Xsum = (G.T @ X).tocsr()  # groups x genes

    obs = gdf.drop(columns="__grp__").copy()
    obs = obs.groupby(groups, as_index=False).first()
    obs.index = uniq

    return ad.AnnData(X=Xsum, obs=obs, var=adata_in.var.copy())

def run_de_fraction(pb, frac, out_csv):
    sub = pb[pb.obs["fraction"].astype(str) == frac].copy()
    sub.obs["grp"] = sub.obs["myeloid_state"].astype(str).astype("category")

    if grp_hi not in sub.obs["grp"].cat.categories or grp_ot not in sub.obs["grp"].cat.categories:
        print(f"[SKIP] {frac}: missing one of the groups.")
        return None

    # Require at least 2 pseudobulk replicates per group for any reasonable DE
    vc = sub.obs["grp"].value_counts()
    if vc.get(grp_hi, 0) < 2 or vc.get(grp_ot, 0) < 2:
        print(f"[SKIP] {frac}: too few pseudobulk replicates per group:\n{vc}")
        return None

    # Normalize logCPM
    sub.layers["counts"] = sub.X.copy()
    sc.pp.normalize_total(sub, target_sum=1e6)
    sc.pp.log1p(sub)
    sub.raw = sub

    # Prefer Wilcoxon for robustness
    sc.tl.rank_genes_groups(
        sub,
        groupby="grp",
        groups=[grp_hi],
        reference=grp_ot,
        method="wilcoxon",
        use_raw=True,
        key_added=f"rg_{frac}"
    )

    de = sc.get.rank_genes_groups_df(sub, group=grp_hi, key=f"rg_{frac}").copy()
    de.rename(columns={"names": "gene"}, inplace=True)
    de["fraction"] = frac
    de.to_csv(out_csv, index=False)
    return de


pb = pseudobulk_sum(adata_tam, ["patient","fraction","myeloid_state"])
print("Pseudobulk shape:", pb)
print("Pseudobulk rows by fraction × state:\n",
      pb.obs.groupby(["fraction","myeloid_state"]).size())

pb.obs["fraction"] = pb.obs["fraction"].astype(str)
pb.obs["grp"] = pb.obs["myeloid_state"].astype(str)

print(pb.obs.groupby(["fraction","grp"]).size().unstack(fill_value=0))

# Run DE for each fraction of interest
fractions_to_test = ["Tumor", "Involved", "Distal", "Benign"]

de_all = []

for frac in sorted(pb.obs["fraction"].unique()):
    pb_f = pb[pb.obs["fraction"].astype(str) == frac].copy()
    pb_f.obs["grp"] = pb_f.obs["myeloid_state"].astype(str)
    pb_f.obs["grp"] = pb_f.obs["grp"].astype("category")

    # Need both groups and enough replicates
    vc = pb_f.obs["grp"].value_counts()
    if (r"SPP1$^{hi}$ TAM" not in vc.index) or ("Other TAM" not in vc.index):
        print("Skip", frac, "missing group")
        continue
    if (vc[r"SPP1$^{hi}$ TAM"] < 3) or (vc["Other TAM"] < 3):
        print("Skip", frac, "too few replicates:", vc.to_dict())
        continue

    sc.tl.rank_genes_groups(
        pb_f,
        groupby="grp",
        groups=[r"SPP1$^{hi}$ TAM"],
        reference="Other TAM",
        method="t-test_overestim_var",   # OK for pseudobulk, but see note below
        key_added="rg"
    )
    d = sc.get.rank_genes_groups_df(pb_f, group=r"SPP1$^{hi}$ TAM", key="rg").copy()
    d.rename(columns={"names":"gene"}, inplace=True)
    d["fraction"] = frac
    de_all.append(d)

de_all = pd.concat(de_all, ignore_index=True)


# -----------------------------
# Fisher per fraction using ADJ p-values
# -----------------------------
MetScoreGenes = list(set(PositiveGenes + NegativeGenes))

def fisher_enrichment(up_genes, background_genes, met_genes):
    bg = {g.upper() for g in background_genes}
    up = {g.upper() for g in up_genes} & bg
    met = {g.upper() for g in met_genes} & bg

    a = len(up & met)
    b = len(up - met)
    c = len(met - up)
    d = len(bg - up - met)

    # avoid OR NaN when a=b=0
    if (a + b) == 0:
        return 0.0, 1.0, dict(a=a,b=b,c=c,d=d, up=len(up), met=len(met), bg=len(bg))

    table = np.array([[a, b],[c, d]], dtype=int)
    OR, p = stats.fisher_exact(table, alternative="greater")
    return OR, p, dict(a=a,b=b,c=c,d=d, up=len(up), met=len(met), bg=len(bg))

fisher_rows = []
TOP_N = 200

for frac in sorted(de_all["fraction"].unique()):
    d = de_all[de_all["fraction"] == frac].copy()

    # define “up genes” = top N positive logFC (no p-value threshold)
    d = d.dropna(subset=["gene","logfoldchanges"]).copy()
    d = d.sort_values("logfoldchanges", ascending=False)
    up = d.loc[d["logfoldchanges"] > 0, "gene"].head(TOP_N).astype(str).tolist()

    bg = d["gene"].astype(str).tolist()
    OR, p, counts = fisher_enrichment(up, bg, MetScoreGenes)

    fisher_rows.append({"fraction": frac, "TOP_N": TOP_N, "OR": OR, "p_fisher": p, **counts})

fisher_df = pd.DataFrame(fisher_rows)
print(fisher_df[["fraction","TOP_N","up","met","a","OR","p_fisher"]])

# -----------------------------
# GSEA prerank
# -----------------------------
MetScoreGenes = list(set(PositiveGenes + NegativeGenes))
met_set = {"Met-Score": sorted({g.upper() for g in MetScoreGenes})}

def make_rank_series(df):
    df = df.dropna(subset=["gene","logfoldchanges","pvals"]).copy()
    p = df["pvals"].astype(float).replace(0, np.nextafter(0,1))
    lfc = df["logfoldchanges"].astype(float)
    r = np.sign(lfc) * (-np.log10(p))
    s = pd.Series(r.values, index=df["gene"].astype(str).values)

    # collapse duplicate genes by max abs rank
    s = s.groupby(s.index).apply(lambda x: x.iloc[np.argmax(np.abs(x.values))])

    # break ties slightly (prevents gseapy warnings / crashes)
    rng = np.random.default_rng(7)
    s = s + rng.normal(0, 1e-8, size=len(s))

    return s.sort_values(ascending=False)

for frac in sorted(de_all["fraction"].unique()):
    d = de_all[de_all["fraction"] == frac].copy()
    r = make_rank_series(d)
    rnk = pd.DataFrame({"gene": [g.upper() for g in r.index], "score": r.values})

    outdir = OUTS_DIR / f"GSE143791_GSEA_{frac}_SPP1hi_vs_Other"
    outdir.mkdir(parents=True, exist_ok=True)

    g = gp.prerank(
        rnk=rnk,
        gene_sets=met_set,
        threads=4,              # gseapy now prefers threads
        permutation_num=1000,
        seed=7,
        outdir=str(outdir),
        min_size=5,
        max_size=5000,
        verbose=False
    )
    g.res2d.to_csv(outdir / "gsea_res2d.csv")


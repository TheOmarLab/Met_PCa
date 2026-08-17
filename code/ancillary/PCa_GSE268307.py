import os
import tarfile
import gzip
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmread
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad

# --------------------------------
# 0) Paths + settings
# --------------------------------
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

BASE_DIR = ROOT / "data" / "GSE268307"
SUPP_DIR = BASE_DIR / "supp"
RAW_DIR  = BASE_DIR / "RAW"      # where tar will be extracted
OUTS_DIR = ROOT / "outs"
FIG_DIR  = ROOT / "figures" / "scRNAseq" / "GSE268307"

for d in [SUPP_DIR, RAW_DIR, OUTS_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir = str(FIG_DIR)
sc.set_figure_params(dpi=150, dpi_save=400, format="png", transparent=False)

# --------------------------------
# 1) Download GSE268307_RAW.tar
# --------------------------------
def download_file(url: str, out_path: Path, chunk_size: int = 1 << 20):
    import requests
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)

tar_path = SUPP_DIR / "GSE268307_RAW.tar"

if not tar_path.exists():
    # GEO standard FTP layout
    # If this ever changes, just paste the (http) link from GEO here.
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268307/suppl/GSE268307_RAW.tar"
    print("Downloading:", url)
    download_file(url, tar_path)
else:
    print("Found existing:", tar_path)

# --------------------------------
# 2) Untar into RAW_DIR
# --------------------------------
# (If you re-run, this is idempotent-ish; you can delete RAW_DIR if you want a clean re-extract.)
print("Extracting tar ->", RAW_DIR)
with tarfile.open(tar_path, "r") as tf:
    tf.extractall(path=RAW_DIR)

# --------------------------------
# 3) Find 10x triplets inside extracted files
# --------------------------------
# Supports:
#   - matrix.mtx / matrix.mtx.gz
#   - features.tsv or genes.tsv (gz or not)
#   - barcodes.tsv (gz or not)
def _read_text_table(path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return pd.read_csv(f, sep="\t", header=None, dtype=str)
    else:
        return pd.read_csv(path, sep="\t", header=None, dtype=str)

def _read_mtx(path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return mmread(f).tocsr()
    else:
        return mmread(path).tocsr()

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

def discover_10x_runs_geo_style(root: Path):
    """
    Handles GEO raw naming like:
      GSMxxxx_matrix.mtx.gz
      GSMxxxx_features.tsv.gz  (or genes.tsv.gz)
      GSMxxxx_barcodes.tsv.gz
    Also still supports canonical 10x folder naming if present.
    """
    runs = []

    # --- Case A: Canonical 10x folders (matrix.mtx(.gz) etc.)
    mtx_files = list(root.rglob("matrix.mtx")) + list(root.rglob("matrix.mtx.gz"))
    for mtx in mtx_files:
        d = mtx.parent
        feat = None
        for fn in ["features.tsv.gz","features.tsv","genes.tsv.gz","genes.tsv"]:
            p = d / fn
            if p.exists():
                feat = p
                break
        bc = None
        for fn in ["barcodes.tsv.gz","barcodes.tsv"]:
            p = d / fn
            if p.exists():
                bc = p
                break
        if feat is None or bc is None:
            continue
        parts = [x for x in mtx.parts if x.startswith("GSM")]
        gsm = parts[0] if parts else d.name
        runs.append({"gsm": gsm, "dir": d, "mtx": mtx, "features": feat, "barcodes": bc})

    # --- Case B: GEO prefix files anywhere under root
    # Find any *_matrix.mtx or *_matrix.mtx.gz
    geo_mtx = list(root.rglob("*_matrix.mtx")) + list(root.rglob("*_matrix.mtx.gz"))
    for mtx in geo_mtx:
        stem = mtx.name
        # remove suffixes to get prefix key
        prefix = stem.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
        d = mtx.parent

        # features can be *_features.tsv(.gz) or *_genes.tsv(.gz)
        feat = None
        for fn in [f"{prefix}_features.tsv.gz", f"{prefix}_features.tsv",
                   f"{prefix}_genes.tsv.gz",    f"{prefix}_genes.tsv"]:
            p = d / fn
            if p.exists():
                feat = p
                break

        bc = None
        for fn in [f"{prefix}_barcodes.tsv.gz", f"{prefix}_barcodes.tsv"]:
            p = d / fn
            if p.exists():
                bc = p
                break

        if feat is None or bc is None:
            continue

        gsm = prefix.split("_")[0] if prefix.startswith("GSM") else prefix
        runs.append({"gsm": gsm, "dir": d, "mtx": mtx, "features": feat, "barcodes": bc, "prefix": prefix})

    df = pd.DataFrame(runs)
    if df.empty:
        return df

    # De-duplicate by the exact mtx path
    df = df.drop_duplicates(subset=["mtx"]).reset_index(drop=True)
    return df

runs = discover_10x_runs_geo_style(RAW_DIR)
print("Discovered runs:", runs.shape[0])
print(runs[["gsm","mtx","features","barcodes"]].head())

if runs.shape[0] == 0:
    raise RuntimeError("No 10x runs found in extracted GSE268307 RAW tar. Check RAW_DIR contents.")

# --------------------------------
# 4) Read all runs into AnnData + concat
# --------------------------------

def read_one_run(prefix: str, mtx: Path, features: Path, barcodes: Path) -> ad.AnnData:
    # features
    with gzip.open(features, "rt") as f:
        feats = pd.read_csv(f, sep="\t", header=None, dtype=str)
    if feats.shape[1] < 2:
        raise ValueError(f"Features file {features} has <2 cols.")
    gene_ids  = feats.iloc[:, 0].astype(str).tolist()
    gene_syms = feats.iloc[:, 1].astype(str).tolist()
    gene_syms = make_unique(gene_syms)

    # matrix (genes x cells)
    with gzip.open(mtx, "rb") as f:
        Xgxc = mmread(f).tocsr()

    # barcodes
    with gzip.open(barcodes, "rt") as f:
        bcs = pd.read_csv(f, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()

    if Xgxc.shape[0] != len(gene_syms):
        raise ValueError(f"{prefix}: matrix rows != features rows")
    if Xgxc.shape[1] != len(bcs):
        raise ValueError(f"{prefix}: matrix cols != barcodes")

    # AnnData is cells x genes
    X = Xgxc.T.tocsr()
    obs_names = [f"{prefix}_{bc}" for bc in bcs]

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame({"gene_id": gene_ids}, index=pd.Index(gene_syms, name="gene_symbol"))
    )
    adata.obs["run_prefix"] = prefix
    adata.obs["GSM"] = prefix.split("_")[0]
    return adata

adata_list = []
for _, r in runs.iterrows():
    adata_list.append(read_one_run(r["prefix"], r["mtx"], r["features"], r["barcodes"]))

adata = ad.concat(adata_list, join="outer", axis=0)
print(adata)

# --------------------------------
# 5) Clinical group annotation (GSE268307 has 4 GSMs)
# --------------------------------
# From GEO page:
# GSM8289374, GSM8289375: localized
# GSM8289376, GSM8289377: mHNPC (metastatic hormone-naïve)
localized = {"GSM8289374","GSM8289375"}
mHNPC     = {"GSM8289376","GSM8289377"}

def assign_group(gsm):
    if gsm in localized:
        return "Localized"
    if gsm in mHNPC:
        return "mHNPC"
    return np.nan

adata.obs["disease_group"] = adata.obs["GSM"].astype(str).map(assign_group)
adata.obs["disease_group"] = pd.Categorical(
    adata.obs["disease_group"],
    categories=["Localized","mHNPC"],
    ordered=True
)

adata.obs["disease_group"].value_counts()

# --------------------------------
# 6) QC + filtering
# --------------------------------
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# tune thresholds as needed
min_genes = 200
max_genes = 6000
max_mt = 20.0

sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :].copy()
adata = adata[adata.obs["pct_counts_mt"] <= max_mt, :].copy()

sc.pp.filter_genes(adata, min_cells=3)

print("After QC:", adata)

# store counts layer BEFORE normalization
adata.layers["counts"] = adata.X.copy()

# --------------------------------
# 7) Normalize/log + set raw + HVG/PCA/UMAP/Leiden
# --------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Keep *full* log1p matrix for scoring later
adata.raw = adata

# Use HVGs for embedding (but raw keeps all genes)
sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=False, flavor="seurat")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)

# Basic UMAP QC
plt.figure(figsize=(6,5))
sc.pl.umap(adata, color=["disease_group","leiden"], wspace=0.4, frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_disease_group_leiden.png", dpi=600, bbox_inches="tight")
plt.close()

# ---------------------------------------
# Save processed object for downstream
# ---------------------------------------

out_file = str(OUTS_DIR / "GSE268307_human_processed.h5ad")
os.makedirs(os.path.dirname(out_file), exist_ok=True)
adata.write_h5ad(out_file)
print(f"\nSaved processed AnnData to: {out_file}")


# ---------------------------------------
# load the processed adata file
# ---------------------------------------

adata = sc.read_h5ad(out_file)
adata


# --------------------------------
# 8)  Annotation
# --------------------------------

# =========================
# Marker sets
# =========================

# Tumor epithelium core
TUMOR_CORE = {"EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "TACSTD2"}

# Luminal tumor
TUMOR_LUMINAL = TUMOR_CORE | {
    "KLK3", "KLK2", "MSMB", "NKX3-1", "AR", "ACPP", "FOLH1", "TMPRSS2", "AMACR"
}

# Basal tumor
TUMOR_BASAL = TUMOR_CORE | {
    "KRT5", "KRT14", "KRT15", "TP63"
}

# Neuroendocrine-like tumor
TUMOR_NE = TUMOR_CORE | {
    "CHGA", "CHGB", "SYP", "ENO2", "NCAM1", "ASCL1", "INSM1"
}

# Cycling / proliferative epithelium
TUMOR_PROLIF = TUMOR_CORE | {
    "MKI67","TOP2A","BIRC5","UBE2C","CCNB1","CCNA2","HMGB2","PCNA","CKS2"
}

# Stromal & vascular
FB_FIBRO = {"COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","FBLN1","VIM","THY1","CXCL12"}
FB_MYOFIBRO = {"ACTA2","TAGLN","MYL9","FAP","POSTN","COL1A1","COL1A2","SPARC"}
SMC = {"MYH11","CNN1","TAGLN","ACTA2","DES"}
PERICYTE = {"RGS5","CSPG4","PDGFRB","KCNJ8","ABCC9","MCAM","DES","NOTCH3"}
ENDO = {"PECAM1","CDH5","KDR","FLT1","VWF","ESAM","CLDN5","KLF2","KLF4", "PROX1","PDPN","LYVE1","FLT4","CCL21"}

####
# Immune lineages
## Core pan-immune / myeloid scaffolding
MYELOID_CORE = {"LYZ","LST1","TYROBP","FCER1G","AIF1","CTSS","SPI1"}
MONO_CLASSICAL = {"S100A8","S100A9","FCN1","VCAN","CTSS","LGALS3","CCR2"}
MONO_NONCLASSICAL = {"FCGR3A","MS4A7","LST1","IFITM3","LGALS3","CTSS"}
# Macrophage/TAM general + subtypes
TAM_C1Q = {"C1QA","C1QB","C1QC","APOE","TREM2","TYROBP","LST1","CTSD","CTSB","LGALS3","SPP1"}  # broad TAM-ish
TAM_RESIDENT = {"C1QA","C1QB","C1QC","MRC1","MSR1","CSF1R","MARCO","SIGLEC1","FCER1G","LST1","APOE"}
# SPP1hi TAM signature-like: high SPP1 with lipid/TREM2/APOE axis; negative vs FCN1/S100A8/9
SPP1_HI_TAM_POS = {"SPP1","TREM2","APOE","LGALS3","CTSB","CTSD","LST1","TYROBP","FCER1G","IFI30","GPNMB"}
SPP1_HI_TAM_NEG = {"FCN1","S100A8","S100A9","VCAN","CCR2"}  # “mono-like” suppressors

# DC states (include LILRA4/IL3RA for pDC, and LAMP3/mreg DC if present)
DC_CDC1 = {"CLEC9A","XCR1","BATF3","IRF8","FCER1A"}  # FCER1A not cDC1-specific but ok
DC_CDC2 = {"CD1C","FCER1A","CLEC10A","CD1E"}
DC_PDC  = {"GZMB","IRF7","TCF4","IL3RA","SPIB","LILRA4"}
#MAST = {"KIT","TPSAB1","TPSB2","CPA3","HDC"}

# Neutrophils
NEU = {"S100A8","S100A9","CSF3R","FCGR3B","MPO","ELANE","LTF"}

# Other immune
T_CORE = {"CD3D","CD3E","CD3G","CD2","TRAC","LCK"}
T_CD4 = {"CD4","IL7R","CCR7","LEF1","TCF7","ICOS"}
T_CD8 = {"CD8A","CD8B","PRF1","NKG7","GZMB","GZMK"}
NK = {"KLRD1","FCGR3A","NKG7","GNLY","PRF1","GZMB","XCL1","XCL2"}
B_NAIVE = {"MS4A1","CD79A","CD79B","BANK1","CD74","HLA-DRA"}
B_PLASMA = {"MZB1","SDC1","XBP1","JCHAIN","PRDM1","IGKC","IGHG1"}
T_TREG = {"FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"}
T_GDT = {"TRDC","TRGC1","TRGC2","TRGV9","TRDV1","TRDV2"}
ERY = {"HBB","HBA1","HBA2","ALAS2"}


MARKER_SETS = {
    # Tumor / epithelium
    "Tumor luminal": TUMOR_LUMINAL,
    "Tumor basal": TUMOR_BASAL,
    "Tumor NE-like": TUMOR_NE,
    "Tumor proliferative": TUMOR_PROLIF,

    # Stroma / vascular
    "Fibroblasts": FB_FIBRO,
    "Myofibroblasts": FB_MYOFIBRO,
    "Smooth muscle": SMC,
    "Pericytes": PERICYTE,
    "Endothelial": ENDO,

    # Immune
    "CD4+ T": T_CORE | T_CD4,
    "CD8+ T": T_CORE | T_CD8,
    "Treg": T_CORE | T_TREG,
    "γδ T": T_CORE | T_GDT,
    "NK": NK,
    "B": B_NAIVE,
    "Plasma": B_PLASMA,
    # Other immune
    #"Mast": MAST,
    "Erythroid": ERY,
    # Immune (myeloid)
    "Neutrophils": NEU,
    "FCN1/S100A8/A9+ Monocytes": MONO_CLASSICAL | MYELOID_CORE,
    "FCGR3A/MS4A7+ Monocytes": MONO_NONCLASSICAL | MYELOID_CORE,
    "APOE/TREM2+ TAMs": TAM_C1Q | MYELOID_CORE,
    "TAMs (resident-like)": TAM_RESIDENT | MYELOID_CORE,
    "cDC1": DC_CDC1,
    "cDC2": DC_CDC2,
    "pDC": DC_PDC,
}

LINEAGE = {
    "Tumor": {
        "Tumor luminal","Tumor basal","Tumor NE-like","Tumor proliferative",
    },
    "Stroma": {
        "Fibroblasts","Myofibroblasts","Pericytes","Endothelial", "Smooth muscle",
    },
    "Immune": {
        "CD4+ T","CD8+ T","Treg","γδ T","NK","B","Plasma",
        #"Mast",
        "Erythroid",
        "Neutrophils","FCN1/S100A8/A9+ Monocytes","FCGR3A/MS4A7+ Monocytes",
        "APOE/TREM2+ TAMs","TAMs (resident-like)",
        "cDC1","cDC2","pDC"
    },
}

def which_lineage(label: str) -> str:
    for lin, labels in LINEAGE.items():
        if label in labels:
            return lin
    return "Other"

def to_compartment(lbl: str) -> str:
    if lbl in LINEAGE["Tumor"]:
        return "Tumor"
    if lbl in LINEAGE["Stroma"]:
        return "Stroma"
    if lbl in LINEAGE["Immune"]:
        return "Immune"
    return "Other"



GROUP_KEY = "leiden"

# Compute DEGs per cluster (Wilcoxon)
sc.tl.rank_genes_groups(
    adata,
    groupby=GROUP_KEY,
    method="t-test_overestim_var",
    use_raw=True
)

# Convert to tidy DataFrame
mt = sc.get.rank_genes_groups_df(adata, group=None)
mt.rename(columns={"names": "SYMBOL"}, inplace=True)

# Add per-group rank (1 = top marker by score/pval, etc.)
mt = (
    mt.sort_values(["group", "scores"], ascending=[True, False])
      .groupby("group", as_index=False)
      .apply(lambda df: df.assign(rank=np.arange(1, len(df) + 1)))
      .reset_index(drop=True)
)

mt.head()


MIN_ACCEPT = 0.1      # minimum score to accept any label (tune)
SECOND_DELTA = 0.25   # how close second label can be (relative) to be marked as "second_if_close"

def weighted_score(dfc: pd.DataFrame, genes):
    dfc = dfc.copy()
    gs  = set(genes)
    hits = dfc[dfc["SYMBOL"].isin(gs)]
    if hits.empty:
        return 0.0, []
    # NaN logFC → 0
    lfc = dfc["logfoldchanges"].astype(float).fillna(0).clip(lower=0)
    hits = hits.assign(logfoldchanges=lfc.loc[hits.index])
    w = (1.0 / hits["rank"].astype(float)) * np.log1p(np.exp(hits["logfoldchanges"]))
    hits = hits.assign(w=w)
    keep = hits.sort_values("w", ascending=False)[["SYMBOL", "rank", "logfoldchanges", "w"]]
    return float(w.sum()), [tuple(x) for x in keep.values]

rows = []

for clust, dfc in mt.groupby("group", sort=True):
    dfc = dfc.sort_values("rank").reset_index(drop=True)

    label_scores = {}
    label_hits   = {}

    for label, genes in MARKER_SETS.items():
        sc_val, hl = weighted_score(dfc, genes)
        label_scores[label] = sc_val
        label_hits[label]   = hl

    # lineage scores
    per_lin = {"Tumor": 0.0, "Stroma": 0.0, "Immune": 0.0, "Other": 0.0}
    for lbl, sc_val in label_scores.items():
        per_lin[which_lineage(lbl)] += sc_val

    # cross-lineage penalty
    for lbl in list(label_scores.keys()):
        lin = which_lineage(lbl)
        other_lin = max(
            [k for k in per_lin.keys() if k != lin],
            key=lambda k: per_lin[k]
        )
        label_scores[lbl] = max(0.0, label_scores[lbl] - 0.1 * per_lin[other_lin])

    best = sorted(label_scores.items(), key=lambda kv: kv[1], reverse=True)
    top_label,    top_score    = best[0]
    second_label, second_score = best[1] if len(best) > 1 else (None, 0.0)

    # **Always accept top label unless no markers at all**
    final_label = top_label if top_score > 0 else "Unknown"
    second = (
        second_label
        if (top_score > 0 and (top_score - second_score)/max(top_score, 1e-9) <= SECOND_DELTA)
        else None
    )

    rows.append({
        "cluster": clust,
        "top_label": top_label,
        "top_score": round(top_score, 4),
        "second_label": second_label,
        "second_score": round(second_score, 4),
        "final_label": final_label,
        "second_if_close": second,
        "support_markers": ", ".join(
            [
                f"{g}(r{int(r)},lfc{(0 if pd.isna(lfc) else lfc):.1f})"
                for g, r, lfc, _w in label_hits.get(top_label, [])[:6]
            ]
        ),
    })

annot = pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)
annot

annot.to_csv(OUTS_DIR / "GSE268307_cluster_annotations_PCa.csv", index=False)
print(annot[["cluster","final_label","top_score","second_if_close","second_score"]])

GROUP_KEY = "leiden"

adata_annot = adata.copy()
adata_annot.obs[GROUP_KEY] = adata_annot.obs[GROUP_KEY].astype(str)

cluster2label  = dict(zip(annot["cluster"].astype(str), annot["final_label"]))
cluster2second = dict(zip(annot["cluster"].astype(str), annot["second_if_close"].fillna("").astype(str)))
cluster2comp   = {k: to_compartment(v) for k, v in cluster2label.items()}

adata.obs["celltype"]     = adata.obs[GROUP_KEY].map(cluster2label).astype("category")
adata.obs["celltype_alt"] = (
    adata.obs[GROUP_KEY].map(cluster2second)
    .replace({"": np.nan})
    .astype("category")
)
adata.obs["compartment"]  = adata.obs[GROUP_KEY].map(cluster2comp).astype("category")

print(adata.obs["celltype"].value_counts())
print(adata.obs["compartment"].value_counts())


# ---------------------------------------
# Save annot object for downstream
# ---------------------------------------

out_file_annot = str(OUTS_DIR / "GSE268307_human_annot.h5ad")
os.makedirs(os.path.dirname(out_file_annot), exist_ok=True)
adata.write_h5ad(out_file_annot)
print(f"\nSaved Annotated AnnData to: {out_file_annot}")


# ---------------------------------------
# load the processed adata file
# ---------------------------------------

adata = sc.read_h5ad(out_file_annot)
adata


#####################################
# Figures
####################################
# Make sure scanpy uses that dir too (for sc.pl.* saving)
sc.settings.figdir = str(FIG_DIR)
sc.set_figure_params(dpi=150, dpi_save=300, format="png", transparent=False)

print(adata)
print("Obs columns:", adata.obs.columns.tolist())
print("Var shape:", adata.var.shape)


# Disease groups
plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["disease_group"],
    size=5,
    legend_fontsize=8,
    legend_loc="upper right",
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_disease_group.png", dpi=400, bbox_inches="tight")
plt.close()


# -------------------------------------------------------------------
# Panel SxA: UMAP colored by celltype (all cells)
# -------------------------------------------------------------------
plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color="celltype",
    size=5,
    legend_fontsize=6,
    legend_loc="on data",
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "SxA_umap_celltype2.png", dpi=400, bbox_inches="tight")
plt.close()

# -------------------------------------------------------------------
# Panel SxA: UMAP colored by celltype and leiden
# -------------------------------------------------------------------
plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["leiden", "celltype"],
    size=5,
    legend_fontsize=8,
    legend_loc="on data",  # remove if too busy and re-run with 'right margin'
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "SxA_umap_leiden_celltype.png", dpi=400, bbox_inches="tight")
plt.close()



# --------------------------------
# Met-Score Enrichment
# met_score below is the directional 45-gene module score computed per cell as
# met_score_pos - met_score_neg (scanpy score_genes over the positive and negative
# signature gene lists). It is not the frozen 41-feature ridge-logistic Met-Score
# classifier: the bulk logistic probability is not defined for a single cell, so
# the directional module score is used as the single-cell surrogate.
# --------------------------------
pos_path = OUTS_DIR / "MetScore_pos_genes.txt"
neg_path = OUTS_DIR / "MetScore_neg_genes.txt"

PositiveGenes = pd.read_csv(pos_path, header=None, sep="\t").iloc[:, 0].astype(str).str.upper().tolist()
NegativeGenes = pd.read_csv(neg_path, header=None, sep="\t").iloc[:, 0].astype(str).str.upper().tolist()

MetScore = PositiveGenes + NegativeGenes

raw_syms_upper = pd.Index([g.upper() for g in adata.raw.var_names])

pos_in = [adata.raw.var_names[i] for i in np.where(raw_syms_upper.isin(set(PositiveGenes)))[0]]
neg_in = [adata.raw.var_names[i] for i in np.where(raw_syms_upper.isin(set(NegativeGenes)))[0]]

print("MetScore overlap:", len(pos_in), "pos,", len(neg_in), "neg")
if len(pos_in) < 5 or len(neg_in) < 5:
    print("[WARN] Low overlap with MetScore genes; interpretation will be noisier.")

sc.tl.score_genes(adata, pos_in, score_name="met_score_pos", use_raw=True)
sc.tl.score_genes(adata, neg_in, score_name="met_score_neg", use_raw=True)
adata.obs["met_score"] = adata.obs["met_score_pos"] - adata.obs["met_score_neg"]

print(adata.obs[["met_score_pos", "met_score_neg", "met_score"]].head())


# -------------------------------------------------------------------
# UMAPs of met_score + disease_group
# -------------------------------------------------------------------
plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["met_score"],
    color_map="magma",
    colorbar_loc = "bottom",
    size=5,
    wspace=0.4,
    #vmin =0,
    #vmax =3,
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_met_score.png", dpi=400, bbox_inches="tight")
plt.close()


plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["met_score_pos"],
    color_map="magma",
    colorbar_loc = "bottom",
    size=5,
    wspace=0.4,
    #vmin =0,
    #vmax =3,
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_met_score_pos.png", dpi=400, bbox_inches="tight")
plt.close()


plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["met_score_neg"],
    color_map="magma",
    colorbar_loc = "bottom",
    size=5,
    wspace=0.4,
    #vmin =0,
    #vmax =3,
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_met_score_neg.png", dpi=400, bbox_inches="tight")
plt.close()


# -----------------------------
# Met-score genes fraction per cell type
# -----------------------------
var_syms = adata.raw.var_names.str.upper()
met_genes = sorted(set(MetScore) & set(var_syms))

rows = []
for gene in MetScore:
    x = adata.raw[:, [gene]].X
    x = x.toarray().ravel() if hasattr(x, "toarray") else np.array(x).ravel()

    df = pd.DataFrame({"gene": gene, "celltype": adata.obs["celltype"].astype(str).values, "expr": x})
    out = df.groupby(["gene","celltype"]).apply(
        lambda d: pd.Series({
            "pct_expr": 100.0 * (d["expr"] > 0).mean(),
            "mean_expr": d["expr"].mean(),
            "median_expr": d["expr"].median(),
            "n_cells": d.shape[0]
        })
    ).reset_index()
    rows.append(out)

out_all = pd.concat(rows, ignore_index=True)
out_all.to_csv(OUTS_DIR / "GSE268307_met_score_gene_celltype_expression_WITH_GENE.csv", index=False)


# --------------------------------
# tumor-compartment MetScore by clinical group
# --------------------------------
tumor = adata[adata.obs["compartment"].astype(str) == "Tumor", :].copy()
print("Tumor-like cells:", tumor.n_obs)

# Per-sample median within tumor-like cells (GSM is effectively patient here; n=4)
tumor_df = tumor.obs[["GSM","disease_group","met_score"]].dropna().copy()
# Make sure met_score is numeric and group keys are plain strings (avoids categorical edge-cases)
tumor_df = tumor_df.copy()
tumor_df["GSM"] = tumor_df["GSM"].astype(str)
tumor_df["disease_group"] = tumor_df["disease_group"].astype(str)
tumor_df["met_score"] = pd.to_numeric(tumor_df["met_score"], errors="coerce")

tumor_df = tumor_df.dropna(subset=["GSM", "disease_group", "met_score"])

per_gsm = (
    tumor_df
    .groupby(["GSM", "disease_group"])   # as_index=True default
    .agg(
        n_cells=("met_score", "size"),
        met_median=("met_score", "median"),
        met_mean=("met_score", "mean"),
    )
    .reset_index()
)

print(per_gsm)

from scipy import stats

x = per_gsm.loc[per_gsm["disease_group"] == "Localized", "met_median"].values
y = per_gsm.loc[per_gsm["disease_group"] == "mHNPC", "met_median"].values

U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
print("MWU on per-GSM median Met-Score:", "U=", U, "p=", p)

# Plot
plt.figure(figsize=(4.2,4.0))
sns.boxplot(data=per_gsm, x="disease_group", y="met_median", color="white",
            boxprops={"edgecolor":"black"}, medianprops={"color":"black"})
sns.stripplot(data=per_gsm, x="disease_group", y="met_median", color="black", size=6, jitter=0.08)
plt.xlabel("")
plt.ylabel("Tumor median 45-gene signature module score per sample")
plt.title("GSE268307 (scRNA): tumor 45-gene signature module score")
plt.tight_layout()
plt.savefig(FIG_DIR / "tumor_met_score_by_group_GSE268307.png", dpi=600, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------
# Per-GSM cell-type composition (stacked bars) + tumor subtype composition
# ----------------------------------------------------------------
# ----------------------------
# Helper: stacked bar plot
# ----------------------------
def stacked_bar(df_wide, title, ylabel, out_png, sort_index=True):
    if sort_index:
        df_wide = df_wide.sort_index()
    ax = df_wide.plot(kind="bar", stacked=True, figsize=(7.2, 4.2), width=0.85)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.legend(
        title="",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close()

# ----------------------------
# 2A) Broad compartment composition per GSM
# ----------------------------
df = adata.obs[["GSM", "disease_group", "celltype", "compartment"]].copy()
df["GSM"] = df["GSM"].astype(str)
df["disease_group"] = df["disease_group"].astype(str)
df["celltype"] = df["celltype"].astype(str)
df["compartment"] = df["compartment"].astype(str)


comp_counts = (
    df.groupby(["GSM", "disease_group", "compartment"])
      .size()
      .rename("n_cells")
      .reset_index()
)

comp_frac = comp_counts.copy()
comp_frac["frac"] = comp_frac.groupby(["GSM"])["n_cells"].transform(lambda x: x / x.sum())

# Save
comp_counts.to_csv(OUTS_DIR / "GSE268307_compartment_counts_per_GSM.csv", index=False)
comp_frac.to_csv(OUTS_DIR / "GSE268307_compartment_fraction_per_GSM.csv", index=False)

# Wide for plotting (rows=GSM, cols=compartment)
comp_wide = comp_frac.pivot_table(index="GSM", columns="compartment", values="frac", fill_value=0.0)

# Optional: order GSMs by disease group then GSM
gsm_order = (
    df[["GSM", "disease_group"]].drop_duplicates()
      .sort_values(["disease_group", "GSM"])
)["GSM"].tolist()
comp_wide = comp_wide.loc[gsm_order]

stacked_bar(
    comp_wide,
    title="GSE268307: Cell compartment composition per sample",
    ylabel="Fraction of cells",
    out_png=FIG_DIR / "GSE268307_compartment_composition_stacked.png",
)

# ----------------------------
# 2B) Celltype composition per GSM (all celltypes)
# ----------------------------
ct_counts = (
    df.groupby(["GSM", "disease_group", "celltype"])
      .size()
      .rename("n_cells")
      .reset_index()
)

ct_frac = ct_counts.copy()
ct_frac["frac"] = ct_frac.groupby(["GSM"])["n_cells"].transform(lambda x: x / x.sum())

ct_counts.to_csv(OUTS_DIR / "GSE268307_celltype_counts_per_GSM.csv", index=False)
ct_frac.to_csv(OUTS_DIR / "GSE268307_celltype_fraction_per_GSM.csv", index=False)

ct_wide = ct_frac.pivot_table(index="GSM", columns="celltype", values="frac", fill_value=0.0)
ct_wide = ct_wide.loc[gsm_order]

stacked_bar(
    ct_wide,
    title="GSE268307: Cell type composition per sample",
    ylabel="Fraction of cells",
    out_png=FIG_DIR / "GSE268307_celltype_composition_stacked.png",
)

# ----------------------------
# 2C) Tumor subtype composition per GSM (Tumor compartment only)
# ----------------------------
tumor_mask = df["compartment"].eq("Tumor")
df_tum = df.loc[tumor_mask, ["GSM", "disease_group", "celltype"]].copy()

if df_tum.shape[0] == 0:
    print("[WARN] No Tumor compartment cells found; skipping tumor subtype composition.")
else:
    # Define tumor subtype labels you care about (adjust if your annotation differs)
    tumor_subtypes = ["Tumor luminal", "Tumor basal", "Tumor proliferative", "Tumor NE-like"]

    df_tum["tumor_subtype"] = np.where(df_tum["celltype"].isin(tumor_subtypes), df_tum["celltype"], "Other tumor/unknown")

    tum_counts = (
        df_tum.groupby(["GSM", "disease_group", "tumor_subtype"])
              .size()
              .rename("n_cells")
              .reset_index()
    )
    tum_frac = tum_counts.copy()
    tum_frac["frac"] = tum_frac.groupby(["GSM"])["n_cells"].transform(lambda x: x / x.sum())

    tum_counts.to_csv(OUTS_DIR / "GSE268307_tumor_subtype_counts_per_GSM.csv", index=False)
    tum_frac.to_csv(OUTS_DIR / "GSE268307_tumor_subtype_fraction_per_GSM.csv", index=False)

    tum_wide = tum_frac.pivot_table(index="GSM", columns="tumor_subtype", values="frac", fill_value=0.0)
    tum_wide = tum_wide.loc[[g for g in gsm_order if g in tum_wide.index]]

    stacked_bar(
        tum_wide,
        title="GSE268307: Tumor subtype composition per sample (Tumor cells only)",
        ylabel="Fraction of tumor cells",
        out_png=FIG_DIR / "GSE268307_tumor_subtype_composition_stacked.png",
    )

print("Done: composition outputs saved to OUTS_DIR and FIG_DIR.")

# ----------------------------------------------------------------
# Tumor-only Met-Score by tumor subtype (per GSM; box/points)
# ----------------------------------------------------------------
# Pull needed obs
obs = adata.obs[["GSM", "disease_group", "compartment", "celltype", "met_score"]].copy()
obs["GSM"] = obs["GSM"].astype(str)
obs["disease_group"] = obs["disease_group"].astype(str)
obs["compartment"] = obs["compartment"].astype(str)
obs["celltype"] = obs["celltype"].astype(str)
obs["met_score"] = pd.to_numeric(obs["met_score"], errors="coerce")
obs = obs.dropna(subset=["met_score", "GSM", "disease_group", "compartment", "celltype"])

# Tumor cells only
tum = obs.loc[obs["compartment"].eq("Tumor")].copy()
if tum.shape[0] == 0:
    raise ValueError("No Tumor compartment cells found; cannot run tumor-only subtype Met-Score analysis.")

# Tumor subtype definition (match your annotation)
tumor_subtypes = ["Tumor luminal", "Tumor basal", "Tumor proliferative", "Tumor NE-like"]
tum = tum[tum["celltype"].isin(tumor_subtypes)].copy()

if tum.shape[0] == 0:
    raise ValueError(f"No tumor cells matched {tumor_subtypes}. Check your 'celltype' labels in adata.obs.")

# Per GSM x subtype summaries
per_gsm_sub = (
    tum.groupby(["GSM", "disease_group", "celltype"], as_index=False)
       .agg(
           n_cells=("met_score", "size"),
           met_median=("met_score", "median"),
           met_mean=("met_score", "mean")
       )
)

per_gsm_sub.to_csv(OUTS_DIR / "GSE268307_tumor_MetScore_per_GSM_by_subtype.csv", index=False)
print(per_gsm_sub)

# Per-GSM dot plot
# Order GSMs by disease_group then GSM
gsm_order = (
    obs[["GSM", "disease_group"]].drop_duplicates()
       .sort_values(["disease_group", "GSM"])
)["GSM"].tolist()

sub_order = [s for s in tumor_subtypes if s in per_gsm_sub["celltype"].unique()]

plt.figure(figsize=(7.2, 4.0))
ax = sns.stripplot(
    data=per_gsm_sub,
    x="celltype",
    y="met_median",
    hue="disease_group",
    order=sub_order,
    dodge=True,
    jitter=0.12,
    size=6
)
ax.set_xlabel("")
ax.set_ylabel("Per-sample median 45-gene signature module score (tumor cells)")
ax.set_title("GSE268307: Tumor 45-gene signature module score by tumor subtype (per sample)")
ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(FIG_DIR / "GSE268307_tumor_met_score_by_subtype_per_GSM.png", dpi=600, bbox_inches="tight")
plt.close()

# Box + points across patients
plt.figure(figsize=(7.2, 4.0))
ax = sns.boxplot(
    data=per_gsm_sub,
    x="celltype",
    y="met_median",
    order=sub_order,
    color="white",
    boxprops={"edgecolor": "black"},
    medianprops={"color": "black"},
)
sns.stripplot(
    data=per_gsm_sub,
    x="celltype",
    y="met_median",
    order=sub_order,
    color="black",
    size=6,
    jitter=0.10,
    alpha=0.8,
)
ax.set_xlabel("")
ax.set_ylabel("Per-sample median 45-gene signature module score (tumor cells)")
ax.set_title("GSE268307: Tumor 45-gene signature module score by subtype (n=4 samples)")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(FIG_DIR / "GSE268307_tumor_met_score_by_subtype_box.png", dpi=600, bbox_inches="tight")
plt.close()

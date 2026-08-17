import os
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import io, sparse
from datetime import datetime
import scipy.sparse as sp
import plotly.express as px
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import glob
from scipy.io import mmread
from scipy import sparse
import anndata as ad
from scipy import stats
from scipy import stats

# ----------------------------
# 0. Paths and basic settings
# ----------------------------

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

base_dir = str(ROOT / "data" / "GSE274229")
raw_dir = os.path.join(base_dir, "RAW")
os.makedirs(base_dir, exist_ok=True)

print("RAW dir:", raw_dir)

# Make Scanpy a bit quieter if you like
sc.settings.verbosity = 2  # 0 = errors only, 1 = warnings, 2 = info
sc.settings.figdir = "./figures/GSE274229"
os.makedirs(sc.settings.figdir, exist_ok=True)

# --------------------------------------------------
# 1. Classify libraries (human vs mouse) by features
# --------------------------------------------------

feat_files = glob.glob(os.path.join(raw_dir, "*_GEX_features.tsv.gz"))
print(f"Found {len(feat_files)} feature files")

def classify_species(feat_file):
    """
    Read the GEX_features.tsv.gz and decide if it's human, mouse, or mixed
    based on gene IDs (ENSG vs ENSMUSG).
    """
    # Features file assumed: <library_id>_GEX_features.tsv.gz
    lib_id = os.path.basename(feat_file).replace("_GEX_features.tsv.gz", "")
    print("Classifying species for:", lib_id)

    with gzip.open(feat_file, "rt") as f:
        # first two columns: gene_id, gene_symbol
        feats = pd.read_csv(
            f,
            sep="\t",
            header=None,
            usecols=[0],
            dtype=str
        )
    gene_ids = feats.iloc[:, 0].fillna("")

    prop_human = np.mean(gene_ids.str.startswith("ENSG"))
    prop_mouse = np.mean(gene_ids.str.startswith("ENSMUSG"))

    if prop_human > 0.5:
        species = "human"
    elif prop_mouse > 0.5:
        species = "mouse"
    else:
        species = "mixed/unknown"

    return {
        "features_file": feat_file,
        "library_id": lib_id,
        "prop_human": prop_human,
        "prop_mouse": prop_mouse,
        "species": species,
    }

species_list = [classify_species(f) for f in feat_files]
species_df = pd.DataFrame(species_list)
print("\nSpecies summary:")
print(species_df["species"].value_counts())
print(species_df.head())

# Keep human libraries only
human_libs = species_df.query("species == 'human'").copy()
print(f"\nUsing {len(human_libs)} human libraries:")
print(human_libs[["library_id", "prop_human", "prop_mouse"]])


# -------------------------------------------------
# 2. Helper to read one human library into AnnData
# -------------------------------------------------

def make_unique(names):
    """
    Make a list of strings unique by adding .1, .2, etc. as needed.
    """
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

def read_one_library(lib_id):
    """
    Read one 10X-like library: matrix, features, barcodes -> AnnData.
    Assumed filenames in RAW:
        <lib_id>_GEX_matrix.mtx.gz
        <lib_id>_GEX_features.tsv.gz
        <lib_id>_GEX_barcodes.tsv.gz
    """
    print(f"\nReading library: {lib_id}")

    mtx_file = os.path.join(raw_dir, f"{lib_id}_GEX_matrix.mtx.gz")
    feat_file = os.path.join(raw_dir, f"{lib_id}_GEX_features.tsv.gz")
    bc_file   = os.path.join(raw_dir, f"{lib_id}_GEX_barcodes.tsv.gz")

    # ---- features: gene_id, gene_symbol, ...
    with gzip.open(feat_file, "rt") as f:
        feats = pd.read_csv(
            f,
            sep="\t",
            header=None,
            dtype=str
        )
    if feats.shape[1] < 2:
        raise ValueError(f"Features file {feat_file} has <2 columns, check file")

    gene_ids  = feats.iloc[:, 0].astype(str)
    gene_syms = feats.iloc[:, 1].astype(str)
    gene_syms = pd.Series(make_unique(list(gene_syms)))

    # ---- matrix (genes x cells)
    with gzip.open(mtx_file, "rb") as f:
        mtx = mmread(f).tocsr()

    # ---- barcodes
    with gzip.open(bc_file, "rt") as f:
        barcodes = pd.read_csv(
            f,
            sep="\t",
            header=None,
            dtype=str
        ).iloc[:, 0].tolist()

    if mtx.shape[0] != len(gene_syms):
        raise ValueError(
            f"{lib_id}: matrix rows {mtx.shape[0]} != number of features {len(gene_syms)}"
        )
    if mtx.shape[1] != len(barcodes):
        raise ValueError(
            f"{lib_id}: matrix cols {mtx.shape[1]} != number of barcodes {len(barcodes)}"
        )

    # ---- Build AnnData: cells x genes, so transpose
    X = mtx.T.tocsr()  # cells x genes

    # Make cell IDs globally unique by prefixing library_id
    cell_ids = [f"{lib_id}_{bc}" for bc in barcodes]

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_ids),
        var=pd.DataFrame(
            {"gene_id": list(gene_ids)},
            index=gene_syms.values
        )
    )
    adata.obs["library_id"] = lib_id

    return adata

# ----------------------------------------------
# 3. Read all human libraries and concatenate
# ----------------------------------------------

adata_list = []
for _, row in human_libs.iterrows():
    lib_id = row["library_id"]
    adata_lib = read_one_library(lib_id)
    adata_list.append(adata_lib)

print(f"\nConcatenating {len(adata_list)} human libraries...")

adata = ad.concat(
    adata_list,
    join="outer",
    label="library_id",
    keys=None,
    index_unique=None
)

print(adata)
print("Number of cells:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# ---------------------------------------------------
# 4. Annotate samples: GSM, disease_group, index
# ---------------------------------------------------
idx = adata.obs_names.to_series()

# Split index: [0]=GSM8445680, [1]=S1, [2]=03, [3]=barcode
parts = idx.str.split("_", expand=True)

meta_cell = pd.DataFrame(index=adata.obs_names)
meta_cell["GSM"]      = parts[0]
meta_cell["lib_id"]   = parts[0] + "_" + parts[1] + "_" + parts[2]  # GSM8445680_S1_03

# Attach back to adata.obs
adata.obs["GSM"]        = meta_cell["GSM"].values
adata.obs["library_id"] = meta_cell["lib_id"].values

print(adata.obs[["GSM", "library_id"]].head())
# Unique GSMs
gsm_map = (
    adata.obs[["GSM"]]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Numeric GSM
gsm_map["gsm_num"] = gsm_map["GSM"].str.replace("GSM", "", regex=False).astype(int)

def assign_group(num):
    if 8445680 <= num <= 8445685:
        return "mCRPC"
    elif 8445686 <= num <= 8445710:
        return "mHSPC"
    elif 8445711 <= num <= 8445723:
        return "Localized"
    else:
        return np.nan

gsm_map["disease_group"] = gsm_map["gsm_num"].apply(assign_group)

# 1–44 index like your R code (for all 8445680–8445723)
gsm_map["sample_index"] = np.where(
    (gsm_map["gsm_num"] >= 8445680) & (gsm_map["gsm_num"] <= 8445723),
    gsm_map["gsm_num"] - 8445679,
    np.nan
)

print(gsm_map.head())

gsm_map = gsm_map.set_index("GSM")

adata.obs = adata.obs.join(
    gsm_map[["disease_group", "sample_index"]],
    on="GSM"
)

adata.obs["disease_group"] = pd.Categorical(
    adata.obs["disease_group"],
    categories=["Localized", "mHSPC", "mCRPC"],
    ordered=True,
)

print(
    adata.obs[["GSM", "library_id", "disease_group", "sample_index"]]
    .head()
)

print("\nObs columns now:")
print(adata.obs.columns)
print("\nDisease group distribution in cells:")
print(adata.obs["disease_group"].value_counts(dropna=False))

# -----------------------------------------------
# 5. QC metrics: n_genes, total_counts, pct MT
# -----------------------------------------------

# Mark mitochondrial genes (human: "MT-")
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["mt"].value_counts()

# Calculate QC metrics per cell
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

print("\nQC summary (first few cells):")
print(
    adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
    .head()
)

# ------------------------------------------
# 6. Basic filtering (you can tune thresholds)
# ------------------------------------------


min_genes = 200
max_genes = 6000
max_mt = 20.0

print("\nFiltering cells with:")
print(f"  n_genes_by_counts >= {min_genes}")
print(f"  n_genes_by_counts <= {max_genes}")
print(f"  pct_counts_mt <= {max_mt}")

initial_n = adata.n_obs

sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :]
adata = adata[adata.obs["pct_counts_mt"] <= max_mt, :]

print(f"Cells before filtering: {initial_n}")
print(f"Cells after filtering:  {adata.n_obs}")

# Also drop genes not expressed in any remaining cells
sc.pp.filter_genes(adata, min_cells=3)
print("Genes after filtering:", adata.n_vars)

# -----------------------------------------------------
# 7. Normalization, log-transform, HVGs, PCA, neighbors
# -----------------------------------------------------

# store raw counts
adata.layers["counts"] = adata.X

# Normalize total counts per cell to 1e4 and log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ✅ Store the full log-normalized matrix as raw BEFORE subsetting
adata.raw = adata

# Highly variable genes (subset to HVGs for .X, but .raw keeps ALL genes)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    subset=True,
    flavor="seurat"
)
print("HVGs selected:", adata.n_vars)

# Scale
sc.pp.scale(adata, max_value=10)



# PCA
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

# Neighborhood graph
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)

# UMAP
sc.tl.umap(adata)

# Leiden clustering
sc.tl.leiden(adata, resolution=1)

print("\nFinished basic processing.")
print("Clusters (leiden_0_6):")
print(adata.obs["leiden"].value_counts())

# ---------------------------------------
# 8. Save processed object for downstream
# ---------------------------------------

out_file = str(ROOT / "outs" / "GSE274229_human_processed.h5ad")
os.makedirs(os.path.dirname(out_file), exist_ok=True)
adata.write_h5ad(out_file)
print(f"\nSaved processed AnnData to: {out_file}")

# ---------------------------------------
# load the processed adata file
# ---------------------------------------

adata = sc.read_h5ad(out_file)
adata


# =========================
# Marker sets — PROSTATE
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

annot.to_csv(ROOT / "outs" / "GSE274229_cluster_annotations_PCa.csv", index=False)
print(annot[["cluster","final_label","top_score","second_if_close","second_score"]])

GROUP_KEY = "leiden"   # or whatever clustering column you used

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



##############################
# update annotation
##############################
FIG_DIR = ROOT / "figures" / "scRNAseq"
FIG_DIR.mkdir(parents=True, exist_ok=True)

plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["leiden", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LST1"],
    size=5,
    legend_fontsize=8,
    legend_loc="on data",  # remove if too busy and re-run with 'right margin'
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "SxA_umap_TAMs_markers.png", dpi=400, bbox_inches="tight")
plt.close()

plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["leiden", "MS4A1", "CD79A"],
    size=5,
    legend_fontsize=8,
    legend_loc="on data",  # remove if too busy and re-run with 'right margin'
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "SxA_umap_Bcell_markers.png", dpi=400, bbox_inches="tight")
plt.close()

plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["leiden", "MKI67","TOP2A","BIRC5","UBE2C","CCNB1","CCNA2","HMGB2","PCNA","CKS2"],
    size=5,
    legend_fontsize=8,
    legend_loc="on data",  # remove if too busy and re-run with 'right margin'
    frameon=False,
    show=False
)
plt.title("", fontsize=14)
plt.tight_layout()
plt.savefig(FIG_DIR / "SxA_umap_Prolif_markers.png", dpi=400, bbox_inches="tight")
plt.close()


# -----------------------------------------
# SECOND LAYER: call SPP1hi TAM within macrophage/TAM only
# -----------------------------------------
# score using raw genes (safer)
raw_genes = set(adata.raw.var_names)

pos_in = [g for g in SPP1_HI_TAM_POS if g in raw_genes]
neg_in = [g for g in SPP1_HI_TAM_NEG if g in raw_genes]

sc.tl.score_genes(adata, pos_in, score_name="SPP1hi_pos", use_raw=True)
sc.tl.score_genes(adata, neg_in, score_name="SPP1hi_neg", use_raw=True)
adata.obs["SPP1hi_score"] = adata.obs["SPP1hi_pos"] - adata.obs["SPP1hi_neg"]

tam_mask = adata.obs["celltype"].isin([
    "APOE/TREM2+ TAMs",
    "TAMs (resident-like)"
])

# --- add the new category BEFORE assignment
# new_label = "SPP1hi TAM"
# if new_label not in adata.obs["celltype"].cat.categories:
#     adata.obs["celltype"] = adata.obs["celltype"].cat.add_categories([new_label])

# robust threshold: top quartile among TAMs (adjust to 0.70–0.85 if needed)
# if tam_mask.sum() > 50:
#     thr = float(np.quantile(adata.obs.loc[tam_mask, "SPP1hi_score"].astype(float), 0.75))
#
#     hit = tam_mask & (adata.obs["SPP1hi_score"].astype(float) >= thr)
#     adata.obs.loc[hit, "celltype"] = new_label

print("\nCelltype counts:")
print(adata.obs["celltype"].value_counts())

print("\nCompartment counts:")
print(adata.obs["compartment"].value_counts())

########################################################################
# save
########################################################################
OUTS_DIR = ROOT / "outs"

out_file_annot = str(ROOT / "outs" / "GSE274229_human_annot.h5ad")
os.makedirs(os.path.dirname(out_file_annot), exist_ok=True)
adata.write_h5ad(out_file_annot)

# -------------------------------------------------------------------
# Figure directory
# -------------------------------------------------------------------
FIG_DIR = ROOT / "figures" / "scRNAseq"
FIG_DIR.mkdir(parents=True, exist_ok=True)

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

#################################################################################
# load met-score
# met_score below is the directional 45-gene module score computed per cell as
# met_score_pos - met_score_neg (scanpy score_genes over the positive and negative
# signature gene lists). It is not the frozen 41-feature ridge-logistic Met-Score
# classifier: the bulk logistic probability is not defined for a single cell, so
# the directional module score is used as the single-cell surrogate.
#################################################################################
pos_path = str(ROOT / "outs" / "MetScore_pos_genes.txt")
neg_path = str(ROOT / "outs" / "MetScore_neg_genes.txt")

PositiveGenes = (
    pd.read_csv(pos_path, header=None, sep="\t")
    .iloc[:, 0]
    .astype(str)
    .str.upper()
    .tolist()
)

NegativeGenes = (
    pd.read_csv(neg_path, header=None, sep="\t")
    .iloc[:, 0]
    .astype(str)
    .str.upper()
    .tolist()
)

MetScore = PositiveGenes + NegativeGenes

# intersect with genes actually present
var_syms = adata.raw.var_names.str.upper()

pos_in = list(sorted(set(PositiveGenes) & set(var_syms)))
neg_in = list(sorted(set(NegativeGenes) & set(var_syms)))

print(f"Positive Met-Score genes in data: {len(pos_in)}")
print(f"Negative Met-Score genes in data: {len(neg_in)}")

if len(pos_in) == 0 or len(neg_in) == 0:
    raise ValueError("No overlap between Met-Score genes and adata.var_names — check symbols!")

##########################################
# Compute Met-Score per cell
##########################################
# Score positive and negative modules separately
sc.tl.score_genes(
    adata,
    gene_list=pos_in,
    score_name="met_score_pos",
    use_raw=True
)

sc.tl.score_genes(
    adata,
    gene_list=neg_in,
    score_name="met_score_neg",
    use_raw=True
)

adata.obs["met_score"] = adata.obs["met_score_pos"] - adata.obs["met_score_neg"]

# Z-score Met-Score across all cells (for nicer cross-plotting)
adata.obs["met_score_z"] = (
    (adata.obs["met_score"] - adata.obs["met_score"].mean()) /
    adata.obs["met_score"].std(ddof=0)
)

print(adata.obs[["met_score_pos", "met_score_neg", "met_score", "met_score_z"]].head())


# -------------------------------------------------------------------
# Panel SxB: UMAP met_score + disease_group
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
plt.savefig(FIG_DIR / "SxB_umap_met_score.png", dpi=400, bbox_inches="tight")
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
plt.savefig(FIG_DIR / "SxB_umap_met_score_pos.png", dpi=400, bbox_inches="tight")
plt.close()

plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=["celltype", "met_score_pos", "SPP1"],
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
plt.savefig(FIG_DIR / "umap_metScore_Spp1_celltypes.png", dpi=400, bbox_inches="tight")
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
plt.savefig(FIG_DIR / "SxB_umap_met_score_neg.png", dpi=400, bbox_inches="tight")
plt.close()


plt.figure(figsize=(6, 5))
sc.pl.umap(
    adata,
    color=MetScore,
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
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_visible(False)

# Improve colorbar labeling
cbar = ax.collections[0].colorbar
cbar.set_label("45-gene signature module score", fontsize=11)
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig(
    FIG_DIR / "SxB_umap_met_score_magma_readable.png",
    dpi=600,
    bbox_inches="tight"
)
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
out_all.to_csv("met_score_gene_celltype_expression_WITH_GENE.csv", index=False)


# -------------------------------------------------------------------
# Restrict to tumor compartment
# -------------------------------------------------------------------
tumor_cells = adata[adata.obs["compartment"] == "Tumor"].copy()
print("Tumor cells:", tumor_cells.n_obs)

# Per-sample tumor summaries
tumor_df = tumor_cells.obs[["GSM", "disease_group", "sample_index", "met_score"]].copy()

summary_tumor = (
    tumor_df
    .groupby(["GSM", "disease_group", "sample_index"], dropna=False)
    .agg(
        n_cells=("met_score", "size"),
        met_median=("met_score", "median"),
        met_mean=("met_score", "mean")
    )
    .reset_index()
)

summary_tumor.to_csv(
    OUTS_DIR / "met_score_tumor_by_sample.csv",
    index=False
)

# Clean data (drop NaNs in met_median and group)
summary_tumor_clean = summary_tumor.dropna(subset=["disease_group", "met_median"]).copy()

# Force clinical order (consistent across plots)
group_order = ["Localized", "mHSPC", "mCRPC"]
summary_tumor_clean["disease_group"] = pd.Categorical(
    summary_tumor_clean["disease_group"],
    categories=group_order,
    ordered=True
)

# If any samples have disease_group outside the expected set, drop them
summary_tumor_clean = summary_tumor_clean[
    summary_tumor_clean["disease_group"].isin(group_order)
].copy()

print("Per-group sample counts:")
print(summary_tumor_clean["disease_group"].value_counts(dropna=False))

# Global test (for text/legend)
def kruskal_by_group(df, value_col, group_col="disease_group"):
    df = df[[group_col, value_col]].dropna()
    groups = [g for g in df[group_col].cat.categories if g in df[group_col].unique()]
    data = []
    kept = []
    for g in groups:
        vals = df.loc[df[group_col] == g, value_col].astype(float).dropna().values
        if len(vals) >= 2:
            data.append(vals)
            kept.append(g)
    if len(data) < 2:
        return np.nan, np.nan, kept
    H, p = stats.kruskal(*data)
    return H, p, kept

H_kw, p_kw, groups_used = kruskal_by_group(summary_tumor_clean, "met_median")
print(f"Kruskal–Wallis (met_median ~ group): H={H_kw:.3f}, p={p_kw:.3e}, groups={groups_used}")

# Pairwise Mann–Whitney U + BH adjusted p-values
def pairwise_mwu(df, value_col, group_col="disease_group", order=None, min_n=3):
    df = df[[group_col, value_col]].dropna()
    if order is None:
        order = df[group_col].unique().tolist()

    rows = []
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            g1, g2 = order[i], order[j]
            x = df.loc[df[group_col] == g1, value_col].astype(float).dropna().values
            y = df.loc[df[group_col] == g2, value_col].astype(float).dropna().values
            if (len(x) >= min_n) and (len(y) >= min_n):
                U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
                rows.append({"g1": g1, "g2": g2, "U": U, "p_raw": p, "n1": len(x), "n2": len(y)})

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    res["p_adj"] = multipletests(res["p_raw"], method="fdr_bh")[1]
    return res

pairwise_stats = pairwise_mwu(summary_tumor_clean, "met_median", order=group_order, min_n=3)
print("\nPairwise MWU (BH-adjusted):")
print(pairwise_stats)

# Plot helper: clean brackets with “adj. p = …”
def add_pval_bracket(ax, x1, x2, y, h, text, lw=1.2, fs=11):
    """
    Draw a bracket between x1 and x2 at height y, with bracket height h,
    and place text centered above.
    """
    ax.plot([x1, x1, x2, x2],
            [y, y + h, y + h, y],
            lw=lw, c="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + h + 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            text, ha="center", va="bottom", fontsize=fs)


# Make the figure
sns.set_theme(context="paper", style="whitegrid", font_scale=1.25)

fig, ax = plt.subplots(figsize=(4.8, 4.8))

# Boxplot (white fill, black edges)
sns.boxplot(
    data=summary_tumor_clean,
    x="disease_group",
    y="met_median",
    order=group_order,
    color="white",
    showcaps=True,
    boxprops={"edgecolor": "black", "linewidth": 1.8},
    medianprops={"color": "black", "linewidth": 2.0},
    whiskerprops={"color": "black", "linewidth": 1.6},
    capprops={"color": "black", "linewidth": 1.6},
    ax=ax
)

# Points
sns.stripplot(
    data=summary_tumor_clean,
    x="disease_group",
    y="met_median",
    order=group_order,
    color="black",
    alpha=0.65,
    size=4.5,
    jitter=0.18,
    ax=ax
)

ax.set_xlabel("")
ax.set_ylabel("Tumor median 45-gene signature module score per sample")
ax.set_title("Tumor 45-gene signature module score by clinical group")  # <- NO global p-value
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Add significant adjusted p-values
group_to_x = {g: i for i, g in enumerate(group_order)}

# Determine baseline height above data
y_max = float(summary_tumor_clean["met_median"].max())
y_min = float(summary_tumor_clean["met_median"].min())
yr = y_max - y_min if y_max > y_min else 1.0

y_start = y_max + 0.01 * yr
h = 0.03 * yr
step = 0.1 * yr

# Keep only significant comparisons
sig = pairwise_stats.loc[pairwise_stats["p_adj"] < 0.05].copy()
# Sort so shorter brackets go on top first (optional aesthetic)
# (helps avoid long bracket hiding short one)
if not sig.empty:
    sig["span"] = sig.apply(lambda r: abs(group_to_x[r["g2"]] - group_to_x[r["g1"]]), axis=1)
    sig = sig.sort_values(["span", "p_adj"], ascending=[True, True]).reset_index(drop=True)

for i, row in sig.iterrows():
    x1 = group_to_x[row["g1"]]
    x2 = group_to_x[row["g2"]]
    y = y_start + i * step
    label = f"adj. p = {row['p_adj']:.3g}"
    add_pval_bracket(ax, x1, x2, y=y, h=h, text=label, lw=1.2, fs=11)

# Expand ylim to make room for brackets
if not sig.empty:
    ax.set_ylim(y_min - 0.05 * yr, y_start + len(sig) * step + 0.10 * yr)

plt.tight_layout()

# Save
FIG_DIR.mkdir(parents=True, exist_ok=True)
out_path = FIG_DIR / "SxC_tumor_met_score_by_group_with_adj_pvals.png"
plt.savefig(out_path, dpi=600, bbox_inches="tight")
plt.close()

print(f"\nSaved: {out_path}")



# -------------------------------------------------------------------
# Panel SxD: Heatmap of mean tumor Met-Score by tumor subtype & group
# -------------------------------------------------------------------
obs_sub = adata.obs[["celltype", "compartment", "disease_group", "met_score"]].dropna()

tumor_ct = obs_sub[obs_sub["compartment"] == "Tumor"].copy()
# Restrict to the three tumor subtypes if you want a clean 3×3 grid
tumor_ct = tumor_ct[tumor_ct["celltype"].isin(
    ["Tumor luminal", "Tumor basal", "Tumor proliferative"]
)]

celltype_summary = (
    tumor_ct
    .groupby(["celltype", "disease_group"], dropna=False)
    .agg(
        n_cells=("met_score", "size"),
        met_mean=("met_score", "mean"),
        met_median=("met_score", "median")
    )
    .reset_index()
)

celltype_summary.to_csv(
    OUTS_DIR / "met_score_by_tumor_celltype_and_group.csv",
    index=False
)

heat_df = celltype_summary.pivot_table(
    index="celltype",
    columns="disease_group",
    values="met_mean"
)

plt.figure(figsize=(6, 4))
sns.heatmap(
    heat_df,
    annot=True,
    fmt=".02f",
    cmap="magma",
    cbar_kws={"label": "Mean 45-gene signature module score"}
)
plt.ylabel("Tumor cell type")
plt.xlabel("Clinical group")
plt.title("")
plt.tight_layout()
plt.savefig(FIG_DIR / "SxD_heatmap_mean_tumor_met_by_group.png", dpi=600, bbox_inches="tight")
plt.close()





# ============================================================
# Myeloid / SPP1hi TAM-like / Met-Score downstream analyses
# ============================================================


# ----------------------------
# Config / dirs
# ----------------------------
OUTS_DIR = ROOT / "outs"
FIG_DIR  = ROOT / "figures" / "scRNAseq" / "GSE274229_myeloid_SPP1hi"
OUTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.set_figure_params(dpi=150, dpi_save=400, format="png", transparent=False)

# ----------------------------
# Helpers
# ----------------------------
def _present_genes(adata_obj, genes, use_raw=True):
    genes_up = {g.upper() for g in genes}
    if use_raw and adata_obj.raw is not None:
        var_up = pd.Index([g.upper() for g in adata_obj.raw.var_names])
        var_orig = pd.Index(adata_obj.raw.var_names)
    else:
        var_up = pd.Index([g.upper() for g in adata_obj.var_names])
        var_orig = pd.Index(adata_obj.var_names)
    m = var_up.isin(list(genes_up))
    return var_orig[m].tolist()

def _bh(pvals):
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
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
    table = np.array([[a, b], [c, d]], dtype=int)
    OR, p = stats.fisher_exact(table, alternative="greater")
    return OR, p, table, dict(a=a, b=b, c=c, d=d, up=len(up), met=len(met), bg=len(bg))

def make_rank_series(de_df, method="signed_logp", padj_col="pvals_adj", lfc_col="logfoldchanges"):
    df = de_df.copy()
    if lfc_col not in df.columns and "logfoldchange" in df.columns:
        lfc_col = "logfoldchange"
    df = df.dropna(subset=["gene", padj_col, lfc_col]).copy()

    p = df[padj_col].astype(float).replace(0, np.nextafter(0, 1))
    lfc = df[lfc_col].astype(float)

    if method == "signed_logp":
        r = np.sign(lfc) * (-np.log10(p))
    elif method == "lfc":
        r = lfc
    elif method == "score":
        if "scores" not in df.columns:
            raise ValueError("No 'scores' column for method='score'")
        r = df["scores"].astype(float)
    else:
        raise ValueError("Unknown rank method")

    s = pd.Series(r.values, index=df["gene"].astype(str).values)
    # resolve duplicate genes
    s = s.groupby(s.index).apply(lambda x: x.iloc[np.argmax(np.abs(x.values))])
    s = s.sort_values(ascending=False)
    return s

def prep_up_genes(de_df, alpha=0.05):
    df = de_df.dropna(subset=["gene", "pvals_adj", "logfoldchanges"]).copy()
    up = df[(df["pvals_adj"] < alpha) & (df["logfoldchanges"] > 0)]["gene"].astype(str).tolist()
    bg = df["gene"].astype(str).tolist()
    return up, bg

# ============================================================
# 1) Subset myeloid using your existing celltype labels
# ============================================================
adata.obs["celltype"] = adata.obs["celltype"].astype("category")
adata_raw = adata.raw.to_adata().copy()
adata_raw


MYELOID_CELLTYPES = [
    "FCN1/S100A8/A9+ Monocytes",
    "FCGR3A/MS4A7+ Monocytes",
    "APOE/TREM2+ TAMs",
    "TAMs (resident-like)",
    "cDC1", "cDC2", "pDC"
]


INCLUDE_NEUTROPHILS = False
if INCLUDE_NEUTROPHILS and ("Neutrophils" in adata_raw.obs["celltype"].cat.categories):
    MYELOID_CELLTYPES += ["Neutrophils"]

myeloid_mask = adata_raw.obs["celltype"].astype(str).isin(MYELOID_CELLTYPES)
adata_my = adata_raw[myeloid_mask].copy()
print("Myeloid cells:", adata_my.n_obs)
print(adata_my.obs["celltype"].value_counts())

# # reset X to COUNTS (so we are not using scaled/HVG-subset expression)
# if "counts" not in adata_my.layers:
#     raise ValueError("layers['counts'] not found in adata_my. (It exists in your parent adata.)")
#
# adata_my.X = adata_my.layers["counts"].copy()
#
# # normalize + log1p inside the myeloid subset
# sc.pp.normalize_total(adata_my, target_sum=1e4)
# sc.pp.log1p(adata_my)

# keep full-gene log1p data for scoring/DE later
adata_my.raw = adata_my.copy()

# sanity check: should not have NaNs now
X = adata_my.X
if hasattr(X, "A"):  # sparse
    X = X.A
print("Any NaNs in adata_my.X?", np.isnan(X).any())

adata_my.X.min()
adata_my.X.max()

adata_my.raw.X.min()
adata_my.raw.X.max()

# ============================================================
# Re-run HVGs/PCA/neighbors/UMAP/Leiden on myeloid subset
# ============================================================

#sc.pp.highly_variable_genes(adata_my, n_top_genes=4000, flavor="seurat", subset=True)
sc.pp.scale(adata_my, max_value=10)
sc.tl.pca(adata_my, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata_my, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata_my)
sc.tl.leiden(adata_my, resolution=0.8, key_added="leiden_myeloid")

plt.figure(figsize=(6,5))
sc.pl.umap(adata_my, color=["celltype","leiden_myeloid"], legend_loc="on data", frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "myeloid_umap_celltype_leiden.png", dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Compute SPP1hi score and Define SPP1hi TAM-like (BOOLEAN), but do NOT relabel celltype
#    Strategy:
#      - restrict to TAMs only
#      - compute SPP1hi_score (you already did on full adata)
#      - define SPP1hi-like as top quartile within TAMs (robust)
# ============================================================
raw_genes = set(adata_my.raw.var_names)
SPP1_HI_TAM_POS = {
    "SPP1","TREM2","APOE","GPNMB","MARCO","MSR1","CSTA",
    "CTSB","CTSD","LPL","SLC40A1","LGALS3"
}
SPP1_HI_TAM_NEG = {"FCN1","S100A8","S100A9","VCAN","CCR2"}

pos_in = [g for g in SPP1_HI_TAM_POS if g in raw_genes]
neg_in = [g for g in SPP1_HI_TAM_NEG if g in raw_genes]

sc.tl.score_genes(adata_my, pos_in, score_name="SPP1hi_pos", use_raw=True)
sc.tl.score_genes(adata_my, neg_in, score_name="SPP1hi_neg", use_raw=True)
adata_my.obs["SPP1hi_score"] = adata_my.obs["SPP1hi_pos"] - adata_my.obs["SPP1hi_neg"]

tam_mask_my = adata_my.obs["celltype"].astype(str).isin(["APOE/TREM2+ TAMs", "TAMs (resident-like)"])

adata_my.obs["SPP1hi_like"] = False
adata_my.obs["SPP1hi_thr_TAM_q75"] = np.nan

if tam_mask_my.sum() >= 50:
    # z-score the module within TAMs (more stable than absolute)
    x = adata_my.obs.loc[tam_mask_my, "SPP1hi_score"].astype(float)
    z = (x - x.mean()) / x.std(ddof=0)
    adata_my.obs["SPP1hi_score_TAMz"] = np.nan
    adata_my.obs.loc[tam_mask_my, "SPP1hi_score_TAMz"] = z

    # require BOTH: module high AND SPP1 high-ish (anti-false positives)
    thr = float(np.quantile(adata_my.obs.loc[tam_mask_my, "SPP1hi_score_TAMz"], 0.75))
    adata_my.obs["SPP1hi_thr_TAM_q75"] = thr

    # optional: sg_SPP1 median gate
    if "sg_SPP1" in adata_my.obs.columns:
        spp1_thr = float(np.quantile(adata_my.obs.loc[tam_mask_my, "sg_SPP1"].astype(float), 0.60))
        hit = tam_mask_my & (adata_my.obs["SPP1hi_score_TAMz"] >= thr) & (adata_my.obs["sg_SPP1"] >= spp1_thr)
    else:
        hit = tam_mask_my & (adata_my.obs["SPP1hi_score_TAMz"] >= thr)

    adata_my.obs.loc[hit, "SPP1hi_like"] = True



print("TAM cells in myeloid:", int(tam_mask_my.sum()))
print("SPP1hi-like TAM cells:", int(adata_my.obs["SPP1hi_like"].sum()))

# create a myeloid_state label for downstream comparisons (no touching celltype)
def _myeloid_state(ct, spp1hi):
    ct = str(ct)
    if ct in ["APOE/TREM2+ TAMs", "TAMs (resident-like)"]:
        return "SPP1hi TAMs" if bool(spp1hi) else "Other TAMs"
    if ct.startswith("FCN1") or "Monocytes" in ct:
        return "Monocyte"
    if ct in ["cDC1","cDC2","pDC"]:
        return ct
    if ct == "Neutrophils":
        return "Neutrophil"
    return "Other myeloid"

adata_my.obs["myeloid_state"] = [
    _myeloid_state(ct, spp1) for ct, spp1 in zip(adata_my.obs["celltype"], adata_my.obs["SPP1hi_like"])
]
adata_my.obs["myeloid_state"] = adata_my.obs["myeloid_state"].astype("category")
adata_my.obs["myeloid_state"].value_counts()

plt.figure(figsize=(6,5))
sc.pl.umap(adata_my, color=["myeloid_state", "SPP1", "SPP1hi_pos", "SPP1hi_score","met_score"], wspace=0.4, frameon=False, show=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "myeloid_umap_state_spp1hi_met.png", dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Met-score stats across myeloid_state (cell-level)
# ============================================================
my_df = adata_my.obs[["myeloid_state","met_score"]].dropna().copy()
groups = [g for g in my_df["myeloid_state"].cat.categories if g in my_df["myeloid_state"].unique()]

kw_data = [my_df.loc[my_df["myeloid_state"] == g, "met_score"].astype(float).values for g in groups]
H_kw, p_kw = stats.kruskal(*kw_data) if len(kw_data) >= 2 else (np.nan, np.nan)
print(f"[Myeloid states] Kruskal–Wallis met_score: H={H_kw:.3f}, p={p_kw:.3e}")

pairs = []
for i in range(len(groups)):
    for j in range(i+1, len(groups)):
        g1, g2 = groups[i], groups[j]
        x = my_df.loc[my_df["myeloid_state"] == g1, "met_score"].astype(float).values
        y = my_df.loc[my_df["myeloid_state"] == g2, "met_score"].astype(float).values
        if len(x) >= 20 and len(y) >= 20:
            U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            pairs.append((g1, g2, len(x), len(y), U, p))

pair_df = pd.DataFrame(pairs, columns=["g1","g2","n1","n2","U","p_raw"])
if not pair_df.empty:
    pair_df["p_adj"] = _bh(pair_df["p_raw"].values)

pair_df.to_csv(OUTS_DIR / "myeloid_met_score_pairwise_MWU.csv", index=False)

plt.figure(figsize=(9,4.5))
order = [g for g in ["SPP1hi TAMs","Other TAMs","Monocyte","cDC1","cDC2","pDC","Neutrophil"] if g in groups]
sns.violinplot(data=my_df, x="myeloid_state", y="met_score", order=order, cut=0, inner=None)
sns.boxplot(data=my_df, x="myeloid_state", y="met_score", order=order, width=0.25,
            boxprops={"facecolor":"white"}, medianprops={"color":"black"})
plt.xticks(rotation=25, ha="right")
plt.title("Per-cell 45-gene signature module score across myeloid states")
plt.tight_layout()
plt.savefig(FIG_DIR / "met_score_by_myeloid_state_violin.png", dpi=600, bbox_inches="tight")
plt.close()


def cliffs_delta(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    # fast-ish approximation using ranks if large
    import scipy.stats as ss
    allv = np.concatenate([x, y])
    r = ss.rankdata(allv)
    rx = r[:len(x)].sum()
    n1, n2 = len(x), len(y)
    U = rx - n1*(n1+1)/2
    delta = (2*U)/(n1*n2) - 1
    return float(delta)

x = my_df.loc[my_df["myeloid_state"]=="SPP1hi TAMs", "met_score"].astype(float).values
y = my_df.loc[my_df["myeloid_state"]=="Other TAMs", "met_score"].astype(float).values
print("Cliff's delta (SPP1hi-like vs Other TAM):", cliffs_delta(x, y))

# ============================================================
# 5) DE: SPP1hi-like TAM vs Other TAM  (contrast A)
# ============================================================
# -- Contrast A: within TAMs only
mask_A = adata_my.obs["myeloid_state"].astype(str).isin(["SPP1hi TAMs","Other TAMs"])
adata_A = adata_my[mask_A].copy()
adata_A.obs["grp_A"] = adata_A.obs["myeloid_state"].astype(str).astype("category")
adata_A.obs["grp_A"].value_counts()

sc.tl.rank_genes_groups(
    adata_A,
    groupby="grp_A",
    groups=["SPP1hi TAMs"],
    reference="Other TAMs",
    method="wilcoxon",
    use_raw=True,
    key_added="rg_A"
)
de_A = sc.get.rank_genes_groups_df(adata_A, group="SPP1hi TAMs", key="rg_A").copy()
de_A.rename(columns={"names":"gene"}, inplace=True)
de_A.to_csv(OUTS_DIR / "DE_SPP1hiLikeTAM_vs_OtherTAM_wilcoxon.csv", index=False)


# ============================================================
# 6) Fisher + GSEA (Met genes)
# ============================================================
# we already loaded MetScore gene list earlier as `MetScore`
MetScoreGenes = list(MetScore)

upA, bgA = prep_up_genes(de_A)
OR_A, p_A, tab_A, counts_A = fisher_enrichment(upA, bgA, MetScoreGenes)

enrich_df = pd.DataFrame([
    {"contrast":"SPP1hiLike vs OtherTAM", "OR":OR_A, "p_fisher":p_A, **counts_A}])

enrich_df["neglog10p"] = -np.log10(enrich_df["p_fisher"].replace(0, np.nextafter(0,1)))
enrich_df.to_csv(OUTS_DIR / "MetGenes_Fisher_enrichment_SPP1hiLikeTAM.csv", index=False)
print(enrich_df)

# --- GSEA prerank for Met-Score set only
met_set_dict = {"Met-Score": sorted({g.upper() for g in MetScoreGenes})}

def run_prerank(de_df, outdir):
    r = make_rank_series(de_df, method="signed_logp")
    r_df = pd.DataFrame({"gene": [g.upper() for g in r.index], "score": r.values})
    out = gp.prerank(
        rnk=r_df,
        gene_sets=met_set_dict,
        processes=4,
        permutation_num=1000,
        seed=7,
        outdir=str(outdir),
        min_size=5,
        max_size=5000,
        verbose=False
    )
    return out.res2d.copy()

gsea_dirA = OUTS_DIR / "GSEA_SPP1hiLike_vs_OtherTAM"
gsea_dirA.mkdir(exist_ok=True, parents=True)

gsea_A = run_prerank(de_A, gsea_dirA)

gsea_A["contrast"] = "SPP1hiLike vs OtherTAM"
gsea_all = gsea_A.reset_index().rename(columns={"index":"Term"})
gsea_all.to_csv(OUTS_DIR / "GSEA_MetScore_prerank_SPP1hiLikeTAM.csv", index=False)

# ============================================================
# 7) Plots:
#   (i) heatmap of Met genes across myeloid Leiden clusters
#   (ii) enrichment bar
# ============================================================
met_present = _present_genes(adata_my, MetScoreGenes, use_raw=True)
print("Met genes present in myeloid:", len(met_present))

# keep heatmap readable: pick top 40 variable met genes in myeloid (raw)
Xmet = adata_my.raw[:, met_present].X
if sp.issparse(Xmet):
    mean = np.asarray(Xmet.mean(axis=0)).ravel()
    mean_sq = np.asarray(Xmet.multiply(Xmet).mean(axis=0)).ravel()
    var = mean_sq - mean**2
else:
    var = np.var(Xmet, axis=0)

met_var = pd.Series(var, index=met_present).sort_values(ascending=False)
TOP_HEAT = min(45, len(met_var))
met_heat_genes = met_var.head(TOP_HEAT).index.tolist()

plt.figure(figsize=(9, 10))
sc.pl.heatmap(
    adata_my,
    var_names=met_heat_genes,
    groupby="leiden_myeloid",
    use_raw=True,
    swap_axes=True,
    show=False
)
plt.tight_layout()
plt.savefig(FIG_DIR / "heatmap_MetGenes_by_myeloid_leiden.png", dpi=600, bbox_inches="tight")
plt.close()

# enrichment summary bar
plot_enrich = enrich_df[["contrast","neglog10p"]].copy()

# try to pull NES if present
def _get_nes(res2d):
    cols = {c.lower(): c for c in res2d.columns}
    nes_col = cols.get("NES", None)
    if nes_col is None:
        return np.nan
    if "Met-Score" in res2d.index:
        return float(res2d.loc["Met-Score", nes_col])
    # fallback
    hit = res2d.index.to_series().astype(str).str.contains("Met", case=False)
    if hit.any():
        return float(res2d.loc[hit.idxmax(), nes_col])
    return np.nan

plot_enrich["GSEA_NES"] = gsea_A["NES"]
plot_enrich.to_csv(OUTS_DIR / "Enrichment_summary_Fisher_GSEA_SPP1hiLike.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
sns.barplot(data=plot_enrich, x="contrast", y="neglog10p", ax=axes[0])
axes[0].set_ylabel("-log10(Fisher p)")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=20)

sns.barplot(data=plot_enrich, x="contrast", y="GSEA_NES", ax=axes[1])
axes[1].set_ylabel("GSEA NES (45-gene signature)")
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=20)

plt.tight_layout()
plt.savefig(FIG_DIR / "enrichment_bar_Fisher_and_GSEA.png", dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Met-score stats across myeloid_state (sample-level)
# ============================================================

OUTS_DIR = ROOT / "outs"
FIG_DIR  = ROOT / "figures" / "scRNAseq" / "GSE274229_myeloid_SPP1hi"
OUTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

df = adata_my.obs[["GSM", "myeloid_state", "met_score"]].copy()
df["GSM"] = df["GSM"].astype(str)
df["myeloid_state"] = df["myeloid_state"].astype(str)
df["met_score"] = pd.to_numeric(df["met_score"], errors="coerce")
df = df.dropna(subset=["GSM", "myeloid_state", "met_score"])

# restrict to the two TAM groups
df_tam = df[df["myeloid_state"].isin(["SPP1hi TAMs", "Other TAMs"])].copy()

# per-sample median (and n cells)
per_sample = (
    df_tam
    .groupby(["GSM", "myeloid_state"], as_index=False)
    .agg(
        n_cells=("met_score", "size"),
        med_met=("met_score", "median"),
        mean_met=("met_score", "mean")
    )
)

# wide format to get paired medians per GSM
wide = per_sample.pivot(index="GSM", columns="myeloid_state", values="med_met")
wide_n = per_sample.pivot(index="GSM", columns="myeloid_state", values="n_cells")

# keep only GSMs that have both TAM groups
wide = wide.dropna(subset=["SPP1hi TAMs", "Other TAMs"]).copy()
wide_n = wide_n.loc[wide.index].copy()

# delta per sample (paired)
wide["delta_med"] = wide["SPP1hi TAMs"] - wide["Other TAMs"]

# Save tables
per_sample.to_csv(OUTS_DIR / "GSE274229_myeloid_TAM_met_score_per_sample.csv", index=False)
wide.reset_index().to_csv(OUTS_DIR / "GSE274229_myeloid_TAM_met_score_per_sample_WIDE.csv", index=False)
wide_n.reset_index().to_csv(OUTS_DIR / "GSE274229_myeloid_TAM_counts_per_sample_WIDE.csv", index=False)

print("GSMs with both TAM groups:", wide.shape[0])
print(wide[["SPP1hi TAMs", "Other TAMs", "delta_med"]].describe())

# Paired test: Wilcoxon signed-rank on per-sample deltas
d = wide["delta_med"].astype(float).values
# Remove exact zeros if needed (wilcoxon can complain if all zeros)
d_nonzero = d[d != 0]

if len(d_nonzero) < 5:
    print("[WARN] Too few non-zero paired deltas for a stable signed-rank test.")
    W, p = np.nan, np.nan
else:
    W, p = stats.wilcoxon(d_nonzero, alternative="two-sided")  # or "greater" if you pre-specify direction

print(f"Paired Wilcoxon signed-rank on per-sample delta_med: W={W}, p={p:.3e}")
print(f"Median delta_med across samples: {np.median(d):.4f}")


plot_df = wide.reset_index()[["GSM", "SPP1hi TAMs", "Other TAMs", "delta_med"]].copy()

# Ensure numeric
plot_df["SPP1hi TAMs"] = plot_df["SPP1hi TAMs"].astype(float)
plot_df["Other TAMs"] = plot_df["Other TAMs"].astype(float)

# Overlay boxplots of per-sample medians
plot_long = plot_df.melt(
    id_vars=["GSM"],
    value_vars=["Other TAMs", "SPP1hi TAMs"],
    var_name="myeloid_state",
    value_name="med_met"
)

label_other = "Other TAMs"
label_spp1hi = r"SPP1$^{hi}$ TAMs"

plot_long["myeloid_state"] = plot_long["myeloid_state"].replace({
    "SPP1hi TAMs": r"SPP1$^{hi}$ TAMs",
    "Other TAMs": "Other TAMs"
})
plot_long["myeloid_state"].value_counts()

order = [label_other, label_spp1hi]


fig, ax = plt.subplots(figsize=(5.2, 4.6))

sns.boxplot(
    data=plot_long,
    x="myeloid_state",
    y="med_met",
    order=order,
    width=0.35,
    showcaps=True,
    boxprops={"facecolor":"white", "edgecolor":"black", "linewidth":1.5},
    medianprops={"color":"black", "linewidth":1.8},
    whiskerprops={"color":"black", "linewidth":1.2},
    capprops={"color":"black", "linewidth":1.2},
    ax=ax
)

# Overlay points
sns.stripplot(
    data=plot_long,
    x="myeloid_state",
    y="med_met",
    order=order,
    color="black",
    alpha=0.6,
    size=3.5,
    jitter=0.12,
    ax=ax
)

ax.set_xlabel("")
ax.set_ylabel("Per-sample median 45-gene signature module score (TAMs)")
ax.set_title("TAM 45-gene signature module score per sample")

# Add p-value annotation
y_max = plot_long["med_met"].max()
y_min = plot_long["med_met"].min()
yr = y_max - y_min

y_bar = y_max + 0.08 * yr
h = 0.03 * yr

# bracket
ax.plot([0, 0, 1, 1], [y_bar, y_bar + h, y_bar + h, y_bar],
        lw=1.5, c="black", clip_on=False)

ax.text(
    0.5,
    y_bar + h + 0.02 * yr,
    f"paired Wilcoxon p = {p:.2e}",
    ha="center",
    va="bottom",
    fontsize=11
)

ax.set_ylim(y_min - 0.05 * yr, y_bar + h + 0.12 * yr)

sns.despine(ax=ax)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(FIG_DIR / "TAM_met_score_per_sample_paired.png", dpi=600, bbox_inches="tight")
plt.close()


# ---------------------------------------------------------
# SPP1hi-like TAM fraction by disease group
# ---------------------------------------------------------

df = adata_my.obs[["GSM","disease_group","myeloid_state","SPP1hi_like","celltype"]].copy()
df["GSM"] = df["GSM"].astype(str)
df["disease_group"] = df["disease_group"].astype(str)

# Restrict to TAMs only
tam = df[df["myeloid_state"].isin(["SPP1hi TAMs","Other TAMs"])].copy()

# Per GSM counts
frac = (
    tam.groupby(["GSM","disease_group"], as_index=False)
       .agg(n_TAM=("myeloid_state","size"),
            n_SPP1hi=("SPP1hi_like","sum"))
)
frac["frac_SPP1hi_TAM"] = frac["n_SPP1hi"] / frac["n_TAM"]

# Global test across disease_group
groups = [g for g in ["Localized","mHSPC","mCRPC"] if g in frac["disease_group"].unique()]
data = [frac.loc[frac["disease_group"]==g, "frac_SPP1hi_TAM"].values for g in groups]
H, p_kw = stats.kruskal(*data)
print("Kruskal frac_SPP1hi_TAM ~ disease_group:", p_kw)

# Pairwise MWU with BH-adjusted p
pairs=[]
for i in range(len(groups)):
    for j in range(i+1,len(groups)):
        g1,g2=groups[i],groups[j]
        x=frac.loc[frac["disease_group"]==g1,"frac_SPP1hi_TAM"].values
        y=frac.loc[frac["disease_group"]==g2,"frac_SPP1hi_TAM"].values
        U,p=stats.mannwhitneyu(x,y,alternative="two-sided")
        pairs.append((g1,g2,p))
pairs_df=pd.DataFrame(pairs,columns=["g1","g2","p_raw"])
pairs_df["p_adj"] = multipletests(pairs_df["p_raw"], method="fdr_bh")[1]
print(pairs_df)

# ---------------------------------------------------------
# Does Δ median Met-Score (SPP1hi−Other TAM) increase with disease group?
# ---------------------------------------------------------

# I already created wide with delta_med; just add disease_group
wide = wide.reset_index()
wide["disease_group"] = wide["GSM"].map(
    df.drop_duplicates("GSM").set_index("GSM")["disease_group"].to_dict()
)

# Global KW across disease groups
data = [wide.loc[wide["disease_group"]==g,"delta_med"].values for g in groups]
H, p_kw = stats.kruskal(*data)
print("Kruskal delta_med ~ disease_group:", p_kw)

# Pairwise MWU on delta_med
pairs=[]
for i in range(len(groups)):
    for j in range(i+1,len(groups)):
        g1,g2=groups[i],groups[j]
        x=wide.loc[wide["disease_group"]==g1,"delta_med"].values
        y=wide.loc[wide["disease_group"]==g2,"delta_med"].values
        U,p=stats.mannwhitneyu(x,y,alternative="two-sided")
        pairs.append((g1,g2,p))
pairs_df=pd.DataFrame(pairs,columns=["g1","g2","p_raw"])
pairs_df["p_adj"] = multipletests(pairs_df["p_raw"], method="fdr_bh")[1]
print(pairs_df)
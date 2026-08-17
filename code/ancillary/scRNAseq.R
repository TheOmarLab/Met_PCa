#######################################################################
## Single-cell: GSE274229 – Met-Score module in myeloid cells
#######################################################################

rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(MetaIntegrator)
library(GEOquery)
library(purrr)
library(readr)
library(Matrix)

## Load Met-Score objects
load("./outs/PP_filter_MetaScore.rda")          
load("./outs/filtersiggenes_MetaScore.rda")     

PositiveGenes <- filter$posGeneNames
NegativeGenes <- filter$negGeneNames


#######################################################################
## 0. Process the single cell data from GSE274229
#######################################################################

# Paths to the 10X directories
base_dir <- "./data/GSE274229"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Download all supplementary files for GSE274229
#    This will include:
#    - GSE274229_PCADT1_* (mtx / features / barcodes)
#    - GSE274229_PCADT2_*
#    - GSE274229_PCPBS1_*
#    - GSE274229_PCPBS2_*
#    - GSE274229_RAW.tar (all raw files)
getGEOSuppFiles("GSE274229", baseDir = base_dir)


# GEOquery createed a subfolder named "GSE274229"
list.files(base_dir)

supp_dir <- file.path(base_dir, "GSE274229")
list.files(supp_dir)

# Untar the big RAW archive into a "RAW" folder
raw_tar <- file.path(supp_dir, "GSE274229_RAW.tar")

raw_outdir <- file.path(base_dir, "RAW")
dir.create(raw_outdir, showWarnings = FALSE)

untar(raw_tar, exdir = raw_outdir)

# Inspect what was extracted
list.files(raw_outdir)

## Base paths
base_dir <- "./data/GSE274229"
raw_dir  <- file.path(base_dir, "RAW")

## All GEX feature files
feat_files <- list.files(
  raw_dir,
  pattern = "_GEX_features.tsv.gz$",
  full.names = TRUE
)

## Helper: classify species based on gene IDs
classify_species <- function(feat_file) {
  feats <- read_tsv(feat_file, col_names = FALSE, show_col_types = FALSE)
  gene_ids <- feats[[1]]
  
  prop_human <- mean(grepl("^ENSG", gene_ids))
  prop_mouse <- mean(grepl("^ENSMUSG", gene_ids))
  
  species <- dplyr::case_when(
    prop_human > 0.5 ~ "human",
    prop_mouse > 0.5 ~ "mouse",
    TRUE             ~ "mixed/unknown"
  )
  
  tibble(
    features_file = feat_file,
    library_id    = sub("_GEX_features.tsv.gz$", "", basename(feat_file)),
    prop_human    = prop_human,
    prop_mouse    = prop_mouse,
    species       = species
  )
}

species_df <- map_dfr(feat_files, classify_species)
species_df

## Keep human libraries only
human_libs <- species_df %>%
  filter(species == "human")

human_libs


## All GEX matrix files
mat_files <- list.files(
  raw_dir,
  pattern = "_GEX_matrix\\.mtx\\.gz$",
  full.names = TRUE
)


## Helper: read one 10x-like triple (matrix, features, barcodes)
read_one_library <- function(mtx_file) {
  lib_id <- basename(mtx_file)
  lib_id <- sub("_GEX_matrix.mtx.gz$", "", lib_id)
  
  message("Reading library: ", lib_id)
  
  feat_file    <- file.path(raw_dir, paste0(lib_id, "_GEX_features.tsv.gz"))
  barcode_file <- file.path(raw_dir, paste0(lib_id, "_GEX_barcodes.tsv.gz"))
  
  ## --- 1) Read features robustly (no quoting) ---
  feats <- read.delim(
    gzfile(feat_file),
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    sep = "\t"
  )
  if (ncol(feats) < 2) {
    stop("Features file ", feat_file, " has <2 columns after parsing. Check download.")
  }
  
  gene_ids  <- feats[[1]]
  gene_syms <- feats[[2]]
  gene_syms <- make.unique(gene_syms)
  
  ## --- 2) Read matrix ---
  mtx <- readMM(gzfile(mtx_file))
  
  ## --- 3) Read barcodes ---
  barcodes <- read.delim(
    gzfile(barcode_file),
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    sep = "\t"
  )[[1]]
  
  ## --- 4) Check dimensions ---
  if (nrow(mtx) != length(gene_syms)) {
    stop(
      "For ", lib_id,
      ": matrix has ", nrow(mtx), " rows but features has ", length(gene_syms),
      " rows. Features file is likely corrupted or mis-downloaded."
    )
  }
  if (ncol(mtx) != length(barcodes)) {
    stop(
      "For ", lib_id,
      ": matrix has ", ncol(mtx), " columns but barcodes has ", length(barcodes),
      " entries."
    )
  }
  
  rownames(mtx) <- gene_syms      # use SYMBOLS for your Met-Score
  colnames(mtx) <- paste0(lib_id, "_", barcodes)
  
  obj <- CreateSeuratObject(
    counts      = mtx,
    project     = "GSE274229",
    min.cells   = 3,
    min.features= 200
  )
  obj$library_id <- lib_id
  
  return(obj)
}

## Read all libraries into a list
seu_list <- lapply(mat_files, read_one_library)
names(seu_list) <- sub("_GEX_matrix\\.mtx\\.gz$", "", basename(mat_files))

## Merge all libraries into one Seurat object
cat("Merging", length(seu_list), "libraries...\n")

sce274 <- Reduce(function(x, y) {
  merge(x, y, add.cell.ids = c(x$library_id[1], y$library_id[1]))
}, seu_list)

sce274
table(sce274$library_id)

## Save raw merged object
saveRDS(sce274, file = "./outs/GSE274229_seurat_raw_allLibs.rds")

#######################################################################
## annotate samples 
#######################################################################

## load
sce274 <- readRDS("./outs/GSE274229_seurat_raw_allLibs.rds")

# annotate samples 
## Extract GSM ID from library_id
meta_lib <- data.frame(
  library_id = sce274$library_id,
  stringsAsFactors = FALSE
) %>%
  distinct() %>%
  mutate(
    GSM = sub("_.*$", "", library_id)
  )

meta_lib

## Add disease group based on GSM numeric range
meta_lib <- meta_lib %>%
  mutate(
    gsm_num = as.integer(sub("GSM", "", GSM)),
    disease_group = case_when(
      gsm_num >= 8445680 & gsm_num <= 8445685 ~ "mCRPC",
      gsm_num >= 8445686 & gsm_num <= 8445710 ~ "mHSPC",
      gsm_num >= 8445711 & gsm_num <= 8445723 ~ "Localized",
      TRUE ~ NA_character_
    ),
    ## Optional: an index 1–44 (matches “sample 1 … sample 44”)
    sample_index = case_when(
      gsm_num >= 8445680 & gsm_num <= 8445723 ~ gsm_num - 8445679L,
      TRUE ~ NA_integer_
    )
  )

table(meta_lib$disease_group, useNA = "ifany")

## Merge per-library metadata onto each cell
## NOTE: dplyr::left_join() drops the data.frame rownames (the cell barcodes).
## Assigning a row-name-less data.frame back into @meta.data silently renames the
## object's cells to "1..N", which later makes JoinLayers()/[[<- fail with
## "Cannot add new cells with [[<-" because the assay still carries the real
## barcodes. Preserve the cell barcodes across the join.
md <- sce274@meta.data
cell_barcodes <- rownames(md)

md <- md %>%
  left_join(meta_lib, by = "library_id")

## left_join returns row order-preserving output; restore the cell barcodes.
stopifnot(nrow(md) == length(cell_barcodes))
rownames(md) <- cell_barcodes

# Overwrite meta.data and set convenient factors
sce274@meta.data <- md

sce274$GSM           <- sce274@meta.data$GSM
sce274$disease_group <- factor(
  sce274@meta.data$disease_group,
  levels = c("Localized", "mHSPC", "mCRPC")
)
sce274$sample_index  <- sce274@meta.data$sample_index

## Quick checks
table(sce274$disease_group)
length(unique(sce274$GSM))


#######################################################################
## 1. QC, normalization, dimensionality reduction, clustering
#######################################################################

set.seed(1234)

DefaultAssay(sce274) <- "RNA"

# 1) Join all "counts" layers into a single "counts" layer.
# (Seurat v5 uses JoinLayers(); AggregateLayers is not a Seurat function. Calling
# JoinLayers on the assay with no explicit `layers`/`new` collapses every per-library
# counts.* layer into one unified "counts" layer, which is the intent here. Passing an
# explicit list of 44 layers errors because `new` would then need 44 names.)
sce274 <- JoinLayers(
  object = sce274,
  assay  = "RNA"
)


## 1) Pick a sane counts layer (use the last one)
all_layers    <- Layers(sce274[["RNA"]])
counts_layer  <- tail(all_layers, 1)
counts_layer
rna_counts <- LayerData(sce274[["RNA"]], layer = counts_layer)
total_counts <- Matrix::colSums(rna_counts)
stopifnot(length(total_counts) == ncol(sce274))
names(total_counts) <- colnames(sce274)


mt_genes  <- grep("^MT-", rownames(rna_counts), value = TRUE)
length(mt_genes)
mt_counts <- Matrix::colSums(rna_counts[mt_genes, , drop = FALSE])

sce274$nCount_RNA <- total_counts
sce274$percent.mt <- (mt_counts / total_counts) * 100
summary(sce274$percent.mt)



# Basic QC plots
# VlnPlot(
#   sce274,
#   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#   group.by = "disease_group",
#   ncol = 3, pt.size = 0.1
# )

## Filter low-quality cells
sce274 <- subset(
  sce274,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 20
)

cat("Remaining cells after QC:\n")
print(dim(sce274))

## 1.3 SCTransform normalization
sce274 <- SCTransform(
  sce274,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  verbose = TRUE
)

## 1.4 PCA + UMAP + neighbors + clustering
sce274 <- RunPCA(sce274, verbose = FALSE)
ElbowPlot(sce274)

n_pcs <- 30  # adjust based on ElbowPlot if needed

sce274 <- RunUMAP(sce274, dims = 1:n_pcs)
sce274 <- FindNeighbors(sce274, dims = 1:n_pcs)
sce274 <- FindClusters(sce274, resolution = 0.5)

# UMAP by cluster and disease group
DimPlot(sce274, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("GSE274229 – clusters")

DimPlot(sce274, reduction = "umap", group.by = "disease_group") +
  ggtitle("GSE274229 – disease groups")

#######################################################################
## 2. Coarse cell-type annotation using marker gene module scores
#######################################################################

DefaultAssay(sce274) <- "SCT"

## 2.1 Define broad marker gene sets (adjust if needed)
tcell_genes <- c("CD3D", "CD3E", "CD2", "CD4", "CD8A", "CD8B", "TRAC")
bcell_genes <- c("MS4A1", "CD79A", "CD79B", "CD74", "MZB1")
nk_genes    <- c("NKG7", "GNLY", "PRF1", "KLRD1")
myeloid_genes <- c("LYZ", "S100A8", "S100A9", "LST1", "CTSS",
                   "CSF1R", "MS4A7", "FCGR3A", "ITGAM", "ITGAX")
epithelial_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14")
endothelial_genes <- c("PECAM1", "VWF", "KDR", "ENG", "RGS5")
fibro_genes <- c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA")
plasma_genes <- c("MZB1", "XBP1", "SDC1", "IGHG1", "JCHAIN")

gene_sets <- list(
  Tcell       = tcell_genes,
  Bcell       = bcell_genes,
  NK          = nk_genes,
  Myeloid     = myeloid_genes,
  Epithelial  = epithelial_genes,
  Endothelial = endothelial_genes,
  Fibroblast  = fibro_genes,
  Plasma      = plasma_genes
)

## 2.2 Module scores for each lineage
for (nm in names(gene_sets)) {
  genes_use <- intersect(gene_sets[[nm]], rownames(sce274))
  if (length(genes_use) < 3) {
    warning("Skipping ", nm, " – fewer than 3 marker genes found in dataset.")
    next
  }
  sce274 <- AddModuleScore(
    sce274,
    features = list(genes_use),
    name     = nm,
    assay    = "SCT"
  )
  # Rename the last added column to something clean (e.g. "TcellScore")
  new_col <- tail(colnames(sce274@meta.data), 1)
  colnames(sce274@meta.data)[colnames(sce274@meta.data) == new_col] <- paste0(nm, "Score")
}

## 2.3 Assign coarse cell type as argmax of scores
score_cols <- paste0(names(gene_sets), "Score")
score_cols <- intersect(score_cols, colnames(sce274@meta.data))

score_mat <- as.matrix(sce274@meta.data[, score_cols, drop = FALSE])

max_idx <- apply(
  score_mat,
  1,
  function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.max(x)
  }
)

celltype_labels <- ifelse(
  is.na(max_idx),
  "Unknown",
  names(gene_sets)[max_idx]
)

sce274$celltype_coarse <- factor(celltype_labels)

table(sce274$celltype_coarse, useNA = "ifany")

DimPlot(
  sce274,
  reduction = "umap",
  group.by = "celltype_coarse",
  label = TRUE
) + ggtitle("Coarse cell types (marker-based)")

#######################################################################
## 3. Focus on myeloid cells and (optional) subclustering
#######################################################################

## 3.1 Subset myeloid compartment
sce_myeloid <- subset(sce274, subset = celltype_coarse == "Myeloid")

cat("Myeloid subset dimensions:\n")
print(dim(sce_myeloid))
table(sce_myeloid$disease_group)

## 3.2 Re-run normalization / clustering within myeloid (optional but recommended)
DefaultAssay(sce_myeloid) <- "SCT"

sce_myeloid <- RunPCA(sce_myeloid, verbose = FALSE)
ElbowPlot(sce_myeloid)

n_pcs_myeloid <- 20

sce_myeloid <- RunUMAP(sce_myeloid, dims = 1:n_pcs_myeloid)
sce_myeloid <- FindNeighbors(sce_myeloid, dims = 1:n_pcs_myeloid)
sce_myeloid <- FindClusters(sce_myeloid, resolution = 0.4)

DimPlot(
  sce_myeloid,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE
) + ggtitle("Myeloid subclusters")

## 3.3 (Optional) inspect canonical myeloid markers across clusters
myeloid_markers_viz <- c("LYZ", "S100A8", "S100A9", "LST1", "CTSS", "CSF1R",
                         "MS4A7", "FCGR3A")

DotPlot(
  sce_myeloid,
  features = myeloid_markers_viz,
  group.by = "seurat_clusters"
) + RotatedAxis() + ggtitle("Myeloid marker expression per cluster")

#######################################################################
## 4. Compute the 45-gene signature module score in myeloid cells
##    (PositiveGenes / NegativeGenes)
##
## Naming note: the per-cell score built below with AddModuleScore
## (mean of positive genes minus mean of negative genes) is a directional
## 45-gene signature module score. It is NOT the frozen 41-feature
## ridge-logistic Met-Score classifier, which returns a bulk logistic-
## regression probability and is not computable per single cell. Keep the
## two names distinct in every downstream label and comment.
#######################################################################

## 4.1 Make sure the signature gene lists intersect the data
DefaultAssay(sce_myeloid) <- "SCT"

pos_in_data <- intersect(PositiveGenes, rownames(sce_myeloid))
neg_in_data <- intersect(NegativeGenes, rownames(sce_myeloid))

cat("Positive genes in data:", length(pos_in_data), "out of", length(PositiveGenes), "\n")
cat("Negative genes in data:", length(neg_in_data), "out of", length(NegativeGenes), "\n")

if (length(pos_in_data) < 3) {
  warning("Very few positive Met-Score genes found in data.")
}
if (length(neg_in_data) < 3) {
  warning("Very few negative Met-Score genes found in data.")
}

## 4.2 Compute module scores for positive and negative sets
sce_myeloid <- AddModuleScore(
  sce_myeloid,
  features = list(pos_in_data),
  name     = "MetPos",
  assay    = "SCT"
)
sce_myeloid <- AddModuleScore(
  sce_myeloid,
  features = list(neg_in_data),
  name     = "MetNeg",
  assay    = "SCT"
)

# AddModuleScore created MetPos1 and MetNeg1
colnames(sce_myeloid@meta.data)[
  colnames(sce_myeloid@meta.data) == "MetPos1"
] <- "MetPosScore"

colnames(sce_myeloid@meta.data)[
  colnames(sce_myeloid@meta.data) == "MetNeg1"
] <- "MetNegScore"

## 4.3 Define the 45-gene signature module score as positive minus negative
## (directional signature module score, not the frozen 41-feature classifier)
sce_myeloid$MetScore <- sce_myeloid$MetPosScore - sce_myeloid$MetNegScore

summary(sce_myeloid$MetScore)

#######################################################################
## 5. Visualization of Met-Score in myeloid cells
#######################################################################

## 5.1 UMAP colored by Met-Score
FeaturePlot(
  sce_myeloid,
  features = "MetScore",
  cols = c("navy", "white", "firebrick"),
  min.cutoff = "q05",
  max.cutoff = "q95"
) + ggtitle("45-gene signature module score in myeloid cells (per cell)")

## 5.2 Violin plot by disease group
VlnPlot(
  sce_myeloid,
  features = "MetScore",
  group.by = "disease_group",
  pt.size  = 0.1
) + ggtitle("45-gene signature module score in myeloid cells by disease group")

## 5.3 Violin plot by myeloid subcluster
VlnPlot(
  sce_myeloid,
  features = "MetScore",
  group.by = "seurat_clusters",
  pt.size  = 0.1
) + ggtitle("45-gene signature module score in myeloid subclusters")

#######################################################################
## 6. Sample-level summaries of Met-Score in myeloid cells
#######################################################################

library(tidyr)

## 6.1 Extract per-cell Met-Score + metadata
meta_myeloid <- sce_myeloid@meta.data %>%
  mutate(
    MetScore = sce_myeloid$MetScore,
    cell_id  = rownames(sce_myeloid@meta.data)
  )

## 6.2 Aggregate to sample (GSM) x disease_group
met_sample <- meta_myeloid %>%
  group_by(GSM, disease_group) %>%
  summarise(
    mean_MetScore = mean(MetScore, na.rm = TRUE),
    median_MetScore = median(MetScore, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

head(met_sample)

## 6.3 Boxplot of sample-level mean Met-Score by disease group
ggplot(met_sample, aes(x = disease_group, y = mean_MetScore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Disease group",
    y = "Mean 45-gene signature module score (myeloid, per GSM)",
    title = "Sample-level 45-gene signature module score in myeloid cells"
  )

#######################################################################
## 7. Heatmap of Met-Score per sample (ComplexHeatmap)
#######################################################################

## 7.1 Prepare matrix: rows = samples (GSM), cols = "MetScore"
met_mat <- met_sample %>%
  dplyr::select(GSM, mean_MetScore, disease_group) %>%
  distinct() %>%
  arrange(disease_group, GSM)

row_annot_df <- met_mat %>%
  dplyr::select(GSM, disease_group) %>%
  distinct() %>%
  column_to_rownames("GSM")

mat <- as.matrix(
  met_mat %>%
    dplyr::select(GSM, mean_MetScore) %>%
    column_to_rownames("GSM")
)

## 7.2 Row annotation for disease group
disease_cols <- c(
  "Localized" = "#1b9e77",
  "mHSPC"     = "#7570b3",
  "mCRPC"     = "#d95f02"
)

row_ha <- rowAnnotation(
  DiseaseGroup = row_annot_df$disease_group,
  col = list(DiseaseGroup = disease_cols)
)

## 7.3 Draw heatmap
Heatmap(
  mat,
  name           = "Mean module score",
  cluster_rows   = TRUE,
  cluster_columns= FALSE,
  right_annotation = row_ha,
  col            = colorRamp2(
    c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
    c("navy", "white", "firebrick")
  ),
  column_title = "Myeloid 45-gene signature module score (per sample)",
  row_names_gp = gpar(fontsize = 8)
)

#######################################################################
## 8. Save processed objects for reuse
#######################################################################

saveRDS(sce274,     file = "./outs/GSE274229_seurat_processed_allcells.rds")
saveRDS(sce_myeloid, file = "./outs/GSE274229_seurat_myeloid_MetScore.rds")





#######################################################################
## 2. Compute the 45-gene signature module score per myeloid cell
##
## Operates on the myeloid subset created above (sce_myeloid). This
## recomputes the directional 45-gene signature module score on the RNA
## assay under the name MetScore_sc (pos minus neg). As in section 4, this
## is the signature module score, not the frozen 41-feature Met-Score
## classifier.
#######################################################################

DefaultAssay(sce_myeloid) <- "RNA"

## Ensure signature genes are present in the scRNA-seq object
genes_sc <- rownames(sce_myeloid)
pos_sc   <- intersect(PositiveGenes, genes_sc)
neg_sc   <- intersect(NegativeGenes, genes_sc)

length(pos_sc)
length(neg_sc)

## Add module scores separately for pos and neg sets
sce_myeloid <- AddModuleScore(
  object   = sce_myeloid,
  features = list(pos_sc),
  name     = "MetScore_pos"
)
sce_myeloid <- AddModuleScore(
  object   = sce_myeloid,
  features = list(neg_sc),
  name     = "MetScore_neg"
)

## AddModuleScore creates columns "MetScore_pos1", "MetScore_neg1"
sce_myeloid$MetScore_sc <- sce_myeloid$MetScore_pos1 - sce_myeloid$MetScore_neg1

## Save for downstream
saveRDS(sce_myeloid, "./outs/GSE274229_myeloid_with_MetScore_sc.rds")


#######################################################################
## 3. Visualize the signature module score across myeloid subclusters
#######################################################################

## Violin plot of the signature module score across myeloid subclusters.
## seurat_clusters is the only per-cell myeloid grouping defined in this
## script (from FindClusters on the myeloid subset above); this script does
## not annotate the subclusters into named myeloid states.
p_violin_myeloid <- VlnPlot(
  sce_myeloid,
  features = "MetScore_sc",
  group.by = "seurat_clusters",
  pt.size  = 0.1
) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.title  = element_text(size = 14, face = "bold")
  ) +
  ylab("45-gene signature module score (pos minus neg)")

ggsave("./figures/GSE274229_myeloid_MetScore_sc_violin.pdf",
       plot = p_violin_myeloid, width = 6.5, height = 4.5, useDingbats = FALSE)

## UMAP colored by the signature module score
if ("umap_1" %in% colnames(Embeddings(sce_myeloid, "umap"))) {
  p_umap_myeloid <- FeaturePlot(
    sce_myeloid,
    features = "MetScore_sc",
    reduction = "umap"
  ) & theme_classic(base_size = 14)

  ggsave("./figures/GSE274229_myeloid_MetScore_sc_UMAP.pdf",
         plot = p_umap_myeloid, width = 5, height = 5, useDingbats = FALSE)
}


#######################################################################
## 4. Heatmap: average expression of signature genes per myeloid subcluster
#######################################################################

DefaultAssay(sce_myeloid) <- "RNA"

## Signature genes present in the data: positive then negative sets.
## pos_sc / neg_sc were defined in section 2 as the signature genes that
## intersect the myeloid object.
genes_heat <- c(pos_sc, neg_sc)
sce_myeloid <- ScaleData(sce_myeloid, features = genes_heat, verbose = FALSE)

avg_expr <- AverageExpression(
  sce_myeloid,
  features  = genes_heat,
  group.by  = "seurat_clusters",
  assays    = "RNA",
  slot      = "scale.data",
  verbose   = FALSE
)

mat <- avg_expr$RNA[genes_heat, , drop = FALSE]   ## genes x states

## Order: pos then neg
mat <- mat[c(pos_sc, neg_sc), , drop = FALSE]

row_split <- factor(
  c(rep("Up in Mets (pos)", length(pos_sc)),
    rep("Down in Mets (neg)", length(neg_sc))),
  levels = c("Up in Mets (pos)", "Down in Mets (neg)")
)

mx <- max(abs(mat), na.rm = TRUE)
col_fun <- circlize::colorRamp2(c(-mx, 0, mx),
                                c("#2166ac", "white", "#b2182b"))

ht_myeloid <- Heatmap(
  mat,
  name             = "Z-score",
  col              = col_fun,
  row_split        = row_split,
  cluster_rows     = FALSE,
  cluster_columns  = TRUE,
  show_row_names   = FALSE,
  column_title     = "Myeloid subclusters (GSE274229)",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp     = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

pdf("./figures/GSE274229_myeloid_MetScore_genes_heatmap.pdf",
    width = 5.5,
    height = 0.16 * nrow(mat) + 2)
draw(ht_myeloid, heatmap_legend_side = "right")
dev.off()


#######################################################################
## 5. Enrichment: signature positive genes in SPP1hi TAM upregulated genes
#######################################################################

## Use the defined myeloid subcluster identity. seurat_clusters is the only
## per-cell myeloid grouping produced in this script (FindClusters above).
Idents(sce_myeloid) <- "seurat_clusters"
table(Idents(sce_myeloid))

## BLOCKED: the SPP1hi TAM state is not defined anywhere in this script.
## This script annotates cells only at the coarse-lineage level (celltype_coarse)
## and subclusters the myeloid compartment into numeric seurat_clusters; it never
## derives an SPP1hi TAM state. Identifying which subcluster (if any) is the
## SPP1hi TAM requires marker-based TAM annotation followed by an SPP1-module
## call, none of which exists here. The companion script
## code/ancillary/PCa_GSE274229.py derives this state (SPP1 module z-score, top
## quartile within TAM-annotated cells). Set spp1_cluster_name to the correct
## seurat_clusters level after that annotation is added; it is left unchanged
## rather than thresholding SPP1 here, which this script does not define.
spp1_cluster_name <- "SPP1hi_TAM"  ## BLOCKED: set to the annotated SPP1hi TAM subcluster

if (!(spp1_cluster_name %in% levels(Idents(sce_myeloid)))) {
  stop("Cluster name ", spp1_cluster_name, " not found in Idents(sce_myeloid).  Adjust 'spp1_cluster_name'.")
}

## Differential expression: SPP1hi TAM vs all other myeloid subclusters
markers_spp1 <- FindMarkers(
  sce_myeloid,
  ident.1         = spp1_cluster_name,
  ident.2         = setdiff(levels(Idents(sce_myeloid)), spp1_cluster_name),
  logfc.threshold = 0.25,
  min.pct         = 0.10
)

## Upregulated genes in SPP1hi TAM (BH FDR < 0.05, log2FC > 0)
markers_spp1$gene <- rownames(markers_spp1)
up_spp1 <- markers_spp1 %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  pull(gene)

length(up_spp1)

## Universe: all genes expressed in myeloid cells
gene_universe <- rownames(sce_myeloid)

met_risk <- intersect(PositiveGenes, gene_universe)

## Overlap counts
a <- length(intersect(up_spp1, met_risk))                      ## in both
b <- length(setdiff(up_spp1, met_risk))                        ## in up_spp1 only
c <- length(setdiff(met_risk, up_spp1))                        ## in Met-Score pos only
d <- length(setdiff(gene_universe, union(up_spp1, met_risk)))  ## in neither

enrich_mat <- matrix(c(a, b, c, d), nrow = 2,
                     dimnames = list(
                       MetScore_pos = c("In_set", "Not_in_set"),
                       SPP1_up      = c("In_SPP1_up", "Not_SPP1_up")
                     ))

enrich_mat

fisher_res <- fisher.test(enrich_mat)
fisher_res

## Save stats
enrich_df <- data.frame(
  cluster          = spp1_cluster_name,
  a_in_both        = a,
  b_SPP1_only      = b,
  c_MetScore_only  = c,
  d_neither        = d,
  odds_ratio       = fisher_res$estimate,
  p_value          = fisher_res$p.value
)

write.csv(enrich_df,
          "./outs/GSE274229_SPP1hiTAM_MetScorePos_enrichment.csv",
          row.names = FALSE)


#######################################################################
## 6. Simple bar plot summarizing enrichment
#######################################################################

enrich_df$log10_p <- -log10(enrich_df$p_value)

p_enrich <- ggplot(enrich_df,
                   aes(x = cluster, y = log10_p)) +
  geom_col(fill = "#2b8cbe") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    x = "",
    y = expression(-log[10]("P value")),
    title = "Enrichment of signature positive genes in SPP1hi TAM up-genes"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title  = element_text(size = 14, face = "bold"),
    plot.title  = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave("./figures/GSE274229_SPP1hiTAM_MetScorePos_enrichment_bar.pdf",
       plot = p_enrich, width = 4.6, height = 4.6, useDingbats = FALSE)








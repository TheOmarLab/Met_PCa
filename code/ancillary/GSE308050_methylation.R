# ============================================================
# GSE308050: methylation DMRs & progression-linked promoter priming
#   - promoter = +/- 2000 bp around TSS
#   - ORA/Fisher enrichment for Met-Score pos/neg
#   - Concordance: pos genes hypo; neg genes hyper (in aggressive)
#
# EXPLORATORY epigenetic analysis. It uses the 45-gene Met-Score signature
# pos/neg lists (27 positive / 18 negative genes) as gene sets for the
# methylation enrichment/concordance tests. It does not use or refit the frozen
# 41-feature ridge-logistic Met-Score classifier.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(limma)
})


# ----------------------------
# 0) Inputs / outputs
# ----------------------------
out_dir <- "./outs/GSE308050_methylation_met_score"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Met-Score genes (use your exact lists)
met_pos <- c("TMSB10","ENSA","ASPN","YWHAZ","HES6","STC2","ASNS","HAVCR2","ARL6IP1","F5",
             "RFTN1","SOX4","PTPN9","ALDH1A1","MRPL11","GABRD","RC3H2","CST2","CXCR4",
             "SEM1","FOXH1","KIF7","BARD1","CADPS","RNF19A","CAMK2N1","GPR37")

met_neg <- c("KCTD14","AZGP1","PART1","CHRNA2","DPT","EDN3","KIAA1210","LTF","SIDT1","CBLL1",
             "PTN","CCK","UFM1","CPA3","CDC42EP5","TMEM121B","AKAP7","KLF4")

met_pos <- unique(toupper(met_pos))
met_neg <- unique(toupper(met_neg))
met_all <- unique(c(met_pos, met_neg))

# ----------------------------
# Download + read GEO DMR table
# ----------------------------
dmr_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308050/suppl/GSE308050_All_merged_DMRs_sig_methy_group_clinics_infor.csv.gz"
dmr_gz  <- file.path(out_dir, "GSE308050_All_merged_DMRs_sig_methy_group_clinics_infor.csv.gz")

if (!file.exists(dmr_gz)) {
  message("Downloading: ", dmr_url)
  download.file(dmr_url, destfile = dmr_gz, mode = "wb", quiet = TRUE)
} else {
  message("Found existing: ", dmr_gz)
}

dmr <- fread(dmr_gz)

stopifnot(ncol(dmr) > 2)
stopifnot("variables_and_DMR_methy_valus" %in% names(dmr))

id_col <- "variables_and_DMR_methy_valus"
sample_ids <- setdiff(names(dmr), id_col)


# ----------------------------
# Identify which rows are clinical vs DMR
# ----------------------------
# DMR rows look like: chr5_171306377_171318934
is_dmr <- grepl("^chr[0-9XYM]+_\\d+_\\d+$", dmr[[id_col]])
cat("Rows total:", nrow(dmr), "\n")
cat("DMR rows:", sum(is_dmr), "\n")
cat("Clinical rows:", sum(!is_dmr), "\n")

dmr_clin <- dmr[!is_dmr, ]
dmr_dmr  <- dmr[ is_dmr, ]

# ----------------------------
# Build phenotype table from clinical rows
# ----------------------------
# Convert clinical rows (variable x samples) to pheno (samples x variables)
pheno <- as.data.frame(dmr_clin)
rownames(pheno) <- pheno[[id_col]]
pheno[[id_col]] <- NULL
pheno <- as.data.frame(t(pheno), stringsAsFactors = FALSE)
pheno$sample_id <- rownames(pheno)

# Key outcome variable per GEO page: Bad_outcome (0/1)
if (!"Bad_outcome" %in% colnames(pheno)) {
  stop("Bad_outcome row not found. Available clinical rows:\n",
       paste(colnames(pheno), collapse = ", "))
}

# Coerce Bad_outcome into 0/1 numeric
pheno$Bad_outcome <- as.numeric(pheno$Bad_outcome)
stopifnot(all(pheno$Bad_outcome %in% c(0,1)))

table(pheno$Bad_outcome)

# Sensitivity covariates
covars <- c("surg_age", "surg_pgleasnt", "surg_pstagt")
covars <- covars[covars %in% colnames(pheno)]

# Clean numeric covars where possible
for (cv in covars) {
  pheno[[cv]] <- suppressWarnings(as.numeric(pheno[[cv]]))
}

# ----------------------------
# Build DMR methylation matrix (DMRs x samples)
# ----------------------------
dmr_ids <- dmr_dmr[[id_col]]
dmr_mat <- as.matrix(dmr_dmr[, ..sample_ids])
rownames(dmr_mat) <- dmr_ids
mode(dmr_mat) <- "numeric"

# Methylation values should be in [0,1] typically
cat("DMR matrix range:", range(dmr_mat, na.rm = TRUE), "\n")

# Make sure sample order matches between pheno and matrix
pheno <- pheno[match(sample_ids, pheno$sample_id), ]
stopifnot(all(pheno$sample_id == sample_ids))

# ----------------------------
# Differential methylation vs Bad_outcome using limma
# ----------------------------
# Helper: safe numeric coercion
as_num01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!all(is.na(x) | x %in% c(0, 1))) warning("Non {0,1} values detected after coercion.")
  x
}

# Choose outcomes
outcomes <- c("recurrence", "Bad_outcome")

# Prepare covariates
# NOTE: clinicalt may be categorical; treat as factor if non-numeric.
prep_stage <- function(x) {
  if (all(is.na(x) | grepl("^[0-9.]+$", x))) return(suppressWarnings(as.numeric(x)))
  return(factor(x))
}

# Ensure outcomes are numeric 0/1
for (y in outcomes) {
  stopifnot(y %in% colnames(pheno))
  pheno[[y]] <- as_num01(pheno[[y]])
}

# Covariates
pheno$surg_age <- suppressWarnings(as.numeric(pheno$surg_age))
pheno$surg_pgleasnt <- suppressWarnings(as.numeric(pheno$surg_pgleasnt))
pheno$surg_lnmets <- as_num01(pheno$surg_lnmets)

pheno$clinicalt <- prep_stage(pheno$clinicalt)
pheno$clinicalt_collapsed <- as.character(pheno$clinicalt)
pheno$clinicalt_collapsed[pheno$clinicalt_collapsed %in% c("2b","2c")] <- "2b2c"
pheno$clinicalt_collapsed <- factor(pheno$clinicalt_collapsed, levels=c("1c","2a","2b2c"))
table(pheno$clinicalt_collapsed)


pheno$recurrence <- factor(pheno$recurrence, levels=c(0,1))
pheno$Bad_outcome <- factor(pheno$Bad_outcome, levels=c(0,1))

beta_to_m <- function(b, offset=1e-6) {
  b <- pmin(pmax(b, offset), 1 - offset)
  log2(b/(1-b))
}
dmr_m <- beta_to_m(dmr_mat)

# Function to fit DMR models
fit_dmr_limma <- function(yvar, include_lnmets = TRUE) {
  covars <- c("surg_age", 
              "surg_pgleasnt",
              "clinicalt_collapsed"
              )
  if (include_lnmets) covars <- c(covars, "surg_lnmets")
  
  # keep only covars that exist and are not all NA
  covars <- covars[covars %in% colnames(pheno)]
  covars <- covars[sapply(covars, function(v) !all(is.na(pheno[[v]])))]
  
  fml <- if (length(covars) == 0) {
    as.formula(paste("~", yvar))
  } else {
    as.formula(paste("~", yvar, "+", paste(covars, collapse = " + ")))
  }
  design <- model.matrix(fml, data = pheno)
  
  # Use complete cases for the design to avoid limma dropping rows inconsistently
  cc <- complete.cases(design)
  design_cc <- design[cc, , drop = FALSE]
  y_cc <- pheno[[yvar]][cc]
  if (is.factor(y_cc) && nlevels(droplevels(y_cc)) < 2) {
    stop("After complete-case filtering, ", yvar, " has <2 levels.")
  }
  
  # Subset methylation matrix to same samples (columns already match pheno order)
  #dmr_mat_cc <- dmr_mat[, cc, drop = FALSE]
  dmr_mat_cc <- dmr_m[, cc, drop=FALSE]
  
  # limma expects features x samples
  fit <- lmFit(dmr_mat_cc, design_cc)
  fit <- eBayes(fit)
  
  # Identify the coefficient corresponding to yvar (factor or numeric)
  #coef_candidates <- colnames(design_cc)[grepl(paste0("^", yvar), colnames(design_cc))]
  coef_candidates <- colnames(design_cc)[colnames(design_cc) == yvar | startsWith(colnames(design_cc), paste0(yvar))]
  # If y is numeric 0/1, it will usually be exactly yvar
  # If y is factor, it will be yvar<level> (e.g., recurrence1)
  if (yvar %in% colnames(design_cc)) {
    coef_name <- yvar
  } else if (length(coef_candidates) == 1) {
    coef_name <- coef_candidates
  } else if (length(coef_candidates) > 1) {
    # pick the non-baseline level if possible (often ends with "1")
    # otherwise just take the last (usually the "highest" level)
    coef_name <- tail(coef_candidates, 1)
  } else {
    stop("Could not find coefficient for yvar='", yvar,
         "'. Design columns are: ", paste(colnames(design_cc), collapse=", "))
  }
  message("Outcome: ", yvar)
  message("Design cols: ", paste(colnames(design_cc), collapse = ", "))
  message("Using coefficient: ", coef_name)
  
  # Pull the coefficient for yvar
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$dmr_id <- rownames(tt)
  
  # Keep key fields
  out <- data.frame(
    dmr_id = tt$dmr_id,
    beta = tt$logFC,                 
    t = tt$t,
    p = tt$P.Value,
    fdr = tt$adj.P.Val,
    outcome = yvar,
    include_lnmets = include_lnmets,
    stringsAsFactors = FALSE
  )
  out
}

# --- Primary: recurrence WITH LN mets adjustment
res_recur <- fit_dmr_limma("recurrence", include_lnmets = FALSE)

res_recur %>%
  dplyr::summarize(
    n_fdr05 = sum(fdr <= 0.05, na.rm=TRUE),
    n_fdr05_dbeta10 = sum(fdr <= 0.05 & abs(beta) >= 0.10, na.rm=TRUE),
    n_fdr05_dbeta05 = sum(fdr <= 0.05 & abs(beta) >= 0.05, na.rm=TRUE)
  ) %>% print()


# --- Secondary: Bad_outcome WITHOUT LN mets (to avoid adjusting for part of endpoint)
res_bad   <- fit_dmr_limma("Bad_outcome", include_lnmets = FALSE)

res_bad %>%
  dplyr::summarize(
    n_fdr05 = sum(fdr <= 0.05, na.rm=TRUE),
    n_fdr05_dbeta10 = sum(fdr <= 0.05 & abs(beta) >= 0.10, na.rm=TRUE),
    n_fdr05_dbeta05 = sum(fdr <= 0.05 & abs(beta) >= 0.05, na.rm=TRUE)
  ) %>% print()


# Save
fwrite(res_recur, "./outs/GSE308050_DMR_results_recurrence_adjLN.csv")
fwrite(res_bad,   "./outs/GSE308050_DMR_results_BadOutcome_noLN.csv")

# ============================================================
# Promoter mapping + ORA + concordance
#   Run for a chosen outcome result table (res)
# ============================================================

canonical_chr <- function(chr) chr %in% paste0("chr", c(1:22,"X","Y","M"))
 
# parse_dmr_coords <- function(df) {
#   df %>%
#     mutate(
#       chr = str_extract(dmr_id, "^chr[0-9XYM]+"),
#       start = as.integer(str_match(dmr_id, "^chr[0-9XYM]+_(\\d+)_")[,2]),
#       end   = as.integer(str_match(dmr_id, "^chr[0-9XYM]+_\\d+_(\\d+)$")[,2])
#     ) %>%
#     filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
#     mutate(start = pmin(start, end), end = pmax(start, end))
# }

parse_coords <- function(res_df){
  res_df %>%
    mutate(
      chr = str_extract(dmr_id, "^chr[0-9XYM]+"),
      start = as.integer(str_match(dmr_id, "^chr[0-9XYM]+_(\\d+)_")[,2]),
      end   = as.integer(str_match(dmr_id, "^chr[0-9XYM]+_\\d+_(\\d+)$")[,2])
    ) %>%
    filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
    mutate(
      start2 = pmin(start, end),
      end2   = pmax(start, end),
      start = start2, end = end2
    ) %>%
    dplyr::select(-start2, -end2) %>%
    filter(canonical_chr(chr))
}


# promoter mapping helper: returns map table (dmr_id -> SYMBOL)
# map_dmrs_to_promoters <- function(coords_df, upstream=2000, downstream=2000) {
#   if (nrow(coords_df) == 0) {
#     return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character()))
#   }
#   
#   # Build GRanges for DMRs
#   coords_df <- coords_df %>% dplyr::filter(chr %in% paste0("chr", c(1:22,"X","Y","M")))
#   gr <- GRanges(seqnames = coords_df$chr, ranges = IRanges(start=coords_df$start, end=coords_df$end))
#   
#   # Drop non-standard chromosomes (gets rid of *_alt etc.)
#   gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
#   
#   # Trim out-of-bound (defensive)
#   gr <- trim(gr)
#   
#   if (length(gr) == 0) {
#     return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character()))
#   }
#   
#   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#   prom <- promoters(genes(txdb), upstream = upstream, downstream = downstream)
#   
#   hits <- findOverlaps(gr, prom, ignore.strand = TRUE)
#   if (length(hits) == 0) {
#     return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character()))
#   }
#   
#   q_idx <- queryHits(hits)
#   s_idx <- subjectHits(hits)
#   
#   entrez <- names(prom)[s_idx]
#   entrez <- entrez[!is.na(entrez) & entrez != ""]
#   
#   if (length(entrez) == 0) {
#     return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character()))
#   }
#   
#   sym <- mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
#                 column = "SYMBOL", multiVals = "first")
#   
#   out <- tibble(
#     dmr_id = coords_df$dmr_id[q_idx],
#     beta   = coords_df$beta[q_idx],
#     fdr    = coords_df$fdr[q_idx],
#     ENTREZID = names(sym),
#     SYMBOL = toupper(as.character(sym))
#   ) %>%
#     filter(!is.na(SYMBOL), SYMBOL != "")
#   
#   out
# }

canonical_levels <- paste0("chr", c(1:22,"X","Y","M"))

keep_canonical_seq <- function(gr) {
  gr[as.character(seqnames(gr)) %in% canonical_levels]
}


map_to_promoter <- function(dmr_df, upstream=2000, downstream=2000){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  prom <- promoters(genes(txdb), upstream=upstream, downstream=downstream)
  prom <- keep_canonical_seq(prom)
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keep_canonical_seq(gr)
  
  hits <- findOverlaps(gr, prom, ignore.strand=FALSE)
  if (length(hits) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character()))
  }
  
  entrez <- names(prom)[subjectHits(hits)]
  sym <- mapIds(org.Hs.eg.db, keys=entrez, keytype="ENTREZID", column="SYMBOL", multiVals="first")
  
  tibble(
    dmr_id = dmr_df$dmr_id[queryHits(hits)],
    beta   = dmr_df$beta[queryHits(hits)],
    fdr    = dmr_df$fdr[queryHits(hits)],
    ENTREZID = entrez,
    SYMBOL = toupper(as.character(sym))
  ) %>% filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(SYMBOL, dmr_id, .keep_all=TRUE)
}


# Gene-level Nearest TSS mapping (distance-based)
map_to_nearest_tss <- function(dmr_df, max_dist=100000){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene_gr <- genes(txdb)
  gene_gr <- keep_canonical_seq(gene_gr)
  tss <- resize(gene_gr, width=1, fix="start")
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keep_canonical_seq(gr)
  
  nearest <- distanceToNearest(gr, tss, ignore.strand=TRUE)
  if (length(nearest) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character(), dist=integer()))
  }
  
  q <- queryHits(nearest)
  s <- subjectHits(nearest)
  dist <- mcols(nearest)$distance
  
  keep <- dist <= max_dist
  q <- q[keep]; s <- s[keep]; dist <- dist[keep]
  if (length(q) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(), ENTREZID=character(), SYMBOL=character(), dist=integer()))
  }
  
  entrez <- names(tss)[s]
  sym <- mapIds(org.Hs.eg.db, keys=entrez, keytype="ENTREZID", column="SYMBOL", multiVals="first")
  
  tibble(
    dmr_id = dmr_df$dmr_id[q],
    beta   = dmr_df$beta[q],
    fdr    = dmr_df$fdr[q],
    ENTREZID = entrez,
    SYMBOL = toupper(as.character(sym)),
    dist = as.integer(dist)
  ) %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(SYMBOL, dmr_id, .keep_all=TRUE)
}

# Transcript-level TSS mapping
map_to_nearest_tx_tss <- function(dmr_df, max_dist=100000){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  tx_gr <- transcripts(txdb, columns=c("tx_id","gene_id"))
  tx_gr <- keepStandardChromosomes(tx_gr, pruning.mode="coarse")
  tx_tss <- resize(tx_gr, width=1, fix="start")   # strand-aware
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  
  # ---- critical: align seqlevels BEFORE setting seqinfo
  common <- intersect(seqlevels(gr), seqlevels(tx_tss))
  gr <- keepSeqlevels(gr, common, pruning.mode="coarse")
  tx_tss2 <- keepSeqlevels(tx_tss, common, pruning.mode="coarse")
  
  # now safe to propagate seqinfo and trim
  seqinfo(gr) <- seqinfo(tx_tss2)
  gr <- trim(gr)
  
  nearest <- distanceToNearest(gr, tx_tss2, ignore.strand=FALSE)
  if (length(nearest) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(),
                  ENTREZID=character(), SYMBOL=character(), dist=integer()))
  }
  
  q <- queryHits(nearest)
  s <- subjectHits(nearest)
  dist <- mcols(nearest)$distance
  
  keep <- dist <= max_dist
  q <- q[keep]; s <- s[keep]; dist <- dist[keep]
  if (length(q) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(),
                  ENTREZID=character(), SYMBOL=character(), dist=integer()))
  }
  
  entrez <- as.character(mcols(tx_tss2)$gene_id[s])
  sym <- mapIds(org.Hs.eg.db, keys=entrez, keytype="ENTREZID",
                column="SYMBOL", multiVals="first")
  
  tibble(
    dmr_id = dmr_df$dmr_id[q],
    beta   = dmr_df$beta[q],
    fdr    = dmr_df$fdr[q],
    ENTREZID = entrez,
    SYMBOL = toupper(as.character(sym)),
    dist = as.integer(dist)
  ) %>% filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(SYMBOL, dmr_id, .keep_all=TRUE)
}


# Promoter mapping around transcript TSS (±2000)
map_to_tx_promoter <- function(dmr_df, upstream=2000, downstream=2000){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  tx_gr <- transcripts(txdb, columns=c("tx_id","gene_id"))
  tx_gr <- keepStandardChromosomes(tx_gr, pruning.mode="coarse")
  tx_prom <- promoters(tx_gr, upstream=upstream, downstream=downstream)
  tx_prom <- keepStandardChromosomes(tx_prom, pruning.mode="coarse")
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  
  common <- intersect(seqlevels(gr), seqlevels(tx_prom))
  gr <- keepSeqlevels(gr, common, pruning.mode="coarse")
  tx_prom2 <- keepSeqlevels(tx_prom, common, pruning.mode="coarse")
  
  seqinfo(gr) <- seqinfo(tx_prom2)
  gr <- trim(gr)
  
  hits <- findOverlaps(gr, tx_prom2, ignore.strand=FALSE)
  if (length(hits) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), fdr=numeric(),
                  ENTREZID=character(), SYMBOL=character()))
  }
  
  entrez <- as.character(mcols(tx_prom2)$gene_id[subjectHits(hits)])
  sym <- mapIds(org.Hs.eg.db, keys=entrez, keytype="ENTREZID",
                column="SYMBOL", multiVals="first")
  
  tibble(
    dmr_id = dmr_df$dmr_id[queryHits(hits)],
    beta   = dmr_df$beta[queryHits(hits)],
    t      = dmr_df$t[queryHits(hits)],
    fdr    = dmr_df$fdr[queryHits(hits)],
    ENTREZID = entrez,
    SYMBOL = toupper(as.character(sym))
  ) %>% filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(SYMBOL, dmr_id, .keep_all=TRUE)
}


dmr_to_nearest_tx_tss_dist <- function(dmr_df){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  tx_gr <- transcripts(txdb, columns=c("tx_id","gene_id"))
  tx_gr <- keepStandardChromosomes(tx_gr, pruning.mode="coarse")
  tx_tss <- resize(tx_gr, width=1, fix="start")
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  
  common <- intersect(seqlevels(gr), seqlevels(tx_tss))
  gr <- keepSeqlevels(gr, common, pruning.mode="coarse")
  tx_tss2 <- keepSeqlevels(tx_tss, common, pruning.mode="coarse")
  
  seqinfo(gr) <- seqinfo(tx_tss2)
  gr <- trim(gr)
  
  nearest <- distanceToNearest(gr, tx_tss2, ignore.strand=FALSE)
  if (length(nearest) == 0) return(tibble(dmr_id=character(), dist_to_tx_tss=integer()))
  
  tibble(
    dmr_id = dmr_df$dmr_id[queryHits(nearest)],
    dist_to_tx_tss = as.integer(mcols(nearest)$distance)
  )
}


map_to_gene_body <- function(dmr_df){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene_gr <- genes(txdb)
  gene_gr <- keepStandardChromosomes(gene_gr, pruning.mode="coarse")
  
  gr <- GRanges(dmr_df$chr, IRanges(dmr_df$start, dmr_df$end))
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  
  hits <- findOverlaps(gr, gene_gr, ignore.strand=TRUE)
  if (length(hits) == 0) {
    return(tibble(dmr_id=character(), beta=numeric(), t=numeric(), fdr=numeric(),
                  ENTREZID=character(), SYMBOL=character()))
  }
  
  entrez <- names(gene_gr)[subjectHits(hits)]
  sym <- mapIds(org.Hs.eg.db, keys=entrez, keytype="ENTREZID", column="SYMBOL", multiVals="first")
  
  tibble(
    dmr_id = dmr_df$dmr_id[queryHits(hits)],
    beta   = dmr_df$beta[queryHits(hits)],
    t      = dmr_df$t[queryHits(hits)],
    fdr    = dmr_df$fdr[queryHits(hits)],
    ENTREZID = entrez,
    SYMBOL = toupper(as.character(sym))
  ) %>% filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(SYMBOL, dmr_id, .keep_all=TRUE)
}

# gene-level collapse: strongest abs(beta)
collapse_to_gene <- function(map_tbl) {
  if (nrow(map_tbl) == 0) {
    return(tibble(SYMBOL=character(), n_dmr=integer(), beta=numeric(), meth_dir=character()))
  }
  map_tbl %>%
    group_by(SYMBOL) %>%
    dplyr::summarize(
      n_dmr = n(),
      #beta = beta[which.max(abs(beta))],
      beta_med = median(beta, na.rm=TRUE),
      beta = beta_med,
      meth_dir = ifelse(beta > 0, "hyper_in_aggressive", "hypo_in_aggressive"),
      .groups = "drop"
    )
}


fisher_ora <- function(target_set, prog_genes, bg_genes, alternative="greater"){
  target_set <- intersect(toupper(target_set), bg_genes)
  prog_genes <- intersect(toupper(prog_genes), bg_genes)
  
  a <- length(intersect(target_set, prog_genes))
  b <- length(setdiff(target_set, prog_genes))
  c <- length(setdiff(prog_genes, target_set))
  d <- length(setdiff(bg_genes, union(target_set, prog_genes)))
  
  mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  ft <- fisher.test(mat, alternative=alternative)
  
  # finite OR with 0.5 correction (report alongside Fisher p)
  or_haldane <- ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))
  
  tibble(a=a,b=b,c=c,d=d,
         OR=unname(ft$estimate), OR_haldane=or_haldane,
         p=ft$p.value, n_target=length(target_set), n_prog=length(prog_genes), n_bg=length(bg_genes))
}


concordance_met <- function(gene_tbl, met_pos, met_neg){
  met_pos <- toupper(met_pos); met_neg <- toupper(met_neg)
  gene_tbl %>%
    mutate(
      in_pos = SYMBOL %in% met_pos,
      in_neg = SYMBOL %in% met_neg
    ) %>%
    filter(in_pos | in_neg) %>%
    mutate(
      expr_dir = ifelse(in_pos, "up_in_mets", "down_in_mets"),
      expected_meth_dir = ifelse(expr_dir=="up_in_mets", "hypo_in_aggressive", "hyper_in_aggressive"),
      concordant = (meth_dir == expected_meth_dir)
    )
}


run_met_priming <- function(res_df, label, met_pos, met_neg, FDR_CUT=0.05, DBETA_CUT=0.10,
                            promoter_up=10000, promoter_down=10000, tss_dist=100000){
  
  coords_all <- parse_coords(res_df)
  sig <- coords_all %>% filter(fdr <= FDR_CUT, abs(beta) >= DBETA_CUT)
  message(label, " | sig DMRs = ", nrow(sig))
  
  ## ---- PROMOTER mapping
  map_all_prom <- map_to_tx_promoter(coords_all, upstream=promoter_up, downstream=promoter_down)
  map_sig_prom <- map_to_tx_promoter(sig,        upstream=promoter_up, downstream=promoter_down)
  
  gene_sig_prom <- collapse_to_gene(map_sig_prom)
  
  bg_prom   <- sort(unique(map_all_prom$SYMBOL))          
  prog_prom <- sort(unique(gene_sig_prom$SYMBOL))         # progression-linked genes
  
  ora_prom <- bind_rows(
    fisher_ora(met_pos, prog_prom, bg_prom) %>% mutate(set="Met-Score POS", mapping=paste0("prom_",promoter_up,"_",promoter_down)),
    fisher_ora(met_neg, prog_prom, bg_prom) %>% mutate(set="Met-Score NEG", mapping=paste0("prom_",promoter_up,"_",promoter_down))
  )
  
  met_prom <- concordance_met(gene_sig_prom, met_pos, met_neg)  # concordance only for prog-linked met genes
  concord_prom <- if (nrow(met_prom)>0) {
    k <- sum(met_prom$concordant); n <- nrow(met_prom)
    tibble(mapping=paste0("prom_",promoter_up,"_",promoter_down), k=k, n=n,
           binom_p=binom.test(k,n,p=0.5,alternative="greater")$p.value)
  } else tibble(mapping=paste0("prom_",promoter_up,"_",promoter_down), k=0, n=0, binom_p=NA_real_)
  
  ## ---- NEAREST TSS mapping
  map_all_tss <- map_to_nearest_tx_tss(coords_all, max_dist=tss_dist)
  map_sig_tss <- map_to_nearest_tx_tss(sig,        max_dist=tss_dist)
  
  gene_sig_tss <- collapse_to_gene(map_sig_tss)
  
  bg_tss   <- sort(unique(map_all_tss$SYMBOL))            
  prog_tss <- sort(unique(gene_sig_tss$SYMBOL))
  
  ora_tss <- bind_rows(
    fisher_ora(met_pos, prog_tss, bg_tss) %>% mutate(set="Met-Score POS", mapping=paste0("tss_",tss_dist)),
    fisher_ora(met_neg, prog_tss, bg_tss) %>% mutate(set="Met-Score NEG", mapping=paste0("tss_",tss_dist))
  )
  
  met_tss <- concordance_met(gene_sig_tss, met_pos, met_neg)
  concord_tss <- if (nrow(met_tss)>0) {
    k <- sum(met_tss$concordant); n <- nrow(met_tss)
    tibble(mapping=paste0("tss_",tss_dist), k=k, n=n,
           binom_p=binom.test(k,n,p=0.5,alternative="greater")$p.value)
  } else tibble(mapping=paste0("tss_",tss_dist), k=0, n=0, binom_p=NA_real_)
  
  list(
    sig = sig,
    promoter = list(map_all=map_all_prom, map_sig=map_sig_prom, gene=gene_sig_prom, ora=ora_prom, met=met_prom, concord=concord_prom,
                    bg_n=length(bg_prom), prog_n=length(prog_prom)),
    tss = list(map_all=map_all_tss, map_sig=map_sig_tss, gene=gene_sig_tss, ora=ora_tss, met=met_tss, concord=concord_tss,
               bg_n=length(bg_tss), prog_n=length(prog_tss))
  )
}


run_mapping_suite <- function(res_df, label, met_pos, met_neg,
                              FDR_CUT=0.05, DBETA_CUT=0.10,
                              prom_up=2000, prom_down=2000, tss_dist=100000){
  
  coords_all <- parse_coords(res_df)
  sig <- coords_all %>% filter(fdr <= FDR_CUT, abs(beta) >= DBETA_CUT)
  message(label, " | sig DMRs: ", nrow(sig))
  
  # --- Promoter
  map_all_prom <- map_to_tx_promoter(coords_all, upstream=prom_up, downstream=prom_down)
  map_sig_prom <- map_to_tx_promoter(sig, upstream=prom_up, downstream=prom_down)
  
  # --- Gene body
  map_all_gb <- map_to_gene_body(coords_all)
  map_sig_gb <- map_to_gene_body(sig)
  
  # --- Nearest TSS
  map_all_tss <- map_to_nearest_tx_tss(coords_all, max_dist=tss_dist)
  map_sig_tss <- map_to_nearest_tx_tss(sig, max_dist=tss_dist)
  
  # Collapse
  gene_sig_prom <- collapse_to_gene(map_sig_prom)
  gene_sig_gb   <- collapse_to_gene(map_sig_gb)
  gene_sig_tss  <- collapse_to_gene(map_sig_tss)
  
  mk_ora <- function(map_all, gene_sig, mapping_name){
    bg <- sort(unique(map_all$SYMBOL))
    prog <- sort(unique(gene_sig$SYMBOL))
    bind_rows(
      fisher_ora(met_pos, prog, bg) %>% mutate(set="Met-Score POS"),
      fisher_ora(met_neg, prog, bg) %>% mutate(set="Met-Score NEG")
    ) %>% mutate(mapping=mapping_name, outcome=label, bg_n=length(bg), prog_n=length(prog))
  }
  
  ora_df <- bind_rows(
    mk_ora(map_all_prom, gene_sig_prom, paste0("prom_",prom_up,"_",prom_down)),
    mk_ora(map_all_gb,   gene_sig_gb,   "gene_body"),
    mk_ora(map_all_tss,  gene_sig_tss,  paste0("nearest_tss_",tss_dist))
  )
  
  list(sig=sig,
       ora=ora_df,
       gene_sig=list(prom=gene_sig_prom, gene_body=gene_sig_gb, tss=gene_sig_tss))
}


report_met_promoter_priming <- function(res_df, met_pos, met_neg,
                                        FDR_CUT=0.15, DBETA_CUT=0.05,
                                        prom_up=2000, prom_down=2000){
  
  coords_all <- parse_coords(res_df)
  sig <- coords_all %>% filter(fdr <= FDR_CUT, abs(beta) >= DBETA_CUT)
  
  # transcript-promoter mapping
  map_all <- map_to_tx_promoter(coords_all, upstream=prom_up, downstream=prom_down)
  map_sig <- map_to_tx_promoter(sig,        upstream=prom_up, downstream=prom_down)
  
  # attach distance-to-nearest transcript TSS for each DMR (for interpretability)
  dist_all <- dmr_to_nearest_tx_tss_dist(coords_all)
  dist_sig <- dmr_to_nearest_tx_tss_dist(sig)
  
  map_sig2 <- map_sig %>%
    left_join(dist_sig, by="dmr_id") %>%
    mutate(
      meth_dir = ifelse(beta > 0, "hyper_in_outcome1", "hypo_in_outcome1"),
      in_pos = SYMBOL %in% toupper(met_pos),
      in_neg = SYMBOL %in% toupper(met_neg),
      set = case_when(in_pos ~ "MetScore_POS", in_neg ~ "MetScore_NEG", TRUE ~ "other")
    ) %>%
    filter(in_pos | in_neg)
  
  map_sig2 <- map_sig2 %>%
    distinct(set, SYMBOL, dmr_id, beta, fdr, dist_to_tx_tss, meth_dir)
  
  # eligibility = present in promoter background
  bg <- unique(map_all$SYMBOL)
  met_pos_bg <- intersect(toupper(met_pos), bg)
  met_neg_bg <- intersect(toupper(met_neg), bg)
  
  # which met genes are "hit" by sig promoter DMRs
  hit_pos <- sort(unique(map_sig2$SYMBOL[map_sig2$set=="MetScore_POS"]))
  hit_neg <- sort(unique(map_sig2$SYMBOL[map_sig2$set=="MetScore_NEG"]))
  
  # compact summaries
  summary_tbl <- tibble(
    set = c("MetScore_POS","MetScore_NEG"),
    n_in_background = c(length(met_pos_bg), length(met_neg_bg)),
    n_hit = c(length(hit_pos), length(hit_neg)),
    hit_genes = c(paste(hit_pos, collapse=";"), paste(hit_neg, collapse=";"))
  )
  
  # detailed per-DMR view for met genes
  detail_tbl <- map_sig2 %>%
    dplyr::select(set, SYMBOL, dmr_id, beta, fdr, dist_to_tx_tss, meth_dir) %>%
    arrange(set, SYMBOL, fdr)
  
  list(
    n_sig_dmrs = nrow(sig),
    background_genes_n = length(bg),
    summary = summary_tbl,
    details = detail_tbl
  )
}

# ----------------------------
# Quantify what fraction of significant DMRs are promoter-proximal under different windows?
# ----------------------------
# diag_bad  <- diagnose_promoter_overlap(res_bad,  FDR_CUT=0.05, DBETA_CUT=0.10)
# diag_recur<- diagnose_promoter_overlap(res_recur, FDR_CUT=0.05, DBETA_CUT=0.10)
# 
# print(diag_bad)
# print(diag_recur)
########
## report promoter priming with transcript-promoter mapping + hit genes + β + distance to TSS
rep <- report_met_promoter_priming(
  res_bad,
  met_pos=met_pos, met_neg=met_neg,
  FDR_CUT=0.25, DBETA_CUT=0.05,
  prom_up=5000, prom_down=5000
)

rep$summary
rep$details

rep$details %>%
  distinct(SYMBOL, set, beta, meth_dir) %>%
  mutate(expected = ifelse(set=="MetScore_POS", "hypo_in_outcome1", "hyper_in_outcome1"),
         concordant = (meth_dir == expected))



# how many unique genes have ANY promoter overlap
length(unique(map_to_tx_promoter(parse_coords(res_bad), 5000, 5000)$SYMBOL))

# how many Met genes are theoretically mappable
intersect(met_pos, map_to_tx_promoter(parse_coords(res_bad), 5000, 5000)$SYMBOL)
intersect(met_neg, map_to_tx_promoter(parse_coords(res_bad), 5000, 5000)$SYMBOL)


# ----------------------------
# Run primary + secondary analyses
# ----------------------------
# ana_recur <- run_promoter_priming_analysis(
#   res_recur,
#   outcome_label = "recurrence_adjLN",
#   upstream=2000, downstream=2000,
#   FDR_CUT=0.15, DBETA_CUT=0.10
# )


ana_bad <- run_met_priming(
  res_bad, "Bad_outcome_noLN",
  met_pos=toupper(met_pos), met_neg=toupper(met_neg),
  FDR_CUT=0.25, DBETA_CUT=0.10,
  promoter_up=2000, promoter_down=2000,  
  tss_dist=100000
)

# quick readouts
bind_rows(ana_bad$promoter$ora, ana_bad$tss$ora) %>% print(n=Inf)
bind_rows(ana_bad$promoter$concord, ana_bad$tss$concord) %>% print(n=Inf)

# how many Met-Score genes mapped?
cat("PROM background genes:", ana_bad$promoter$bg_n, " progression genes:", ana_bad$promoter$prog_n, "\n")
cat("TSS  background genes:", ana_bad$tss$bg_n,      " progression genes:", ana_bad$tss$prog_n, "\n")
cat("Met genes with progression-linked promoter/nearTSS DMRs:",
    nrow(ana_bad$promoter$met), "/", nrow(ana_bad$tss$met), "\n")

met_gene <- ana_bad$promoter$met %>%
  group_by(SYMBOL) %>%
  dplyr::summarize(beta_gene = median(beta, na.rm=TRUE),
            meth_dir = ifelse(beta_gene > 0, "hyper_in_aggressive", "hypo_in_aggressive"),
            .groups="drop")

#####
ana_bad3 <- run_met_priming(res_bad, "Bad_outcome_noLN",
                            met_pos=toupper(met_pos), met_neg=toupper(met_neg),
                            FDR_CUT=0.15, DBETA_CUT=0.10,
                            promoter_up=2000, promoter_down=2000,
                            tss_dist=250000
)
# quick readouts
bind_rows(ana_bad3$promoter$ora, ana_bad3$tss$ora) %>% print(n=Inf)
bind_rows(ana_bad3$promoter$concord, ana_bad3$tss$concord) %>% print(n=Inf)

# how many Met-Score genes mapped?
cat("PROM background genes:", ana_bad3$promoter$bg_n, " progression genes:", ana_bad3$promoter$prog_n, "\n")
cat("TSS  background genes:", ana_bad3$tss$bg_n,      " progression genes:", ana_bad3$tss$prog_n, "\n")
cat("Met genes with progression-linked promoter/nearTSS DMRs:",
    nrow(ana_bad3$promoter$met), "/", nrow(ana_bad3$tss$met), "\n")

#####
ana_bad4 <- run_met_priming(res_bad, "Bad_outcome_noLN",
                            met_pos=toupper(met_pos), met_neg=toupper(met_neg),
                            FDR_CUT=0.25, DBETA_CUT=0.05,
                            promoter_up=2000, promoter_down=2000,
                            tss_dist=250000
)
# quick readouts
bind_rows(ana_bad4$promoter$ora, ana_bad4$tss$ora) %>% print(n=Inf)
bind_rows(ana_bad4$promoter$concord, ana_bad4$tss$concord) %>% print(n=Inf)

# how many Met-Score genes mapped?
cat("PROM background genes:", ana_bad4$promoter$bg_n, " progression genes:", ana_bad4$promoter$prog_n, "\n")
cat("TSS  background genes:", ana_bad4$tss$bg_n,      " progression genes:", ana_bad4$tss$prog_n, "\n")
cat("Met genes with progression-linked promoter/nearTSS DMRs:",
    nrow(ana_bad4$promoter$met), "/", nrow(ana_bad4$tss$met), "\n")



#####
##############
# BAD_OUTCOME
###############
suite_bad <- run_mapping_suite(
  res_df   = res_bad,
  label    = "Bad_outcome_noLN",
  met_pos  = met_pos,
  met_neg  = met_neg,
  FDR_CUT  = 0.05,
  DBETA_CUT= 0.10,
  prom_up  = 2000, prom_down = 2000,
  tss_dist = 100000
)

# More permissive sensitivity
suite_bad_sens <- run_mapping_suite(
  res_df   = res_bad,
  label    = "Bad_outcome_noLN_FDR0.15_DB0.05",
  met_pos  = met_pos,
  met_neg  = met_neg,
  FDR_CUT  = 0.15,
  DBETA_CUT= 0.05,
  prom_up  = 2000, prom_down = 2000,
  tss_dist = 250000
)

##############
# RECURRENCE
###############
suite_recur <- run_mapping_suite(
  res_df   = res_recur,
  label    = "recurrence_adjLN",
  met_pos  = met_pos,
  met_neg  = met_neg,
  FDR_CUT  = 0.15,    
  DBETA_CUT= 0.10,
  prom_up  = 2000, prom_down = 2000,
  tss_dist = 100000
)


# Inspect ORA tables
print(suite_bad$ora)
print(suite_bad_sens$ora)

# Save ORA
fwrite(as.data.table(suite_bad$ora), file.path(out_dir, "BadOutcome_ORA_mapping_suite.csv"))
fwrite(as.data.table(suite_bad_sens$ora), file.path(out_dir, "BadOutcome_ORA_mapping_suite_sensitivity.csv"))



###################################################
# Rank-based gene scoring + GSEA
###################################################
library(fgsea)

gene_score_from_map <- function(map_tbl){
  if (nrow(map_tbl) == 0) return(tibble(SYMBOL=character(), score=numeric()))
  map_tbl %>%
    group_by(SYMBOL) %>%
    dplyr::summarize(
      beta = beta[which.max(abs(beta))],
      fdr  = min(fdr, na.rm=TRUE),
      score = sign(beta) * (-log10(pmax(fdr, 1e-300))),
      .groups="drop"
    )
}

gene_rank_from_map_t <- function(map_tbl){
  map_tbl %>%
    group_by(SYMBOL) %>%
    dplyr::summarize(t_gene = t[which.max(abs(t))], .groups="drop")
}

# Bad_outcome using gene-body mapping:
coords_all <- parse_coords(res_bad)
map_all_gb <- map_to_gene_body(coords_all)
map_all_prom <- map_to_tx_promoter(coords_all, upstream=2000, downstream=2000)

gs_t <- gene_rank_from_map_t(map_all_gb) %>%
  filter(!is.na(SYMBOL), SYMBOL != "", is.finite(t_gene)) %>%
  distinct(SYMBOL, .keep_all=TRUE)

# gs_t <- map_all_prom %>%
#   group_by(SYMBOL) %>%
#   dplyr::summarize(t_gene = median(t, na.rm=TRUE), .groups="drop") %>%
#   filter(is.finite(t_gene), !is.na(SYMBOL), SYMBOL!="") %>%
#   distinct(SYMBOL, .keep_all=TRUE)


ranks <- gs_t$t_gene
names(ranks) <- gs_t$SYMBOL
ranks <- ranks[!duplicated(names(ranks))]

# break ties deterministically to avoid arbitrary ordering warnings
set.seed(7)
ranks_jitter <- ranks + rnorm(length(ranks), mean = 0, sd = 1e-8)
# gene sets
met_pos_u <- unique(toupper(met_pos))
met_neg_u <- unique(toupper(met_neg))

universe <- names(ranks_jitter)
cat("Universe size:", length(universe), "\n")
cat("Met POS overlap:", length(intersect(met_pos_u, universe)), "\n")
cat("Met NEG overlap:", length(intersect(met_neg_u, universe)), "\n")



pathways <- list(
  "MetScore_POS" = met_pos_u,
  "MetScore_NEG" = met_neg_u
)

# run fgsea
fg <- fgseaMultilevel(
  pathways = pathways,
  stats    = ranks_jitter,
  minSize  = 2,
  maxSize  = 5000,
  eps      = 1e-50
) %>%
  arrange(padj, desc(abs(NES)))

fg

fg_write <- fg %>%
  mutate(
    leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  )

fg_write

# Save
fg_out <- file.path(out_dir, "BadOutcome_GSEA_fgsea_geneBody_ranked.csv")
write.csv(fg_write, fg_out, row.names = FALSE)

# Also save the ranked list (useful for debugging / reproducibility)
rank_out <- file.path(out_dir, "BadOutcome_geneBody_ranked_scores.csv")
write.csv(gs_t, rank_out, row.names = FALSE)
cat("Wrote:", rank_out, "\n")

# quick bar plot
p <- ggplot(fg_write, aes(x = pathway, y = NES)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(
    title = "GSE308050 (Bad_outcome): GSEA on gene-body mapped methylation scores",
    x = "", y = "NES"
  ) +
  coord_flip()

ggsave(
  filename = file.path(out_dir, "BadOutcome_GSEA_fgseaMultilevel_geneBody_ranked_barplot.png"),
  plot = p, width = 6, height = 3.5, dpi = 300
)


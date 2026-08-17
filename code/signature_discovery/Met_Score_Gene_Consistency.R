################################################################################
# Gene-level directional consistency across the 6 discovery cohorts: per-cohort effect sizes for Table S2 
# Plus a leave-one-cohort-out to test whether any of the 45 signature genes change
# direction or lose significance when a single discovery cohort is excluded from the meta-analysis.
################################################################################
# Clean the working directory
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

# Named thresholds used throughout (avoids repeating bare 0.05/50 literals across ~6 places).
FDR_THRESH      <- 0.05  # gene-level significance cutoff (matches filterGenes() upstream)
HET_PVAL_THRESH <- 0.05  # Cochran's Q p-value cutoff for "significant heterogeneity"
I2_HIGH_THRESH  <- 50    # I-squared (%) cutoff for "high heterogeneity"

compute_I_squared <- function(Q, df) {
  raw <- (Q - df) / Q
  i2  <- ifelse(is.finite(raw), pmax(0, raw) * 100, 0)
  i2
}

## ---- Figure S5 panels a/b (heterogeneity + true LOCO reselection) ----------
## In this mode the script writes only outs/FigureS5/ and never touches the
## canonical Table S2 or any other output. Panel b performs genuine
## reselection: for each held-out cohort it reruns the meta-analysis and gene
## filter on the other five cohorts, unlike the fixed-45 loss-of-significance
## diagnostic in the default path below.
if ("--figure-s5-only" %in% commandArgs(trailingOnly = TRUE)) {
  suppressPackageStartupMessages({ library(MetaIntegrator); library(pROC) })
  grDevices::pdf(NULL)   # discard any stray MetaIntegrator plot; avoids an Rplots.pdf
  outdir <- "./outs/FigureS5"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  load("./outs/meta_analysis_results.rda")   # meta_analysis_results
  load("./outs/PP_filter_MetaScore.rda")     # filter
  pos45 <- filter$posGeneNames; neg45 <- filter$negGeneNames; genes45 <- c(pos45, neg45)
  stopifnot(length(genes45) == 45L)
  sig_cfg <- read.csv("config/metscore_signature_v1.csv", stringsAsFactors = FALSE)
  deployed41 <- sig_cfg$gene[as.logical(sig_cfg$used_in_locked_model)]
  stopifnot(length(deployed41) == 41L, all(deployed41 %in% genes45))

  # panel a: pooled Hedges g heterogeneity + 95% prediction interval
  # PI = g +/- t(0.975, k-2) * sqrt(tau^2 + SE^2); the 45-gene biological
  # signature is distinguished from the 41 deployed model features.
  pr <- as.data.frame(meta_analysis_results$metaAnalysis$pooledResults)[genes45, , drop = FALSE]
  k <- pr$numStudies
  i2 <- compute_I_squared(pr$cochranesQ, k - 1)
  pi_half <- qt(0.975, df = pmax(k - 2, 1)) * sqrt(pr$tauSquared + pr$effectSizeStandardError^2)
  panelA <- data.frame(
    gene = genes45, direction = ifelse(genes45 %in% pos45, "POS", "NEG"),
    used_in_locked_model = genes45 %in% deployed41,
    pooled_hedges_g = round(pr$effectSize, 7), se = round(pr$effectSizeStandardError, 7),
    pooled_p = signif(pr$effectSizePval, 7), pooled_fdr = signif(pr$effectSizeFDR, 7),
    tau2 = round(pr$tauSquared, 7), i2_pct = round(i2, 4),
    cochran_q = round(pr$cochranesQ, 6), q_p = signif(pr$heterogeneityPval, 7),
    n_cohorts = as.integer(k),
    pi_lo = round(pr$effectSize - pi_half, 7), pi_hi = round(pr$effectSize + pi_half, 7),
    pi_crosses_zero = (pr$effectSize - pi_half < 0) & (pr$effectSize + pi_half > 0),
    stringsAsFactors = FALSE)
  write.csv(panelA, file.path(outdir, "FigureS5_panelA_heterogeneity.csv"), row.names = FALSE)

  # panel b: true leave-one-cohort-out reselection on the other five cohorts
  od <- meta_analysis_results$originalData; cohorts <- names(od)
  memb <- matrix(0L, length(genes45), length(cohorts),
                 dimnames = list(genes45, paste0("omit_", cohorts)))
  fold <- list(); loso <- list()
  for (om in cohorts) {
    ma <- runMetaAnalysis(list(originalData = od[setdiff(cohorts, om)]),
                          runLeaveOneOutAnalysis = FALSE, maxCores = 1)
    ma <- filterGenes(ma, isLeaveOneOut = FALSE,
                      effectSizeThresh = 0.2, FDRThresh = 0.05, numberStudiesThresh = 4)
    fl <- ma$filterResults[[1]]; resel <- c(fl$posGeneNames, fl$negGeneNames)
    memb[intersect(resel, genes45), paste0("omit_", om)] <- 1L
    fold[[length(fold) + 1]] <- data.frame(
      omitted_cohort = om, n_reselected = length(resel),
      n_overlap_with_45 = length(intersect(resel, genes45)),
      n_added = length(setdiff(resel, genes45)), n_lost = length(setdiff(genes45, resel)),
      added_genes = paste(sort(setdiff(resel, genes45)), collapse = ";"),
      lost_genes = paste(sort(setdiff(genes45, resel)), collapse = ";"),
      stringsAsFactors = FALSE)
    # LOSO reselected-signature score: apply the fold-specific filter to the
    # complete omitted cohort with the canonical MetaIntegrator score, then
    # binary held-out AUC with a DeLong CI.
    stopifnot(length(resel) > 0L)
    test_ds <- od[[om]]; test_ds$keys <- rownames(test_ds$expr)
    y <- as.integer(test_ds$class)
    score <- as.numeric(calculateScore(filterObject = fl, datasetObject = test_ds))
    stopifnot(length(score) == length(y), all(is.finite(score)))
    r <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    loso[[length(loso) + 1]] <- data.frame(
      cohort = om, n = length(y), cases = sum(y == 1L), controls = sum(y == 0L),
      n_pos = length(fl$posGeneNames), n_neg = length(fl$negGeneNames), n_total = length(resel),
      auc = round(as.numeric(pROC::auc(r)), 7), ci_low = round(ci[1], 7), ci_high = round(ci[3], 7),
      stringsAsFactors = FALSE)
  }
  write.csv(do.call(rbind, fold), file.path(outdir, "FigureS5_panelB_fold_summary.csv"), row.names = FALSE)
  write.csv(data.frame(gene = genes45, direction = ifelse(genes45 %in% pos45, "POS", "NEG"),
                       memb, check.names = FALSE),
            file.path(outdir, "FigureS5_panelB_membership.csv"), row.names = FALSE)
  loso_df <- do.call(rbind, loso)
  write.csv(loso_df, file.path(outdir, "FigureS5_panelC_LOSO_signature_auc.csv"), row.names = FALSE)
  # compare against the historical GitHub LOSO result (checks, not fitting targets)
  hist_path <- "./outs/LOSO_internal_validation_AUC.csv"
  if (file.exists(hist_path)) {
    h <- read.csv(hist_path, stringsAsFactors = FALSE)
    cmp <- merge(loso_df[, c("cohort", "auc", "n_pos", "n_neg")], h[, c("heldout", "auc", "n_pos", "n_neg")],
                 by.x = "cohort", by.y = "heldout", suffixes = c("", "_hist"))
    cmp$auc_abs_diff <- abs(cmp$auc - cmp$auc_hist)
    cat("LOSO reselected-signature AUC vs historical GitHub:\n"); print(cmp, row.names = FALSE)
    cat(sprintf("max |AUC - historical| = %.3e\n", max(cmp$auc_abs_diff)))
  }
  cat("Figure S5 panels a/b/c-LOSO written to", outdir, "\n")
  grDevices::dev.off()
  quit(save = "no", status = 0)
}

dir.create("./outs", showWarnings = FALSE, recursive = TRUE)
dir.create("./figures", showWarnings = FALSE, recursive = TRUE)

# ── Load previous analysis  ──────────────────────────────────────────────────
load("./outs/filtersiggenes_MetaScore.rda") # Filter_SignatureGenes (45 genes)
load("./outs/PP_filter_MetaScore.rda")      # filter (posGeneNames/negGeneNames)

# Rename
sig_filter <- filter
rm(filter)

sig_genes <- Filter_SignatureGenes
cat(sprintf("Loaded saved signature: %d genes (%d up, %d down)\n",
            length(sig_genes), length(sig_filter$posGeneNames), length(sig_filter$negGeneNames)))

# ── Table S2 + Heterogeneity  ─────────────────────────────────────────────────

n_mismatches <- NA_integer_
het          <- NULL

if (file.exists("./outs/meta_analysis_results.rda")) {
  load("./outs/meta_analysis_results.rda")
  
  es_mat <- meta_analysis_results$metaAnalysis$datasetEffectSizes
  # Matrix indexing (es_mat) errors loudly on a missing gene, but data.frame indexing 
  #(pooled_df below) silently returns NA rows
  pooled_df <- as.data.frame(meta_analysis_results$metaAnalysis$pooledResults)
  missing_genes <- setdiff(sig_genes, rownames(pooled_df))
  if (length(missing_genes) > 0) {
    warning(sprintf(
      "%d signature gene(s) not found in meta_analysis_results$metaAnalysis$pooledResults (likely a stale/mismatched .rda pair): %s. These rows will be NA in Table S2.",
      length(missing_genes), paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  
  es_sig <- as.data.frame(es_mat[sig_genes, , drop = FALSE])
  es_sig$gene              <- rownames(es_sig)
  pooled                   <- pooled_df[sig_genes, ]
  es_sig$pooled_effect_size <- pooled$effectSize
  es_sig$pooled_FDR        <- pooled$effectSizeFDR
  es_sig$direction         <- ifelse(es_sig$gene %in% sig_filter$posGeneNames, "Up in Mets", "Down in Mets")
  
  cohort_cols   <- colnames(es_mat)
  expected_sign <- ifelse(es_sig$direction == "Up in Mets", 1, -1)
  sign_flags    <- sapply(cohort_cols, function(cn) {
    cohort_sign <- sign(es_sig[[cn]])
    ifelse(is.na(cohort_sign), NA, cohort_sign != expected_sign)
  })
  colnames(sign_flags) <- paste0(cohort_cols, "_sign_mismatch")
  
  TableS2 <- cbind(es_sig[, c("gene", "direction", "pooled_effect_size", "pooled_FDR", cohort_cols)],
                   sign_flags)
  
  n_mismatches <- sum(sign_flags, na.rm = TRUE)
  cat(sprintf("\nTable S2: %d genes x %d cohorts; %d sign mismatches.\n",
              nrow(TableS2), length(cohort_cols), n_mismatches))
  if (n_mismatches > 0) {
    mismatch_rows <- which(rowSums(sign_flags, na.rm = TRUE) > 0)
    cat("Genes with at least one cohort-level sign mismatch:\n")
    print(TableS2[mismatch_rows, c("gene", "direction", cohort_cols)], row.names = FALSE)
  }
  
  het           <- pooled[, c("cochranesQ", "heterogeneityPval", "tauSquared", "numStudies")]
  het$gene      <- sig_genes
  het$df        <- het$numStudies - 1
  het$I_squared <- compute_I_squared(het$cochranesQ, het$df)
  
  TableS2 <- merge(TableS2,
                   het[, c("gene", "cochranesQ", "heterogeneityPval", "tauSquared", "I_squared")],
                   by = "gene")
  TableS2 <- TableS2[order(TableS2$direction, -abs(TableS2$pooled_effect_size)), ]
  write.csv(TableS2, "./outs/TableS2_per_cohort_effect_sizes.csv", row.names = FALSE)
  
  cat(sprintf("Median I-squared: %.1f%% (IQR %.1f-%.1f%%)\n",
              median(het$I_squared), quantile(het$I_squared, 0.25), quantile(het$I_squared, 0.75)))
  cat(sprintf("Genes with significant heterogeneity (p<%.2f): %d / %d\n",
              HET_PVAL_THRESH, sum(het$heterogeneityPval < HET_PVAL_THRESH, na.rm = TRUE), nrow(het)))
  cat(sprintf("Genes with I-squared > %d%%: %d / %d\n",
              I2_HIGH_THRESH, sum(het$I_squared > I2_HIGH_THRESH, na.rm = TRUE), nrow(het)))
  write.csv(het[order(-het$I_squared),
                c("gene", "cochranesQ", "df", "heterogeneityPval", "tauSquared", "I_squared")],
            "./outs/Heterogeneity_by_gene.csv", row.names = FALSE)
  
} else {
  message("meta_analysis_results.rda not found — skipping Table S2 regeneration.\n",
          "Run MET_PCa_MetaIntegrator.R to create it. Loading existing CSVs if available.")
  if (file.exists("./outs/TableS2_per_cohort_effect_sizes.csv")) {
    TableS2 <- read.csv("./outs/TableS2_per_cohort_effect_sizes.csv", stringsAsFactors = FALSE)
    mismatch_cols <- grep("_sign_mismatch$", colnames(TableS2), value = TRUE)
    if (length(mismatch_cols) > 0)
      n_mismatches <- as.integer(sum(TableS2[, mismatch_cols], na.rm = TRUE))
    cat(sprintf("Loaded existing TableS2 (%d rows).\n", nrow(TableS2)))
  }
  if (file.exists("./outs/Heterogeneity_by_gene.csv")) {
    het <- read.csv("./outs/Heterogeneity_by_gene.csv", stringsAsFactors = FALSE)
    cat(sprintf("Loaded existing Heterogeneity_by_gene.csv (%d rows).\n", nrow(het)))
  }
}

# ── Leave-one-cohort-out gene stability  ──────────────────────────────────────

# MetaIntegrator already did the leave-one-out reruns automatically
# flag whether each one still agrees with the original finding (direction + significance), 
# and calculate the one heterogeneity statistic (I²)

compute_loo_summary <- function(ma_results, sig_genes, filter_obj) {
  loo       <- ma_results$leaveOneOutAnalysis
  loo_names <- names(loo)
  dplyr::bind_rows(lapply(loo_names, function(nm) {
    pr_loo <- as.data.frame(loo[[nm]]$pooledResults)
    missing_genes <- setdiff(sig_genes, rownames(pr_loo))
    if (length(missing_genes) > 0) {
      warning(sprintf(
        "%d signature gene(s) not found in LOO fold '%s' pooledResults: %s. These rows will be NA.",
        length(missing_genes), nm, paste(missing_genes, collapse = ", ")), call. = FALSE)
    }
    pr_loo <- pr_loo[sig_genes, , drop = FALSE]
    pr_loo$gene           <- sig_genes
    pr_loo$removed_cohort <- sub("^removed_", "", nm)
    pr_loo$direction_full <- ifelse(pr_loo$gene %in% filter_obj$posGeneNames,
                                    "Up in Mets", "Down in Mets")
    pr_loo$loo_sign_match <- sign(pr_loo$effectSize) ==
      ifelse(pr_loo$direction_full == "Up in Mets", 1, -1)
    pr_loo$loo_still_sig  <- pr_loo$effectSizeFDR <= FDR_THRESH
    gc <- function(df, nm) if (nm %in% colnames(df)) df[[nm]] else NA_real_
    q  <- gc(pr_loo, "cochranesQ"); n <- as.numeric(gc(pr_loo, "numStudies"))
    pr_loo$I_squared <- round(compute_I_squared(q, n - 1), 1)
    keep <- c("removed_cohort","gene","direction_full",
              "effectSize","effectSizeSD","effectSizePval","effectSizeFDR",
              "tauSquared","numStudies","cochranesQ","heterogeneityPval","I_squared",
              "loo_sign_match","loo_still_sig")
    pr_loo[, intersect(keep, colnames(pr_loo))]
  }))
}

# Tracks whether loo_summary was recomputed this run, so a stale cached LOO_by_cohort.csv isn't silently reused.
loo_summary_freshly_computed <- FALSE

if (file.exists("./outs/LOO_gene_consistency.csv")) {
  loo_summary <- read.csv("./outs/LOO_gene_consistency.csv", stringsAsFactors = FALSE)
  cat(sprintf("Loaded existing LOO_gene_consistency.csv (%d rows, %d cols).\n",
              nrow(loo_summary), ncol(loo_summary)))
  if (!"cochranesQ" %in% colnames(loo_summary) && exists("meta_analysis_results")) {
    cat("Extended LOO columns absent — recomputing from meta_analysis_results.\n")
    loo_summary <- compute_loo_summary(meta_analysis_results, sig_genes, sig_filter)
    loo_summary_freshly_computed <- TRUE
    write.csv(loo_summary, "./outs/LOO_gene_consistency.csv", row.names = FALSE)
  }
} else if (exists("meta_analysis_results")) {
  loo_summary <- compute_loo_summary(meta_analysis_results, sig_genes, sig_filter)
  loo_summary_freshly_computed <- TRUE
  write.csv(loo_summary, "./outs/LOO_gene_consistency.csv", row.names = FALSE)
} else {
  stop("Neither LOO_gene_consistency.csv nor meta_analysis_results.rda found.\n",
       "Run MET_PCa_MetaIntegrator.R first.")
}

n_sign_flips_loo <- sum(!loo_summary$loo_sign_match)
n_lost_sig_loo   <- sum(loo_summary$loo_sign_match & !loo_summary$loo_still_sig)
cat(sprintf("\nLOO gene consistency (%d checks): %d sign flips, %d pairs lose FDR<0.05.\n",
            nrow(loo_summary), n_sign_flips_loo, n_lost_sig_loo))
if (n_sign_flips_loo > 0) {
  cat("Genes/cohorts with a LOO sign flip:\n")
  print(loo_summary[!loo_summary$loo_sign_match,
                    c("removed_cohort", "gene", "direction_full", "effectSize")],
        row.names = FALSE)
}

# Reuse the cached cohort summary only when LOO results were not recomputed.
if (!loo_summary_freshly_computed && file.exists("./outs/LOO_by_cohort.csv")) {
  loo_by_cohort <- read.csv("./outs/LOO_by_cohort.csv", stringsAsFactors = FALSE)
  cat("Loaded existing LOO_by_cohort.csv.\n")
} else {
  loo_by_cohort <- loo_summary %>%
    dplyr::group_by(removed_cohort) %>%
    dplyr::summarise(
      n_sign_flips = sum(!loo_sign_match),
      n_lost_sig   = sum(loo_sign_match & !loo_still_sig),
      n_still_sig  = sum(loo_sign_match & loo_still_sig),
      pct_lost_sig = round(100 * n_lost_sig / dplyr::n(), 1),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_lost_sig))
  write.csv(loo_by_cohort, "./outs/LOO_by_cohort.csv", row.names = FALSE)
}
cat("\n=== Lost significance by removed cohort (out of 45 genes each) ===\n")
print(as.data.frame(loo_by_cohort))

p_loo <- ggplot(loo_by_cohort, aes(x = reorder(removed_cohort, n_lost_sig), y = n_lost_sig)) +
  geom_col(fill = "#07768f") +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n_lost_sig, pct_lost_sig)),
            hjust = -0.1, size = 4) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 45), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Cohort removed", y = "Signature genes losing FDR < 0.05 significance",
       title = "Leave-one-cohort-out: loss of gene-level significance") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(face = "bold"),
        plot.margin = margin(10, 40, 10, 10))
ggsave("./figures/LOO_lost_significance_by_cohort.pdf", p_loo, width = 7.5, height = 4.5)

# ──  Gene recovery table + cohort metadata  ───────────────────────────────────
# "Recovery" = genes from the full 45-gene panel still FDR<0.05 significant when that cohort is held out (direction never flips — see above).

gene_recovery_df <- loo_by_cohort %>%
  dplyr::mutate(
    n_full_panel  = length(sig_genes),
    n_recovered   = n_still_sig,
    pct_recovered = round(100 * n_recovered / n_full_panel, 1),
    jaccard_index = round(n_recovered / n_full_panel, 3)
  ) %>%
  dplyr::rename(held_out = removed_cohort) %>%
  dplyr::select(held_out, n_full_panel, n_recovered, pct_recovered,
                jaccard_index, n_sign_flips, n_lost_sig) %>%
  dplyr::arrange(desc(pct_recovered))

cohort_meta_gene <- data.frame(
  held_out = c("GSE116918", "GSE55935", "GSE51066",
               "GSE46691", "GSE41408", "GSE70769"),
  platform = c("Affymetrix HGU133Plus2", "Affymetrix HGU133A",
               "Affymetrix HGU133Plus2", "Affymetrix HGU133Plus2",
               "Affymetrix Human Gene 1.0 ST", "Affymetrix HGU133Plus2"),
  stringsAsFactors = FALSE
)

gene_recovery_df <- gene_recovery_df %>%
  dplyr::left_join(cohort_meta_gene, by = "held_out")

loco_thr <- tryCatch(
  read.csv("./outs/LOCO_threshold_stability.csv", stringsAsFactors = FALSE),
  error = function(e) { message("LOCO_threshold_stability.csv not found; skipping merge."); NULL }
)

if (!is.null(loco_thr)) {
  thr_cols <- c("held_out_cohort", "n_held_out", "n_mets_held", "threshold_fold", "lambda_fold")
  thr_cols_present <- intersect(thr_cols, colnames(loco_thr))
  loco_thr_sub <- loco_thr[, thr_cols_present, drop = FALSE]
  loco_thr_sub <- dplyr::rename(loco_thr_sub, held_out = held_out_cohort)
  if ("n_held_out" %in% colnames(loco_thr_sub) && "n_mets_held" %in% colnames(loco_thr_sub)) {
    loco_thr_sub$event_rate_pct <- round(100 * loco_thr_sub$n_mets_held / loco_thr_sub$n_held_out, 1)
  }
  gene_recovery_df <- gene_recovery_df %>%
    dplyr::left_join(loco_thr_sub, by = "held_out")
}

write.csv(gene_recovery_df, "./outs/LOCO_gene_stability_with_metadata.csv", row.names = FALSE)
cat("\nGene recovery + cohort metadata: ./outs/LOCO_gene_stability_with_metadata.csv\n")
print(gene_recovery_df[, c("held_out", "n_full_panel", "n_recovered", "pct_recovered",
                           "jaccard_index", "n_sign_flips", "platform")],
      row.names = FALSE)

# ──  GConsolidated Excel supplementary table: Table S2  ───────────────────────

if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  
  # Reads canonical Table S2 values (owned by generate_supplementary_tables.py; read-only here).
  # The LOO-augmented workbook is written to a distinct file (HET_S2_PATH below) so this
  # script never overwrites the canonical TableS2_metscore_genes.xlsx.
  CANONICAL_S2_PATH <- "./outs/TableS2_metscore_genes.xlsx"
  HET_S2_PATH       <- "./outs/TableS2_metscore_gene_heterogeneity.xlsx"
  
  load_canonical_s2 <- function(path = CANONICAL_S2_PATH) {
    if (!file.exists(path)) {
      stop(sprintf(
        "Canonical Table S2 file not found at '%s'. Run generate_supplementary_tables.py first, ",
        path),
        "or update CANONICAL_S2_PATH to point at its actual location.")
    }
    raw <- readxl::read_excel(path, sheet = "Met-Score signature genes")
    # Rename Python's manuscript-facing headers to the plain internal names 
    data.frame(
      gene       = raw$Gene,
      effectSize = raw$`Pooled effect size (Hedges' g)`,
      se         = raw$`Standard error`,
      pval       = raw$`Pooled p-value`,
      fdr        = raw$`Pooled FDR (BH)`,
      tau2       = raw$`Tau²`,
      numStudies = raw$`N studies contributing`,
      stringsAsFactors = FALSE
    )
  }
  
  canonical_s2 <- load_canonical_s2()
  cat(sprintf("Loaded canonical Table S2 values from %s (%d genes)\n",
              CANONICAL_S2_PATH, nrow(canonical_s2)))
  
  # Wraps the canonical_s2 gene lookup with an explicit warning,  
  # match() would otherwise silently return all-NA rows for any unmatched gene.
  lookup_canonical <- function(genes) {
    idx <- match(genes, canonical_s2$gene)
    if (anyNA(idx)) {
      warning(sprintf(
        "%d gene(s) not found in the canonical Table S2 file (%s): %s. Corresponding Excel rows will be NA -- canonical_s2 is likely out of date relative to the current signature.",
        sum(is.na(idx)), CANONICAL_S2_PATH, paste(genes[is.na(idx)], collapse = ", ")), call. = FALSE)
    }
    canonical_s2[idx, ]
  }
  
  # Build one Table S2-format data.frame from a pooledResults matrix. 
  # Falls back to canonical_s2 lookup if effectSizeSD/effectSizePval are absent
  build_s2_sheet <- function(pooled_mat, g_genes, filter_obj) {
    pr      <- as.data.frame(pooled_mat)[g_genes, , drop = FALSE]
    pr$gene <- g_genes
    
    get_col <- function(df, nm) if (nm %in% colnames(df)) df[[nm]] else NA_real_
    
    se_vec   <- get_col(pr, "effectSizeStandardError")
    pval_vec <- get_col(pr, "effectSizePval")
    
    if (all(is.na(se_vec)) || all(is.na(pval_vec))) {
      lkp <- lookup_canonical(g_genes)
      if (all(is.na(se_vec)))   se_vec   <- lkp$se
      if (all(is.na(pval_vec))) pval_vec <- lkp$pval
    }
    
    q_vals <- get_col(pr, "cochranesQ")
    n_vals <- as.numeric(get_col(pr, "numStudies"))
    i2     <- compute_I_squared(q_vals, n_vals - 1)
    
    df <- data.frame(
      Gene                             = pr$gene,
      Direction                        = ifelse(pr$gene %in% filter_obj$posGeneNames,
                                                "Up in metastasis", "Down in metastasis"),
      "Pooled effect size (Hedges' g)" = round(get_col(pr, "effectSize"), 6),
      "Standard error"                 = round(se_vec, 6),
      "Pooled p-value"                 = pval_vec,
      "Pooled FDR (BH)"               = get_col(pr, "effectSizeFDR"),
      "Tau²"                           = round(get_col(pr, "tauSquared"), 6),
      "N studies contributing"         = as.integer(get_col(pr, "numStudies")),
      "Q statistic (Cochran)"          = round(q_vals, 4),
      "Q p-value"                      = get_col(pr, "heterogeneityPval"),
      "I-squared (%)"                  = round(i2, 1),
      check.names = FALSE
    )
    df[order(df$Direction, -abs(df[["Pooled effect size (Hedges' g)"]])), ]
  }

  build_s2_from_csv <- function(tbl) {
    lkp <- lookup_canonical(tbl$gene)
    df <- data.frame(
      Gene                             = tbl$gene,
      Direction                        = tbl$direction,
      "Pooled effect size (Hedges' g)" = round(lkp$effectSize, 6),
      "Standard error"                 = round(lkp$se, 6),
      "Pooled p-value"                 = lkp$pval,
      "Pooled FDR (BH)"               = lkp$fdr,
      "Tau²"                           = round(lkp$tau2, 6),
      "N studies contributing"         = as.integer(lkp$numStudies),
      "Q statistic (Cochran)"          = round(tbl$cochranesQ, 4),
      "Q p-value"                      = tbl$heterogeneityPval,
      "I-squared (%)"                  = round(tbl$I_squared, 1),
      check.names = FALSE
    )
    df[order(df$Direction, -abs(df[["Pooled effect size (Hedges' g)"]])), ]
  }
  
  # Build LOO sheet from a pooledResults matrix.
  build_loo_sheet <- function(pooled_mat, g_genes, filter_obj) {
    pr     <- as.data.frame(pooled_mat)[g_genes, , drop = FALSE]
    gc_col <- function(df, nm) if (nm %in% colnames(df)) df[[nm]] else NA_real_
    df <- data.frame(
      Gene                             = g_genes,
      Direction                        = ifelse(g_genes %in% filter_obj$posGeneNames,
                                                "Up in metastasis", "Down in metastasis"),
      "Pooled effect size (Hedges' g)" = round(pr$effectSize, 6),
      "Standard error"                 = round(gc_col(pr, "effectSizeStandardError"), 6),
      "Pooled FDR (BH)"               = pr$effectSizeFDR,
      "N studies contributing"         = as.integer(gc_col(pr, "numStudies")),
      check.names = FALSE
    )
    df[order(df$Direction, -abs(df[["Pooled effect size (Hedges' g)"]])), ]
  }
  
  # Fallback for LOO sheets using LOO_gene_consistency.csv.
  build_loo_from_csv <- function(loo_df, cid) {
    sub <- loo_df[loo_df$removed_cohort == cid, ]
    data.frame(
      Gene                             = sub$gene,
      Direction                        = sub$direction_full,
      "Pooled effect size (Hedges' g)" = round(sub$effectSize, 6),
      "Standard error"                 = NA_real_,
      "Pooled FDR (BH)"               = sub$effectSizeFDR,
      "N studies contributing"         = 5L,
      check.names = FALSE
    )
  }
  
  hdr_style <- createStyle(fontSize = 11, fontColour = "white", fgFill = "#305496",
                           halign = "CENTER", textDecoration = "bold",
                           wrapText = TRUE, border = "TopBottomLeftRight",
                           borderColour = "#BFBFBF")
  num_style <- createStyle(fontSize = 10, halign = "RIGHT",
                           border = "TopBottomLeftRight", borderColour = "#BFBFBF")
  txt_style <- createStyle(fontSize = 10,
                           border = "TopBottomLeftRight", borderColour = "#BFBFBF")
  
  write_styled_sheet <- function(wb, sheet_name, df) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = hdr_style)
    n_rows   <- nrow(df)
    num_cols <- which(sapply(df, is.numeric))
    txt_cols <- setdiff(seq_len(ncol(df)), num_cols)
    if (length(num_cols) > 0)
      addStyle(wb, sheet_name, num_style,
               rows = 2:(n_rows + 1), cols = num_cols, gridExpand = TRUE)
    if (length(txt_cols) > 0)
      addStyle(wb, sheet_name, txt_style,
               rows = 2:(n_rows + 1), cols = txt_cols, gridExpand = TRUE)
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  }
  
  wb <- createWorkbook()
  
  # Sheet 1: full 6-cohort meta-analysis
  if (exists("meta_analysis_results")) {
    s1 <- build_s2_sheet(meta_analysis_results$metaAnalysis$pooledResults,
                         sig_genes, sig_filter)
  } else {
    tbl <- if (exists("TableS2")) TableS2 else
      read.csv("./outs/TableS2_per_cohort_effect_sizes.csv", stringsAsFactors = FALSE)
    s1 <- build_s2_from_csv(tbl)
  }
  write_styled_sheet(wb, "Met-Score signature genes", s1)
  
  # Sheets 2-7: one per held-out cohort, same 10-column structure as Sheet 1
  cohort_ids <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE70769")
  for (cid in cohort_ids) {
    loo_key <- paste0("removed_", cid)
    if (exists("meta_analysis_results") &&
        loo_key %in% names(meta_analysis_results$leaveOneOutAnalysis)) {
      sn <- build_loo_sheet(
        meta_analysis_results$leaveOneOutAnalysis[[loo_key]]$pooledResults,
        sig_genes, sig_filter)
    } else {
      sn <- build_loo_from_csv(loo_summary, cid)
    }
    write_styled_sheet(wb, paste0("LOO ", cid), sn)
  }
  
  saveWorkbook(wb, HET_S2_PATH, overwrite = TRUE)
  cat(sprintf("\nConsolidated Excel (7 sheets): %s\n", HET_S2_PATH))
} else {
  message("openxlsx not installed — install with install.packages('openxlsx') for Excel output.")
}
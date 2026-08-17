############################################################################
# Figure 4: cross-sectional localized vs de novo mHNPC discrimination in the
# diagnostic biopsy cohorts (GSE268308, GSE268309), using the cohort-
# standardized 45-gene MetaIntegrator signature score. Statistics are computed
# once, written to outs/Figure4/, and re-read so the figure uses the saved values.
#
# Inputs : data/GSE268308_GSE268309.rda, config/metscore_signature_v1.csv
# Outputs: outs/Figure4/{biopsy_validation_stats,biopsy_gene_coverage,
#          biopsy_plot_data,biopsy_roc_coordinates}.csv
#          figures/Figure4_Biopsy_Validation.{pdf,tiff}
############################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
  library(patchwork)
})

DATA_RDA <- "./data/GSE268308_GSE268309.rda"
SIG_CSV  <- "./config/metscore_signature_v1.csv"
OUT_DIR  <- "./outs/Figure4"
FIG_DIR  <- "./figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE)

# ---- inputs -------------------------------------------------------------
e <- new.env()
load(DATA_RDA, envir = e)
stopifnot(all(c("expr_GSE268308", "group_GSE268308", "mapped_symbols_GSE268308",
                "expr_GSE268309", "group_GSE268309", "mapped_symbols_GSE268309") %in% ls(e)))

sig <- read.csv(SIG_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "direction") %in% names(sig)))
pos_genes <- sig$gene[grepl("^POS", sig$direction)]
neg_genes <- sig$gene[grepl("^NEG", sig$direction)]
# exactly 45 unique configured genes, 27 positive / 18 negative, valid directions
stopifnot(length(unique(sig$gene)) == 45L,
          nrow(sig) == 45L,
          length(pos_genes) == 27L,
          length(neg_genes) == 18L,
          length(intersect(pos_genes, neg_genes)) == 0L,
          all(grepl("^(POS|NEG)", sig$direction)))

# ---- MetaIntegrator 45-gene geometric-mean score (exact reproduction) ----
signature_score <- function(expr, pos, neg) {
  expr <- as.matrix(expr)
  stopifnot(is.numeric(expr), !anyDuplicated(rownames(expr)))
  m <- min(expr, na.rm = TRUE)
  if (m < 0) expr <- expr + (abs(m) + 1)          # matches calculateScore shift
  n <- ncol(expr)
  geo <- function(genes) {
    mat <- matrix(1, nrow = length(genes), ncol = n)   # |genes| x n, missing = 1
    present <- intersect(genes, rownames(expr))
    if (length(present)) mat[match(present, genes), ] <- expr[present, , drop = FALSE]
    mat[!is.finite(mat) | mat <= 0] <- 1
    exp(colMeans(log(mat)))                            # per-sample geometric mean
  }
  s <- as.numeric(scale(geo(pos) - geo(neg)))          # within-cohort z
  stopifnot(all(is.finite(s)))
  s
}

# coverage provenance: genes mapped before filtering vs retained/removed by the
# expression filter (removed genes are filtered from the matrix, not unmapped)
gene_coverage <- function(expr, mapped, cohort) {
  rn <- rownames(as.matrix(expr))                 # retained after filtering
  mp_pos <- intersect(pos_genes, mapped); mp_neg <- intersect(neg_genes, mapped)
  rt_pos <- intersect(pos_genes, rn);     rt_neg <- intersect(neg_genes, rn)
  rm_pos <- setdiff(mp_pos, rt_pos);      rm_neg <- setdiff(mp_neg, rt_neg)
  data.frame(cohort = cohort,
             pos_total = length(pos_genes), neg_total = length(neg_genes),
             pos_mapped = length(mp_pos), neg_mapped = length(mp_neg),
             pos_retained = length(rt_pos), neg_retained = length(rt_neg),
             pos_removed = length(rm_pos), neg_removed = length(rm_neg),
             pos_removed_by_filter = paste(rm_pos, collapse = ";"),
             neg_removed_by_filter = paste(rm_neg, collapse = ";"),
             stringsAsFactors = FALSE)
}

# ---- per-cohort statistics ---------------------------------------------
analyze <- function(expr, group, cohort, exp_loc, exp_met) {
  expr <- as.matrix(expr)
  # fail-closed alignment: group is named and identical to expr columns
  stopifnot(!is.null(names(group)), identical(names(group), colnames(expr)))
  src <- as.character(group)
  # only No_Mets / Mets may appear; recode explicitly (no silent default)
  stopifnot(setequal(unique(src), c("No_Mets", "Mets")))
  grp <- factor(c(No_Mets = "Localized", Mets = "mHNPC")[src],
                levels = c("Localized", "mHNPC"))
  stopifnot(!any(is.na(grp)))
  s <- signature_score(expr, pos_genes, neg_genes)
  n_loc <- sum(grp == "Localized"); n_met <- sum(grp == "mHNPC")
  stopifnot(n_loc == exp_loc, n_met == exp_met)         # fail closed on group counts
  # primary: two-sided Wilcoxon rank-sum, normal approximation with continuity correction
  w_asym <- wilcox.test(s ~ grp, exact = FALSE, correct = TRUE)$p.value
  # sensitivity: exact Wilcoxon, only when there are no ties
  has_ties <- any(duplicated(s))
  w_exact <- if (!has_ties) wilcox.test(s ~ grp, exact = TRUE)$p.value else NA_real_
  # descriptive effect size (Cohen's d, mHNPC minus Localized)
  x1 <- s[grp == "mHNPC"]; x0 <- s[grp == "Localized"]
  sp <- sqrt(((length(x1) - 1) * var(x1) + (length(x0) - 1) * var(x0)) / (length(s) - 2))
  d  <- (mean(x1) - mean(x0)) / sp
  # oriented ROC: Localized = control (0), mHNPC = case (1), higher score = case
  y  <- ifelse(grp == "mHNPC", 1L, 0L)
  ro <- pROC::roc(response = y, predictor = s, levels = c(0, 1),
                  direction = "<", quiet = TRUE)
  au <- as.numeric(pROC::auc(ro))
  ci <- as.numeric(pROC::ci.auc(ro, method = "delong"))
  roc_df <- data.frame(cohort = cohort,
                       fpr = 1 - ro$specificities,
                       tpr = ro$sensitivities,
                       stringsAsFactors = FALSE)
  roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]
  list(cohort = cohort, n = length(s), n_localized = n_loc, n_mHNPC = n_met,
       wilcox_p_asymptotic = w_asym, wilcox_p_exact_sensitivity = w_exact,
       has_ties = has_ties, cohens_d = d,
       auc = au, ci_low = ci[1], ci_high = ci[3],
       score = s, group = grp, roc = roc_df)
}

A <- analyze(e$expr_GSE268308, e$group_GSE268308, "GSE268308", 47L, 31L)
B <- analyze(e$expr_GSE268309, e$group_GSE268309, "GSE268309", 17L, 15L)

# BH adjustment across the two primary (asymptotic) p-values
bh <- p.adjust(c(A$wilcox_p_asymptotic, B$wilcox_p_asymptotic), method = "BH")

# ---- canonical outputs --------------------------------------------------
stats_df <- data.frame(
  cohort = c(A$cohort, B$cohort),
  n = c(A$n, B$n),
  n_localized = c(A$n_localized, B$n_localized),
  n_mHNPC = c(A$n_mHNPC, B$n_mHNPC),
  primary_test = "two-sided Wilcoxon rank-sum, normal approximation with continuity correction (exact=FALSE, correct=TRUE)",
  wilcox_p_asymptotic = c(A$wilcox_p_asymptotic, B$wilcox_p_asymptotic),
  bh_q = bh,
  wilcox_p_exact_sensitivity = c(A$wilcox_p_exact_sensitivity, B$wilcox_p_exact_sensitivity),
  has_ties = c(A$has_ties, B$has_ties),
  cohens_d = c(A$cohens_d, B$cohens_d),
  auc = c(A$auc, B$auc),
  auc_ci_low = c(A$ci_low, B$ci_low),
  auc_ci_high = c(A$ci_high, B$ci_high),
  auc_ci_method = "DeLong",
  stringsAsFactors = FALSE)
write.csv(stats_df, file.path(OUT_DIR, "biopsy_validation_stats.csv"), row.names = FALSE)

cov_df <- rbind(gene_coverage(e$expr_GSE268308, e$mapped_symbols_GSE268308, "GSE268308"),
                gene_coverage(e$expr_GSE268309, e$mapped_symbols_GSE268309, "GSE268309"))
# all 45 configured genes map before filtering in both cohorts
stopifnot(all(cov_df$pos_mapped == 27L), all(cov_df$neg_mapped == 18L))
write.csv(cov_df, file.path(OUT_DIR, "biopsy_gene_coverage.csv"), row.names = FALSE)

# plot data: no sample identifiers (cohort, sequential index, group, score)
plot_df <- rbind(
  data.frame(cohort = "GSE268308", sample_index = seq_along(A$score),
             group = as.character(A$group), score = A$score, stringsAsFactors = FALSE),
  data.frame(cohort = "GSE268309", sample_index = seq_along(B$score),
             group = as.character(B$group), score = B$score, stringsAsFactors = FALSE))
write.csv(plot_df, file.path(OUT_DIR, "biopsy_plot_data.csv"), row.names = FALSE)

roc_out <- rbind(A$roc, B$roc)
write.csv(roc_out, file.path(OUT_DIR, "biopsy_roc_coordinates.csv"), row.names = FALSE)

# ---- figure: consume the saved primary values ---------------------------
S <- read.csv(file.path(OUT_DIR, "biopsy_validation_stats.csv"), stringsAsFactors = FALSE)
P <- read.csv(file.path(OUT_DIR, "biopsy_plot_data.csv"), stringsAsFactors = FALSE)
R <- read.csv(file.path(OUT_DIR, "biopsy_roc_coordinates.csv"), stringsAsFactors = FALSE)
P$group <- factor(P$group, levels = c("Localized", "mHNPC"))

COL <- c("Localized" = "#2b2eb5", "mHNPC" = "#ed6905")
# shared y-limits that fully contain the untrimmed violin density of every group,
# so the violins taper smoothly and are not cut off at either end
vext <- range(unlist(lapply(split(P$score, list(P$cohort, P$group)),
                            function(x) range(stats::density(x)$x))))
sp   <- diff(vext)
YLIM <- c(vext[1] - 0.06 * sp, vext[2] + 0.14 * sp)   # headroom at top for the p label
fmt_p <- function(p) formatC(p, format = "e", digits = 1)   # ASCII, e.g. 1.6e-12

dist_panel <- function(coh, tag) {
  d  <- P[P$cohort == coh, ]
  st <- S[S$cohort == coh, ]
  labs_x <- c(sprintf("Localized (n=%d)", st$n_localized),
              sprintf("mHNPC (n=%d)", st$n_mHNPC))
  plab <- sprintf("Wilcoxon p = %s", fmt_p(st$wilcox_p_asymptotic))
  ggplot(d, aes(group, score, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.55, width = 0.9, colour = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.4) +
    geom_jitter(width = 0.08, size = 1.5, alpha = 0.75, colour = "grey20") +
    annotate("text", x = 1.5, y = YLIM[2], label = plab, size = 3.4, vjust = 1) +
    scale_fill_manual(values = COL) +
    scale_x_discrete(labels = labs_x) +
    coord_cartesian(ylim = YLIM) +
    labs(tag = tag, x = NULL, y = "Met-Score (cohort z)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text = element_text(colour = "black"))
}

roc_panel <- function(coh, tag) {
  r  <- R[R$cohort == coh, ]
  st <- S[S$cohort == coh, ]
  lab <- sprintf("AUC = %.3f\n95%% CI %.3f - %.3f", st$auc, st$auc_ci_low, st$auc_ci_high)
  ggplot(r, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 1.1, colour = "#2b2eb5") +
    annotate("text", x = 0.50, y = 0.12, label = lab, size = 3.4, hjust = 0) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(tag = tag, x = "1 - specificity", y = "sensitivity") +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(colour = "black"), plot.margin = margin(5, 8, 5, 5))
}

# cohort name shown once per row (centred over its two panels)
cohort_banner <- function(txt)
  ggplot() + labs(title = txt) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13,
                                    margin = margin(t = 4, b = 2)))

fig <- cohort_banner("GSE268308") /
       (dist_panel("GSE268308", "a") | roc_panel("GSE268308", "b")) /
       cohort_banner("GSE268309") /
       (dist_panel("GSE268309", "c") | roc_panel("GSE268309", "d")) +
  patchwork::plot_layout(heights = c(0.08, 1, 0.08, 1)) &
  theme(plot.tag = element_text(size = 15, face = "bold"))

ggsave(file.path(FIG_DIR, "Figure4_Biopsy_Validation.pdf"), fig,
       width = 8.0, height = 8.0, useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "Figure4_Biopsy_Validation.tiff"), fig,
       width = 8.0, height = 8.0, dpi = 400, compression = "lzw")

cat("\n---- Figure 4 biopsy validation ----\n")
print(stats_df[, c("cohort", "n", "n_localized", "n_mHNPC",
                    "wilcox_p_asymptotic", "bh_q", "wilcox_p_exact_sensitivity",
                    "cohens_d", "auc", "auc_ci_low", "auc_ci_high")])
print(cov_df[, c("cohort", "pos_mapped", "pos_retained", "pos_removed",
                 "neg_mapped", "neg_retained", "neg_removed")])
cat("Wrote outs/Figure4/*.csv and figures/Figure4_Biopsy_Validation.{pdf,tiff}\n")
cat("=== DONE ===\n")

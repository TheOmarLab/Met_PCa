################################################################################
# Purity_ESTIMATE.R
# ESTIMATE-based tumor purity sensitivity analysis for JHU and Durham cohorts.
#
# Motivation: since MetScore integrates TME signals, tumor purity variation
# could in principle confound the score in bulk RNA-seq / microarray cohorts.
#
# Approach:
#   1. Run ESTIMATE on JHU (Affymetrix microarray) and Durham (microarray)
#      to derive StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity.
#   2. Compute Spearman correlation: MetScore prob vs ESTIMATE TumorPurity.
#   3. Figure S1a: per-cohort-SD Met-Score HR before vs after adding ESTIMATE
#      TumorPurity, on one common complete-case set (JHU Lin-Ying case-cohort;
#      Durham complete-cohort cause-specific Cox).
#   4. Descriptive scatter / violin figures (JHU + Durham).
#
# Outputs (all under outs/purity/ and figures/purity/):
#   outs/purity/JHU_ESTIMATE_scores.csv
#   outs/purity/Durham_ESTIMATE_scores.csv
#   outs/purity/JHU_purity_correlation.txt
#   outs/purity/Durham_purity_correlation.txt
#   outs/purity/JHU_purity_adjusted_Cox_models.rda   (before/after primary models)
#   outs/purity/Durham_purity_adjusted_Cox_models.rda (before/after primary models)
#   outs/purity/FigureS1a_purity_robustness.csv
#   outs/purity/Durham_platform_provenance.txt
#   figures/purity/JHU_ESTIMATE_purity_vs_MetScore.pdf/.tiff
#   figures/purity/Durham_ESTIMATE_purity_vs_MetScore.pdf/.tiff
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(patchwork)
})

# --- 0. Load ESTIMATE (dependency; not installed here) ---------------------
# ESTIMATE (R-Forge) is a required dependency documented in README/setup; this
# analysis script never installs packages. Fail closed if it is unavailable.
if (!requireNamespace("estimate", quietly = TRUE))
  stop("Package 'estimate' is required (see README/setup for installation); it is not installed.")
library(estimate)

# --- 1. Output directories -------------------------------------------------
dir.create("./outs/purity",            showWarnings = FALSE, recursive = TRUE)
dir.create("./figures/purity",         showWarnings = FALSE, recursive = TRUE)
tmp_dir <- "./outs/purity/tmp_estimate"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Expression matrix writer for ESTIMATE input -----------------------
# filterCommonGenes() reads a plain tab-delimited file (row.names=1, header=TRUE)
# WITHOUT GCT prefix lines; do not include "#1.2" or dim header.
write_for_estimate <- function(expr_mat, out_path) {
  # Collapse duplicate gene symbols to the highest-variance row rather than the
  # first occurrence, so the retained probe is the most informative one.
  if (any(duplicated(rownames(expr_mat)))) {
    row_var  <- apply(expr_mat, 1, stats::var, na.rm = TRUE)
    expr_mat <- expr_mat[order(row_var, decreasing = TRUE), , drop = FALSE]
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  }
  df <- data.frame(
    Name        = rownames(expr_mat),
    Description = rownames(expr_mat),
    as.data.frame(expr_mat, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  write.table(df, file = out_path,
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  invisible(out_path)
}

# --- 3. GCT parser helper --------------------------------------------------
# Parses ESTIMATE output GCT → data.frame with sample_id + score columns.
# ESTIMATE output has header "NAME  Description  Description.1  sample1  ..."
parse_estimate_gct <- function(gct_path) {
  lns <- readLines(gct_path)
  # Lines 1-2 are "#1.2" and dim; line 3 onward is the table
  dat <- read.table(
    text = paste(lns[-(1:2)], collapse = "\n"),
    sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  # Use first column (NAME / Name) as row identifiers for scores
  rownames(dat) <- dat[[1]]
  # Drop all non-sample columns: first column + any Description* variants
  skip <- union(1L, grep("^Description", colnames(dat)))
  dat  <- dat[, -skip, drop = FALSE]
  # Transpose so rows = samples, cols = score types
  out <- as.data.frame(t(dat), stringsAsFactors = FALSE)
  out[] <- lapply(out, as.numeric)
  out$sample_id <- rownames(out)
  rownames(out) <- NULL
  return(out)
}

# Event colors used consistently across all scatter plots
EVENT_COLORS <- c("No metastasis" = "#9CA3AF", "Metastasis" = "#C0392B")

# --- 4. Scatter plot helper ------------------------------------------------
# color_col: optional column in df to map to point color (e.g. "MetEvent").
#            When NULL, all points use pt_color. When provided, uses EVENT_COLORS.
make_purity_scatter <- function(df, x_col, y_col, rho, pval, title,
                                pt_color, color_col = NULL) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]]))

  if (!is.null(color_col)) {
    p <- p + geom_point(aes(color = .data[[color_col]]), alpha = 0.65, size = 1.8)
  } else {
    p <- p + geom_point(alpha = 0.65, size = 1.8, color = pt_color)
  }

  p <- p +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
    labs(x = "ESTIMATE Tumor Purity", y = "Met-Score Probability", title = title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(face = "bold"),
      legend.title = element_blank()
    ) +
    annotate("text",
      x     = min(df[[x_col]], na.rm = TRUE) + 0.02,
      y     = max(df[[y_col]], na.rm = TRUE) - 0.02,
      hjust = 0, size = 4.5,
      label = paste0(
        "Spearman rho = ", sprintf("%.2f", rho),
        "\nP = ", format(pval, digits = 2, scientific = TRUE)
      )
    )

  if (!is.null(color_col)) {
    p <- p + scale_color_manual(values = EVENT_COLORS)
  }
  p
}

# --- 5. Load JHU data -------------------------------------------------------
message("\n=== JHU ===")
load("./data/Dataset7.rda")        # provides Dataset7 (MetaIntegrator object)
load("./outs/coxdata.rda")         # provides CoxData_jhu
load("./outs/jhu_pheno_filter_MetaScorer_Zscore.rda")  # provides pheno_jhu

# Dataset7$expr: full gene-symbol × sample matrix, log2-transformed Affymetrix
jhu_expr_full <- Dataset7$expr

# Identify common samples between expression matrix and CoxData_jhu
jhu_ids_cox  <- rownames(CoxData_jhu)
jhu_ids_expr <- colnames(jhu_expr_full)
common_jhu   <- intersect(jhu_ids_expr, jhu_ids_cox)
message(sprintf("JHU: expr cols = %d | CoxData rows = %d | common = %d",
                length(jhu_ids_expr), length(jhu_ids_cox), length(common_jhu)))
stopifnot(length(common_jhu) > 0)
jhu_expr_aligned <- jhu_expr_full[, common_jhu, drop = FALSE]

# Run ESTIMATE
jhu_raw_gct    <- file.path(tmp_dir, "jhu_full.gct")
jhu_common_gct <- file.path(tmp_dir, "jhu_common.gct")
jhu_scores_gct <- file.path(tmp_dir, "jhu_estimate_scores.gct")

write_for_estimate(jhu_expr_aligned, jhu_raw_gct)
filterCommonGenes(input.f = jhu_raw_gct, output.f = jhu_common_gct, id = "GeneSymbol")
# JHU is Affymetrix Human Exon 1.0 ST (GPL5188; ESTIMATE affymetrix)
estimateScore(input.ds = jhu_common_gct, output.ds = jhu_scores_gct, platform = "affymetrix")

jhu_est <- parse_estimate_gct(jhu_scores_gct)
message(sprintf("ESTIMATE scores computed for %d JHU samples", nrow(jhu_est)))
write.csv(jhu_est, "./outs/purity/JHU_ESTIMATE_scores.csv", row.names = FALSE)

# Merge ESTIMATE scores with MetScore + survival data.
# rownames(CoxData_jhu) and colnames(Dataset7$expr) are both CEL file names
# (e.g. GDX.4612.JHU063.CEL); jhu_est$sample_id is also keyed on CEL names.
# Canonical pathological Grade Group from total pathological Gleason sum and
# pathological primary grade. GG2 (3+4) is the reference used across the repo.
gg_from <- function(total_gs, primary) {
  gg <- rep(NA_character_, length(total_gs))
  gg[total_gs <= 6]                <- "GG1"
  gg[total_gs == 7 & primary == 3] <- "GG2"
  gg[total_gs == 7 & primary == 4] <- "GG3"
  gg[total_gs == 8]                <- "GG4"
  gg[total_gs >= 9]                <- "GG5"
  gg
}

# Pathological T stage collapsed to T2 / T3 / T4 (T2 reference), matching the
# canonical multivariable comparator.
pt_collapse <- function(x) {
  x   <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^T2", x)] <- "T2"
  out[grepl("^T3", x)] <- "T3"
  out[grepl("^T4", x)] <- "T4"
  out
}

# PRIMARY purity comparison (Figure S1a). The Met-Score enters as the continuous
# per-cohort-SD z of the locked classifier probability (used, never refit),
# adjusted for the common comparator (pathological Grade Group + log2(PSA+1) +
# pathological pT). Fit twice on one common complete-case set: BEFORE purity and
# AFTER adding ESTIMATE TumorPurity. The cause-specific metastasis clock censors a
# death without prior metastasis at its death time. JHU uses a Lin-Ying case-cohort
# Cox (cohort.size 745, robust, verified design strata); Durham uses complete-cohort
# cause-specific Cox. JHU fails rather than falling back to ordinary Cox.
fit_purity_perSD <- function(df, cohort, estimator, met_col, mettime_col,
                             death_col, deathtime_col) {
  m  <- as.integer(df[[met_col]]);   mt <- as.numeric(df[[mettime_col]])
  de <- as.integer(df[[death_col]]); dt <- as.numeric(df[[deathtime_col]])
  status <- ifelse(m == 1L, 1L, ifelse(de == 1L & m == 0L & !is.na(dt) & dt <= mt, 2L, 0L))
  df$analysis_time <- ifelse(status == 2L, dt, mt)
  df$status1       <- as.integer(status == 1L)
  need <- c("GG", "log2PSA", "pT", "MetScoreProb", "TumorPurity", "analysis_time", "status1")
  if (estimator == "cch") need <- c(need, "insub", "id")
  cc <- df[stats::complete.cases(df[, need]) & is.finite(df$analysis_time) & df$analysis_time > 0, , drop = FALSE]
  cc$GG   <- relevel(factor(as.character(cc$GG)), ref = "GG2")
  cc$pT   <- relevel(factor(as.character(cc$pT)), ref = "T2")
  cc$ms_z <- as.numeric(scale(cc$MetScoreProb))          # per common-set SD; HR-unit transform only
  f_before <- Surv(analysis_time, status1) ~ ms_z + GG + log2PSA + pT
  f_after  <- update(f_before, . ~ . + TumorPurity)
  if (estimator == "cch") {
    fitb <- survival::cch(f_before, data = cc, subcoh = ~ insub, id = ~ id,
                          cohort.size = 745L, method = "LinYing", robust = TRUE)
    fita <- survival::cch(f_after,  data = cc, subcoh = ~ insub, id = ~ id,
                          cohort.size = 745L, method = "LinYing", robust = TRUE)
    hr <- function(fit) { cf <- summary(fit)$coefficients; b <- cf["ms_z", "Value"]; se <- cf["ms_z", "SE"]
      c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = cf["ms_z", "p"]) }
    est_label  <- "Lin-Ying case-cohort Cox (cohort.size=745, robust)"
    strata_txt <- paste(sprintf("%s=%d", names(table(cc$cch)), as.integer(table(cc$cch))), collapse = "; ")
    basis      <- sprintf("complete-case n=%d of 239; design strata %s", nrow(cc), strata_txt)
  } else {
    fitb <- survival::coxph(f_before, data = cc); fita <- survival::coxph(f_after, data = cc)
    hr <- function(fit) { s <- summary(fit)
      c(hr = s$conf.int["ms_z", "exp(coef)"], lo = s$conf.int["ms_z", "lower .95"],
        hi = s$conf.int["ms_z", "upper .95"], p = s$coefficients["ms_z", "Pr(>|z|)"]) }
    est_label <- "complete-cohort cause-specific Cox"
    basis     <- sprintf("complete-case n=%d of 555", nrow(cc))
  }
  b <- hr(fitb); a <- hr(fita)
  mk <- function(model, x) data.frame(
    cohort = cohort, model = model, n = nrow(cc), events = sum(cc$status1),
    HR_perSD = round(unname(x["hr"]), 4), CI_lo = round(unname(x["lo"]), 4),
    CI_hi = round(unname(x["hi"]), 4), p_value = signif(unname(x["p"]), 4),
    estimator = est_label,
    covariates = if (model == "before") "GG + log2(PSA+1) + pT + MetScore_perSD"
                 else "GG + log2(PSA+1) + pT + MetScore_perSD + TumorPurity",
    event_clock = "cause-specific metastasis; death without metastasis censored at death time",
    platform = if (cohort == "JHU") "Affymetrix Human Exon 1.0 ST (GPL5188; ESTIMATE affymetrix)" else DUR_PLATFORM,
    design_basis = basis, stringsAsFactors = FALSE)
  list(cohort = cohort, before = fitb, after = fita, rows = rbind(mk("before", b), mk("after", a)))
}

CoxData_jhu$cel_file <- rownames(CoxData_jhu)
jhu_merged <- CoxData_jhu %>%
  dplyr::rename(MetScoreProb = `Met-Score prob`) %>%
  left_join(jhu_est %>% dplyr::rename(cel_file = sample_id), by = "cel_file") %>%
  filter(!is.na(TumorPurity), !is.na(MetScoreProb))
jhu_merged$MetEvent <- factor(jhu_merged$met,
  levels = c(0, 1), labels = c("No metastasis", "Metastasis"))

# Pathological Grade Group (GG2 reference) for the purity-adjusted models
jhu_merged$GG <- relevel(
  factor(gg_from(as.numeric(as.character(jhu_merged$`Pathological GS`)),
                 as.numeric(as.character(jhu_merged$pathgs_p)))),
  ref = "GG2"
)

# Common-comparator covariates for the PRIMARY purity model:
#   log2PSA = log2(preop PSA + 1), keeps PSA == 0 finite
#   pT      = pathological T stage collapsed to T2/T3/T4
#   ms_z    = per-SD z of the locked classifier probability (used, not refit)
jhu_merged$log2PSA <- log2(as.numeric(jhu_merged$preop_psa) + 1)
jhu_merged$pT      <- pt_collapse(jhu_merged$pstage)
jhu_merged$ms_z    <- as.numeric(scale(jhu_merged$MetScoreProb))
message(sprintf("JHU merged N = %d", nrow(jhu_merged)))

# Spearman correlation
jhu_cor <- cor.test(jhu_merged$TumorPurity, jhu_merged$MetScoreProb, method = "spearman")
message(sprintf("JHU Spearman rho = %.4f  p = %.4g",
                unname(jhu_cor$estimate), jhu_cor$p.value))
cat("JHU Spearman rho =", unname(jhu_cor$estimate),
    " p =", jhu_cor$p.value, "\n",
    file = "./outs/purity/JHU_purity_correlation.txt")

# PRIMARY purity-adjusted model (MFS endpoint): continuous per-SD Met-Score
# with the common comparator (Grade Group + log2PSA + pT), fit BEFORE and AFTER
# adding ESTIMATE TumorPurity on one complete-case set. This specification
# reproduces the per-SD Met-Score anchors.
# Full 239-person case-cohort ledger, verified before complete-case subsetting.
stopifnot(nrow(CoxData_jhu) == 239L, sum(as.integer(CoxData_jhu$met) == 1L) == 93L)
.jhu_cch <- as.character(CoxData_jhu[["post_rp_patients_cchdef"]])
stopifnot(sum(.jhu_cch == "Sub-cohort cases") == 28L,
          sum(.jhu_cch == "Sub-cohort controls") == 146L,
          sum(.jhu_cch == "cases") == 65L)
message("JHU full ledger: 239 (28 subcohort cases / 146 subcohort controls / 65 outside cases), 93 metastases")
jhu_merged$insub <- as.integer(as.character(jhu_merged$post_rp_patients_cchdef) %in%
                               c("Sub-cohort cases", "Sub-cohort controls"))
jhu_merged$id    <- seq_len(nrow(jhu_merged))
jhu_merged$cch   <- as.character(jhu_merged$post_rp_patients_cchdef)
jhu_primary    <- fit_purity_perSD(jhu_merged, "JHU", "cch", "met", "met_time", "os", "os_time")
jhu_cox_before <- jhu_primary$before
jhu_cox_after  <- jhu_primary$after
cat("\n-- JHU primary purity model: per-SD Met-Score, after adjustment --\n")
print(summary(jhu_cox_after))

save(jhu_cox_before, jhu_cox_after,
     file = "./outs/purity/JHU_purity_adjusted_Cox_models.rda")

# Scatter plot
p_jhu <- make_purity_scatter(
  df        = jhu_merged,
  x_col     = "TumorPurity",
  y_col     = "MetScoreProb",
  rho       = unname(jhu_cor$estimate),
  pval      = jhu_cor$p.value,
  title     = "JHU: Met-Score vs ESTIMATE Tumor Purity",
  pt_color  = "#0072B2",
  color_col = "MetEvent"
)
ggsave("./figures/purity/JHU_ESTIMATE_purity_vs_MetScore.pdf",
       p_jhu, width = 6.2, height = 4.8, useDingbats = FALSE)
ggsave("./figures/purity/JHU_ESTIMATE_purity_vs_MetScore.tiff",
       p_jhu, width = 6.2, height = 4.8, dpi = 450, compression = "lzw")

# --- 6. Load Durham data ----------------------------------------------------
message("\n=== Durham ===")
# Durham assay-platform provenance: Sheet2 of the cohort workbook lists the raw
# array identifiers. Assert the Affymetrix Human Exon 1.0 ST-family evidence
# that fixes the ESTIMATE platform to "affymetrix".
DUR_XLSX <- "data/Durham_cohort_and_GRID_cohort/Durham_cohort_011526.xlsx"
.prov <- suppressMessages(read_excel(DUR_XLSX, sheet = "Sheet2"))
.cel  <- unlist(lapply(.prov, as.character)); .cel <- .cel[!is.na(.cel) & grepl("\\.CEL$", .cel, ignore.case = TRUE)]
.huex <- .cel[grepl("HuEx-1_0-st-v2", .cel)]
stopifnot(length(.cel) == 887L, length(.huex) == 121L)
DUR_PLATFORM <- "Affymetrix Human Exon 1.0 ST-v2 (887 CEL identifiers, 121 HuEx-1_0-st-v2; ESTIMATE affymetrix)"
writeLines(c(sprintf("Durham platform provenance from %s Sheet2", DUR_XLSX),
             sprintf("total .CEL identifiers: %d", length(.cel)),
             sprintf("HuEx-1_0-st-v2 identifiers: %d", length(.huex)),
             "platform = affymetrix (Human Exon 1.0 ST-family)"),
           "./outs/purity/Durham_platform_provenance.txt")
message(sprintf("Durham platform: %d .CEL (%d HuEx-1_0-st-v2) -> affymetrix", length(.cel), length(.huex)))
durham_expr_raw <- read_excel(DUR_XLSX, sheet = "eset_gene_filtered")
expr_mat_dur    <- as.data.frame(durham_expr_raw)
gene_syms_dur   <- expr_mat_dur$Symbol
expr_mat_dur    <- expr_mat_dur[, -c(1, 2), drop = FALSE]
expr_mat_dur    <- apply(expr_mat_dur, 2, as.numeric)
rownames(expr_mat_dur) <- gene_syms_dur

load("./output/Durham/durham_metscore_batchcorrected.rda")  # clin_valid

# Drop sentinel rows (pogl == 0) if present
if ("pogl" %in% colnames(clin_valid)) {
  clin_valid <- clin_valid[!is.na(clin_valid$pogl) & clin_valid$pogl != 0, ]
}
clin_valid <- clin_valid[!is.na(clin_valid$MetScore_prob), ]

common_dur <- intersect(colnames(expr_mat_dur), clin_valid$sample_id)
message(sprintf("Durham: expr cols = %d | clin_valid rows = %d | common = %d",
                ncol(expr_mat_dur), nrow(clin_valid), length(common_dur)))
stopifnot(length(common_dur) > 0)
dur_expr_aligned <- expr_mat_dur[, common_dur, drop = FALSE]

# Run ESTIMATE
dur_raw_gct    <- file.path(tmp_dir, "durham_full.gct")
dur_common_gct <- file.path(tmp_dir, "durham_common.gct")
dur_scores_gct <- file.path(tmp_dir, "durham_estimate_scores.gct")

write_for_estimate(dur_expr_aligned, dur_raw_gct)
filterCommonGenes(input.f = dur_raw_gct, output.f = dur_common_gct, id = "GeneSymbol")
# Durham is an Affymetrix Human Exon 1.0 ST-family dataset (provenance asserted
# above from Sheet2); ESTIMATE platform = "affymetrix".
estimateScore(input.ds = dur_common_gct, output.ds = dur_scores_gct, platform = "affymetrix")

dur_est <- parse_estimate_gct(dur_scores_gct)
message(sprintf("ESTIMATE scores computed for %d Durham samples", nrow(dur_est)))
write.csv(dur_est, "./outs/purity/Durham_ESTIMATE_scores.csv", row.names = FALSE)

# Merge with clin_valid
dur_merged <- clin_valid %>%
  left_join(dur_est, by = "sample_id") %>%
  filter(!is.na(TumorPurity), !is.na(MetScore_prob))
dur_merged$MetEvent <- factor(dur_merged$mets,
  levels = c(0, 1), labels = c("No metastasis", "Metastasis"))

# Pathological Grade Group (GG2 reference) for the purity-adjusted models
dur_merged$GG <- relevel(
  factor(gg_from(as.numeric(as.character(dur_merged$PathGleason)),
                 as.numeric(as.character(dur_merged$pogl1)))),
  ref = "GG2"
)

# Common-comparator covariates for the PRIMARY purity model (see JHU above).
dur_merged$log2PSA <- log2(as.numeric(dur_merged$psapresurg) + 1)
dur_merged$pT      <- pt_collapse(dur_merged$stg)
dur_merged$ms_z    <- as.numeric(scale(dur_merged$MetScore_prob))
message(sprintf("Durham merged N = %d", nrow(dur_merged)))

dur_cor <- cor.test(dur_merged$TumorPurity, dur_merged$MetScore_prob, method = "spearman")
message(sprintf("Durham Spearman rho = %.4f  p = %.4g",
                unname(dur_cor$estimate), dur_cor$p.value))
cat("Durham Spearman rho =", unname(dur_cor$estimate),
    " p =", dur_cor$p.value, "\n",
    file = "./outs/purity/Durham_purity_correlation.txt")

# PRIMARY purity-adjusted model (MFS endpoint: surgmets / mets): continuous
# per-SD Met-Score with the common comparator (Grade Group + log2PSA + pT),
# before and after adding ESTIMATE TumorPurity on one complete-case set.
dur_merged$MetScoreProb <- dur_merged$MetScore_prob
dur_primary    <- fit_purity_perSD(dur_merged, "Durham", "cox", "mets", "surgmets", "dead", "limbo")
dur_cox_before <- dur_primary$before
dur_cox_after  <- dur_primary$after
cat("\n-- Durham primary purity model: per-SD Met-Score, after adjustment --\n")
print(summary(dur_cox_after))

save(dur_cox_before, dur_cox_after,
     file = "./outs/purity/Durham_purity_adjusted_Cox_models.rda")

# Scatter plot
p_dur <- make_purity_scatter(
  df        = dur_merged,
  x_col     = "TumorPurity",
  y_col     = "MetScore_prob",
  rho       = unname(dur_cor$estimate),
  pval      = dur_cor$p.value,
  title     = "Durham: Met-Score vs ESTIMATE Tumor Purity",
  pt_color  = "#009E73",
  color_col = "MetEvent"
)
ggsave("./figures/purity/Durham_ESTIMATE_purity_vs_MetScore.pdf",
       p_dur, width = 6.2, height = 4.8, useDingbats = FALSE)
ggsave("./figures/purity/Durham_ESTIMATE_purity_vs_MetScore.tiff",
       p_dur, width = 6.2, height = 4.8, dpi = 450, compression = "lzw")

# --- 8. Figure S1a canonical output ----------------------------------------
message("\n=== Figure S1a summary ===")

# Per-cohort-SD Met-Score HR, before vs after adding ESTIMATE TumorPurity, on one
# common complete-case set per cohort. Sole purity-robustness output.
s1a <- rbind(jhu_primary$rows, dur_primary$rows)

# Fail-closed checks: finite four-row output; before/after share complete-case n
# and events within a cohort; cause-specific clock and complete-case ledgers hold.
stopifnot(
  nrow(s1a) == 4L,
  all(is.finite(c(s1a$HR_perSD, s1a$CI_lo, s1a$CI_hi, s1a$p_value))),
  all(s1a$event_clock == "cause-specific metastasis; death without metastasis censored at death time")
)
for (co in unique(s1a$cohort)) {
  rr <- s1a[s1a$cohort == co, ]
  stopifnot(length(unique(rr$n)) == 1L, length(unique(rr$events)) == 1L)
}
stopifnot(
  unique(s1a$n[s1a$cohort == "JHU"])    == 235L, unique(s1a$events[s1a$cohort == "JHU"])    == 93L,
  unique(s1a$n[s1a$cohort == "Durham"]) == 555L, unique(s1a$events[s1a$cohort == "Durham"]) == 40L
)

cat("\nFigure S1a: per-cohort-SD Met-Score HR, before vs after ESTIMATE purity:\n")
print(s1a[, c("cohort", "model", "n", "events", "HR_perSD", "CI_lo", "CI_hi", "p_value", "estimator")],
      row.names = FALSE)
write.csv(s1a, "./outs/purity/FigureS1a_purity_robustness.csv", row.names = FALSE)

# --- Figure S1a microenvironment correlations (Met-Score vs ESTIMATE) -------
# Met-Score rank correlation with tumour purity, stromal score, immune score.
# JHU uses phase-two design weights (controls 745/265) with a conditional
# design-stratified bootstrap; Durham uses complete-cohort Spearman with a
# patient bootstrap. B fixed-seed replicates; percentile 95% CIs.
ALPHA_S1A <- 265 / 745
S1A_B <- as.integer(Sys.getenv("METPCA_PURITY_BOOTSTRAPS", "2000"))
S1A_SEED <- 20260814L
wt_spearman <- function(x, y, w) { rx <- rank(x); ry <- rank(y)
  mx <- sum(w * rx) / sum(w); my <- sum(w * ry) / sum(w)
  num <- sum(w * (rx - mx) * (ry - my)); dx <- sum(w * (rx - mx)^2); dy <- sum(w * (ry - my)^2)
  num / sqrt(dx * dy) }
me_scores <- c(TumorPurity = "Tumour purity", StromalScore = "Stromal score", ImmuneScore = "Immune score")
jw <- ifelse(as.character(jhu_merged$cch) == "Sub-cohort controls", 1 / ALPHA_S1A, 1)
s1a_me <- list()
for (sc in names(me_scores)) {
  ok <- is.finite(jhu_merged[[sc]]) & is.finite(jhu_merged$MetScoreProb)
  xj <- jhu_merged[[sc]][ok]; yj <- jhu_merged$MetScoreProb[ok]; wj <- jw[ok]
  strj <- split(seq_along(xj), as.character(jhu_merged$cch)[ok])
  set.seed(S1A_SEED); bj <- numeric(0); ju <- 0L; jf <- 0L
  for (b in seq_len(S1A_B)) { ix <- unlist(lapply(strj, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
    v <- tryCatch(wt_spearman(xj[ix], yj[ix], wj[ix]), error = function(e) NA_real_)
    if (is.finite(v)) { bj <- c(bj, v); ju <- ju + 1L } else jf <- jf + 1L }
  qj <- quantile(bj, c(.025, .975))
  s1a_me[[length(s1a_me) + 1L]] <- data.frame(cohort = "JHU", score = unname(me_scores[[sc]]), n = sum(ok),
    rho = round(wt_spearman(xj, yj, wj), 5), ci_lo = round(unname(qj[1]), 5), ci_hi = round(unname(qj[2]), 5),
    estimator = "design-weighted Spearman", resampling = "conditional design-stratified bootstrap",
    boot_seed = S1A_SEED, boot_used = ju, boot_failed = jf, stringsAsFactors = FALSE)
  okd <- is.finite(dur_merged[[sc]]) & is.finite(dur_merged$MetScore_prob)
  xd <- dur_merged[[sc]][okd]; yd <- dur_merged$MetScore_prob[okd]; nd <- length(xd)
  set.seed(S1A_SEED); bd <- numeric(0); du <- 0L; dfl <- 0L
  for (b in seq_len(S1A_B)) { ix <- sample.int(nd, nd, replace = TRUE)
    v <- tryCatch(suppressWarnings(cor(xd[ix], yd[ix], method = "spearman")), error = function(e) NA_real_)
    if (is.finite(v)) { bd <- c(bd, v); du <- du + 1L } else dfl <- dfl + 1L }
  qd <- quantile(bd, c(.025, .975))
  s1a_me[[length(s1a_me) + 1L]] <- data.frame(cohort = "Durham", score = unname(me_scores[[sc]]), n = sum(okd),
    rho = round(suppressWarnings(cor(xd, yd, method = "spearman")), 5), ci_lo = round(unname(qd[1]), 5), ci_hi = round(unname(qd[2]), 5),
    estimator = "Spearman", resampling = "patient bootstrap",
    boot_seed = S1A_SEED, boot_used = du, boot_failed = dfl, stringsAsFactors = FALSE)
}
s1a_me <- do.call(rbind, s1a_me)
write.csv(s1a_me, "./outs/purity/FigureS1a_microenvironment_correlations.csv", row.names = FALSE)
cat("\nFigure S1a microenvironment correlations (Met-Score vs ESTIMATE):\n")
print(s1a_me[, c("cohort", "score", "n", "rho", "ci_lo", "ci_hi")], row.names = FALSE)

# --- 9. Option A: Purity distribution by Met-Score class -------------------
# If MetScore were just tracking purity, High-risk samples would have
# systematically lower purity. This violin tests that directly.
message("\n=== Option A: Purity by MetScore class ===")

df_purity_class <- bind_rows(
  jhu_merged %>%
    transmute(Cohort = "JHU\n(MFS)", TumorPurity,
              MetScoreClass = as.character(MetScoreClass)),
  dur_merged %>%
    transmute(Cohort = "Durham\n(MFS)", TumorPurity,
              MetScoreClass = as.character(MetScoreClass))
) %>%
  filter(!is.na(MetScoreClass)) %>%
  mutate(
    MetScoreClass = factor(MetScoreClass, levels = c("Low risk", "High risk")),
    Cohort        = factor(Cohort, levels = c("JHU\n(MFS)", "Durham\n(MFS)"))
  )

risk_colors <- c("Low risk" = "#6baed6", "High risk" = "#e63946")

p_purity_class <- ggplot(df_purity_class,
                          aes(x = MetScoreClass, y = TumorPurity, fill = MetScoreClass)) +
  geom_violin(trim = FALSE, alpha = 0.55, width = 0.8, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     size = 3.5, label.y.npc = 0.96) +
  scale_fill_manual(values = risk_colors) +
  facet_wrap(~Cohort, nrow = 1) +
  labs(
    x     = NULL,
    y     = "Tumor Purity (ESTIMATE)",
    title = "Tumor Purity by Met-Score Risk Class"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text  = element_text(face = "bold", size = 11),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    axis.text.x = element_text(size = 10)
  )

ggsave("./figures/purity/Purity_by_MetScoreClass.pdf",
       p_purity_class, width = 10, height = 4.5, useDingbats = FALSE)
ggsave("./figures/purity/Purity_by_MetScoreClass.tiff",
       p_purity_class, width = 10, height = 4.5, dpi = 450, compression = "lzw")

# --- 10. JHU + Durham purity vs Met-Score scatter (descriptive) ------------
message("\n=== JHU + Durham scatter ===")

df_jhu_2 <- jhu_merged %>%
  transmute(TumorPurity, MetScoreProb, MetEvent, Cohort = "JHU")

df_dur_2 <- dur_merged %>%
  transmute(TumorPurity, MetScoreProb = MetScore_prob, MetEvent, Cohort = "VA Durham")

df_2cohort <- bind_rows(df_jhu_2, df_dur_2) %>%
  mutate(Cohort = factor(Cohort, levels = c("JHU", "VA Durham")))

annot_2 <- data.frame(
  Cohort    = factor(c("JHU", "VA Durham"), levels = c("JHU", "VA Durham")),
  rho_label = c(
    sprintf("rho = %.2f\np = %s", unname(jhu_cor$estimate),
            format(jhu_cor$p.value, digits = 2, scientific = TRUE)),
    sprintf("rho = %.2f\np = %s", unname(dur_cor$estimate),
            format(dur_cor$p.value, digits = 2, scientific = TRUE))
  ),
  x = c(0.35, 0.35),
  y = c(0.94, 0.94),
  stringsAsFactors = FALSE
)

p_scatter_2 <- ggplot(df_2cohort,
                       aes(x = TumorPurity, y = MetScoreProb, color = MetEvent)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE,
              aes(x = TumorPurity, y = MetScoreProb)) +
  scale_color_manual(values = EVENT_COLORS,
                     name   = "Metastasis status") +
  geom_text(data = annot_2,
            aes(x = x, y = y, label = rho_label),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = 3.6, color = "black") +
  facet_wrap(~Cohort, nrow = 1, scales = "free_y") +
  labs(
    x     = "Tumor Purity (ESTIMATE)",
    y     = "Met-Score Probability",
    title = "Tumor Purity vs Met-Score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 10),
    strip.text         = element_text(face = "bold", size = 12),
    strip.background   = element_blank(),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    axis.title         = element_text(face = "bold"),
    axis.text          = element_text(face = "bold")
  )

ggsave("./figures/purity/JHU_VA_purity_MetScore_scatter.pdf",
       p_scatter_2, width = 9, height = 4.5, useDingbats = FALSE)
ggsave("./figures/purity/JHU_VA_purity_MetScore_scatter.tiff",
       p_scatter_2, width = 9, height = 4.5, dpi = 450, compression = "lzw")

cat("\n=== Purity_ESTIMATE.R complete ===\n")
cat("Outputs:\n")
cat("  outs/purity/JHU_ESTIMATE_scores.csv\n")
cat("  outs/purity/Durham_ESTIMATE_scores.csv\n")
cat("  outs/purity/JHU_purity_correlation.txt\n")
cat("  outs/purity/Durham_purity_correlation.txt\n")
cat("  outs/purity/JHU_purity_adjusted_Cox_models.rda\n")
cat("  outs/purity/Durham_purity_adjusted_Cox_models.rda\n")
cat("  outs/purity/FigureS1a_purity_robustness.csv\n")
cat("  outs/purity/Durham_platform_provenance.txt\n")
cat("  figures/purity/JHU_ESTIMATE_purity_vs_MetScore.pdf/.tiff\n")
cat("  figures/purity/Durham_ESTIMATE_purity_vs_MetScore.pdf/.tiff\n")
cat("  figures/purity/Purity_by_MetScoreClass.pdf/.tiff\n")
cat("  figures/purity/JHU_VA_purity_MetScore_scatter.pdf/.tiff\n")

# Remove temporary ESTIMATE intermediates after successful execution.
unlink(tmp_dir, recursive = TRUE)
cat("Removed temporary ESTIMATE intermediates:", tmp_dir, "\n")

############################################################################
# Met-Score External Validation - Durham Cohort (Batch-Corrected)
#
# Two scores are computed for external validation WITHOUT data leakage:
#
# 1. ACCEPTED DEPLOYED SCORE - the frozen locked Met-Score classifier:
#    - Quantile-normalize Durham signature genes to the training reference
#    - Training quantiles from training data only (no Durham information)
#    - Apply the frozen tracked logistic-regression contract (no refit)
#    - Dichotomize at the frozen training-derived threshold
#    This is the deployed Met-Score behind the accepted validation results.
#
# 2. SUPPORTING / DESCRIPTIVE - a model-free gene-set direction score:
#    - 27 pos / 18 neg gene directions from MetaIntegrator
#    - Per-gene rank Z-score within Durham (no outcome used)
#    - Score = mean(Z of pos genes) - mean(Z of neg genes), median split
#    Shown for platform-agnostic robustness only; not the deployed classifier.
#
# Classification & survival analyses match JHU approach exactly:
#   - ROC, Confusion Matrix, MCC (classification endpoint: metastasis)
#   - KM curves with matched styling
#   - Univariate Cox (MetScoreClass) + a restricted-subset, Gleason-only
#     adjusted Cox that is DESCRIPTIVE ONLY (see note below)
#   - C-index comparison (Gleason vs Met-Score vs Combined)
#   - Violin plots with Wilcoxon, Cliff's delta, AUC
#   - Time-dependent AUC at 5y/10y
#
# PRIMARY ADJUSTED DURHAM HR IS OWNED ELSEWHERE:
#   The primary/multivariable adjusted Durham association (continuous
#   per-SD Met-Score adjusted for the common comparator: pathological
#   Gleason grade group + log2(PSA+1) + pT stage) is fitted by the
#   canonical script
#     code/survival_analysis/Met_PCa_Survival_Multivariate.R
#   which consumes the clin_valid object saved at the end of THIS script
#   and reports the Durham adjusted per-SD Met-Score HR anchor
#   (1.60, 95% CI 1.19-2.14). The adjusted Cox in section 12 below adjusts
#   only for Gleason category (binary MetScoreClass + PathGleason) on the
#   GS 7-9 subset; it does NOT include log2(PSA+1), pT, or the continuous
#   Met-Score, so it is a restricted-subset DESCRIPTIVE analysis and is NOT
#   the contract common comparator. Do not treat its HRs as the primary
#   adjusted Durham effect; defer to Met_PCa_Survival_Multivariate.R for
#   that. The section-12 table/forest are retained only as the descriptive
#   Gleason-adjusted OS/PCSM rows already cited in the manuscript.
############################################################################
rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pROC)
  library(caret)
  library(mltools)
  library(survival)
  library(prodlim)
  library(pec)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  library(timeROC)
  library(survcomp)
  library(ggpubr)
  library(limma)
  library(effsize)
  library(tidyverse)
  library(forcats)
})

out_dir <- "./output/Durham"
fig_dir <- "./figures/Durham"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("./figures/Durham/time_dependent", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Met-Score model, gene definitions, and TRAINING data
# ============================================================
## Load the frozen Met-Score coefficients and decision threshold from config/.
source("./code/utils/locked_metscore.R")
.locked_model    <- load_locked_metscore()
LAMBDA_LOCKED    <- .locked_model$lambda           # reported only; scoring uses beta
LOCKED_THRESHOLD <- .locked_model$threshold

## ---- JHU-frozen Gleason + Met-Score combined predictors --------------------
## Combined = raw pathological Gleason sum + design-weighted-standardized
## Met-Score, developed only in the JHU case-cohort (Ross et al.) and frozen for
## unchanged application to Durham. MFS/PCSM use a phase-two-weighted Fine-Gray
## subdistribution model (competing event: death before metastasis, non-PC
## death); BCR/OS use inverse-probability-weighted Cox. Same design weights as
## the DCA producer (controls 745/265, cases 1). No Durham fitting.
jhu_frozen_combined <- function() {
  je <- new.env(); load("./outs/coxdata.rda", envir = je); j <- get("CoxData_jhu", envir = je)
  ALPHA <- 265 / 745
  cch <- as.character(j[["post_rp_patients_cchdef"]])
  w <- ifelse(cch == "Sub-cohort controls", 1 / ALPHA, 1)
  stopifnot(nrow(j) == 239L, sum(cch == "Sub-cohort cases") == 28L,
            sum(cch == "Sub-cohort controls") == 146L, sum(cch == "cases") == 65L)
  gl <- as.numeric(as.character(j[["Pathological GS"]]))
  ms <- as.numeric(j[["Met-Score prob"]])
  mean_JHU <- sum(w * ms) / sum(w); sd_JHU <- sqrt(sum(w * (ms - mean_JHU)^2) / sum(w))
  ms_dev <- (ms - mean_JHU) / sd_JHU
  id <- seq_len(nrow(j))
  met <- as.integer(j$met); os <- as.integer(j$os)
  mt <- as.numeric(j$met_time); ot <- as.numeric(j$os_time)
  mfs_evt  <- ifelse(met == 1L, 1L, ifelse(os == 1L & met == 0L & !is.na(ot) & ot <= mt, 2L, 0L))
  mfs_t    <- ifelse(mfs_evt == 2L, ot, mt)
  pcsm <- as.integer(j$pcsm)
  pcsm_evt <- ifelse(pcsm == 1L, 1L, ifelse(os == 1L & pcsm == 0L, 2L, 0L))
  bcr <- as.integer(j$bcr); bt <- as.numeric(j$bcr_time)
  D <- data.frame(gl = gl, ms_dev = ms_dev, w = w, id = id,
                  mfs_evt = mfs_evt, mfs_t = mfs_t, pcsm_evt = pcsm_evt, pcsm_t = ot,
                  bcr = bcr, bt = bt, os = os, ot = ot, stringsAsFactors = FALSE)
  Dcc <- D[!is.na(D$gl) & is.finite(D$ms_dev), , drop = FALSE]  # one missing Gleason
  evtf <- function(e) factor(e, levels = c(0, 1, 2), labels = c("cens", "event", "comp"))
  fit_fg <- function(df, time, evt) {
    d <- data.frame(time = time, evtf = evtf(evt), gl = df$gl, ms_dev = df$ms_dev, w = df$w, id = df$id)
    d <- d[d$time > 0, ]
    fg <- survival::finegray(Surv(time, evtf) ~ ., data = d, etype = "event", weights = w, id = id)
    coef(survival::coxph(Surv(fgstart, fgstop, fgstatus) ~ gl + ms_dev, data = fg,
                         weights = fgwt, cluster = id, robust = TRUE, ties = "breslow"))
  }
  fit_cox <- function(df, time, ev) {
    d <- data.frame(time = time, ev = ev, gl = df$gl, ms_dev = df$ms_dev, w = df$w, id = df$id)
    d <- d[d$time > 0, ]
    coef(survival::coxph(Surv(time, ev) ~ gl + ms_dev, data = d, weights = w,
                         cluster = id, robust = TRUE, ties = "breslow"))
  }
  list(beta = list(MFS = fit_fg(Dcc, Dcc$mfs_t, Dcc$mfs_evt),
                   PCSM = fit_fg(Dcc, Dcc$pcsm_t, Dcc$pcsm_evt),
                   BCR = fit_cox(Dcc, Dcc$bt, Dcc$bcr),
                   OS = fit_cox(Dcc, Dcc$ot, Dcc$os)),
       mean_JHU = mean_JHU, sd_JHU = sd_JHU,
       n_model = nrow(Dcc), n_cohort = nrow(D), n_missing_gs = sum(is.na(D$gl)))
}

## ---- Figure S2 aggregate producer (shared by the full run and --figure-s2-only) ----
## Writes five identifier-free CSVs: PCSM Aalen-Johansen CIF curves by locked
## class (cmprsk::cuminc gives the AJ CIF and Gray's test), a PCSM summary with
## the fixed 10-year high-minus-low CIF difference and a full-size bootstrap CI,
## high-vs-low secondary-endpoint Cox summaries (model-based + Lin-Wei robust),
## fixed-marker MFS time-dependent AUC, and fixed-marker endpoint C-index.
write_figure_s2_aggregates <- function(cv, out_dir, threshold,
                                       boot_B = 2000L, boot_seed = 20260812L) {
  req <- c("limbo", "fu", "dead", "deadofpc", "recurrence", "MetScoreClass", "MetScore_prob")
  stopifnot(all(req %in% colnames(cv)), nrow(cv) == 555L,
            all(is.finite(cv$MetScore_prob)),
            all(is.finite(cv$limbo)), all(cv$limbo > 0),
            all(is.finite(cv$fu)),    all(cv$fu > 0),
            all(cv$dead %in% c(0, 1)), all(cv$deadofpc %in% c(0, 1)),
            all(cv$recurrence %in% c(0, 1, 2)),
            !any(cv$deadofpc == 1 & cv$dead != 1))
  cls <- ifelse(cv$MetScore_prob >= threshold, "High risk", "Low risk")
  stopifnot(identical(cls, as.character(cv$MetScoreClass)))
  grp <- factor(cv$MetScoreClass, levels = c("Low risk", "High risk"))
  stopifnot(sum(grp == "Low risk") == 223L, sum(grp == "High risk") == 332L)
  # competing-risk status: 1 = PC death, 2 = non-PC death, 0 = censored
  st <- ifelse(cv$deadofpc == 1, 1L, ifelse(cv$dead == 1 & cv$deadofpc == 0, 2L, 0L))
  stopifnot(sum(st == 1) == 18L, sum(st == 2) == 175L, all(st %in% c(0, 1, 2)))
  gch <- as.character(grp); tt <- cv$limbo

  step_at <- function(el, t) { i <- findInterval(t, el$time); ifelse(i == 0, 0, el$est[pmax(i, 1)]) }
  var_at  <- function(el, t) { i <- findInterval(t, el$time); ifelse(i == 0, 0, el$var[pmax(i, 1)]) }

  ci <- cmprsk::cuminc(tt, st, gch, cencode = 0)
  gray_p <- ci$Tests["1", "pv"]
  months <- 0:120
  curve <- do.call(rbind, lapply(c("Low risk", "High risk"), function(cl) {
    el <- ci[[paste(cl, "1")]]; cif <- step_at(el, months); v <- var_at(el, months)
    data.frame(class = cl, month = months, cif = round(cif, 7),
               cif_lo = round(pmax(0, cif - 1.96 * sqrt(v)), 7),
               cif_hi = round(pmin(1, cif + 1.96 * sqrt(v)), 7),
               n_risk = vapply(months, function(m) sum(tt[gch == cl] >= m), integer(1)),
               stringsAsFactors = FALSE)
  }))
  write.csv(curve, file.path(out_dir, "FigureS2_Durham_PCSM_CIF_curve.csv"), row.names = FALSE)

  cif10 <- c(Low = step_at(ci[["Low risk 1"]], 120), High = step_at(ci[["High risk 1"]], 120))
  cif5  <- c(Low = step_at(ci[["Low risk 1"]], 60),  High = step_at(ci[["High risk 1"]], 60))
  diff10 <- unname(cif10["High"] - cif10["Low"])

  # full-size patient bootstrap for the fixed 10-year high-minus-low CIF difference;
  # a present zero-PC-death group is valid (CIF 0), only genuine failures are dropped
  set.seed(boot_seed)
  n <- nrow(cv); diffs <- numeric(0)
  att <- 0L; f_absent <- 0L; f_fit <- 0L; f_nonfin <- 0L
  for (b in seq_len(boot_B)) {
    att <- att + 1L
    idx <- sample.int(n, n, replace = TRUE); g <- gch[idx]
    if (length(unique(g)) < 2L) { f_absent <- f_absent + 1L; next }
    r <- tryCatch(cmprsk::cuminc(tt[idx], st[idx], g, cencode = 0), error = function(e) NULL)
    if (is.null(r) || is.null(r[["High risk 1"]]) || is.null(r[["Low risk 1"]])) { f_fit <- f_fit + 1L; next }
    d <- step_at(r[["High risk 1"]], 120) - step_at(r[["Low risk 1"]], 120)
    if (!is.finite(d)) { f_nonfin <- f_nonfin + 1L; next }
    diffs <- c(diffs, d)
  }
  used <- length(diffs)
  qb <- if (used >= 50L) unname(quantile(diffs, c(0.025, 0.975))) else c(NA_real_, NA_real_)

  summ <- data.frame(
    metric = c("n_total", "pc_deaths", "competing_deaths", "n_low", "n_high", "gray_p",
               "cif_5y_low", "cif_5y_high", "cif_10y_low", "cif_10y_high",
               "diff_10y_high_minus_low", "diff_10y_boot_lo", "diff_10y_boot_hi",
               "boot_seed", "boot_attempted", "boot_used", "boot_failed",
               "boot_fail_absent_group", "boot_fail_fit_error", "boot_fail_nonfinite"),
    value = c(n, sum(st == 1), sum(st == 2), sum(grp == "Low risk"), sum(grp == "High risk"),
              signif(gray_p, 7), signif(cif5["Low"], 7), signif(cif5["High"], 7),
              signif(cif10["Low"], 7), signif(cif10["High"], 7), signif(diff10, 7),
              signif(qb[1], 7), signif(qb[2], 7), boot_seed, att, used, att - used,
              f_absent, f_fit, f_nonfin),
    stringsAsFactors = FALSE)
  write.csv(summ, file.path(out_dir, "FigureS2_Durham_PCSM_summary.csv"), row.names = FALSE)

  # secondary endpoints: high vs low Cox, model-based + Lin-Wei robust + cox.zph
  zc <- stats::qnorm(0.975)
  cox_row <- function(time, event, endpoint, tf, ef, estimand) {
    f  <- survival::coxph(survival::Surv(time, event) ~ grp, robust = TRUE)
    co <- summary(f)$coefficients
    beta <- co["grpHigh risk", "coef"]; se_m <- co["grpHigh risk", "se(coef)"]
    se_r <- co["grpHigh risk", "robust se"]
    zph <- tryCatch(survival::cox.zph(f)$table["GLOBAL", "p"], error = function(e) NA_real_)
    data.frame(endpoint = endpoint, time_field = tf, event_field = ef, estimand = estimand,
               n = f$n, events = f$nevent, HR = round(exp(beta), 5),
               se_model = round(se_m, 6), ci_lo_model = round(exp(beta - zc * se_m), 5),
               ci_hi_model = round(exp(beta + zc * se_m), 5), p_model = signif(2 * pnorm(-abs(beta / se_m)), 5),
               se_robust = round(se_r, 6), ci_lo_robust = round(exp(beta - zc * se_r), 5),
               ci_hi_robust = round(exp(beta + zc * se_r), 5), p_robust = signif(2 * pnorm(-abs(beta / se_r)), 5),
               cox_zph_p = signif(zph, 5), stringsAsFactors = FALSE)
  }
  hr <- rbind(
    cox_row(cv$fu,    as.integer(cv$recurrence >= 1), "BCR",  "fu",    "recurrence>=1", "cause-specific hazard, high vs low"),
    cox_row(cv$limbo, cv$dead,                        "OS",   "limbo", "dead",          "all-cause hazard, high vs low"),
    cox_row(cv$limbo, cv$deadofpc,                    "PCSM", "limbo", "deadofpc",      "cause-specific hazard, non-PC death censored, high vs low"))
  write.csv(hr, file.path(out_dir, "FigureS2_Durham_secondary_HR.csv"), row.names = FALSE)

  # JHU-frozen Gleason+Met-Score combined predictors, applied unchanged to Durham
  JC <- jhu_frozen_combined()
  dgl <- as.numeric(cv$PathGleason)
  dms_dev <- (as.numeric(cv$MetScore_prob) - JC$mean_JHU) / JC$sd_JHU
  LP <- lapply(JC$beta, function(b) b[["gl"]] * dgl + b[["ms_dev"]] * dms_dev)
  cat(sprintf("JHU-frozen combined: model-complete n=%d (%d missing Gleason of %d cohort); ms mean=%.5f sd=%.5f\n",
              JC$n_model, JC$n_missing_gs, JC$n_cohort, JC$mean_JHU, JC$sd_JHU))

  # ---- panel c: fixed-marker MFS time-dependent AUC (def-2 timeROC, cause 1) ----
  # MFS competing-event table matches the main Figure-2 producer: metastasis =
  # surgmets/mets, death at limbo before metastasis is the competing event.
  mfs_status <- ifelse(cv$mets == 1L, 1L,
                       ifelse(cv$dead == 1L & is.finite(cv$limbo) & cv$limbo <= cv$surgmets, 2L, 0L))
  mfs_time <- ifelse(mfs_status == 2L, cv$limbo, cv$surgmets)
  stopifnot(sum(mfs_status == 1L) == 40L, sum(mfs_status == 2L) == 167L,
            all(is.finite(mfs_time)), all(mfs_time > 0))
  td_rows <- function(marker, mname) {
    keep <- is.finite(marker) & is.finite(mfs_time) & mfs_time > 0 & !is.na(mfs_status)
    m <- marker[keep]; T <- mfs_time[keep]; D <- mfs_status[keep]
    do.call(rbind, lapply(c(60, 120), function(t0) {
      tr <- timeROC::timeROC(T = T, delta = D, marker = m, cause = 1, times = t0, iid = TRUE)
      nm <- paste0("t=", t0)
      if (!is.null(tr$AUC_2)) { aucv <- tr$AUC_2; cim <- confint(tr)$CI_AUC_2 }
      else                    { aucv <- tr$AUC;   cim <- confint(tr)$CI_AUC }
      cir <- if (!is.null(rownames(cim)) && nm %in% rownames(cim)) cim[nm, ] else cim[nrow(cim), ]
      cd <- sum(D == 2L & T <= t0)
      data.frame(marker = mname, horizon_months = t0, n = length(m),
                 cases = sum(D == 1L & T <= t0), competing_deaths = cd,
                 auc = round(as.numeric(aucv[nm]), 7),
                 ci_lo = round(as.numeric(cir[1]) / 100, 7), ci_hi = round(as.numeric(cir[2]) / 100, 7),
                 control_def = if (cd > 0) "definition 2 (competing death as control)" else "standard",
                 n_missing = sum(!keep), stringsAsFactors = FALSE)
    }))
  }
  # three-series panel: Gleason and frozen Met-Score alone, plus the JHU-frozen
  # MFS combined predictor (Gleason + Met-Score) applied unchanged to Durham
  mfsauc <- rbind(td_rows(cv$PathGleason, "Gleason"),
                  td_rows(cv$MetScore_prob, "Met-Score"),
                  td_rows(LP$MFS, "Combined"))
  write.csv(mfsauc, file.path(out_dir, "FigureS2_Durham_MFS_timeAUC.csv"), row.names = FALSE)

  # ---- panel d: fixed-marker 10-year time-truncated IPCW concordance (pec::cindex) ----
  # No model fitting: raw PathGleason and frozen MetScore_prob rankings only, marginal
  # IPCW censoring, eval.times = pred.times = 120. Competing-risk endpoints use the
  # cause-1 cumulative-incidence interface (higher marker = higher risk); right-censored
  # endpoints use the survival-probability interface, so the marker is negated to keep
  # higher marker = higher risk. Full-size paired bootstrap recomputes IPCW per resample.
  gl <- as.numeric(cv$PathGleason); ms <- as.numeric(cv$MetScore_prob)
  eps <- list(
    MFS  = list(mode = "cr",   time = mfs_time, status = mfs_status, comb = LP$MFS,
                basis = "competing risk cause=1 (metastasis; death before metastasis event 2); surgmets/mets, limbo"),
    BCR  = list(mode = "surv", time = cv$fu,    status = as.integer(cv$recurrence >= 1), comb = LP$BCR,
                basis = "right-censored fu/recurrence>=1"),
    OS   = list(mode = "surv", time = cv$limbo, status = as.integer(cv$dead), comb = LP$OS,
                basis = "right-censored limbo/dead"),
    PCSM = list(mode = "cr",   time = cv$limbo, status = st, comb = LP$PCSM,
                basis = "competing risk cause=1 (PC death; non-PC death event 2); limbo"))
  # one call returns all three markers on the same data (paired), recomputing IPCW
  # weights; the JHU-frozen combined predictor enters as a fixed per-patient marker
  cindex_tri <- function(mode, time, status, glm, msm, cbm) {
    d <- data.frame(time = time, status = as.integer(status))
    sgn <- if (mode == "surv") -1 else 1
    obj <- list(g = matrix(sgn * glm, ncol = 1), m = matrix(sgn * msm, ncol = 1), c = matrix(sgn * cbm, ncol = 1))
    fm <- if (mode == "cr") Hist(time, status) ~ 1 else Surv(time, status) ~ 1
    a <- list(object = obj, formula = fm, data = d, eval.times = 120, pred.times = 120,
              cens.model = "marginal", verbose = FALSE)
    if (mode == "cr") a$cause <- 1
    r <- do.call(pec::cindex, a)
    c(as.numeric(r$AppCindex$g), as.numeric(r$AppCindex$m), as.numeric(r$AppCindex$c))
  }
  cidx_seed <- 20260813L; nn <- nrow(cv)
  set.seed(cidx_seed)
  boot_idx <- lapply(seq_len(boot_B), function(b) sample.int(nn, nn, replace = TRUE))
  crows <- list()
  for (en in names(eps)) {
    ep <- eps[[en]]
    pt <- cindex_tri(ep$mode, ep$time, ep$status, gl, ms, ep$comb)
    bg <- numeric(0); bm <- numeric(0); bc <- numeric(0); a2 <- 0L; f_fit <- 0L; f_nonfin <- 0L
    for (b in seq_len(boot_B)) {
      a2 <- a2 + 1L; ix <- boot_idx[[b]]
      v <- tryCatch(cindex_tri(ep$mode, ep$time[ix], ep$status[ix], gl[ix], ms[ix], ep$comb[ix]),
                    error = function(e) NULL)
      if (is.null(v)) { f_fit <- f_fit + 1L; next }
      if (!all(is.finite(v))) { f_nonfin <- f_nonfin + 1L; next }
      bg <- c(bg, v[1]); bm <- c(bm, v[2]); bc <- c(bc, v[3])
    }
    ub <- length(bg)
    q <- function(x) if (ub >= 50L) unname(quantile(x, c(0.025, 0.975))) else c(NA_real_, NA_real_)
    qg <- q(bg); qm <- q(bm); qc <- q(bc)
    mkrow <- function(marker, cpt, qq) data.frame(
      endpoint = en, marker = marker, horizon_months = 120L, endpoint_basis = ep$basis,
      estimator = "pec::cindex IPCW 10y time-truncated concordance", censoring_model = "marginal",
      c_index = round(cpt, 7), ci_lo = round(qq[1], 7), ci_hi = round(qq[2], 7),
      boot_seed = cidx_seed, boot_attempted = a2, boot_used = ub, boot_failed = a2 - ub,
      boot_fail_fit_error = f_fit, boot_fail_nonfinite = f_nonfin, stringsAsFactors = FALSE)
    crows[[paste(en, "Gleason")]]   <- mkrow("Gleason",   pt[1], qg)
    crows[[paste(en, "Met-Score")]] <- mkrow("Met-Score", pt[2], qm)
    crows[[paste(en, "Combined")]]  <- mkrow("Combined",  pt[3], qc)
  }
  cidx_tab <- do.call(rbind, crows)
  write.csv(cidx_tab, file.path(out_dir, "FigureS2_Durham_fixed_marker_Cindex.csv"), row.names = FALSE)

  cat(sprintf("Figure S2: PCSM 10y CIF low=%.5f high=%.5f Gray p=%.5f; boot used %d/%d; HR BCR/OS/PCSM=%.3f/%.3f/%.3f\n",
              cif10["Low"], cif10["High"], gray_p, used, att, hr$HR[1], hr$HR[2], hr$HR[3]))
  cat(sprintf("Figure S2 panel c AUC Gleason/Met-Score/Combined: 5y=%.5f/%.5f/%.5f 10y=%.5f/%.5f/%.5f\n",
              mfsauc$auc[1], mfsauc$auc[3], mfsauc$auc[5], mfsauc$auc[2], mfsauc$auc[4], mfsauc$auc[6]))
  cat("Figure S2 panel d IPCW 10y concordance Gleason/Met-Score/Combined:\n")
  for (en in names(eps)) { rr <- cidx_tab[cidx_tab$endpoint == en, ]
    cat(sprintf("  %-4s %.5f/%.5f/%.5f\n", en, rr$c_index[1], rr$c_index[2], rr$c_index[3])) }
  invisible(TRUE)
}

## narrow mode: emit only the Figure-S2 aggregates from the saved full-cohort RDA
if ("--figure-s2-only" %in% commandArgs(TRUE)) {
  rda <- file.path(out_dir, "durham_metscore_batchcorrected.rda")
  stopifnot(file.exists(rda))
  .s2 <- new.env(); load(rda, envir = .s2)
  stopifnot("clin_valid" %in% ls(.s2))
  write_figure_s2_aggregates(.s2$clin_valid, out_dir, LOCKED_THRESHOLD)
  cat("--figure-s2-only: wrote Figure-S2 aggregates; no other output regenerated\n")
  quit(save = "no", status = 0)
}

load("./outs/PP_filter_MetaScore.rda")
load("./outs/filtersiggenes_MetaScore.rda")
load("./outs/MetastasisData_JHUOut.rda")  # trainMat, testMat, trainGroup, testGroup

pos_genes <- filter$posGeneNames   # 27 positive genes
neg_genes <- filter$negGeneNames   # 18 negative genes
all_sig_genes <- Filter_SignatureGenes
# The 41 deployed model genes come from the locked contract (frozen order).
model_genes <- .locked_model$feature_names

cat("=== Met-Score Gene Signature ===\n")
cat("Positive genes:", length(pos_genes), "\n")
cat("Negative genes:", length(neg_genes), "\n")
cat("Model genes:", length(model_genes), "\n")
cat(sprintf("Locked LR: ridge cv.glmnet, lambda = %.6f, threshold = %.6f\n\n",
            LAMBDA_LOCKED, LOCKED_THRESHOLD))

# Quantile normalize training data (as in original pipeline)
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")

# Build training data frame (same as original)
Data_train <- as.data.frame(cbind(t(usedTrainMat), trainGroup))
Data_train$trainGroup <- as.factor(Data_train$trainGroup)
levels(Data_train$trainGroup) <- c(0, 1)
colnames(Data_train)[colnames(Data_train) == "trainGroup"] <- "label"

# Score the complete training matrix; the helper selects the model genes.
train_probs <- as.numeric(locked_metscore_score(
  t(usedTrainMat), .locked_model)$prob)
train_thr <- LOCKED_THRESHOLD

cat("=== Training Data Distribution ===\n")
cat("Training prob range:", round(range(train_probs), 4), "\n")
cat("Training prob median:", round(median(train_probs), 4), "\n")
cat("Training % above threshold:", round(mean(train_probs >= train_thr) * 100, 1), "%\n\n")

# ============================================================
# 2. Load Durham expression data
# ============================================================
cat("Loading Durham expression data (this may take a few minutes)...\n")
durham_expr_raw <- read_excel(
  "data/Durham_cohort_and_GRID_cohort/Durham_cohort_011526.xlsx",
  sheet = "eset_gene_filtered"
)
cat("Loaded:", nrow(durham_expr_raw), "genes x", ncol(durham_expr_raw) - 2, "samples\n")

expr_mat <- as.data.frame(durham_expr_raw)
gene_symbols <- expr_mat$Symbol
expr_mat <- expr_mat[, -c(1, 2)]
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- gene_symbols

# ============================================================
# 3. Load clinical data & link to expression
# ============================================================
cat("Loading clinical data...\n")
clin <- read_excel("data/Durham_cohort_and_GRID_cohort/Durham_cohort_clinical_data_022526.xlsx",
                   sheet = "clin")
sample_map <- read_excel("data/Durham_cohort_and_GRID_cohort/Durham_cohort_011526.xlsx",
                         sheet = "Sheet2")
colnames(sample_map) <- c("cell_file_name", "sample_id")

clin_mapped <- merge(clin, sample_map, by = "cell_file_name")
clin_valid <- clin_mapped[!is.na(clin_mapped$mets), ]
cat("Samples with valid clinical data:", nrow(clin_valid), "\n\n")

# ----- Remove placeholder/sentinel clinical records --------------------------
# In the raw Durham clinical Excel, 3 rows (External_IDs GDX371, GDX909,
# GDX922) have pogl == 0 alongside EVERY other clinical field set to 0/0.0
# (race=0, age=0, dead=0, ece=0, recurrence=0, fu=0, surgmets=0, mets=0,
# BMI=0, weight=0, height=0, psapresurg=0). These rows are placeholder /
# sentinel records — `0` here means "unknown / not entered" rather than
# literal values, and they do NOT represent valid Gleason or survival
# information. Removing them here is a data-quality clean-up of placeholder
# records, NOT a Gleason-based biological exclusion. They survive the
# !is.na(mets) filter above because 0 is non-NA.
n_before_sentinel <- nrow(clin_valid)
clin_valid <- clin_valid[!(is.na(clin_valid$pogl) | clin_valid$pogl == 0), ]
n_dropped_sentinel <- n_before_sentinel - nrow(clin_valid)
cat(sprintf("Removed %d Durham sentinel clinical rows with pogl == 0.\n",
            n_dropped_sentinel))
cat("Cohort size after sentinel removal:", nrow(clin_valid),
    "(was", n_before_sentinel, ")\n\n")

# Derive BCR binary
clin_valid$bcr_binary <- ifelse(clin_valid$recurrence >= 1, 1, 0)

# Map Gleason: pogl (pathological overall gleason).
# Stored as integer here; the GS 7-9 restriction is applied ONLY inside
# the multivariate-Cox block below (where the model needs it). KM curves,
# univariate HRs, time-dependent AUC, C-index, and all violin/ROC panels
# run on the full Durham cohort, so Figure 2 panel d remains the
# headline prognostic display and is internally consistent with the
# published Met-Score curves.
clin_valid$PathGleason <- as.integer(clin_valid$pogl)
cat("Pathological Gleason distribution (full cohort, used everywhere\n",
    "EXCEPT the multivariate Cox in section 12):\n", sep = "")
print(table(clin_valid$PathGleason, useNA = "ifany"))
cat("\n")

# ============================================================
# 4. APPROACH 1: Gene-Set Z-Score (Model-Free, Platform-Agnostic)
# ============================================================
cat("================================================================\n")
cat("=== APPROACH 1: Gene-Set Z-Score (Model-Free) ===\n")
cat("================================================================\n\n")

pos_in_durham <- pos_genes[pos_genes %in% rownames(expr_mat)]
neg_in_durham <- neg_genes[neg_genes %in% rownames(expr_mat)]
cat("Positive genes present:", length(pos_in_durham), "/", length(pos_genes), "\n")
cat("Negative genes present:", length(neg_in_durham), "/", length(neg_genes), "\n")

durham_samples <- clin_valid$sample_id
sig_genes_present <- c(pos_in_durham, neg_in_durham)
expr_sig <- expr_mat[sig_genes_present, durham_samples]

# Rank-based Z-score transformation within Durham (no outcome info used)
expr_z <- t(apply(expr_sig, 1, function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(r))
}))

# Compute gene-set score: mean(pos Z) - mean(neg Z)
pos_z_mean <- colMeans(expr_z[pos_in_durham, , drop = FALSE])
neg_z_mean <- colMeans(expr_z[neg_in_durham, , drop = FALSE])
geneset_score <- pos_z_mean - neg_z_mean

# Median split for classification (unsupervised)
geneset_median <- median(geneset_score)
geneset_class <- ifelse(geneset_score >= geneset_median, "High", "Low")

clin_valid$GeneSet_score <- geneset_score[match(clin_valid$sample_id, names(geneset_score))]
clin_valid$GeneSet_class <- factor(geneset_class[match(clin_valid$sample_id, names(geneset_class))],
                                    levels = c("Low", "High"))
clin_valid$GeneSet_numeric <- ifelse(clin_valid$GeneSet_class == "High", 1, 0)

cat("\nGene-Set Score distribution:\n")
cat("  Range:", round(range(clin_valid$GeneSet_score), 4), "\n")
cat("  Median (threshold):", round(geneset_median, 4), "\n")
cat("  High:", sum(clin_valid$GeneSet_class == "High"),
    "| Low:", sum(clin_valid$GeneSet_class == "Low"), "\n\n")

# ============================================================
# 5. APPROACH 2: Quantile-Normalized Logistic Regression
# ============================================================
cat("================================================================\n")
cat("=== APPROACH 2: Quantile-Normalized Logistic Regression ===\n")
cat("================================================================\n\n")

# Stop if a required gene has duplicate rows in the Durham expression.
.dup_req <- model_genes[vapply(model_genes,
                               function(g) sum(rownames(expr_mat) == g) > 1L,
                               logical(1))]
if (length(.dup_req) > 0)
  stop(sprintf(paste0("Durham validation: required locked gene(s) appear >1x ",
                      "in Durham expression row names: %s"),
               paste(unique(.dup_req), collapse = ", ")))

genes_in_both <- model_genes[model_genes %in% rownames(expr_mat)]
genes_missing <- model_genes[!model_genes %in% rownames(expr_mat)]
cat("Model genes in Durham:", length(genes_in_both), "/", length(model_genes), "\n")
if (length(genes_missing) > 0) cat("Missing:", paste(genes_missing, collapse = ", "), "\n")

durham_model_expr <- expr_mat[genes_in_both, durham_samples]
train_model_expr <- usedTrainMat[genes_in_both, ]

# Per-gene quantile normalization: map Durham to training distribution
cat("Performing per-gene quantile normalization to training reference...\n")

durham_norm <- matrix(NA, nrow = nrow(durham_model_expr), ncol = ncol(durham_model_expr),
                      dimnames = dimnames(durham_model_expr))

for (g in rownames(durham_model_expr)) {
  train_vals <- train_model_expr[g, ]
  durham_vals <- durham_model_expr[g, ]
  durham_ranks <- rank(durham_vals, ties.method = "average") / (length(durham_vals) + 1)
  durham_norm[g, ] <- quantile(train_vals, probs = durham_ranks)
}

expr_df_norm <- as.data.frame(t(durham_norm))
# For any model gene absent from the Durham platform, impute with the
# per-gene mean of the training-reference distribution (usedTrainMat)
# rather than 0. Limma between-arrays QN does NOT zero-centre genes —
# imputing 0 would shift the linear predictor by Σ_missing β_g · (0 −
# x̄_g_train), which is a small but non-zero bias on every Durham
# patient. Imputing the gene mean leaves the linear-predictor expectation
# unchanged for missing genes, which is the calibrated convention.
for (g in genes_missing) expr_df_norm[[g]] <- mean(usedTrainMat[g, ])

# Score the complete bridged matrix; the helper selects the model genes.
durham_prob_norm <- as.numeric(
  locked_metscore_score(as.matrix(expr_df_norm), .locked_model)$prob)
names(durham_prob_norm) <- rownames(expr_df_norm)

cat("Normalized prob range:", round(range(durham_prob_norm), 4), "\n")
cat("Normalized prob median:", round(median(durham_prob_norm), 4), "\n")

# Use training threshold
durham_pred_norm <- ifelse(durham_prob_norm >= train_thr, 1, 0)
cat("With training threshold (", train_thr, "):\n")
cat("  High:", sum(durham_pred_norm == 1), "| Low:", sum(durham_pred_norm == 0), "\n")
cat("  % High:", round(mean(durham_pred_norm) * 100, 1), "%\n\n")

clin_valid$MetScore_prob <- durham_prob_norm[match(clin_valid$sample_id, names(durham_prob_norm))]
clin_valid$MetScore_class <- factor(durham_pred_norm[match(clin_valid$sample_id, names(durham_pred_norm))],
                                     levels = c(0, 1), labels = c("Low", "High"))
clin_valid$MetScore_numeric <- ifelse(clin_valid$MetScore_class == "High", 1, 0)

# Also build MetScoreClass matching JHU naming convention
clin_valid$MetScoreClass <- factor(
  clin_valid$MetScore_class,
  levels = c("Low", "High"),
  labels = c("Low risk", "High risk")
)

# ============================================================
# 6. Also keep raw (uncorrected) predictions for comparison
# ============================================================
expr_raw <- t(expr_mat[genes_in_both, durham_samples])
expr_df_raw <- as.data.frame(expr_raw)
# Same mean-imputation rationale as the QN'd path above (line ~213):
# limma between-arrays QN does not zero-centre genes, so imputing 0
# would shift the linear predictor by a small bias on every patient.
for (g in genes_missing) expr_df_raw[[g]] <- mean(usedTrainMat[g, ])
# Score the complete raw matrix; the helper selects the model genes.
durham_prob_raw <- as.numeric(
  locked_metscore_score(as.matrix(expr_df_raw), .locked_model)$prob)
names(durham_prob_raw) <- rownames(expr_df_raw)
durham_pred_raw <- ifelse(durham_prob_raw >= train_thr, 1, 0)

clin_valid$MetScore_prob_raw <- durham_prob_raw[match(clin_valid$sample_id, names(durham_prob_raw))]
clin_valid$MetScore_class_raw <- factor(durham_pred_raw[match(clin_valid$sample_id, names(durham_pred_raw))],
                                         levels = c(0, 1), labels = c("Low", "High"))

# ============================================================
# 7. CLASSIFICATION: Metastasis endpoint (matching JHU approach)
# ============================================================
cat("\n================================================================\n")
cat("=== CLASSIFICATION: Metastasis (matching JHU) ===\n")
cat("================================================================\n\n")

# --- QN-LogReg (primary) ---
cat("--- QN-LogReg Model ---\n")
ROC_durham_qn <- roc(clin_valid$mets, clin_valid$MetScore_prob,
                      levels = c(0, 1), direction = "<",
                      ci = TRUE, quiet = TRUE)
cat("AUC:", round(as.numeric(ROC_durham_qn$auc), 4),
    "(95% CI:", round(ROC_durham_qn$ci[1], 4), "-", round(ROC_durham_qn$ci[3], 4), ")\n")

pred_f <- factor(clin_valid$MetScore_numeric, levels = c(0, 1))
actual_f <- factor(clin_valid$mets, levels = c(0, 1))
Confusion_durham_qn <- confusionMatrix(pred_f, actual_f, positive = "1", mode = "everything")
cat(capture.output(Confusion_durham_qn), sep = "\n")

MCC_durham_qn <- mltools::mcc(pred = pred_f, actuals = actual_f)
cat("\nMCC:", round(MCC_durham_qn, 4), "\n\n")

# --- Gene-Set Z-Score ---
cat("--- Gene-Set Z-Score ---\n")
ROC_durham_gs <- roc(clin_valid$mets, clin_valid$GeneSet_score,
                      levels = c(0, 1), direction = "<",
                      ci = TRUE, quiet = TRUE)
cat("AUC:", round(as.numeric(ROC_durham_gs$auc), 4),
    "(95% CI:", round(ROC_durham_gs$ci[1], 4), "-", round(ROC_durham_gs$ci[3], 4), ")\n")

pred_gs_f <- factor(clin_valid$GeneSet_numeric, levels = c(0, 1))
Confusion_durham_gs <- confusionMatrix(pred_gs_f, actual_f, positive = "1", mode = "everything")
cat(capture.output(Confusion_durham_gs), sep = "\n")

MCC_durham_gs <- mltools::mcc(pred = pred_gs_f, actuals = actual_f)
cat("\nMCC:", round(MCC_durham_gs, 4), "\n\n")

# --- Raw (uncorrected) for comparison ---
cat("--- Raw LogReg (uncorrected) ---\n")
ROC_durham_raw <- roc(clin_valid$mets, clin_valid$MetScore_prob_raw,
                       levels = c(0, 1), direction = "<",
                       ci = TRUE, quiet = TRUE)
cat("AUC:", round(as.numeric(ROC_durham_raw$auc), 4),
    "(95% CI:", round(ROC_durham_raw$ci[1], 4), "-", round(ROC_durham_raw$ci[3], 4), ")\n")

pred_raw_f <- factor(ifelse(clin_valid$MetScore_class_raw == "High", 1, 0), levels = c(0, 1))
Confusion_durham_raw <- confusionMatrix(pred_raw_f, actual_f, positive = "1", mode = "everything")
MCC_durham_raw <- mltools::mcc(pred = pred_raw_f, actuals = actual_f)
cat("MCC:", round(MCC_durham_raw, 4), "\n\n")

# ============================================================
# 8. ROC for all phenotypes (table format)
# ============================================================
cat("================================================================\n")
cat("=== ROC for All Phenotypes (QN-LogReg + GeneSet) ===\n")
cat("================================================================\n\n")

phenotypes <- list(
  list(name = "Metastasis", var = "mets", time_var = "surgmets"),
  list(name = "BCR", var = "bcr_binary", time_var = "fu"),
  list(name = "All-cause Death", var = "dead", time_var = "limbo"),
  list(name = "PC-specific Death", var = "deadofpc", time_var = "limbo")
)

auc_table <- data.frame()
for (ph in phenotypes) {
  events <- sum(clin_valid[[ph$var]] == 1, na.rm = TRUE)
  for (score_name in c("QN-LogReg", "GeneSet Z-Score", "Raw (uncorrected)")) {
    prob_var <- switch(score_name,
                       "QN-LogReg" = "MetScore_prob",
                       "GeneSet Z-Score" = "GeneSet_score",
                       "Raw (uncorrected)" = "MetScore_prob_raw")
    roc_obj <- tryCatch({
      roc(clin_valid[[ph$var]], clin_valid[[prob_var]],
          levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
    }, error = function(e) NULL)
    if (!is.null(roc_obj)) {
      auc_table <- rbind(auc_table, data.frame(
        Phenotype = ph$name, Events = events, N = nrow(clin_valid),
        Approach = score_name,
        AUC = round(as.numeric(roc_obj$auc), 4),
        CI_lower = round(roc_obj$ci[1], 4),
        CI_upper = round(roc_obj$ci[3], 4),
        stringsAsFactors = FALSE
      ))
    }
  }
}

print(auc_table, row.names = FALSE)
write.csv(auc_table, file.path(out_dir, "MetScore_AUC_BatchCorrected.csv"), row.names = FALSE)

# ============================================================
# 9. ROC Figure (metastasis, all approaches)
# ============================================================
pdf(file.path(fig_dir, "ROC_Metastasis_AllApproaches.pdf"), width = 9, height = 7)
colors_app <- c("#999999", "#E41A1C", "#4DAF4A")
plot(ROC_durham_raw, col = colors_app[1], lwd = 2, lty = 2,
     main = "Met-Score ROC - Metastasis (Durham)\nComparison of Approaches",
     legacy.axes = TRUE)
plot(ROC_durham_qn, col = colors_app[2], lwd = 2, add = TRUE)
plot(ROC_durham_gs, col = colors_app[3], lwd = 2, add = TRUE)
legend("bottomright",
       legend = c(
         paste0("Raw LogReg (AUC=", round(pROC::auc(ROC_durham_raw), 3), ")"),
         paste0("QN-LogReg (AUC=", round(pROC::auc(ROC_durham_qn), 3), ")"),
         paste0("GeneSet Z-Score (AUC=", round(pROC::auc(ROC_durham_gs), 3), ")")
       ),
       col = colors_app, lwd = 2, lty = c(2, 1, 1), cex = 0.9)
dev.off()
cat("\nROC comparison figure saved.\n")

# ============================================================
# 10. Probability distribution comparison figure
# ============================================================
pdf(file.path(fig_dir, "Probability_Distribution_Comparison.pdf"), width = 10, height = 6)
par(mfrow = c(1, 3))
hist(clin_valid$MetScore_prob_raw, breaks = 30, col = "#CCCCCC", border = "grey40",
     main = "Raw LogReg Probabilities", xlab = "Probability", xlim = c(0, 1))
abline(v = train_thr, col = "red", lwd = 2, lty = 2)
legend("topright", "Training Thr", col = "red", lty = 2, lwd = 2, cex = 0.8)

hist(clin_valid$MetScore_prob, breaks = 30, col = "#FFD6D6", border = "#E41A1C",
     main = "QN-Corrected LogReg Probabilities", xlab = "Probability", xlim = c(0, 1))
abline(v = train_thr, col = "red", lwd = 2, lty = 2)
legend("topright", "Training Thr", col = "red", lty = 2, lwd = 2, cex = 0.8)

hist(clin_valid$GeneSet_score, breaks = 30, col = "#D6FFD6", border = "#4DAF4A",
     main = "Gene-Set Z-Score", xlab = "Score")
abline(v = geneset_median, col = "blue", lwd = 2, lty = 3)
legend("topright", "Median", col = "blue", lty = 3, lwd = 2, cex = 0.8)
dev.off()
cat("Distribution comparison figure saved.\n")

# ============================================================
# 11. SURVIVAL ANALYSIS (matching JHU approach exactly)
# ============================================================
cat("\n================================================================\n")
cat("=== SURVIVAL ANALYSIS (JHU-matched) ===\n")
cat("================================================================\n\n")

# --- Consistent KM styling (matching JHU) ---
cols_ms <- c("Low risk" = "#2B2EB5", "High risk" = "#ED6905")

km_theme <- theme_survminer(
  base_size = 22,
  font.main      = c(24, "bold"),
  font.x         = c(22, "plain", "black"),
  font.y         = c(22, "plain", "black"),
  font.tickslab  = c(20, "plain", "black"),
  font.legend    = c(22, "plain", "black")
) +
  theme(
    plot.margin     = margin(15, 25, 15, 15),
    legend.position = "bottom",
    legend.key.size = unit(0.9, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    panel.border = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  )

# --- KM: Metastasis-free survival (QN-LogReg) ---
Fit_durham_mfs <- survfit(Surv(surgmets, mets) ~ MetScoreClass, data = clin_valid)

p_mfs_durham <- ggsurvplot(
  Fit_durham_mfs,
  data = clin_valid,
  pval = TRUE,
  pval.size = 6.5,
  pval.method = TRUE,
  pval.coord = c(max(clin_valid$surgmets, na.rm = TRUE) * 0.02, 0.15),
  palette = cols_ms,
  legend.title = "",
  legend.labs = c("Low risk", "High risk"),
  risk.table = FALSE,
  censor.shape = "|",
  censor.size = 3.2,
  size = 1.8,
  xlab = "Time (months)",
  ylab = "Metastasis-free survival probability",
  ggtheme = km_theme
)

ggsave(file.path(fig_dir, "KM_MFS_QNLogReg.pdf"),
       p_mfs_durham$plot, width = 6.8, height = 6.8, useDingbats = FALSE)
ggsave(file.path(fig_dir, "KM_MFS_QNLogReg.tiff"),
       p_mfs_durham$plot, width = 6.8, height = 6.8,
       dpi = 600, compression = "lzw")

# --- KM: Overall survival (QN-LogReg) ---
Fit_durham_os <- survfit(Surv(limbo, dead) ~ MetScoreClass, data = clin_valid)

p_os_durham <- ggsurvplot(
  Fit_durham_os,
  data = clin_valid,
  pval = TRUE,
  pval.size = 6.5,
  pval.method = TRUE,
  pval.coord = c(max(clin_valid$limbo, na.rm = TRUE) * 0.02, 0.15),
  palette = cols_ms,
  legend.title = "",
  legend.labs = c("Low risk", "High risk"),
  risk.table = FALSE,
  censor.shape = "|",
  censor.size = 3.2,
  size = 1.8,
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  ggtheme = km_theme
)

ggsave(file.path(fig_dir, "KM_OS_QNLogReg.pdf"),
       p_os_durham$plot, width = 6.8, height = 6.8, useDingbats = FALSE)
ggsave(file.path(fig_dir, "KM_OS_QNLogReg.tiff"),
       p_os_durham$plot, width = 6.8, height = 6.8,
       dpi = 600, compression = "lzw")

# --- KM: BCR-free survival (QN-LogReg) ---
Fit_durham_bcr <- survfit(Surv(fu, bcr_binary) ~ MetScoreClass, data = clin_valid)

p_bcr_durham <- ggsurvplot(
  Fit_durham_bcr,
  data = clin_valid,
  pval = TRUE,
  pval.size = 6.5,
  pval.method = TRUE,
  pval.coord = c(max(clin_valid$fu, na.rm = TRUE) * 0.02, 0.15),
  palette = cols_ms,
  legend.title = "",
  legend.labs = c("Low risk", "High risk"),
  risk.table = FALSE,
  censor.shape = "|",
  censor.size = 3.2,
  size = 1.8,
  xlab = "Time (months)",
  ylab = "BCR-free survival probability",
  ggtheme = km_theme
)

ggsave(file.path(fig_dir, "KM_BCR_QNLogReg.pdf"),
       p_bcr_durham$plot, width = 6.8, height = 6.8, useDingbats = FALSE)

# --- KM: PC-specific survival (QN-LogReg) ---
Fit_durham_pcsm <- survfit(Surv(limbo, deadofpc) ~ MetScoreClass, data = clin_valid)

p_pcsm_durham <- ggsurvplot(
  Fit_durham_pcsm,
  data = clin_valid,
  pval = TRUE,
  pval.size = 6.5,
  pval.method = TRUE,
  pval.coord = c(max(clin_valid$limbo, na.rm = TRUE) * 0.02, 0.15),
  palette = cols_ms,
  legend.title = "",
  legend.labs = c("Low risk", "High risk"),
  risk.table = FALSE,
  censor.shape = "|",
  censor.size = 3.2,
  size = 1.8,
  xlab = "Time (months)",
  ylab = "PC-specific survival probability",
  ggtheme = km_theme
)

ggsave(file.path(fig_dir, "KM_PCSM_QNLogReg.pdf"),
       p_pcsm_durham$plot, width = 6.8, height = 6.8, useDingbats = FALSE)

# --- Also generate GeneSet Z-Score KMs ---
clin_valid$GeneSet_MetScoreClass <- factor(
  clin_valid$GeneSet_class,
  levels = c("Low", "High"),
  labels = c("Low risk", "High risk")
)

for (ph in phenotypes) {
  surv_obj <- Surv(clin_valid[[ph$time_var]], clin_valid[[ph$var]])
  fit <- survfit(surv_obj ~ GeneSet_MetScoreClass, data = clin_valid)
  fname <- paste0("KM_GeneSetZScore_", gsub("[- ]", "_", ph$name))

  p_gs <- ggsurvplot(
    fit, data = clin_valid,
    pval = TRUE, pval.size = 6.5, pval.method = TRUE,
    pval.coord = c(max(clin_valid[[ph$time_var]], na.rm = TRUE) * 0.02, 0.15),
    palette = cols_ms,
    legend.title = "", legend.labs = c("Low risk", "High risk"),
    risk.table = FALSE, censor.shape = "|", censor.size = 3.2, size = 1.8,
    xlab = "Time (months)",
    ylab = paste0(ph$name, "-free survival probability"),
    ggtheme = km_theme
  )
  ggsave(file.path(fig_dir, paste0(fname, ".pdf")),
         p_gs$plot, width = 6.8, height = 6.8, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(fname, ".tiff")),
         p_gs$plot, width = 6.8, height = 6.8, dpi = 600, compression = "lzw")
}
cat("All KM figures saved.\n")

# ============================================================
# 12. Cox regression: univariate (MetScoreClass, primary here) +
#     RESTRICTED-SUBSET, GLEASON-ONLY ADJUSTED Cox (descriptive only)
#
# The univariate MetScoreClass Cox, KM, C-index, and time-AUC are the
# associations this script legitimately owns. The "multivariate" block
# below is a Gleason-only adjusted, GS 7-9 restricted DESCRIPTIVE model
# (binary MetScoreClass + PathGleason). It intentionally omits log2(PSA+1),
# pT, and the continuous per-SD Met-Score, so it does NOT reproduce the
# Durham per-SD adjusted HR anchor (1.60, 95% CI 1.19-2.14) and is NOT the
# contract common comparator. The primary/multivariable adjusted Durham HR
# is fitted canonically in
#   code/survival_analysis/Met_PCa_Survival_Multivariate.R
# (GG + log2PSA + pT + continuous ms_z), which consumes clin_valid saved at
# the end of this script. The Table1/forest emitted here are the descriptive
# Gleason-adjusted OS/PCSM rows referenced in the manuscript fact-check, not
# the primary adjusted effect.
# ============================================================
cat("\n=== Cox Regression: univariate (MetScoreClass) + restricted-subset",
    "Gleason-only adjusted (descriptive; primary adjusted HR is in",
    "Met_PCa_Survival_Multivariate.R) ===\n\n")

# Helper: build forest data frame (matching JHU)
make_forest_df <- function(fit, cohort, endpoint) {
  sm <- summary(fit)
  coef_tbl <- as.data.frame(sm$coef)
  ci_tbl <- as.data.frame(sm$conf.int)[, c("lower .95", "upper .95")]
  df <- coef_tbl %>%
    transmute(
      Term  = rownames(coef_tbl),
      HR    = `exp(coef)`,
      lower = ci_tbl[, 1],
      upper = ci_tbl[, 2],
      p     = `Pr(>|z|)`
    ) %>%
    mutate(
      Cohort   = cohort,
      Endpoint = endpoint,
      label_hr = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper),
      label_p  = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
      label    = paste0("HR ", label_hr, "; p=", label_p)
    )
  df
}

# --- Univariate Cox for each endpoint ---
surv_endpoints <- list(
  list(name = "MFS", time = "surgmets", event = "mets"),
  list(name = "BCR", time = "fu", event = "bcr_binary"),
  list(name = "OS",  time = "limbo",   event = "dead"),
  list(name = "PCSM", time = "limbo",  event = "deadofpc")
)

surv_table <- data.frame()

for (ep in surv_endpoints) {
  surv_obj <- Surv(clin_valid[[ep$time]], clin_valid[[ep$event]])

  # Univariate Cox by classes
  cox_uni <- tryCatch(
    coxph(surv_obj ~ MetScoreClass, data = clin_valid),
    error = function(e) NULL,
    warning = function(w) tryCatch(coxph(surv_obj ~ MetScoreClass, data = clin_valid), error = function(e) NULL)
  )

  # C-index
  ci_obj <- tryCatch(
    concordance.index(clin_valid$MetScore_prob, clin_valid[[ep$time]],
                      clin_valid[[ep$event]], method = "noether"),
    error = function(e) NULL
  )

  if (!is.null(cox_uni)) {
    sm <- summary(cox_uni)
    hr <- exp(coef(cox_uni))[1]
    hr_ci <- exp(confint(cox_uni))
    cox_p <- sm$coefficients[1, 5]
    lr <- survdiff(surv_obj ~ MetScoreClass, data = clin_valid)
    lr_p <- 1 - pchisq(lr$chisq, 1)

    surv_table <- rbind(surv_table, data.frame(
      Endpoint = ep$name, Model = "Met-Score (univariate)",
      HR = round(hr, 3),
      HR_CI = paste0(round(hr_ci[1], 3), "-", round(hr_ci[2], 3)),
      Cox_p = format.pval(cox_p, digits = 3),
      LogRank_p = format.pval(lr_p, digits = 3),
      C_index = if (!is.null(ci_obj)) round(ci_obj$c.index, 4) else NA,
      C_CI = if (!is.null(ci_obj)) paste0(round(ci_obj$lower, 4), "-", round(ci_obj$upper, 4)) else NA,
      stringsAsFactors = FALSE
    ))
  }

  # Restricted-subset, Gleason-only ADJUSTED Cox (DESCRIPTIVE). This
  # adjusts the binary MetScoreClass for Gleason category only; it does
  # NOT include log2(PSA+1), pT, or the continuous per-SD Met-Score, so it
  # is not the contract common comparator and does not reproduce the
  # Durham per-SD adjusted HR anchor (1.60) -- that primary adjusted model
  # lives in code/survival_analysis/Met_PCa_Survival_Multivariate.R.
  # Restrict to GS 7-9 ONLY for this descriptive analysis: pogl=0/5/6 has
  # too few patients (3+7+61) to support stable category-level estimates
  # and produces a singular design matrix in the unpenalised Cox model.
  # GS 7 is set as the reference so the table aligns with the JHU rows.
  multi_df <- clin_valid[!is.na(clin_valid$PathGleason) &
                          clin_valid$PathGleason %in% c(7, 8, 9), ]
  multi_df$PathGleason <- factor(multi_df$PathGleason, levels = c(7, 8, 9))
  surv_obj_multi <- Surv(multi_df[[ep$time]], multi_df[[ep$event]])
  cox_multi <- tryCatch(
    coxph(surv_obj_multi ~ MetScoreClass + PathGleason, data = multi_df),
    error = function(e) NULL,
    warning = function(w) tryCatch(
      coxph(surv_obj_multi ~ MetScoreClass + PathGleason, data = multi_df),
      error = function(e) NULL)
  )

  if (!is.null(cox_multi)) {
    sm <- summary(cox_multi)
    cat("--- Restricted-subset Gleason-only adjusted Cox (descriptive;",
        "not the primary comparator):", ep$name, "---\n")
    print(sm)
    cat("\n")

    # The descriptive restricted-GS7-9 Gleason-adjusted rows are no longer
    # exported as a table; the console summary above is kept for diagnostics
    # and the per-endpoint forest below is retained. Durham secondary-endpoint
    # table rows are produced by Met_PCa_Survival_Multivariate.R.

    # Forest plot (matching JHU)
    tryCatch({
      fname <- paste0("Forest_Durham_", ep$name, "_MetScore_Gleason")
      tiff(file.path(fig_dir, paste0(fname, ".tiff")),
           width = 6500, height = 7500, res = 550)
      print(ggforest(cox_multi, fontsize = 1.7, main = "HR"))
      dev.off()
      pdf(file.path(fig_dir, paste0(fname, ".pdf")),
          width = 11.8, height = 13.6)
      print(ggforest(cox_multi, fontsize = 1.7, main = "HR"))
      dev.off()
    }, error = function(e) cat("  Forest plot skipped:", e$message, "\n"))
  }
}

cat("\n=== SURVIVAL SUMMARY TABLE ===\n")
print(surv_table, row.names = FALSE)
write.csv(surv_table, file.path(out_dir, "MetScore_Survival_BatchCorrected.csv"), row.names = FALSE)

# Durham secondary-endpoint table rows are produced canonically by
# Met_PCa_Survival_Multivariate.R (outs/TableS3_Durham_secondary_endpoints.csv);
# this script no longer emits a descriptive endpoint-association table.

# ============================================================
# 14. Violin plots (matching JHU approach)
# ============================================================
cat("\n=== VIOLIN PLOTS ===\n\n")

make_violin <- function(df, score_col, event_col, title = "") {
  if (event_col == "mets") {
    labels_vec <- c("No metastasis", "Metastasis")
  } else if (event_col == "bcr_binary") {
    labels_vec <- c("No recurrence", "Recurrence")
  } else if (event_col == "dead") {
    labels_vec <- c("Alive", "Dead")
  } else if (event_col == "deadofpc") {
    labels_vec <- c("Alive", "Dead of PCa")
  } else {
    labels_vec <- c("No event", "Event")
  }

  df[[event_col]] <- factor(df[[event_col]], levels = c(0, 1), labels = labels_vec)
  p_res <- wilcox.test(df[[score_col]] ~ df[[event_col]])
  p_txt <- paste0("p = ", format(p_res$p.value, digits = 3, scientific = TRUE))
  y_max <- max(df[[score_col]], na.rm = TRUE)
  y_min <- min(df[[score_col]], na.rm = TRUE)
  y_range <- y_max - y_min
  label_y <- y_max + 0.12 * y_range

  fill_vals <- c("#2b2eb5", "#ed6905")
  names(fill_vals) <- labels_vec

  ggplot(df, aes(x = .data[[event_col]], y = .data[[score_col]],
                 fill = .data[[event_col]])) +
    geom_violin(trim = FALSE, alpha = 0.8, width = 0.9, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1, color = "black") +
    geom_jitter(width = 0.08, size = 1.8, alpha = 0.7) +
    annotate("text", x = 1.5, y = label_y, label = p_txt, size = 7, fontface = "bold") +
    scale_fill_manual(values = fill_vals) +
    labs(x = "", y = "Met-Score (Probability)", title = title) +
    theme_classic(base_size = 22) +
    theme(
      axis.title.y = element_text(size = 22, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold"),
      axis.text.y = element_text(size = 20, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      legend.position = "none",
      plot.margin = margin(15, 25, 15, 15)
    ) +
    expand_limits(y = label_y + 0.1 * y_range)
}

compute_violin_stats <- function(df, score_col, event_col) {
  score <- df[[score_col]]
  event_factor <- df[[event_col]]
  keep <- !(is.na(score) | is.na(event_factor))
  score <- score[keep]; event_factor <- event_factor[keep]
  event_num <- as.numeric(event_factor) - 1

  p_val <- wilcox.test(score ~ event_num)$p.value
  cd <- effsize::cliff.delta(score ~ event_num)
  delta_raw <- cd$estimate
  median0 <- median(score[event_num == 0], na.rm = TRUE)
  median1 <- median(score[event_num == 1], na.rm = TRUE)
  delta_adj <- ifelse(median1 > median0, abs(delta_raw), -abs(delta_raw))
  roc_obj <- pROC::roc(event_num, score, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))

  list(n_no_event = sum(event_num == 0), n_event = sum(event_num == 1),
       p_value = p_val, cliffs_delta = delta_adj, auc = auc_val)
}

# Violin: Metastasis
clin_violin <- clin_valid  # fresh copy for violin
clin_violin$met_event <- factor(clin_violin$mets, levels = c(0, 1))
stats_durham_met <- compute_violin_stats(clin_violin, "MetScore_prob", "met_event")
cat("Metastasis violin stats:", str(stats_durham_met), "\n")

p_met_violin <- make_violin(clin_valid, "MetScore_prob", "mets",
                             title = "Durham - Met-Score by Metastasis Event")
ggsave(file.path(fig_dir, "Violin_Durham_MetScore_by_met_event.pdf"),
       p_met_violin, width = 7.5, height = 6.5, dpi = 600, useDingbats = FALSE)
ggsave(file.path(fig_dir, "Violin_Durham_MetScore_by_met_event.tiff"),
       p_met_violin, width = 7.5, height = 6.5, dpi = 600, compression = "lzw")

# Violin: BCR
p_bcr_violin <- make_violin(clin_valid, "MetScore_prob", "bcr_binary",
                              title = "Durham - Met-Score by BCR Event")
ggsave(file.path(fig_dir, "Violin_Durham_MetScore_by_BCR_event.pdf"),
       p_bcr_violin, width = 7.5, height = 6.5, dpi = 600, useDingbats = FALSE)

# Violin: OS
p_os_violin <- make_violin(clin_valid, "MetScore_prob", "dead",
                             title = "Durham - Met-Score by Mortality")
ggsave(file.path(fig_dir, "Violin_Durham_MetScore_by_OS_event.pdf"),
       p_os_violin, width = 7.5, height = 6.5, dpi = 600, useDingbats = FALSE)

cat("Violin plots saved.\n")

# ============================================================
# 15.5  GS7 RMST analyses — All-GS7, GG2, GG3 (parallel to JHU)
# ============================================================
# Mirrors the canonical RMST block in code/Met_PCa_Survival.R for JHU.
# Restricts to pathological Gleason 7 with valid Grade Group
# (GG2 = pogl1=3 & pogl2=4; GG3 = pogl1=4 & pogl2=3) and computes
# RMST(High) - RMST(Low) with 95% CI and p via survRM2::rmst2.
# Default tau = min(max time per arm) so both arms are observable.

suppressPackageStartupMessages({
  library(survRM2); library(dplyr)
})

.rmst_diff_only <- function(df, time_col, event_col, group_col, tau = NULL) {
  df <- df %>% dplyr::filter(!is.na(.data[[time_col]]),
                              !is.na(.data[[event_col]]),
                              !is.na(.data[[group_col]]))
  time  <- as.numeric(df[[time_col]])
  event <- as.numeric(as.character(df[[event_col]]))
  g     <- droplevels(factor(df[[group_col]]))
  stopifnot(length(levels(g)) == 2)
  arm   <- ifelse(g == levels(g)[2], 1, 0)   # arm = 1 → 2nd level (High)
  if (is.null(tau)) {
    tau <- min(max(time[arm == 0], na.rm = TRUE),
               max(time[arm == 1], na.rm = TRUE))
  }
  fit <- survRM2::rmst2(time = time, status = event, arm = arm, tau = tau)
  ur  <- as.data.frame(fit$unadjusted.result, stringsAsFactors = FALSE)
  rn  <- rownames(ur)
  idx <- which(grepl("RMST", rn, ignore.case = TRUE) &
                grepl("\\(arm=1\\)-\\(arm=0\\)", rn))
  if (length(idx) != 1) {
    stop("Could not uniquely identify RMST difference row.")
  }
  data.frame(
    tau         = tau,
    RMST_diff   = as.numeric(ur[idx, "Est."]),
    RMST_diff_L = as.numeric(ur[idx, "lower .95"]),
    RMST_diff_U = as.numeric(ur[idx, "upper .95"]),
    p           = as.numeric(ur[idx, "p"]),
    n           = nrow(df),
    events      = sum(event == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Build GS7 subset with Grade Group classification
durham_gs7 <- clin_valid %>%
  dplyr::filter(!is.na(pogl), pogl == 7,
                !is.na(MetScoreClass),
                !is.na(surgmets), !is.na(mets)) %>%
  dplyr::mutate(
    GG = dplyr::case_when(
      as.numeric(pogl1) == 3 & as.numeric(pogl2) == 4 ~ "GG2_3plus4",
      as.numeric(pogl1) == 4 & as.numeric(pogl2) == 3 ~ "GG3_4plus3",
      TRUE ~ NA_character_
    ),
    MetScoreClass = factor(MetScoreClass, levels = c("Low risk", "High risk"))
  ) %>%
  dplyr::filter(!is.na(GG))

cat("\n=== Durham GS7 RMST (All / GG2 / GG3) ===\n")
cat(sprintf("GS7 N = %d (GG2 = %d, GG3 = %d), events = %d\n",
             nrow(durham_gs7),
             sum(durham_gs7$GG == "GG2_3plus4"),
             sum(durham_gs7$GG == "GG3_4plus3"),
             sum(durham_gs7$mets == 1, na.rm = TRUE)))

rmst_durham_all  <- .rmst_diff_only(durham_gs7,                                "surgmets", "mets", "MetScoreClass")
rmst_durham_GG2  <- .rmst_diff_only(durham_gs7 %>% dplyr::filter(GG == "GG2_3plus4"),
                                     "surgmets", "mets", "MetScoreClass")
rmst_durham_GG3  <- .rmst_diff_only(durham_gs7 %>% dplyr::filter(GG == "GG3_4plus3"),
                                     "surgmets", "mets", "MetScoreClass")

durham_rmst <- dplyr::bind_rows(
  data.frame(Stratum = "All GS7",          rmst_durham_all,  stringsAsFactors = FALSE),
  data.frame(Stratum = "GG2 (3+4)",        rmst_durham_GG2,  stringsAsFactors = FALSE),
  data.frame(Stratum = "GG3 (4+3)",        rmst_durham_GG3,  stringsAsFactors = FALSE)
)
print(durham_rmst, row.names = FALSE)
write.csv(durham_rmst,
          file.path(out_dir, "Durham_GS7_RMST.csv"),
          row.names = FALSE)
cat("Wrote -> ", file.path(out_dir, "Durham_GS7_RMST.csv"), "\n", sep = "")

# ============================================================
# 15.6  GS7 DEEP-DIVE (Cox / LRT / stratified / GG-specific / ΔC)
# ============================================================
# Parallel to the JHU GS7 block in code/Met_PCa_Survival.R.
# Uses the same subset (durham_gs7, GS7 + valid GG) constructed above
# for RMST. Writes one canonical CSV used directly by manuscript Section 3.

cat("\n=== Durham GS7 deep-dive (Cox + LRT + ΔC) ===\n")
.dd <- durham_gs7 %>%
  dplyr::mutate(
    GG_factor = factor(GG, levels = c("GG2_3plus4", "GG3_4plus3")),
    MetScoreProb = MetScore_prob
  )

# 1) Univariate Cox (Met-Score binary)
cox_dd_uni <- coxph(Surv(surgmets, mets) ~ MetScoreClass, data = .dd)
sm_uni     <- summary(cox_dd_uni)
ph_uni     <- cox.zph(cox_dd_uni, transform = "rank")
lr_uni     <- survdiff(Surv(surgmets, mets) ~ MetScoreClass, data = .dd)
lr_p_uni   <- 1 - pchisq(lr_uni$chisq, df = length(lr_uni$n) - 1)

# 2) Multivariable Cox: Met-Score + GG (GG2 reference)
cox_dd_full <- coxph(Surv(surgmets, mets) ~ MetScoreClass + GG_factor, data = .dd)
sm_full     <- summary(cox_dd_full)

# 3) Cox: GG-only baseline (for LRT and ΔC)
cox_dd_gg <- coxph(Surv(surgmets, mets) ~ GG_factor, data = .dd)

# 4) LRT — Met-Score class added to GG-only
lrt_class <- anova(cox_dd_gg, cox_dd_full, test = "LRT")

# 5) LRT — Met-Score continuous probability added to GG-only
.dd_p <- .dd %>% dplyr::filter(!is.na(MetScoreProb))
cox_dd_gg_p   <- coxph(Surv(surgmets, mets) ~ GG_factor, data = .dd_p)
cox_dd_full_p <- coxph(Surv(surgmets, mets) ~ MetScoreProb + GG_factor, data = .dd_p)
lrt_prob      <- anova(cox_dd_gg_p, cox_dd_full_p, test = "LRT")

# 6) Stratified Cox (GG as strata)
cox_dd_strat <- coxph(Surv(surgmets, mets) ~ MetScoreClass + strata(GG_factor),
                       data = .dd)
sm_strat <- summary(cox_dd_strat)

# 7) GG-specific univariate Cox
cox_dd_GG2 <- coxph(Surv(surgmets, mets) ~ MetScoreClass,
                    data = dplyr::filter(.dd, GG == "GG2_3plus4"))
cox_dd_GG3 <- coxph(Surv(surgmets, mets) ~ MetScoreClass,
                    data = dplyr::filter(.dd, GG == "GG3_4plus3"))
sm_GG2 <- summary(cox_dd_GG2); sm_GG3 <- summary(cox_dd_GG3)

# 8) C-index (Harrell) via survival::concordance
C_GG       <- as.numeric(survival::concordance(cox_dd_gg)$concordance)
C_GG_class <- as.numeric(survival::concordance(cox_dd_full)$concordance)
C_GG_prob  <- as.numeric(survival::concordance(cox_dd_full_p)$concordance)
C_GG_for_p <- as.numeric(survival::concordance(cox_dd_gg_p)$concordance)

# Assemble a clean tidy CSV: one row per metric, every cell populated.
# Columns: metric | description | estimate | ci_low | ci_high | p_value | n | events
# Cells are blank ("") rather than NA when not applicable to that metric,
# so the CSV reads naturally without scattered NAs.
.fmt   <- function(x, d = 3) ifelse(is.na(x), "", as.character(round(x, d)))
.fmtp  <- function(x)        ifelse(is.na(x), "", signif(x, 4))
durham_deepdive <- rbind(
  data.frame(metric = "n_total",        description = "Patients (GS7 + valid GG)",
             estimate = nrow(.dd), ci_low = "", ci_high = "", p_value = "",
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "n_GG2",          description = "GG2 (3+4) subset",
             estimate = sum(.dd$GG == "GG2_3plus4"), ci_low = "", ci_high = "", p_value = "",
             n = sum(.dd$GG == "GG2_3plus4"),
             events = sum(.dd$mets == 1 & .dd$GG == "GG2_3plus4", na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "n_GG3",          description = "GG3 (4+3) subset",
             estimate = sum(.dd$GG == "GG3_4plus3"), ci_low = "", ci_high = "", p_value = "",
             n = sum(.dd$GG == "GG3_4plus3"),
             events = sum(.dd$mets == 1 & .dd$GG == "GG3_4plus3", na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "univariate_HR",  description = "Met-Score High vs Low (Cox)",
             estimate = .fmt(sm_uni$conf.int[1, "exp(coef)"]),
             ci_low   = .fmt(sm_uni$conf.int[1, "lower .95"]),
             ci_high  = .fmt(sm_uni$conf.int[1, "upper .95"]),
             p_value  = .fmtp(sm_uni$coefficients[1, "Pr(>|z|)"]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "univariate_logrank_p", description = "Log-rank p (univariate KM)",
             estimate = .fmtp(lr_p_uni), ci_low = "", ci_high = "",
             p_value  = .fmtp(lr_p_uni),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "univariate_C_Harrell", description = "Univariate Cox concordance (Harrell)",
             estimate = .fmt(as.numeric(survival::concordance(cox_dd_uni)$concordance)),
             ci_low = "", ci_high = "", p_value = "",
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "univariate_schoenfeld_p", description = "Schoenfeld global PH test (rank-transformed)",
             estimate = .fmtp(ph_uni$table[nrow(ph_uni$table), "p"]),
             ci_low = "", ci_high = "",
             p_value  = .fmtp(ph_uni$table[nrow(ph_uni$table), "p"]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "multivar_HR_msclass", description = "Met-Score (mult. Cox; +GG, GG2 ref)",
             estimate = .fmt(sm_full$conf.int["MetScoreClassHigh risk", "exp(coef)"]),
             ci_low   = .fmt(sm_full$conf.int["MetScoreClassHigh risk", "lower .95"]),
             ci_high  = .fmt(sm_full$conf.int["MetScoreClassHigh risk", "upper .95"]),
             p_value  = .fmtp(sm_full$coefficients["MetScoreClassHigh risk", "Pr(>|z|)"]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "multivar_HR_GG3vGG2", description = "GG3 vs GG2 (mult. Cox; +Met-Score)",
             estimate = .fmt(sm_full$conf.int["GG_factorGG3_4plus3", "exp(coef)"]),
             ci_low   = .fmt(sm_full$conf.int["GG_factorGG3_4plus3", "lower .95"]),
             ci_high  = .fmt(sm_full$conf.int["GG_factorGG3_4plus3", "upper .95"]),
             p_value  = .fmtp(sm_full$coefficients["GG_factorGG3_4plus3", "Pr(>|z|)"]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "LRT_msclass_chi2", description = "LRT: GG -> GG + Met-Score class (chi2, df=1)",
             estimate = .fmt(lrt_class$Chisq[2]),
             ci_low = "", ci_high = "",
             p_value  = .fmtp(lrt_class$`Pr(>|Chi|)`[2]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "LRT_msprob_chi2", description = "LRT: GG -> GG + Met-Score probability (chi2, df=1)",
             estimate = .fmt(lrt_prob$Chisq[2]),
             ci_low = "", ci_high = "",
             p_value  = .fmtp(lrt_prob$`Pr(>|Chi|)`[2]),
             n = nrow(.dd_p), events = sum(.dd_p$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "stratified_HR_msclass", description = "Met-Score (stratified Cox; strata = GG)",
             estimate = .fmt(sm_strat$conf.int[1, "exp(coef)"]),
             ci_low   = .fmt(sm_strat$conf.int[1, "lower .95"]),
             ci_high  = .fmt(sm_strat$conf.int[1, "upper .95"]),
             p_value  = .fmtp(sm_strat$coefficients[1, "Pr(>|z|)"]),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "C_GG_only", description = "C-index: GG only (Harrell)",
             estimate = .fmt(C_GG), ci_low = "", ci_high = "", p_value = "",
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "C_GG_plus_class", description = "C-index: GG + Met-Score class (deltaC vs GG only)",
             estimate = .fmt(C_GG_class), ci_low = "", ci_high = "",
             p_value  = paste0("deltaC=", .fmt(C_GG_class - C_GG)),
             n = nrow(.dd), events = sum(.dd$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "C_GG_plus_prob", description = "C-index: GG + Met-Score probability (deltaC vs GG only)",
             estimate = .fmt(C_GG_prob), ci_low = "", ci_high = "",
             p_value  = paste0("deltaC=", .fmt(C_GG_prob - C_GG_for_p)),
             n = nrow(.dd_p), events = sum(.dd_p$mets == 1, na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "gg_specific_HR_GG2", description = "Met-Score Cox restricted to GG2 (univariate)",
             estimate = .fmt(sm_GG2$conf.int[1, "exp(coef)"]),
             ci_low   = .fmt(sm_GG2$conf.int[1, "lower .95"]),
             ci_high  = .fmt(sm_GG2$conf.int[1, "upper .95"]),
             p_value  = .fmtp(sm_GG2$coefficients[1, "Pr(>|z|)"]),
             n = sum(.dd$GG == "GG2_3plus4"),
             events = sum(.dd$mets == 1 & .dd$GG == "GG2_3plus4", na.rm = TRUE),
             stringsAsFactors = FALSE),
  data.frame(metric = "gg_specific_HR_GG3", description = "Met-Score Cox restricted to GG3 (univariate)",
             estimate = .fmt(sm_GG3$conf.int[1, "exp(coef)"]),
             ci_low   = .fmt(sm_GG3$conf.int[1, "lower .95"]),
             ci_high  = .fmt(sm_GG3$conf.int[1, "upper .95"]),
             p_value  = .fmtp(sm_GG3$coefficients[1, "Pr(>|z|)"]),
             n = sum(.dd$GG == "GG3_4plus3"),
             events = sum(.dd$mets == 1 & .dd$GG == "GG3_4plus3", na.rm = TRUE),
             stringsAsFactors = FALSE)
)
print(durham_deepdive, row.names = FALSE)
write.csv(durham_deepdive,
          file.path(out_dir, "Durham_GS7_deepdive.csv"),
          row.names = FALSE)
cat("Wrote -> ", file.path(out_dir, "Durham_GS7_deepdive.csv"), "\n", sep = "")

# ============================================================
# 15.7  GS7 score-by-event Wilcoxon + classification AUC
#       Canonical CSV for Figure 3a/b: Wilcoxon raw p + ROC AUC + 95% CI
# ============================================================
.dd_w <- .dd %>%
  dplyr::filter(!is.na(MetScoreProb), !is.na(mets)) %>%
  dplyr::mutate(event = as.integer(mets))

.wp <- tryCatch(
  wilcox.test(.dd_w$MetScoreProb ~ .dd_w$event)$p.value,
  error = function(e) NA_real_)

.roc_gs7 <- tryCatch(
  pROC::roc(response = .dd_w$event, predictor = .dd_w$MetScoreProb,
            levels = c(0, 1), direction = "<", ci = TRUE),
  error = function(e) NULL)
if (!is.null(.roc_gs7)) {
  .auc_gs7  <- as.numeric(pROC::auc(.roc_gs7))
  .auc_ci   <- as.numeric(pROC::ci.auc(.roc_gs7))
  .auc_lo   <- .auc_ci[1]
  .auc_hi   <- .auc_ci[3]
} else {
  .auc_gs7 <- .auc_lo <- .auc_hi <- NA_real_
}

durham_gs7_aw <- data.frame(
  Cohort         = "Durham",
  n              = nrow(.dd_w),
  events         = sum(.dd_w$event == 1),
  Wilcoxon_raw_p = signif(.wp, 4),
  AUC            = round(.auc_gs7, 4),
  AUC_lo         = round(.auc_lo,  4),
  AUC_hi         = round(.auc_hi,  4)
)
cat("\n=== Durham GS7 score Wilcoxon + AUC ===\n")
print(durham_gs7_aw, row.names = FALSE)
write.csv(durham_gs7_aw,
          file.path(out_dir, "Durham_GS7_score_Wilcoxon_AUC.csv"),
          row.names = FALSE)
cat("Wrote -> ", file.path(out_dir, "Durham_GS7_score_Wilcoxon_AUC.csv"), "\n", sep = "")

# 15.8  (removed) The Durham GS7 Decipher-surrogate median-split head-to-head
#       was superseded by the canonical Figure S2 producer
#       (code/ancillary/Met_PCa_Survival_DECIPHER.R), which consumes the frozen
#       surrogate probabilities and writes identifier-free aggregates.

# ============================================================
# 16. Calibration metrics (locked model)
# Metrics: Brier score, calibration slope, CITL, scaled Brier, E/O ratio.
# Uses MetScore_prob (QN-ridge LR probability) and mets (binary outcome).
# ============================================================
suppressPackageStartupMessages(library(scales))

.logit_dur <- function(p) log(p / (1 - p))

.calib_metrics_dur <- function(y_obs, p_pred, cohort = "", n_bins = 10) {
  df <- data.frame(y = as.numeric(as.character(y_obs)), p = p_pred)
  df <- df[complete.cases(df), ]

  BS        <- mean((df$y - df$p)^2)
  BS_null   <- mean(df$y) * (1 - mean(df$y))
  BS_scaled <- 1 - BS / BS_null

  lp <- .logit_dur(pmin(pmax(df$p, 1e-6), 1 - 1e-6))
  cf <- glm(y ~ lp, data = df, family = binomial())

  df$bin <- dplyr::ntile(df$p, n_bins)
  bin_df <- df %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(mean_pred = mean(p), mean_obs = mean(y),
                     n = dplyr::n(), .groups = "drop")

  list(cohort = cohort, n = nrow(df), events = sum(df$y),
       prevalence = mean(df$y),
       brier = round(BS, 4), brier_scaled = round(BS_scaled, 4),
       slope = round(coef(cf)[["lp"]], 3),
       intercept = round(coef(cf)[["(Intercept)"]], 3),
       EO_ratio = round(mean(df$p) / mean(df$y), 3),
       bin_df = bin_df)
}

.plot_calib_dur <- function(calib, out_prefix, title_text = "") {
  bdf <- calib$bin_df
  ann <- sprintf("Brier = %.3f\nSlope = %.2f\nCITL = %.2f",
                 calib$brier, calib$slope, calib$intercept)

  # Zoom axes to the observed data range + 5% padding, cap at 1
  ax_max <- min(1, ceiling(max(c(bdf$mean_pred, bdf$mean_obs),
                               na.rm = TRUE) / 0.05) * 0.05 + 0.05)
  ax_lim <- c(0, ax_max)
  ann_x  <- ax_max * 0.03
  ann_y  <- ax_max * 0.97

  p <- ggplot(bdf, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", linewidth = 0.8, color = "grey55") +
    geom_smooth(aes(weight = n), method = "loess", formula = y ~ x,
                se = TRUE, linewidth = 1.1, color = "#ED6905",
                fill = "#ED6905", alpha = 0.15) +
    geom_point(size = 3.5, color = "#2B2EB5", alpha = 0.85) +
    geom_text(aes(label = n), vjust = -0.9, size = 3.2,
              color = "#2B2EB5", fontface = "bold") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = ax_lim, expand = expansion(mult = 0.02)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = ax_lim, expand = expansion(mult = 0.02)) +
    annotate("text", x = ann_x, y = ann_y, label = ann,
             hjust = 0, vjust = 1, size = 4.5, color = "black") +
    labs(x = "Predicted probability", y = "Observed event rate",
         title = title_text) +
    theme_classic(base_size = 14) +
    theme(
      axis.title   = element_text(face = "bold", size = 14),
      axis.text    = element_text(face = "bold", size = 13),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 15),
      legend.position = "none"
    )

  ggsave(paste0(out_prefix, ".pdf"), plot = p, width = 5, height = 5,
         useDingbats = FALSE)
  ggsave(paste0(out_prefix, ".tiff"), plot = p, width = 5, height = 5,
         dpi = 400, compression = "lzw")
  invisible(p)
}

calib_durham <- .calib_metrics_dur(
  y_obs  = clin_valid$mets,
  p_pred = clin_valid$MetScore_prob,
  cohort = "Durham VA"
)
.plot_calib_dur(
  calib_durham,
  out_prefix = file.path(fig_dir, "Durham_calibration_lockedLR"),
  title_text = "Durham VA: Locked model calibration"
)

.dur_calib_tbl <- data.frame(
  Cohort      = calib_durham$cohort,
  N           = calib_durham$n,
  Events      = calib_durham$events,
  Prevalence  = round(calib_durham$prevalence, 3),
  Brier       = calib_durham$brier,
  Brier_scaled= calib_durham$brier_scaled,
  Calib_slope = calib_durham$slope,
  CITL        = calib_durham$intercept,
  EO_ratio    = calib_durham$EO_ratio
)
write.csv(.dur_calib_tbl,
          file.path(out_dir, "CalibrationMetrics_Durham.csv"),
          row.names = FALSE)
cat("\n=== Durham calibration metrics — locked ridge LR ===\n")
print(.dur_calib_tbl, row.names = FALSE)

# ============================================================
# 17. Figure S2 aggregates + save workspace
# ============================================================
write_figure_s2_aggregates(clin_valid, out_dir, LOCKED_THRESHOLD)

save(clin_valid, auc_table, surv_table,
     pos_genes, neg_genes, geneset_median,
     train_thr, Confusion_durham_qn, Confusion_durham_gs,
     MCC_durham_qn, MCC_durham_gs, stats_durham_met,
     durham_rmst, durham_deepdive,
     file = file.path(out_dir, "durham_metscore_batchcorrected.rda"))

cat("\n================================================================\n")
cat("=== ALL RESULTS SAVED ===\n")
cat("Output:", out_dir, "\n")
cat("Figures:", fig_dir, "\n")
cat("================================================================\n")
cat("=== DONE ===\n")

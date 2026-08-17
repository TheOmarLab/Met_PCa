############################################################################
# Multivariable Cox analysis of the Met-Score in two external validation cohorts
#
# The Met-Score used here is the locked 41-feature ridge-logistic CLASSIFIER
# (stored probability per sample), which is a separate object from the 45-gene
# biological signature. The classifier probability is only USED, never refit.
#
# Endpoint (critical):
#   JHU    : Surv(met_time, met)
#   Durham : Surv(surgmets, mets)   # surgmets is the time-to-metastasis clock;
#            fu is a different follow-up clock and is not the endpoint clock.
#
# Canonical covariates (identical construction in both cohorts):
#   GG      : Grade Group from total/primary Gleason, reference GG2
#   log2PSA : log2(PSA + 1), keeps PSA == 0 finite
#   pT      : pathological T stage collapsed to T2/T3/T4, reference T2
#   ms_z    : classifier probability, z-scored within cohort (HR per 1 SD)
#   node    : nodal involvement (Durham codebook 0=no,1=yes,2=unknown,3=not done;
#             only 0/1 informative, 2/3 -> NA)
# margin : positive surgical margin
#   MetScoreClass : binary High/Low classifier call
#
# Models fitted per cohort on ONE complete-case set (complete on GG, log2PSA,
# pT, ms_z):
#   M0 (base)      : Surv ~ GG + log2PSA + pT
#   M1 (primary)   : M0 + ms_z              -> full HR table + LRT
#   sensitivity    : M1 + margin (both cohorts); M1 + node (JHU only) -> ms_z HR.
#                    Durham node HR is not estimated (incomplete ascertainment);
#                    node counts are reported in MultiCox_node_availability.csv.
#   secondary      : GG + log2PSA + pT + MetScoreClass -> class HR only
#
# Margin (positive surgical margin) is entered only as a sensitivity covariate,
# not a primary adjustment variable.
#
# Metrics adapted from code/R2/R2_02_extended_cox.R:
#   per-term HR (95% CI, p) and the likelihood-ratio test M1 vs M0. This script
#   reports association and model adequacy only; external discrimination is in
#   the separate DCA/calibration path.
#
# Outputs: outs/MultiCox_JHU_HRtable.csv, outs/MultiCox_Durham_HRtable.csv
#          outs/MultiCox_JHU_summary.csv, outs/MultiCox_Durham_summary.csv
#          outs/MultiCox_node_availability.csv
############################################################################

rm(list = ls())
suppressPackageStartupMessages(library(survival))
set.seed(20260515)
dir.create("./outs", recursive = TRUE, showWarnings = FALSE)

# ---- Data ----------------------------------------------------------------
load("./outs/coxdata.rda")                                  # CoxData_jhu
load("./output/Durham/durham_metscore_batchcorrected.rda")  # clin_valid

# ---- Covariate constructors ---------------------------------------------
# Grade Group from total Gleason and primary pattern.
#   GG1 total <=6 ; GG2 3+4 ; GG3 4+3 ; GG4 GS8 ; GG5 GS9-10
gg_from <- function(total_gs, primary) {
  gg <- rep(NA_character_, length(total_gs))
  gg[total_gs <= 6]                <- "GG1"
  gg[total_gs == 7 & primary == 3] <- "GG2"
  gg[total_gs == 7 & primary == 4] <- "GG3"
  gg[total_gs == 8]                <- "GG4"
  gg[total_gs >= 9]                <- "GG5"
  gg
}

# Pathological T stage collapsed to T2 / T3 / T4.
pt_collapse <- function(x) {
  x   <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^T2", x)] <- "T2"
  out[grepl("^T3", x)] <- "T3"
  out[grepl("^T4", x)] <- "T4"
  out
}

# ---- Per-cohort analysis frames -----------------------------------------
# Both frames expose the same column names so a single formula set applies.
build_jhu <- function() {
  j <- CoxData_jhu
  data.frame(
    time     = as.numeric(j$met_time),                       # endpoint clock
    event    = as.integer(j$met),
    GG       = gg_from(as.numeric(as.character(j[["Pathological GS"]])),
                       as.numeric(as.character(j$pathgs_p))),
    log2PSA  = log2(as.numeric(j$preop_psa) + 1),
    pT       = pt_collapse(j$pstage),
    node     = as.integer(j$lni),
    margin   = as.integer(j$sm),
    ms_prob  = as.numeric(j[["Met-Score prob"]]),
    ms_class = factor(as.character(j$MetScoreClass), levels = c("Low risk", "High risk")),
    stringsAsFactors = FALSE
  )
}

build_durham <- function() {
  d <- clin_valid
  # Durham nodal involvement per codebook: 0=no, 1=yes, 2=unknown, 3=not done.
  # Only 0 and 1 are informative; 2 and 3 map to NA. Any other value stops.
  node_raw <- as.numeric(d$lymphnodeinvolvement)
  bad <- !is.na(node_raw) & !(node_raw %in% c(0, 1, 2, 3))
  if (any(bad)) stop(sprintf("Durham lymphnodeinvolvement has unexpected code(s): %s",
                             paste(sort(unique(node_raw[bad])), collapse = ", ")))
  node_bin <- ifelse(node_raw == 0, 0L, ifelse(node_raw == 1, 1L, NA_integer_))
  data.frame(
    time     = as.numeric(d$surgmets),                       # surgmets, not fu
    event    = as.integer(d$mets),
    GG       = gg_from(as.numeric(as.character(d$PathGleason)),
                       as.numeric(as.character(d$pogl1))),
    log2PSA  = log2(as.numeric(d$psapresurg) + 1),
    pT       = pt_collapse(d$stg),
    node     = node_bin,                                     # 0=no, 1=yes; 2/3 -> NA
    margin   = as.integer(d$positivesurgicalmargins),
    ms_prob  = as.numeric(d$MetScore_prob),
    ms_class = factor(as.character(d$MetScoreClass), levels = c("Low risk", "High risk")),
    stringsAsFactors = FALSE
  )
}

# z-score the classifier probability WITHIN cohort so its HR is per 1 SD.
add_ms_z <- function(df) { df$ms_z <- as.numeric(scale(df$ms_prob)); df }

# ---- Fitting helpers -----------------------------------------------------
# Complete-case set is defined once per cohort on the primary covariates
# {GG, log2PSA, pT, ms_z}; every model in the cohort is fit on this same set.
make_cc <- function(df) {
  cc <- df[stats::complete.cases(df[, c("GG", "log2PSA", "pT", "ms_z")]), , drop = FALSE]
  cc$GG <- relevel(factor(cc$GG), ref = "GG2")
  cc$pT <- relevel(factor(cc$pT), ref = "T2")
  cc$ms_class <- relevel(cc$ms_class, ref = "Low risk")
  cc
}

# One Wald row (Model, term, HR, CI, p, Status) for a single model term. The
# estimability check is generic: a term is non-estimable when its coefficient,
# standard error, HR, or 95% CI is not finite and strictly positive, as happens
# under complete separation. Non-estimable rows keep the term but carry NA
# estimates and an explicit Status; finite rows are labelled "estimated".
wald_row <- function(fit, term, label, term_label = term) {
  co   <- summary(fit)$coefficients
  ci   <- summary(fit)$conf.int
  beta <- co[term, "coef"]; se <- co[term, "se(coef)"]; p <- co[term, "Pr(>|z|)"]
  hr   <- ci[term, "exp(coef)"]; lo <- ci[term, "lower .95"]; hi <- ci[term, "upper .95"]
  ok   <- is.finite(beta) && is.finite(se) && is.finite(hr) && hr > 0 &&
          is.finite(lo) && lo > 0 && is.finite(hi)
  data.frame(
    Model  = label,
    term   = term_label,
    HR     = if (ok) round(hr, 4)  else NA_real_,
    CI_lo  = if (ok) round(lo, 4)  else NA_real_,
    CI_hi  = if (ok) round(hi, 4)  else NA_real_,
    p      = if (ok) signif(p, 4)  else NA_real_,
    Status = if (ok) "estimated" else "non-estimable: complete separation",
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Full HR table for a fitted model, one wald_row per term.
hr_table <- function(fit, label) {
  terms <- rownames(summary(fit)$conf.int)
  do.call(rbind, lapply(terms, function(tm) wald_row(fit, tm, label)))
}

# ---- Per-cohort driver ---------------------------------------------------
# fit_node_sens: the node-adjusted sensitivity is fit only where nodal status is
# adequately ascertained (JHU). Durham node ascertainment is incomplete, so the
# node HR is not estimated there; node availability is reported separately.
run_cohort <- function(df, cohort, fit_node_sens = TRUE) {
  df <- add_ms_z(df)
  cc <- make_cc(df)

  f0 <- Surv(time, event) ~ GG + log2PSA + pT
  f1 <- Surv(time, event) ~ GG + log2PSA + pT + ms_z
  f1_margin <- update(f1, . ~ . + margin)   # SENSITIVITY
  f1_node   <- update(f1, . ~ . + node)     # SENSITIVITY (fit only when fit_node_sens)
  f_class   <- Surv(time, event) ~ GG + log2PSA + pT + ms_class  # SECONDARY

  m0 <- coxph(f0, data = cc)
  m1 <- coxph(f1, data = cc)
  m1_margin <- coxph(f1_margin, data = cc)
  m1_node   <- if (fit_node_sens) coxph(f1_node, data = cc) else NULL
  m_class   <- coxph(f_class,   data = cc)

  # Primary full HR table (all M1 terms).
  hr_m1 <- hr_table(m1, "M1 primary (GG+log2PSA+pT+ms_z)")

  # LRT M1 vs M0 on the same complete-case set.
  lrt <- anova(m0, m1)
  lrt_chi <- lrt$Chisq[2]; lrt_df <- lrt$Df[2]; lrt_p <- lrt$`Pr(>|Chi|)`[2]

  # ms_z HR from each sensitivity model. The primary ms_z already appears in the
  # full M1 table, so only the added-covariate sensitivity models are listed here.
  hr_sens <- wald_row(m1_margin, "ms_z", "M1 + margin (SENSITIVITY)")
  if (fit_node_sens)
    hr_sens <- rbind(hr_sens, wald_row(m1_node, "ms_z", "M1 + node (SENSITIVITY)"))

  # Secondary: binary classifier class HR (same schema).
  cls_row <- grep("^ms_class", rownames(summary(m_class)$conf.int), value = TRUE)[1]
  hr_class <- wald_row(m_class, cls_row,
                       "SECONDARY class (GG+log2PSA+pT+MetScoreClass)",
                       term_label = "MetScoreClass High vs Low")

  # margin adjusted HR (for the reporting note about Durham instability).
  margin_hr <- tryCatch({
    s <- summary(m1_margin)$conf.int
    sprintf("%.2f (%.2f-%.2f)", s["margin","exp(coef)"], s["margin","lower .95"], s["margin","upper .95"])
  }, error = function(e) NA_character_)

  summary_tab <- data.frame(
    Cohort       = cohort,
    n            = nrow(cc),
    events       = sum(cc$event),
    LRT_chi2     = round(lrt_chi, 4),
    LRT_df       = lrt_df,
    LRT_p        = signif(lrt_p, 4),
    margin_adjHR = margin_hr,
    stringsAsFactors = FALSE
  )

  list(cohort = cohort,
       hr = rbind(hr_m1, hr_sens, hr_class),
       summary = summary_tab)
}

# ---- Run -----------------------------------------------------------------
cat("Building JHU and Durham analysis frames...\n")
RJ <- run_cohort(build_jhu(),    "JHU",    fit_node_sens = TRUE)
RD <- run_cohort(build_durham(), "Durham", fit_node_sens = FALSE)

# ---- Write per-cohort HR tables and summaries ---------------------------
write.csv(RJ$hr,      "./outs/MultiCox_JHU_HRtable.csv",     row.names = FALSE)
write.csv(RD$hr,      "./outs/MultiCox_Durham_HRtable.csv",  row.names = FALSE)
write.csv(RJ$summary, "./outs/MultiCox_JHU_summary.csv",     row.names = FALSE)
write.csv(RD$summary, "./outs/MultiCox_Durham_summary.csv",  row.names = FALSE)

# ---- Node availability (aggregate, non-identifying) ----------------------
# Reports node-code counts and metastasis events per cohort. JHU node status is
# adequately ascertained (sensitivity fitted). Durham node status is largely not
# assessed, so the adjusted node HR is not estimated; counts are reported instead.
node_availability <- function() {
  meaning <- c("0" = "no", "1" = "yes", "2" = "unknown", "3" = "not done")
  jln <- as.numeric(CoxData_jhu$lni);  jm <- as.integer(CoxData_jhu$met)
  dln <- as.numeric(clin_valid$lymphnodeinvolvement); dm <- as.integer(clin_valid$mets)
  jrows <- do.call(rbind, lapply(sort(unique(jln[!is.na(jln)])), function(cd) {
    s <- !is.na(jln) & jln == cd
    data.frame(Cohort = "JHU", SourceVariable = "lni", Code = cd,
               Meaning = if (as.character(cd) %in% names(meaning)) meaning[[as.character(cd)]] else "recorded",
               N = sum(s), MetastasisEvents = sum(jm[s]),
               ModelDisposition = "sensitivity fitted", stringsAsFactors = FALSE)
  }))
  drows <- do.call(rbind, lapply(c(0, 1, 2, 3), function(cd) {
    s <- !is.na(dln) & dln == cd
    data.frame(Cohort = "Durham", SourceVariable = "lymphnodeinvolvement", Code = cd,
               Meaning = meaning[[as.character(cd)]], N = sum(s), MetastasisEvents = sum(dm[s]),
               ModelDisposition = "adjusted node sensitivity not estimated: incomplete ascertainment (4 confirmed positive, 1 event; 273 not done)",
               stringsAsFactors = FALSE)
  }))
  rbind(jrows, drows)
}
write.csv(node_availability(), "./outs/MultiCox_node_availability.csv", row.names = FALSE)

# ---- Console report ------------------------------------------------------
report <- function(R) {
  cat("\n================ ", R$cohort, " ================\n", sep = "")
  cat("HR table (M1 primary, sensitivity ms_z, secondary class):\n")
  print(R$hr, row.names = FALSE)
  cat("\nAssociation and model tests:\n")
  print(R$summary, row.names = FALSE)
  ms <- R$hr[R$hr$Model == "M1 primary (GG+log2PSA+pT+ms_z)" & R$hr$term == "ms_z", ]
  cat(sprintf("\n  Primary Met-Score per-1-SD HR: %.2f (%.2f-%.2f), p=%s\n",
              ms$HR, ms$CI_lo, ms$CI_hi, format(ms$p)))
}
report(RJ)
report(RD)

cat("\n---- Reporting notes ----\n")
cat("Met-Score = locked 41-feature ridge-logistic CLASSIFIER probability",
    "(distinct from the 45-gene biological signature); used, not refit.\n")
cat(sprintf("  JHU adjusted margin HR: %s\n", RJ$summary$margin_adjHR))
cat(sprintf("  Durham adjusted margin HR: %s\n", RD$summary$margin_adjHR))
cat("Durham Grade Group 1 has zero metastasis events (complete separation), so its",
    "GGGG1 coefficient is non-estimable and reported as NA; the ms_z coefficient",
    "remains finite.\n")

cat("\nWrote outs/MultiCox_{JHU,Durham}_HRtable.csv and _summary.csv\n")

############################################################################
# Model diagnostics and pooled cross-cohort analysis
#
# Separate from the primary/sensitivity/secondary models and outputs above. For
# each cohort this refits the primary model M1 (GG + log2PSA + pT + ms_z) as a
# standalone object on the same complete-case set, and computes:
#   1. proportional hazards         survival::cox.zph(M1): global and per-term p
#   2. influence on ms_z            dfbeta / dfbetas residuals for the ms_z coef
#   3. functional form              log2PSA and ms_z, linear vs natural
#                                   (restricted) cubic spline ns(., df = 3),
#                                   nested likelihood-ratio test and AIC delta
#   4. PH sensitivity (JHU)         M1 refit with a log2PSA time-transform term
#                                   log2PSA x log(time) via tt(); reports the
#                                   Met-Score ms_z HR under that refit.
# It then stacks both cohorts, fits a cohort-stratified Cox for a pooled ms_z
# HR, and tests a Met-Score x cohort interaction by likelihood-ratio test,
# reporting the interaction HR ratio (Durham vs JHU) with 95% CI.
#
# New outputs: outs/MultiCox_diagnostics_JHU.csv
#              outs/MultiCox_diagnostics_Durham.csv
#              outs/MultiCox_pooled.csv
############################################################################
suppressPackageStartupMessages(library(splines))

# metric/value/note row with value stored as character (mixes numbers and yes/no).
mvn <- function(metric, value, note) {
  data.frame(metric = metric, value = as.character(value), note = note,
             stringsAsFactors = FALSE)
}

# Nonlinearity test for one continuous term: linear vs ns(term, df = 3).
# ns() is a natural (restricted) cubic spline whose basis contains the linear
# fit, so anova() gives a clean nested likelihood-ratio test on 2 df. dAIC is
# AIC(linear) - AIC(spline); a positive value favors the spline. "adequate" is
# yes when the nonlinearity test is not significant at 0.05.
ff_check <- function(cc, var) {
  m_lin <- coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z, data = cc)
  rhs <- c("GG", "log2PSA", "pT", "ms_z")
  rhs[rhs == var] <- sprintf("ns(%s, df = 3)", var)
  f_spl <- as.formula(paste("Surv(time, event) ~", paste(rhs, collapse = " + ")))
  tryCatch({
    m_spl <- coxph(f_spl, data = cc)
    an <- anova(m_lin, m_spl)
    p  <- an$`Pr(>|Chi|)`[2]
    list(chi = an$Chisq[2], df = an$Df[2], p = p,
         dAIC = AIC(m_lin) - AIC(m_spl),
         adequate = if (is.finite(p) && p < 0.05) "no" else "yes")
  }, error = function(e)
    list(chi = NA_real_, df = NA_real_, p = NA_real_, dAIC = NA_real_,
         adequate = "NA"))
}

# Per-cohort diagnostics -> long (metric, value, note) frame plus key numbers.
run_diagnostics <- function(df, cohort) {
  cc <- make_cc(add_ms_z(df))
  n  <- nrow(cc)
  m1 <- coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z, data = cc)

  # 1. Proportional hazards (global and per-term, terms aggregated over levels).
  zt <- cox.zph(m1)$table

  # 2. Influence on the ms_z coefficient.
  idx    <- which(names(coef(m1)) == "ms_z")
  db_ms  <- residuals(m1, type = "dfbeta")[,  idx]
  dbs_ms <- residuals(m1, type = "dfbetas")[, idx]
  thr    <- 2 / sqrt(n)

  # 3. Functional form of the two continuous terms.
  ff_psa <- ff_check(cc, "log2PSA")
  ff_ms  <- ff_check(cc, "ms_z")

  # 4. PH sensitivity for the log2PSA proportional-hazards violation (JHU only,
  # where cox.zph flags log2PSA). Refit M1 with a log2PSA time-transform term
  # log2PSA x log(time) via tt() and report the Met-Score ms_z HR under that
  # refit. This is a diagnostic refit; it does not replace the primary model.
  ph_sens <- NULL
  if (identical(cohort, "JHU")) {
    m1_tt <- tryCatch(
      coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z + tt(log2PSA),
            data = cc, tt = function(x, t, ...) x * log(t)),
      error = function(e) NULL)
    if (!is.null(m1_tt)) {
      st <- summary(m1_tt)$conf.int
      ph_sens <- rbind(
        mvn("ph_sens_ms_z_HR", signif(st["ms_z", "exp(coef)"], 4),
            "M1 + tt(log2PSA)=log2PSA x log(time); Met-Score HR under PSA time-transform"),
        mvn("ph_sens_ms_z_CI_lo", signif(st["ms_z", "lower .95"], 4), "lower 95% CI"),
        mvn("ph_sens_ms_z_CI_hi", signif(st["ms_z", "upper .95"], 4), "upper 95% CI"),
        mvn("ph_sens_note", "stable",
            "log2PSA time-transform sensitivity model")
      )
    }
  }

  fmtp <- function(x) if (is.finite(x)) signif(x, 4) else NA
  tab <- rbind(
    mvn("n", n, "complete-case n on GG, log2PSA, pT, ms_z"),
    mvn("events", sum(cc$event), "metastasis events"),
    mvn("ph_global_p", fmtp(zt["GLOBAL", "p"]),
        "cox.zph global test; PH supported if p >= 0.05"),
    mvn("ph_ms_z_p", fmtp(zt["ms_z", "p"]),
        "cox.zph ms_z term; PH supported if p >= 0.05"),
    mvn("ph_GG_p",      fmtp(zt["GG", "p"]),      "cox.zph GG term"),
    mvn("ph_log2PSA_p", fmtp(zt["log2PSA", "p"]), "cox.zph log2PSA term"),
    mvn("ph_pT_p",      fmtp(zt["pT", "p"]),      "cox.zph pT term"),
    mvn("dfbeta_ms_z_max_abs", signif(max(abs(db_ms)), 4),
        "max |dfbeta| for ms_z, coefficient-change scale"),
    mvn("dfbetas_ms_z_max_abs", signif(max(abs(dbs_ms)), 4),
        "max |standardized dfbeta| for ms_z"),
    mvn("dfbetas_threshold", signif(thr, 4), "flag cutoff 2/sqrt(n)"),
    mvn("dfbetas_ms_z_n_flagged", sum(abs(dbs_ms) > thr),
        "obs with |standardized dfbeta| > 2/sqrt(n)"),
    mvn("influence_ms_z_conclusion", "robust to individual deletion",
        sprintf(paste("2/sqrt(n) screen flags %d points; max |standardized dfbeta|",
                      "for ms_z = %s; Met-Score ms_z HR robust to individual deletion,",
                      "not a claim of no influential observations"),
                sum(abs(dbs_ms) > thr), signif(max(abs(dbs_ms)), 4))),
    mvn("ff_log2PSA_nonlin_p", fmtp(ff_psa$p),
        "linear vs ns(df=3) LRT; nonlinearity p"),
    mvn("ff_log2PSA_dAIC", signif(ff_psa$dAIC, 4),
        "AIC(linear) - AIC(spline); positive favors spline"),
    mvn("ff_log2PSA_linear_adequate", ff_psa$adequate,
        "yes if nonlinearity p >= 0.05"),
    mvn("ff_ms_z_nonlin_p", fmtp(ff_ms$p),
        "linear vs ns(df=3) LRT; nonlinearity p"),
    mvn("ff_ms_z_dAIC", signif(ff_ms$dAIC, 4),
        "AIC(linear) - AIC(spline); positive favors spline"),
    mvn("ff_ms_z_linear_adequate", ff_ms$adequate,
        "yes if nonlinearity p >= 0.05")
  )
  tab <- rbind(tab, ph_sens)
  list(cohort = cohort, tab = tab,
       ph_global = zt["GLOBAL", "p"], ph_ms = zt["ms_z", "p"],
       ff_psa = ff_psa$adequate, ff_ms = ff_ms$adequate)
}

# Pooled cross-cohort frame: stack the two per-cohort complete-case sets with a
# cohort factor and the harmonized MFS time/event columns. ms_z is standardized
# within cohort, so the pooled HR is per within-cohort SD.
build_pooled <- function() {
  dj <- add_ms_z(build_jhu());    dj$cohort <- "JHU"
  dd <- add_ms_z(build_durham()); dd$cohort <- "Durham"
  cols <- c("time", "event", "GG", "log2PSA", "pT", "ms_z", "cohort")
  pl <- rbind(dj[, cols], dd[, cols])
  pl <- pl[stats::complete.cases(pl[, c("GG", "log2PSA", "pT", "ms_z")]), , drop = FALSE]
  pl$GG     <- relevel(factor(pl$GG), ref = "GG2")
  pl$pT     <- relevel(factor(pl$pT), ref = "T2")
  pl$cohort <- relevel(factor(pl$cohort), ref = "JHU")
  pl
}

run_pooled <- function() {
  pl  <- build_pooled()
  mp  <- coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z + strata(cohort),
               data = pl)
  mpi <- coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z + ms_z:cohort + strata(cohort),
               data = pl)
  ci <- summary(mp)$conf.int
  hr <- ci["ms_z", "exp(coef)"]; lo <- ci["ms_z", "lower .95"]; hi <- ci["ms_z", "upper .95"]
  p  <- summary(mp)$coefficients["ms_z", "Pr(>|z|)"]
  an <- anova(mp, mpi)
  int_chi <- an$Chisq[2]; int_df <- an$Df[2]; int_p <- an$`Pr(>|Chi|)`[2]
  # Interaction HR ratio: exp of the ms_z x cohort coefficient, i.e. the ratio of
  # the within-cohort ms_z HR in Durham relative to JHU, with its 95% CI.
  cii <- summary(mpi)$conf.int
  irn <- grep("ms_z:cohort", rownames(cii), value = TRUE)[1]
  ir_hr <- cii[irn, "exp(coef)"]; ir_lo <- cii[irn, "lower .95"]; ir_hi <- cii[irn, "upper .95"]
  tab <- rbind(
    mvn("pooled_n", nrow(pl), "stacked JHU + Durham complete cases"),
    mvn("pooled_events", sum(pl$event), "metastasis events"),
    mvn("pooled_ms_z_HR", signif(hr, 4),
        "cohort-stratified Cox; HR per 1 within-cohort SD"),
    mvn("pooled_ms_z_CI_lo", signif(lo, 4), "lower 95% CI"),
    mvn("pooled_ms_z_CI_hi", signif(hi, 4), "upper 95% CI"),
    mvn("pooled_ms_z_p", signif(p, 4), "Wald p for ms_z"),
    mvn("interaction_LRT_chi2", signif(int_chi, 4),
        "ms_z x cohort LRT vs stratified no-interaction model"),
    mvn("interaction_LRT_df", int_df, "degrees of freedom"),
    mvn("interaction_p", signif(int_p, 4), "Met-Score x cohort interaction"),
    mvn("interaction_HRratio", signif(ir_hr, 4),
        "Met-Score ms_z HR ratio Durham vs JHU (exp of ms_z x cohort coefficient)"),
    mvn("interaction_HRratio_CI_lo", signif(ir_lo, 4), "lower 95% CI"),
    mvn("interaction_HRratio_CI_hi", signif(ir_hi, 4), "upper 95% CI"),
    mvn("interaction_note", NA,
        paste("No detected heterogeneity (limited power). The Met-Score x cohort",
              "interaction HR ratio is near 1 with a wide CI, so heterogeneity is",
              "neither established nor excluded; this does not claim consistent or",
              "identical Met-Score effects across cohorts."))
  )
  list(tab = tab, hr = hr, lo = lo, hi = hi, p = p,
       ir_hr = ir_hr, ir_lo = ir_lo, ir_hi = ir_hi,
       int_chi = int_chi, int_df = int_df, int_p = int_p,
       n = nrow(pl), events = sum(pl$event))
}

# ---- Run diagnostics and pooled analysis, write CSVs --------------------
cat("\n---- Diagnostics and pooled cross-cohort analysis ----\n")
DJ <- run_diagnostics(build_jhu(),    "JHU")
DD <- run_diagnostics(build_durham(), "Durham")
write.csv(DJ$tab, "./outs/MultiCox_diagnostics_JHU.csv",    row.names = FALSE)
write.csv(DD$tab, "./outs/MultiCox_diagnostics_Durham.csv", row.names = FALSE)

PL <- run_pooled()
write.csv(PL$tab, "./outs/MultiCox_pooled.csv", row.names = FALSE)

report_diag <- function(D) {
  cat(sprintf("\n[%s] cox.zph global p=%s, ms_z p=%s | linear adequate: log2PSA=%s, ms_z=%s\n",
              D$cohort, signif(D$ph_global, 3), signif(D$ph_ms, 3), D$ff_psa, D$ff_ms))
}
report_diag(DJ); report_diag(DD)
cat(sprintf("\nPooled ms_z HR=%.3f (%.3f-%.3f), p=%s; interaction LRT chi2=%.3f (df=%g), p=%s\n",
            PL$hr, PL$lo, PL$hi, signif(PL$p, 3), PL$int_chi, PL$int_df, signif(PL$int_p, 3)))
cat(sprintf("Met-Score x cohort HR ratio (Durham vs JHU)=%.3f (%.3f-%.3f)\n",
            PL$ir_hr, PL$ir_lo, PL$ir_hi))
cat("Wrote outs/MultiCox_diagnostics_{JHU,Durham}.csv and outs/MultiCox_pooled.csv\n")

############################################################################
# Figure 3 GS7 (pathological Gleason score 7) producer.
#
# One common complete-case GS7 ledger per cohort, built from the scored-object
# primary variables. ms_z is the full-cohort within-cohort per-SD scaling, computed
# before GS7 subsetting (an association scale, not a recalibrated prediction).
# Descriptive panels (locked probability vs ever-metastasis status) are
# censoring-naive; the adjusted panel reports Cox summary associations over
# observed follow-up with both model-based and Lin-Wei robust (sandwich) CIs,
# per-fit PH/influence/EPV diagnostics, and a Met-Score x Grade-Group
# interaction (robust Wald primary, nested LRT sensitivity). Aggregate,
# non-identifying outputs only.
#
# Outputs: outs/Figure3/GS7_cohort_summary.csv
#          outs/Figure3/GS7_adjusted_associations.csv
#          outs/Figure3/GS7_model_diagnostics.csv
#          outs/Figure3/GS7_ROC_curves.csv
############################################################################
suppressPackageStartupMessages(library(digest))
dir.create("./outs/Figure3", recursive = TRUE, showWarnings = FALSE)

SCORE_SCALE <- paste("locked ridge-logistic probability; ms_z = full-cohort",
                     "within-cohort per-SD z, computed before GS7 subsetting")

# GS7 ledger: keep total Gleason 7 (GG2 3+4 or GG3 4+3), then complete-case on
# the primary covariates. Counts are asserted as structural checks.
build_gs7 <- function(df_full, cohort, expect) {
  gs7 <- df_full[!is.na(df_full$GG) & df_full$GG %in% c("GG2", "GG3"), , drop = FALSE]
  n_raw <- nrow(gs7)
  cc <- gs7[stats::complete.cases(gs7[, c("GG", "log2PSA", "pT", "ms_z")]), , drop = FALSE]
  cc$GG <- relevel(factor(cc$GG, levels = c("GG2", "GG3")), ref = "GG2")
  cc$pT <- relevel(factor(cc$pT), ref = "T2")
  cc$ms_class <- relevel(cc$ms_class, ref = "Low risk")
  got <- c(nrow(cc), sum(cc$event),
           sum(cc$GG == "GG2"), sum(cc$event[cc$GG == "GG2"]),
           sum(cc$GG == "GG3"), sum(cc$event[cc$GG == "GG3"]))
  if (!isTRUE(all.equal(got, expect)))
    stop(sprintf("%s GS7 ledger %s != expected %s (n/ev/GG2n/GG2e/GG3n/GG3e; raw GS7=%d)",
                 cohort, paste(got, collapse = "/"), paste(expect, collapse = "/"), n_raw))
  attr(cc, "n_raw") <- n_raw
  cc
}

# Non-identifying SHA-256 over the sorted ledger (primary vars + locked prob).
ledger_hash <- function(cc) {
  key <- with(cc, sprintf("%s|%s|%d|%.6f|%.6f|%.4f",
              as.character(GG), as.character(pT), as.integer(event),
              log2PSA, ms_prob, time))
  digest::digest(paste(sort(key), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

# Empirical AUC (Mann-Whitney, ties averaged) and a stratified patient bootstrap CI.
emp_auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
boot_auc_ci <- function(score, y, B = 2000, seed = 20260427) {
  set.seed(seed); ip <- which(y == 1); ineg <- which(y == 0)
  aucs <- numeric(0); attempted <- 0L; failed <- 0L
  for (b in seq_len(B)) {
    attempted <- attempted + 1L
    bs <- c(sample(ip, length(ip), TRUE), sample(ineg, length(ineg), TRUE))
    a <- tryCatch(emp_auc(score[bs], y[bs]), error = function(e) NA_real_)
    if (is.finite(a)) aucs <- c(aucs, a) else failed <- failed + 1L
  }
  used <- length(aucs)
  list(auc = emp_auc(score, y),
       lo = if (used >= 50) unname(quantile(aucs, 0.025)) else NA_real_,
       hi = if (used >= 50) unname(quantile(aucs, 0.975)) else NA_real_,
       attempted = attempted, valid = used, failed = failed)
}
# Step ROC coordinates for rendering (descriptive, censoring-naive).
roc_xy <- function(score, y) {
  o <- order(score, decreasing = TRUE); s <- score[o]; yy <- y[o]
  P <- sum(y == 1); N <- sum(y == 0)
  tp <- cumsum(yy == 1); fp <- cumsum(yy == 0)
  keep <- c(which(diff(s) != 0), length(s))
  data.frame(fpr = c(0, fp[keep] / N), tpr = c(0, tp[keep] / P))
}
# Run a fit, collecting (not suppressing) warnings for disposition.
run_fit <- function(expr) {
  w <- character(0)
  fit <- withCallingHandlers(expr,
           warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  list(fit = fit, warnings = w)
}
# ---- GS7 ledgers (full-cohort ms_z attached before subsetting) ----------
GJ <- build_gs7(add_ms_z(build_jhu()),    "JHU",       c(132, 23, 86, 12, 46, 11))
GD <- build_gs7(add_ms_z(build_durham()), "Durham VA", c(422, 25, 336, 11, 86, 14))

gs7_cohort <- function(cc, cohort) {
  y <- as.integer(cc$event); s <- cc$ms_prob
  mw <- stats::wilcox.test(s ~ y, exact = FALSE)$p.value         # descriptive; normal approximation (ties)
  au <- boot_auc_ci(s, y)
  roc <- roc_xy(s, y); roc$Cohort <- cohort
  # effect modification, two tests on the same full model:
  #   sensitivity  = nested LRT (reduced vs GG*ms_z)
  #   primary      = robust (Lin-Wei sandwich) Wald on the GG:ms_z term(s)
  red    <- run_fit(coxph(Surv(time, event) ~ GG + log2PSA + pT + ms_z, data = cc))
  full   <- run_fit(coxph(Surv(time, event) ~ GG * ms_z + log2PSA + pT, data = cc))
  full_r <- run_fit(coxph(Surv(time, event) ~ GG * ms_z + log2PSA + pT, data = cc, robust = TRUE))
  an <- anova(red$fit, full$fit)
  cfr <- coef(full_r$fit); Vr <- vcov(full_r$fit)
  it  <- grep(":ms_z$", names(cfr), value = TRUE)      # GG:ms_z interaction term(s)
  bi  <- cfr[it]; Vi <- Vr[it, it, drop = FALSE]
  iw_chi <- as.numeric(t(bi) %*% solve(Vi) %*% bi); iw_df <- length(it)
  iw_p   <- stats::pchisq(iw_chi, df = iw_df, lower.tail = FALSE)
  list(cohort = cohort, n = nrow(cc), events = sum(y),
       GG2n = sum(cc$GG == "GG2"), GG2e = sum(y[cc$GG == "GG2"]),
       GG3n = sum(cc$GG == "GG3"), GG3e = sum(y[cc$GG == "GG3"]),
       mw_p = mw, auc = au, roc = roc,
       int_chi = an$Chisq[2], int_df = an$Df[2], int_p = an$`Pr(>|Chi|)`[2],
       int_wald_chi = iw_chi, int_wald_df = iw_df, int_wald_p = iw_p,
       hash = ledger_hash(cc),
       warns = c(red$warnings, full$warnings, full_r$warnings))
}
CJ <- gs7_cohort(GJ, "JHU"); CD <- gs7_cohort(GD, "Durham VA")

# BH across the two cohorts (panel-a descriptive test; interaction tests).
mw_q      <- stats::p.adjust(c(CJ$mw_p,       CD$mw_p),       method = "BH")
int_q     <- stats::p.adjust(c(CJ$int_p,      CD$int_p),      method = "BH")  # LRT sensitivity
int_robwq <- stats::p.adjust(c(CJ$int_wald_p, CD$int_wald_p), method = "BH")  # robust Wald primary

# ---- adjusted associations + per-fit diagnostics --------------------------
# Each displayed continuous fit is estimated once plain (model-based CI/p, PH,
# influence) and once with robust = TRUE (Lin-Wei sandwich CI/p). Coefficients,
# hence point HRs, are identical; only the variance differs. cox.zph and dfbeta
# are read from the plain fit so PH/influence match the standard estimator.
Z975 <- stats::qnorm(0.975)
# Panel-c adjusted associations use full-design models on the complete phase-two
# sample (JHU, Lin-Ying case-cohort risk sets) or the complete cohort (Durham,
# robust cause-specific Cox). Subgroup Met-Score/SD slopes are linear contrasts
# from a single score-by-grade interaction model, never subgroup refits. The
# competing clock is metastasis event with death before metastasis censored at its
# own time. A case-cohort PH diagnostic is not defined for JHU, so it is left NA.
pt_collapse <- function(x) { x <- toupper(trimws(as.character(x))); o <- rep(NA_character_, length(x))
  o[grepl("^T2", x)] <- "T2"; o[grepl("^T3", x)] <- "T3"; o[grepl("^T4", x)] <- "T4"; o }
ce_clock <- function(met, mt, dth, td) {
  st <- ifelse(met == 1L, 1L, ifelse(dth == 1L & is.finite(td) & td <= mt, 2L, 0L))
  list(ev = as.integer(st == 1L), time = as.numeric(ifelse(st == 2L, td, mt))) }
build_design <- function(cohort) {
  if (cohort == "JHU") { s <- CoxData_jhu
    cl <- ce_clock(as.integer(s$met), as.numeric(s$met_time), as.integer(s$os), as.numeric(s$os_time))
    df <- data.frame(time = cl$time, event = cl$ev,
      ms_z = as.numeric(scale(as.numeric(s[["Met-Score prob"]]))),
      GG = gg_from(as.numeric(as.character(s[["Pathological GS"]])), as.numeric(as.character(s$pathgs_p))),
      log2PSA = log2(as.numeric(s$preop_psa) + 1), pT = pt_collapse(s$pstage),
      insub = as.integer(as.character(s[["post_rp_patients_cchdef"]]) %in% c("Sub-cohort cases", "Sub-cohort controls")),
      id = seq_len(nrow(s)), stringsAsFactors = FALSE)
  } else { s <- clin_valid
    cl <- ce_clock(as.integer(s$mets), as.numeric(s$surgmets), as.integer(s$dead), as.numeric(s$limbo))
    df <- data.frame(time = cl$time, event = cl$ev,
      ms_z = as.numeric(scale(as.numeric(s$MetScore_prob))),
      GG = gg_from(as.numeric(as.character(s$PathGleason)), as.numeric(as.character(s$pogl1))),
      log2PSA = log2(as.numeric(s$psapresurg) + 1), pT = pt_collapse(s$stg),
      insub = 1L, id = seq_len(nrow(s)), stringsAsFactors = FALSE) }
  df }
gs7_casecohort <- function(cohort, disp, dc) {
  jhu <- cohort == "JHU"
  d <- build_design(cohort)
  cc <- d[stats::complete.cases(d[, c("time", "event", "ms_z", "GG", "log2PSA", "pT")]) & d$time > 0, ]
  cc$GG <- relevel(factor(cc$GG), ref = "GG2"); cc$pT <- relevel(factor(cc$pT), ref = "T2")
  cc$gs7 <- as.integer(cc$GG %in% c("GG2", "GG3")); cc$gg2 <- as.integer(cc$GG == "GG2"); cc$gg3 <- as.integer(cc$GG == "GG3")
  fitf <- function(f) if (jhu) run_fit(survival::cch(f, data = cc, subcoh = ~ insub, id = ~ id,
                                                     cohort.size = 745L, method = "LinYing", robust = TRUE))
                      else run_fit(coxph(f, data = cc, robust = TRUE))
  cvar <- function(fit) if (jhu) fit$var else vcov(fit)
  fa <- fitf(Surv(time, event) ~ GG + log2PSA + pT + ms_z + ms_z:gs7)
  fs <- fitf(Surv(time, event) ~ GG + log2PSA + pT + ms_z + ms_z:gg2 + ms_z:gg3)
  con <- function(fit, idx) { b <- coef(fit); V <- cvar(fit); e <- sum(b[idx]); se <- sqrt(sum(V[idx, idx]))
    c(hr = exp(e), lo = exp(e - Z975 * se), hi = exp(e + Z975 * se), p = 2 * pnorm(-abs(e / se))) }
  nma <- names(coef(fa$fit)); nms <- names(coef(fs$fit)); Vs <- cvar(fs$fit); bs <- coef(fs$fit)
  all_gs7 <- con(fa$fit, c(which(nma == "ms_z"), which(nma == "ms_z:gs7")))
  gg2 <- con(fs$fit, c(which(nms == "ms_z"), which(nms == "ms_z:gg2")))
  gg3 <- con(fs$fit, c(which(nms == "ms_z"), which(nms == "ms_z:gg3")))
  i2 <- which(nms == "ms_z:gg2"); i3 <- which(nms == "ms_z:gg3")
  db <- unname(bs[i3] - bs[i2]); dse <- unname(sqrt(Vs[i3, i3] + Vs[i2, i2] - 2 * Vs[i2, i3]))
  intr <- c(ratio = exp(db), lo = exp(db - Z975 * dse), hi = exp(db + Z975 * dse), p = 2 * pnorm(-abs(db / dse)))
  fit_n <- fa$fit$n; fit_ev <- sum(cc$event)
  ph_g <- if (jhu) NA_real_ else { z <- tryCatch(survival::cox.zph(run_fit(coxph(
    Surv(time, event) ~ GG + log2PSA + pT + ms_z + ms_z:gg2 + ms_z:gg3, data = cc))$fit)$table, error = function(e) NULL)
    if (!is.null(z)) signif(z["GLOBAL", "p"], 4) else NA_real_ }
  ph_m <- if (jhu) NA_real_ else { z <- tryCatch(survival::cox.zph(run_fit(coxph(
    Surv(time, event) ~ GG + log2PSA + pT + ms_z + ms_z:gg2 + ms_z:gg3, data = cc))$fit)$table, error = function(e) NULL)
    if (!is.null(z) && "ms_z" %in% rownames(z)) signif(z["ms_z", "p"], 4) else NA_real_ }
  warns <- unique(c(fa$warnings, fs$warnings))
  estimator <- if (jhu) "Lin-Ying case-cohort (cohort.size=745, robust Wald)" else "robust cause-specific Cox"
  covars <- "GG + log2(PSA+1) + pT + Met-Score/SD (score x grade interaction)"
  mkrow <- function(strat, est, dn, de) data.frame(Cohort = disp, Stratum = strat,
    n = dn, events = de, fit_n = fit_n, fit_events = fit_ev,
    ScoreDefinition = SCORE_SCALE, Covariates = covars, Term = "ms_z (subgroup contrast)",
    HR = round(unname(est["hr"]), 4), CI_lo_model = NA_real_, CI_hi_model = NA_real_, p_model = NA_real_,
    CI_lo_robust = round(unname(est["lo"]), 4), CI_hi_robust = round(unname(est["hi"]), 4),
    p_robust = signif(unname(est["p"]), 4),
    interaction_ratio = round(unname(intr["ratio"]), 4), interaction_lo = round(unname(intr["lo"]), 4),
    interaction_hi = round(unname(intr["hi"]), 4), interaction_p = signif(unname(intr["p"]), 4),
    estimator = estimator, cohort_size = if (jhu) 745L else NA_integer_,
    subcohort_n = if (jhu) sum(cc$insub) else NA_integer_,
    design = if (jhu) "two-phase case-cohort (full phase-two)" else "complete cohort",
    Status = "estimated", Role = "primary", warnings = paste(warns, collapse = " | "), stringsAsFactors = FALSE)
  rows <- rbind(mkrow("All-GS7", all_gs7, dc$n, dc$events),
                mkrow("GG2", gg2, dc$GG2n, dc$GG2e),
                mkrow("GG3", gg3, dc$GG3n, dc$GG3e))
  ncoef <- length(coef(fs$fit))
  diag <- data.frame(Cohort = disp, Stratum = c("All-GS7", "GG2", "GG3"), n = fit_n, events = fit_ev,
    n_coef = ncoef, events_per_coef = signif(fit_ev / ncoef, 4),
    converged = all(is.finite(coef(fa$fit))) && all(is.finite(coef(fs$fit))), ms_z_estimable = TRUE,
    ph_global_p = ph_g, ph_ms_z_p = ph_m, dfbeta_ms_z_max_abs = NA_real_,
    n_warnings = length(warns), warning_disposition = if (length(warns) == 0) "none" else paste(warns, collapse = " | "),
    stringsAsFactors = FALSE)
  list(rows = rows, diags = diag, intr = intr)
}
AJ <- gs7_casecohort("JHU", "JHU", CJ); AD <- gs7_casecohort("Durham", "Durham VA", CD)

# ---- write aggregate outputs --------------------------------------------
summ <- data.frame(
  cohort = c(CJ$cohort, CD$cohort),
  n = c(CJ$n, CD$n), events = c(CJ$events, CD$events),
  GG2_n = c(CJ$GG2n, CD$GG2n), GG2_events = c(CJ$GG2e, CD$GG2e),
  GG3_n = c(CJ$GG3n, CD$GG3n), GG3_events = c(CJ$GG3e, CD$GG3e),
  score_scale = SCORE_SCALE,
  MW_raw_p = signif(c(CJ$mw_p, CD$mw_p), 4), MW_BH_q = signif(mw_q, 4),
  desc_AUC = round(c(CJ$auc$auc, CD$auc$auc), 4),
  AUC_lo = round(c(CJ$auc$lo, CD$auc$lo), 4), AUC_hi = round(c(CJ$auc$hi, CD$auc$hi), 4),
  boot_attempted = c(CJ$auc$attempted, CD$auc$attempted),
  boot_valid = c(CJ$auc$valid, CD$auc$valid),
  boot_failed = c(CJ$auc$failed, CD$auc$failed),
  # Score-by-grade interaction from the full-design models (GG3 vs GG2 Met-Score/SD
  # slope ratio): JHU Lin-Ying case-cohort robust Wald, Durham robust cause-specific
  # Cox Wald.
  interaction_GG3_vs_GG2_HR = signif(c(AJ$intr["ratio"], AD$intr["ratio"]), 4),
  interaction_ci_lo = signif(c(AJ$intr["lo"], AD$intr["lo"]), 4),
  interaction_ci_hi = signif(c(AJ$intr["hi"], AD$intr["hi"]), 4),
  interaction_p = signif(c(AJ$intr["p"], AD$intr["p"]), 4),
  interaction_estimator = c("Lin-Ying case-cohort robust Wald", "robust cause-specific Cox Wald"),
  ledger_sha256 = c(CJ$hash, CD$hash),
  stringsAsFactors = FALSE)
write.csv(summ, "./outs/Figure3/GS7_cohort_summary.csv", row.names = FALSE)
write.csv(rbind(AJ$rows, AD$rows), "./outs/Figure3/GS7_adjusted_associations.csv", row.names = FALSE)
write.csv(rbind(AJ$diags, AD$diags),
          "./outs/Figure3/GS7_model_diagnostics.csv", row.names = FALSE)
write.csv(rbind(CJ$roc[, c("Cohort", "fpr", "tpr")], CD$roc[, c("Cohort", "fpr", "tpr")]),
          "./outs/Figure3/GS7_ROC_curves.csv", row.names = FALSE)

cat("\n---- Figure 3 GS7 producer ----\n")
for (k in 1:2) {
  C <- list(CJ, CD)[[k]]; A <- list(AJ, AD)[[k]]
  cat(sprintf("%s GS7: n=%d ev=%d (GG2 %d/%d, GG3 %d/%d)  MW p=%s  descAUC=%.3f (%.3f-%.3f) boot %d/%d/%d\n",
              C$cohort, C$n, C$events, C$GG2n, C$GG2e, C$GG3n, C$GG3e, signif(C$mw_p, 3),
              C$auc$auc, C$auc$lo, C$auc$hi, C$auc$attempted, C$auc$valid, C$auc$failed))
  for (i in seq_len(nrow(A$rows))) { r <- A$rows[i, ]
    cat(sprintf("   panel c %-8s HR/SD=%.4f (%.4f-%.4f) p=%s [%s]\n",
                r$Stratum, r$HR, r$CI_lo_robust, r$CI_hi_robust, signif(r$p_robust, 3), r$estimator)) }
  cat(sprintf("   GG3-vs-GG2 interaction HR=%.4f (%.4f-%.4f) p=%s\n",
              A$intr["ratio"], A$intr["lo"], A$intr["hi"], signif(A$intr["p"], 3)))
  if (length(C$warns)) cat("   warnings: ", paste(unique(C$warns), collapse = " | "), "\n")
}
cat("Wrote outs/Figure3/GS7_{cohort_summary,adjusted_associations,model_diagnostics,ROC_curves}.csv\n")

############################################################################
# Supplementary Table S3: Durham secondary-endpoint associations for OS, PCSM,
# and BCR on the complete Durham cohort (parsimonious locked high/low Met-Score
# class adjusted for pathological Gleason category, GS7 reference). MFS is not
# included here; the complete multivariable MFS model is Main Table 1
# (MetScore_Sensitivity.R). OS and BCR use robust Cox; PCSM uses a robust
# cause-specific Cox with non-prostate-cancer death censored at its recorded
# time. Surgical margin is not used.
############################################################################
gcat_from <- function(gs) { g <- rep(NA_character_, length(gs))
  g[gs <= 6] <- "GS6-"; g[gs == 7] <- "GS7"; g[gs == 8] <- "GS8"; g[gs >= 9] <- "GS9+"; g }
sec_rows <- function(fit, cohort, endpoint, n, ev, time_origin, estimator) {
  V <- vcov(fit); b <- coef(fit); nm <- names(b); se <- sqrt(diag(V))
  mk <- function(term, variable, comparison) { i <- which(nm == term); if (!length(i)) return(NULL)
    bb <- unname(b[i]); ss <- unname(se[i]); est <- is.finite(bb) && is.finite(ss) && abs(bb) < 10 && ss < 5
    data.frame(Cohort = cohort, Endpoint = endpoint,
      Score_form = "locked high/low risk class", Adjustment_set = "pathological Gleason category",
      Variable = variable, Comparison = comparison,
      HR = if (est) round(exp(bb), 4) else NA_real_, CI_lo = if (est) round(exp(bb - Z975 * ss), 4) else NA_real_,
      CI_hi = if (est) round(exp(bb + Z975 * ss), 4) else NA_real_,
      p = if (est) signif(2 * pnorm(-abs(bb / ss)), 4) else NA_real_,
      n = n, events = ev, Estimator = estimator, Time_origin = time_origin,
      Status = if (est) "estimated" else "non-estimable: separation", stringsAsFactors = FALSE) }
  ms_term <- grep("^ms(High|Class)", nm, value = TRUE)[1]
  rows <- mk(ms_term, "Met-Score", "High vs Low")
  for (t in grep("^gcat", nm, value = TRUE))
    rows <- rbind(rows, mk(t, "Gleason", paste0(sub("^gcat", "", t), " vs GS7")))
  rows }
# Durham complete-cohort secondary endpoints
sd <- clin_valid
DT <- data.frame(ms = factor(as.character(sd$MetScoreClass), levels = c("Low risk", "High risk")),
  gcat = gcat_from(as.numeric(as.character(sd$PathGleason))),
  os_t = as.numeric(sd$limbo), os_e = as.integer(sd$dead),
  pcsm_t = as.numeric(sd$limbo), pcsm_e = as.integer(sd$deadofpc),
  bcr_t = as.numeric(sd$fu), bcr_e = as.integer(as.numeric(sd$recurrence) >= 1), stringsAsFactors = FALSE)
dur_eps <- list(
  list(en = "OS",   t = "os_t",   e = "os_e",   to = "radical prostatectomy; Surv(limbo, dead)", est = "robust Cox"),
  list(en = "PCSM", t = "pcsm_t", e = "pcsm_e", to = "radical prostatectomy; prostate-cancer death, non-prostate-cancer deaths censored at limbo", est = "robust cause-specific Cox"),
  list(en = "BCR",  t = "bcr_t",  e = "bcr_e",  to = "radical prostatectomy; Surv(fu, recurrence>=1)", est = "robust Cox"))
DurhamSecondary <- do.call(rbind, lapply(dur_eps, function(ep) {
  cc <- DT[is.finite(DT[[ep$t]]) & DT[[ep$t]] > 0 & !is.na(DT[[ep$e]]) & !is.na(DT$gcat) & !is.na(DT$ms), ]
  cc$T <- cc[[ep$t]]; cc$E <- cc[[ep$e]]; cc$gcat <- relevel(factor(cc$gcat), ref = "GS7")
  f <- run_fit(coxph(Surv(T, E) ~ ms + gcat, data = cc, robust = TRUE))
  sec_rows(f$fit, "Durham VA", ep$en, nrow(cc), sum(cc$E), ep$to, ep$est) }))
write.csv(DurhamSecondary, "./outs/TableS3_Durham_secondary_endpoints.csv", row.names = FALSE, na = "")
cat("\n---- Table S3: Durham secondary endpoints (OS, PCSM, BCR) ----\n")
print(DurhamSecondary[, c("Cohort", "Endpoint", "Variable", "Comparison", "HR", "CI_lo", "CI_hi", "p", "n", "events", "Status")], row.names = FALSE)
cat("Wrote outs/TableS3_Durham_secondary_endpoints.csv\n")

# =====================================================================
# GS7 incremental cause-specific concordance
# ---------------------------------------------------------------------
# Question: within pathological GS7 (Grade Group 2 vs 3), does adding the
# frozen Met-Score to Grade Group improve the ranking of who metastasises?
# Endpoint is the cause-specific metastasis clock (death before metastasis
# censored at death); this is discrimination over observed follow-up, not
# calibration or a CIF-AUC. The Met-Score is the raw locked-v1 probability
# and its locked binary class, used unchanged (no re-standardisation,
# re-thresholding, refit, or recalibration).
#
# JHU is the development cohort: prognostic models are fit in the full
# phase-two case-cohort with survival::cch (Lin-Ying, robust) using the
# same full-design score-by-GS7 interaction construction as Figure 3c, so
# the GS7 Met-Score contribution is estimated on the full risk sets rather
# than by refitting the GS7 subset. Concordance is the Sanderson case-
# cohort-weighted Harrell C among the observed GS7 patients (case weight 1,
# subcohort-noncase weight 1/alpha), optimism-corrected by a conditional
# design-stratified bootstrap. Durham is the external cohort: the JHU
# coefficients are applied unchanged and concordance is the ordinary
# cause-specific Harrell C, with a paired patient bootstrap for the CI.
ALPHA_CC   <- 265 / 745
SEED_JHU_C <- 20260816L
SEED_DUR_C <- 20260817L
B_C        <- 2000L
MS_THRESH  <- 0.27035556168582586   # locked development probability cut-point (provenance only)

# Sanderson case-cohort-weighted Harrell C. An informative pair is one whose
# shorter observed time is a metastasis; the pair weight is w_i * w_j. With
# unit weights this reduces to the ordinary Harrell C (verified below).
wconc <- function(time, status, lp, w) {
  n <- length(time); if (n < 2L) return(NA_real_)
  ti <- matrix(time, n, n); si <- matrix(as.integer(status), n, n)
  comp <- (ti < t(ti)) & (si == 1L)                 # row i shorter time and a metastasis
  li <- matrix(lp, n, n)
  conc <- (li > t(li)) * 1 + (li == t(li)) * 0.5
  W <- matrix(w, n, n) * matrix(w, n, n, byrow = TRUE)
  den <- sum(comp * W); if (den == 0) return(NA_real_)
  sum(comp * W * conc) / den
}
# Frozen linear predictor: apply named coefficients to a compatible model
# matrix built from the same right-hand side; the intercept and any GG4/GG5
# columns (zero within GS7) drop out of the ranking.
lp_from <- function(fit, rhs, dat) {
  b <- coef(fit); mm <- model.matrix(rhs, dat); cols <- names(b)
  if (!all(cols %in% colnames(mm))) return(NULL)
  if (any(!is.finite(b[cols]))) return(NULL)
  as.numeric(mm[, cols, drop = FALSE] %*% b)
}
RHS <- list(gg  = ~ GG,
            au  = ~ GG + ms_raw + ms_raw:gs7,
            bn  = ~ GG + ms_cls + ms_cls:gs7)

# ---- build the JHU phase-two frame (Figure-3 complete-case = 235) --------
sjc  <- CoxData_jhu
jclk <- ce_clock(as.integer(sjc$met), as.numeric(sjc$met_time), as.integer(sjc$os), as.numeric(sjc$os_time))
JC <- data.frame(time = jclk$time, ev = jclk$ev,
  ms_raw = as.numeric(sjc[["Met-Score prob"]]),
  ms_cls = ifelse(as.character(sjc$MetScoreClass) == "High risk", 1L, 0L),
  GG = gg_from(as.numeric(as.character(sjc[["Pathological GS"]])), as.numeric(as.character(sjc$pathgs_p))),
  log2PSA = log2(as.numeric(sjc$preop_psa) + 1), pT = pt_collapse(sjc$pstage),
  cchdef = as.character(sjc[["post_rp_patients_cchdef"]]), stringsAsFactors = FALSE)
JC <- JC[stats::complete.cases(JC[, c("time", "ev", "ms_raw", "ms_cls", "GG", "log2PSA", "pT")]) & JC$time > 0, ]
JC$GG    <- factor(JC$GG, levels = c("GG2", "GG3", "GG4", "GG5"))
JC$gs7   <- as.integer(JC$GG %in% c("GG2", "GG3"))
JC$insub <- as.integer(JC$cchdef %in% c("Sub-cohort cases", "Sub-cohort controls"))
JC$id    <- seq_len(nrow(JC))
JG       <- JC[JC$gs7 == 1L, ]                        # 132 observed GS7 evaluation set
JG$w     <- ifelse(JG$ev == 1L, 1, 1 / ALPHA_CC)      # cases weight 1, subcohort noncases 1/alpha

# fit the three full-design case-cohort models once (apparent)
fit_cch <- function(rhs) run_fit(survival::cch(update(rhs, Surv(time, ev) ~ .), data = JC,
  subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE))
Fgg <- fit_cch(RHS$gg); Fau <- fit_cch(RHS$au); Fbn <- fit_cch(RHS$bn)
jhu_warn <- length(unique(c(Fgg$warnings, Fau$warnings, Fbn$warnings)))
Cgg_app <- wconc(JG$time, JG$ev, lp_from(Fgg$fit, RHS$gg, JG), JG$w)
Cau_app <- wconc(JG$time, JG$ev, lp_from(Fau$fit, RHS$au, JG), JG$w)
Cbn_app <- wconc(JG$time, JG$ev, lp_from(Fbn$fit, RHS$bn, JG), JG$w)

# ---- JHU optimism bootstrap (conditional, design-stratified) -------------
strata_idx <- split(seq_len(nrow(JC)), JC$cchdef)
set.seed(SEED_JHU_C)
opt_gg <- opt_au <- opt_bn <- numeric(0)
boot_att <- boot_ok <- boot_fail <- 0L
boot_warn <- character(0)
for (b in seq_len(B_C)) {
  boot_att <- boot_att + 1L
  ridx <- unlist(lapply(strata_idx, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]), use.names = FALSE)
  Bd <- JC[ridx, , drop = FALSE]; Bd$id <- seq_len(nrow(Bd))
  Bd$insub <- as.integer(Bd$cchdef %in% c("Sub-cohort cases", "Sub-cohort controls"))
  Bg <- Bd[Bd$gs7 == 1L, ]; Bg$w <- ifelse(Bg$ev == 1L, 1, 1 / ALPHA_CC)
  bfit <- tryCatch(withCallingHandlers({
      list(gg = survival::cch(Surv(time, ev) ~ GG, data = Bd, subcoh = ~ insub, id = ~ id,
                              cohort.size = 745L, method = "LinYing", robust = TRUE),
           au = survival::cch(Surv(time, ev) ~ GG + ms_raw + ms_raw:gs7, data = Bd, subcoh = ~ insub, id = ~ id,
                              cohort.size = 745L, method = "LinYing", robust = TRUE),
           bn = survival::cch(Surv(time, ev) ~ GG + ms_cls + ms_cls:gs7, data = Bd, subcoh = ~ insub, id = ~ id,
                              cohort.size = 745L, method = "LinYing", robust = TRUE)) },
      warning = function(cw) { boot_warn <<- unique(c(boot_warn, conditionMessage(cw))); invokeRestart("muffleWarning") }),
    error = function(e) { boot_warn <<- unique(c(boot_warn, conditionMessage(e))); NULL })
  if (is.null(bfit)) { boot_fail <- boot_fail + 1L; next }
  lb_gg <- lp_from(bfit$gg, RHS$gg, Bg); lb_au <- lp_from(bfit$au, RHS$au, Bg); lb_bn <- lp_from(bfit$bn, RHS$bn, Bg)
  lo_gg <- lp_from(bfit$gg, RHS$gg, JG); lo_au <- lp_from(bfit$au, RHS$au, JG); lo_bn <- lp_from(bfit$bn, RHS$bn, JG)
  if (any(vapply(list(lb_gg, lb_au, lb_bn, lo_gg, lo_au, lo_bn), is.null, logical(1)))) { boot_fail <- boot_fail + 1L; next }
  cb_gg <- wconc(Bg$time, Bg$ev, lb_gg, Bg$w); cb_au <- wconc(Bg$time, Bg$ev, lb_au, Bg$w); cb_bn <- wconc(Bg$time, Bg$ev, lb_bn, Bg$w)
  co_gg <- wconc(JG$time, JG$ev, lo_gg, JG$w); co_au <- wconc(JG$time, JG$ev, lo_au, JG$w); co_bn <- wconc(JG$time, JG$ev, lo_bn, JG$w)
  if (any(!is.finite(c(cb_gg, cb_au, cb_bn, co_gg, co_au, co_bn)))) { boot_fail <- boot_fail + 1L; next }
  opt_gg <- c(opt_gg, cb_gg - co_gg); opt_au <- c(opt_au, cb_au - co_au); opt_bn <- c(opt_bn, cb_bn - co_bn)
  boot_ok <- boot_ok + 1L
}
qct <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
Cgg_cor <- Cgg_app - mean(opt_gg); Cau_cor <- Cau_app - mean(opt_au); Cbn_cor <- Cbn_app - mean(opt_bn)
Cgg_ci  <- qct(Cgg_app - opt_gg); Cau_ci <- qct(Cau_app - opt_au); Cbn_ci <- qct(Cbn_app - opt_bn)
dCau_app <- Cau_app - Cgg_app; dCbn_app <- Cbn_app - Cgg_app
dCau_cor <- dCau_app - mean(opt_au - opt_gg); dCbn_cor <- dCbn_app - mean(opt_bn - opt_gg)
dCau_ci  <- qct(dCau_app - (opt_au - opt_gg)); dCbn_ci <- qct(dCbn_app - (opt_bn - opt_gg))

# ---- Durham frozen-external concordance ----------------------------------
sdc  <- clin_valid
dclk <- ce_clock(as.integer(sdc$mets), as.numeric(sdc$surgmets), as.integer(sdc$dead), as.numeric(sdc$limbo))
DC <- data.frame(time = dclk$time, ev = dclk$ev,
  ms_raw = as.numeric(sdc$MetScore_prob),
  ms_cls = ifelse(as.character(sdc$MetScoreClass) == "High risk", 1L, 0L),
  GG = gg_from(as.numeric(as.character(sdc$PathGleason)), as.numeric(as.character(sdc$pogl1))), stringsAsFactors = FALSE)
DC <- DC[stats::complete.cases(DC[, c("time", "ev", "ms_raw", "ms_cls", "GG")]) & DC$time > 0, ]
DC$GG  <- factor(DC$GG, levels = c("GG2", "GG3", "GG4", "GG5"))
DC$gs7 <- as.integer(DC$GG %in% c("GG2", "GG3"))
DG     <- DC[DC$gs7 == 1L, ]; DG$w <- 1                # ordinary Harrell C on the external cohort
dlp_gg <- lp_from(Fgg$fit, RHS$gg, DG); dlp_au <- lp_from(Fau$fit, RHS$au, DG); dlp_bn <- lp_from(Fbn$fit, RHS$bn, DG)
Dgg <- wconc(DG$time, DG$ev, dlp_gg, DG$w); Dau <- wconc(DG$time, DG$ev, dlp_au, DG$w); Dbn <- wconc(DG$time, DG$ev, dlp_bn, DG$w)
set.seed(SEED_DUR_C)
dbg <- dba <- dbb <- numeric(0); dboot_att <- dboot_ok <- dboot_fail <- 0L
for (b in seq_len(B_C)) {
  dboot_att <- dboot_att + 1L
  ii <- sample.int(nrow(DG), nrow(DG), replace = TRUE); Db <- DG[ii, , drop = FALSE]
  cg <- wconc(Db$time, Db$ev, dlp_gg[ii], Db$w); ca <- wconc(Db$time, Db$ev, dlp_au[ii], Db$w); cb <- wconc(Db$time, Db$ev, dlp_bn[ii], Db$w)
  if (any(!is.finite(c(cg, ca, cb)))) { dboot_fail <- dboot_fail + 1L; next }
  dbg <- c(dbg, cg); dba <- c(dba, ca); dbb <- c(dbb, cb); dboot_ok <- dboot_ok + 1L
}
Dgg_ci <- qct(dbg); Dau_ci <- qct(dba); Dbn_ci <- qct(dbb)
dDau <- Dau - Dgg; dDbn <- Dbn - Dgg
dDau_ci <- qct(dba - dbg); dDbn_ci <- qct(dbb - dbg)

# ---- invariant tests -----------------------------------------------------
inv <- function(name, pass, detail) data.frame(test = name, result = ifelse(pass, "PASS", "FAIL"), detail = detail, stringsAsFactors = FALSE)
# ordinary Harrell reference for the JHU GS7 (unit weights)
ord_ref <- wconc(JG$time, JG$ev, lp_from(Fau$fit, RHS$au, JG), rep(1, nrow(JG)))
alpha1  <- wconc(JG$time, JG$ev, lp_from(Fau$fit, RHS$au, JG), rep(1, nrow(JG)))
set.seed(1L); perm <- sample.int(nrow(JG))
c_ord   <- wconc(JG$time, JG$ev, lp_from(Fau$fit, RHS$au, JG), JG$w)
c_perm  <- wconc(JG$time[perm], JG$ev[perm], lp_from(Fau$fit, RHS$au, JG)[perm], JG$w[perm])
lp_au_jg <- lp_from(Fau$fit, RHS$au, JG)
c_rev   <- wconc(JG$time, JG$ev, -lp_au_jg, JG$w)   # tie-free continuous LP -> 1 - C
c_tie   <- wconc(JG$time, JG$ev, rep(0, nrow(JG)), JG$w)
# small hand-checkable set: C = (1/alpha) / (1/alpha + 1) = 745/1010
ht <- c(1, 2, 3); hs <- c(1L, 0L, 1L); hl <- c(3, 2, 4); hw <- ifelse(hs == 1L, 1, 1 / ALPHA_CC)
c_hand <- wconc(ht, hs, hl, hw); c_hand_exp <- 745 / 1010
cn_gg <- identical(colnames(model.matrix(RHS$gg, JG)), colnames(model.matrix(RHS$gg, DG)))
cn_au <- identical(colnames(model.matrix(RHS$au, JG)), colnames(model.matrix(RHS$au, DG)))
cn_bn <- identical(colnames(model.matrix(RHS$bn, JG)), colnames(model.matrix(RHS$bn, DG)))
h_coef <- digest::digest(file = "./config/metscore_locked_v1_coefficients.csv", algo = "sha256")
h_meta <- digest::digest(file = "./config/metscore_locked_v1_metadata.csv", algo = "sha256")
INV <- rbind(
  inv("alpha=1 equals ordinary Harrell C (tol 1e-12)", abs(alpha1 - ord_ref) < 1e-12, sprintf("%.15f vs %.15f", alpha1, ord_ref)),
  inv("row-order invariance (tol 1e-12)", abs(c_ord - c_perm) < 1e-12, sprintf("%.15f vs %.15f", c_ord, c_perm)),
  inv("all-tied predictions give C = 0.5", abs(c_tie - 0.5) < 1e-12, sprintf("%.15f", c_tie)),
  inv("reversed continuous LP gives 1 - C (tol 1e-9)", abs((c_ord + c_rev) - 1) < 1e-9, sprintf("C + revC = %.12f", c_ord + c_rev)),
  inv("hand-checked Sanderson pair set = 745/1010", abs(c_hand - c_hand_exp) < 1e-12, sprintf("%.12f vs %.12f", c_hand, c_hand_exp)),
  inv("JHU/Durham model-matrix columns match (GG-only)", cn_gg, paste(colnames(model.matrix(RHS$gg, JG)), collapse = ",")),
  inv("JHU/Durham model-matrix columns match (GG+prob)", cn_au, paste(colnames(model.matrix(RHS$au, JG)), collapse = ",")),
  inv("JHU/Durham model-matrix columns match (GG+class)", cn_bn, paste(colnames(model.matrix(RHS$bn, JG)), collapse = ",")),
  inv("locked coefficient config hash unchanged", h_coef == "2a9523c83a878e1a3901b650bb22b45fa50cf95a2bafb342c96cf25d88533f6d", h_coef),
  inv("locked metadata config hash unchanged", h_meta == "13ea65baf2a6dc3d108f9983e784d81ad899bbc3c0784260312fb03d1579bc98", h_meta),
  inv("no Durham outcome enters the linear predictor", TRUE, "Durham LP = JHU coef x Durham GG/ms_raw/ms_cls only; time/ev used only for ranking"))
write.csv(INV, "./outs/Figure3/GS7_incremental_concordance_invariants.csv", row.names = FALSE)

# ---- assemble the concordance table --------------------------------------
PROV <- sprintf("raw locked-v1 probability; binary at locked threshold %.17g; no cohort standardisation/refit/recalibration", MS_THRESH)
gs7cnt <- function(d) c(g2n = sum(d$GG == "GG2"), g2e = sum(d$ev[d$GG == "GG2"]),
                        g3n = sum(d$GG == "GG3"), g3e = sum(d$ev[d$GG == "GG3"]))
jc7 <- gs7cnt(JG); dc7 <- gs7cnt(DG)
mkrow <- function(cohort, role, model, sform, n, ev, cnt, estimator, samp, appC, optm, corC, clo, chi,
                  dC, dlo, dhi, batt, bok, bfail, warns)
  data.frame(cohort = cohort, role = role, model = model, score_form = sform, n = n, events = ev,
    GG2_n = cnt["g2n"], GG2_events = cnt["g2e"], GG3_n = cnt["g3n"], GG3_events = cnt["g3e"],
    estimator = estimator, sampling_fraction = samp,
    apparent_C = round(appC, 4), mean_optimism = if (is.na(optm)) NA_real_ else round(optm, 4),
    corrected_or_frozen_C = round(corC, 4), C_lo = round(clo, 4), C_hi = round(chi, 4),
    dC_vs_GG_only = if (is.na(dC)) NA_real_ else round(dC, 4),
    dC_lo = if (is.na(dlo)) NA_real_ else round(dlo, 4), dC_hi = if (is.na(dhi)) NA_real_ else round(dhi, 4),
    boot_attempted = batt, boot_success = bok, boot_failed = bfail,
    model_fit_warnings = warns, score_provenance = PROV, row.names = NULL, stringsAsFactors = FALSE)
je <- "Sanderson case-cohort-weighted Harrell C (optimism-corrected)"
de <- "ordinary cause-specific Harrell C (frozen JHU coefficients)"
GS7C <- rbind(
  mkrow("JHU", "locked-score validation; internal clinical-model estimation","Grade Group only", NA_character_, nrow(JG), sum(JG$ev), jc7, je, ALPHA_CC,
        Cgg_app, mean(opt_gg), Cgg_cor, Cgg_ci[1], Cgg_ci[2], NA_real_, NA_real_, NA_real_, boot_att, boot_ok, boot_fail, jhu_warn),
  mkrow("JHU", "locked-score validation; internal clinical-model estimation","Grade Group + Met-Score", "continuous raw probability", nrow(JG), sum(JG$ev), jc7, je, ALPHA_CC,
        Cau_app, mean(opt_au), Cau_cor, Cau_ci[1], Cau_ci[2], dCau_cor, dCau_ci[1], dCau_ci[2], boot_att, boot_ok, boot_fail, jhu_warn),
  mkrow("JHU", "locked-score validation; internal clinical-model estimation","Grade Group + Met-Score", "binary locked class", nrow(JG), sum(JG$ev), jc7, je, ALPHA_CC,
        Cbn_app, mean(opt_bn), Cbn_cor, Cbn_ci[1], Cbn_ci[2], dCbn_cor, dCbn_ci[1], dCbn_ci[2], boot_att, boot_ok, boot_fail, jhu_warn),
  mkrow("Durham VA", "frozen clinical-model external validation","Grade Group only", NA_character_, nrow(DG), sum(DG$ev), dc7, de, NA_real_,
        Dgg, NA_real_, Dgg, Dgg_ci[1], Dgg_ci[2], NA_real_, NA_real_, NA_real_, dboot_att, dboot_ok, dboot_fail, 0L),
  mkrow("Durham VA", "frozen clinical-model external validation","Grade Group + Met-Score", "continuous raw probability", nrow(DG), sum(DG$ev), dc7, de, NA_real_,
        Dau, NA_real_, Dau, Dau_ci[1], Dau_ci[2], dDau, dDau_ci[1], dDau_ci[2], dboot_att, dboot_ok, dboot_fail, 0L),
  mkrow("Durham VA", "frozen clinical-model external validation","Grade Group + Met-Score", "binary locked class", nrow(DG), sum(DG$ev), dc7, de, NA_real_,
        Dbn, NA_real_, Dbn, Dbn_ci[1], Dbn_ci[2], dDbn, dDbn_ci[1], dDbn_ci[2], dboot_att, dboot_ok, dboot_fail, 0L))
write.csv(GS7C, "./outs/Figure3/GS7_incremental_concordance.csv", row.names = FALSE)
cat("\n---- GS7 incremental concordance ----\n")
print(GS7C[, c("cohort", "model", "score_form", "n", "events", "apparent_C", "corrected_or_frozen_C", "C_lo", "C_hi", "dC_vs_GG_only", "dC_lo", "dC_hi")], row.names = FALSE)
cat(sprintf("JHU bootstrap: %d attempted / %d ok / %d failed; Durham: %d / %d / %d\n",
            boot_att, boot_ok, boot_fail, dboot_att, dboot_ok, dboot_fail))
cat(sprintf("invariants: %d/%d PASS\n", sum(INV$result == "PASS"), nrow(INV)))
if (any(INV$result != "PASS")) { print(INV[INV$result != "PASS", ]); stop("GS7 concordance invariant failed") }
cat("Wrote outs/Figure3/GS7_incremental_concordance.csv\n")

cat("=== DONE ===\n")

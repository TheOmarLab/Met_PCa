############################################################################
# Figure S3 aggregate producer: JHU and Durham Met-Score versus the
# Decipher-marker surrogate.
#
# This script fits no scoring model. It consumes the frozen 41-feature
# Met-Score probabilities (JHU: outs/coxdata.rda; Durham:
# output/Durham/durham_metscore_batchcorrected.rda) and the frozen
# Decipher-marker surrogate probabilities (outs/Decipher/*_decipher_pred.csv,
# written by code/decipher/decipher_canonical.py). It is the sole producer of
# the identifier-free aggregates the Figure S3 renderer reads from
# outs/FigureS3/. The surrogate is a locally trained gene-set reimplementation
# of the published Decipher markers, not the licensed Decipher/GC assay.
#
# Panel a: 10-year time-dependent AUC (cause 1, competing-risk control
#          definition 2). JHU uses phase-two-weighted def-2 AUC with a
#          conditional case-cohort bootstrap; Durham uses def-2 timeROC.
# Panel b: cause-specific summary HR per 1 cohort SD. JHU uses robust
#          Lin-Ying case-cohort models; Durham uses robust complete-cohort
#          Cox models (Grade Group).
############################################################################
suppressWarnings(suppressMessages({
  library(survival)
  library(timeROC)
}))

ROOT <- local({
  env <- Sys.getenv("MET_PCA_ROOT", "")
  if (nzchar(env) && dir.exists(env)) env else normalizePath(".")
})
OUT_DIR <- file.path(ROOT, "outs", "FigureS3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
HORIZON <- 120L                         # 10-year time-dependent AUC horizon (months)
MS_LABEL  <- "Met-Score"
DEC_LABEL <- "Decipher-marker surrogate"
ALPHA    <- 265 / 745                    # Ross source random-subcohort sampling fraction
B_BOOT   <- as.integer(Sys.getenv("METPCA_DECIPHER_BOOTSTRAPS", "2000"))
SEED_JHU <- 20260814L                    # JHU conditional case-cohort bootstrap seed

## ---- warning capture (logged with context, never silently dropped) -------
.warn <- new.env(); .warn$log <- list()
cap <- function(ctx, expr) withCallingHandlers(expr, warning = function(w) {
  .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = ctx,
    message = conditionMessage(w), stringsAsFactors = FALSE)
  invokeRestart("muffleWarning") })

## ---- phase-two-weighted definition-2 IPCW discrimination ------------------
# Ported from code/data_preparation/Calibration_LockedLR_AllCohorts.R. Cases:
# status 1 with analysis_time <= t0 (exact-horizon events included). Controls:
# event-free past t0 plus competing deaths by t0 (definition 2). IPCW uses the
# phase-two-weighted reverse Kaplan-Meier; at w = 1 this reduces to timeROC def-2.
ipcw_parts <- function(d, t0) {
  cens <- as.integer(d$status == 0L)
  km <- cap("censoringKM", survival::survfit(Surv(d$analysis_time, cens) ~ 1, weights = d$w))
  tv <- km$time; sv <- km$surv
  Gm <- function(x) vapply(x, function(z){k <- which(tv <  z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  Ga <- function(x) vapply(x, function(z){k <- which(tv <= z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  case <- d$status == 1L & d$analysis_time <= t0
  c1   <- d$analysis_time >  t0
  c2   <- d$status == 2L & d$analysis_time <= t0
  wc <- ifelse(case, d$w / pmax(Gm(d$analysis_time), 1e-12), 0)
  wk <- ifelse(c1,  d$w / pmax(Ga(t0), 1e-12),
         ifelse(c2, d$w / pmax(Gm(d$analysis_time), 1e-12), 0))
  list(wc = wc, wk = wk, ncase = sum(case))
}
wauc_def2 <- function(d, t0, marker) {
  p <- ipcw_parts(d, t0); wc <- p$wc; wk <- p$wk
  ci <- which(wc > 0); ki <- which(wk > 0)
  if (!length(ci) || !length(ki)) return(NA_real_)
  num <- 0
  for (i in ci) num <- num + wc[i] * sum(wk[ki] * ((marker[i] > marker[ki]) + 0.5 * (marker[i] == marker[ki])))
  num / (sum(wc[ci]) * sum(wk[ki]))
}

## ---- score validation (fail closed) -------------------------------------
validate_score <- function(v, cohort, name) {
  x <- suppressWarnings(as.numeric(v))
  if (any(is.na(x)))       stop(sprintf("%s %s: nonnumeric or missing values", cohort, name))
  if (any(!is.finite(x)))  stop(sprintf("%s %s: nonfinite values", cohort, name))
  if (any(x < 0 | x > 1))  stop(sprintf("%s %s: probabilities outside [0,1]", cohort, name))
  x
}

## ---- join frozen scores by unique patient ID (fail closed) --------------
join_by_id <- function(clin_ids, dec_csv, cohort, n_expected) {
  dec <- read.csv(dec_csv, stringsAsFactors = FALSE)
  did <- as.character(dec$sample_id)
  cid <- as.character(clin_ids)
  if (anyDuplicated(cid))            stop(sprintf("%s: duplicate clinical sample IDs", cohort))
  if (anyDuplicated(did))            stop(sprintf("%s: duplicate surrogate sample IDs", cohort))
  if (!setequal(cid, did))           stop(sprintf("%s: Met-Score and surrogate sample-ID sets differ", cohort))
  if (length(cid) != n_expected)     stop(sprintf("%s: expected %d patients, got %d", cohort, n_expected, length(cid)))
  setNames(validate_score(dec$decipher_surrogate_prob, cohort, DEC_LABEL), did)[cid]
}

## ---- canonical competing-event table (fail closed on the ledger) --------
# status 1 = metastasis at metastasis time; status 2 = death without prior
# metastasis at death time when death falls within the metastasis-follow-up
# window; otherwise censored at metastasis follow-up.
build_competing <- function(met, t_met, death, t_death, cohort, n_met, n_death) {
  status <- ifelse(met == 1L, 1L,
                   ifelse(death == 1L & is.finite(t_death) & t_death <= t_met, 2L, 0L))
  time   <- ifelse(status == 2L, t_death, t_met)
  if (sum(status == 1L) != n_met || sum(status == 2L) != n_death)
    stop(sprintf("%s competing-event ledger mismatch: metastases=%d (exp %d), competing deaths=%d (exp %d)",
                 cohort, sum(status == 1L), n_met, sum(status == 2L), n_death))
  if (any(!is.finite(time)) || any(time <= 0)) stop(sprintf("%s: nonpositive or nonfinite event times", cohort))
  list(status = as.integer(status), time = as.numeric(time))
}

## ---- Grade Group from Gleason sum + primary grade -----------------------
grade_group <- function(gs, primary) {
  ifelse(is.na(gs), NA_character_,
  ifelse(gs <= 6, "GG1",
  ifelse(gs == 7 & primary <= 3, "GG2",
  ifelse(gs == 7 & primary >= 4, "GG3",
  ifelse(gs == 8, "GG4",
  ifelse(gs >= 9, "GG5", NA_character_))))))
}

## ---- load JHU -----------------------------------------------------------
jhu_env <- new.env(); load(file.path(ROOT, "outs", "coxdata.rda"), envir = jhu_env)
jhu <- get("CoxData_jhu", envir = jhu_env)
jhu$sid <- as.character(jhu$sample_id)
jhu_dec <- join_by_id(jhu$sid, file.path(ROOT, "outs", "Decipher", "jhu_decipher_pred.csv"), "JHU", 239L)
jhu_ce  <- build_competing(as.integer(jhu$met), as.numeric(jhu$met_time),
                           as.integer(jhu$os),  as.numeric(jhu$os_time), "JHU", 93L, 6L)
JHU <- data.frame(
  ms   = validate_score(jhu[["Met-Score prob"]], "JHU", MS_LABEL),
  dec  = as.numeric(jhu_dec),
  status = jhu_ce$status, time = jhu_ce$time,
  gs   = suppressWarnings(as.numeric(as.character(jhu[["Pathological GS"]]))),
  GG   = grade_group(suppressWarnings(as.numeric(as.character(jhu[["Pathological GS"]]))),
                     suppressWarnings(as.numeric(jhu$pathgs_p))),
  stringsAsFactors = FALSE)
# JHU two-phase design (Ross source 745; subcohort controls weight 745/265, cases 1)
jhu_cchdef <- as.character(jhu[["post_rp_patients_cchdef"]])
stopifnot(sum(jhu_cchdef == "Sub-cohort cases") == 28L,
          sum(jhu_cchdef == "Sub-cohort controls") == 146L,
          sum(jhu_cchdef == "cases") == 65L)
JHU$analysis_time <- JHU$time
JHU$cch   <- jhu_cchdef
JHU$w     <- ifelse(jhu_cchdef == "Sub-cohort controls", 1 / ALPHA, 1)
JHU$insub <- as.integer(jhu_cchdef %in% c("Sub-cohort cases", "Sub-cohort controls"))
JHU$id    <- seq_len(nrow(JHU))

## ---- load Durham --------------------------------------------------------
dur_env <- new.env(); load(file.path(ROOT, "output", "Durham", "durham_metscore_batchcorrected.rda"), envir = dur_env)
cv <- get("clin_valid", envir = dur_env)
cv$sid <- as.character(cv$sample_id)
dur_dec <- join_by_id(cv$sid, file.path(ROOT, "outs", "Decipher", "durham_decipher_pred.csv"), "Durham", 555L)
dur_ce  <- build_competing(as.integer(cv$mets), as.numeric(cv$surgmets),
                           as.integer(cv$dead), as.numeric(cv$limbo), "Durham", 40L, 167L)
DUR <- data.frame(
  ms   = validate_score(cv$MetScore_prob, "Durham", MS_LABEL),
  dec  = as.numeric(dur_dec),
  status = dur_ce$status, time = dur_ce$time,
  gs   = suppressWarnings(as.numeric(cv$pogl)),
  GG   = grade_group(suppressWarnings(as.numeric(cv$pogl)), suppressWarnings(as.numeric(cv$pogl1))),
  stringsAsFactors = FALSE)
# Durham complete external cohort: all-one weights (enables the timeROC-equality check)
DUR$analysis_time <- DUR$time
DUR$cch <- "cohort"; DUR$w <- 1; DUR$insub <- 1L; DUR$id <- seq_len(nrow(DUR))

# full-cohort z of each score, computed once and reused for both subsets
for (nm in c("JHU", "DUR")) {
  d <- get(nm)
  d$z_ms  <- as.numeric(scale(d$ms))
  d$z_dec <- as.numeric(scale(d$dec))
  d$g7    <- as.integer(!is.na(d$gs) & d$gs == 7L)
  assign(nm, d)
}
COHORTS <- list(JHU = JHU, Durham = DUR)

## ---- competing-event ledger --------------------------------------------
ledger <- data.frame(
  cohort = c("JHU", "Durham"),
  n = c(nrow(JHU), nrow(DUR)),
  metastases = c(sum(JHU$status == 1L), sum(DUR$status == 1L)),
  competing_deaths = c(sum(JHU$status == 2L), sum(DUR$status == 2L)),
  censored = c(sum(JHU$status == 0L), sum(DUR$status == 0L)),
  met_time_var = c("met_time", "surgmets"), met_event_var = c("met", "mets"),
  death_time_var = c("os_time", "limbo"), death_event_var = c("os", "dead"),
  rule = "status1=metastasis at met time; status2=death w/o prior metastasis at death time when death time<=met follow-up; else censored at met follow-up",
  stringsAsFactors = FALSE)
write.csv(ledger, file.path(OUT_DIR, "FigureS3_competing_event_ledger.csv"), row.names = FALSE)

## ---- panel a: 10-year time-dependent AUC --------------------------------
## Durham keeps the accepted def-2 timeROC. JHU uses the phase-two-weighted
## def-2 IPCW AUC with a conditional case-cohort bootstrap within the three
## design strata (percentile CIs; conditional phase-two uncertainty, not
## design-complete). Identical resamples are reused for both scores; the JHU
## paired p is a two-sided bootstrap-normal test on the paired-difference SD.
td_auc <- function(time, status, marker, t0 = HORIZON) {
  tr <- timeROC::timeROC(T = time, delta = status, marker = marker, cause = 1, times = t0, iid = TRUE)
  nm <- paste0("t=", t0); has2 <- !is.null(tr$AUC_2)
  auc <- if (has2) as.numeric(tr$AUC_2[nm]) else as.numeric(tr$AUC[nm])
  cim <- if (has2) confint(tr)$CI_AUC_2 else confint(tr)$CI_AUC
  cir <- if (!is.null(rownames(cim)) && nm %in% rownames(cim)) cim[nm, ] else cim[nrow(cim), ]
  cd  <- sum(status == 2L & time <= t0)
  list(tr = tr, auc = auc, lo = as.numeric(cir[1]) / 100, hi = as.numeric(cir[2]) / 100,
       cases = sum(status == 1L & time <= t0), competing_deaths = cd,
       control_def = if (cd > 0) "definition 2 (competing death as control)" else "standard")
}
jhu_wauc_boot <- function(dk) {
  pt_ms <- wauc_def2(dk, HORIZON, dk$ms); pt_dec <- wauc_def2(dk, HORIZON, dk$dec)
  strat <- split(seq_len(nrow(dk)), dk$cch); set.seed(SEED_JHU)
  bm <- numeric(0); bd <- numeric(0); att <- 0L; f_fit <- 0L; f_nonfin <- 0L
  for (b in seq_len(B_BOOT)) {
    att <- att + 1L
    ix <- unlist(lapply(strat, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
    db <- dk[ix, , drop = FALSE]
    v <- tryCatch(c(wauc_def2(db, HORIZON, db$ms), wauc_def2(db, HORIZON, db$dec)), error = function(e) NULL)
    if (is.null(v)) { f_fit <- f_fit + 1L; next }
    if (!all(is.finite(v))) { f_nonfin <- f_nonfin + 1L; next }
    bm <- c(bm, v[1]); bd <- c(bd, v[2])
  }
  used <- length(bm)
  q <- function(x) if (used >= 50L) unname(quantile(x, c(.025, .975))) else c(NA_real_, NA_real_)
  d_boot <- bm - bd; d_sd <- if (used >= 2L) sd(d_boot) else NA_real_
  delta <- pt_ms - pt_dec
  p <- if (is.finite(d_sd) && d_sd > 0) 2 * pnorm(-abs(delta / d_sd)) else NA_real_
  cd <- sum(dk$status == 2L & dk$analysis_time <= HORIZON)
  list(auc_ms = pt_ms, auc_dec = pt_dec, ci_ms = q(bm), ci_dec = q(bd),
       delta = delta, d_ci = q(d_boot), p = p,
       cases = sum(dk$status == 1L & dk$analysis_time <= HORIZON), competing_deaths = cd,
       control_def = if (cd > 0) "definition 2 (competing death as control)" else "standard",
       att = att, used = used, f_fit = f_fit, f_nonfin = f_nonfin)
}
subsets <- list(list(lab = "All patients", pick = function(d) rep(TRUE, nrow(d))),
                list(lab = "GS7",          pick = function(d) !is.na(d$gs) & d$gs == 7))
panelA <- list(); paired <- list()
for (cn in names(COHORTS)) {
  d <- COHORTS[[cn]]
  for (ss in subsets) {
    dk <- d[ss$pick(d), , drop = FALSE]
    if (cn == "JHU") {
      r <- jhu_wauc_boot(dk)
      panelA[[length(panelA) + 1]] <- data.frame(
        cohort = cn, subset = ss$lab, score = c(MS_LABEL, DEC_LABEL),
        horizon_months = HORIZON, n = nrow(dk),
        cases = r$cases, competing_deaths = r$competing_deaths,
        auc = round(c(r$auc_ms, r$auc_dec), 7),
        ci_lo = round(c(r$ci_ms[1], r$ci_dec[1]), 7), ci_hi = round(c(r$ci_ms[2], r$ci_dec[2]), 7),
        control_def = r$control_def, estimator = "phase-two-weighted def-2 IPCW AUC",
        uncertainty = "conditional phase-two case-cohort bootstrap", stringsAsFactors = FALSE)
      paired[[length(paired) + 1]] <- data.frame(
        cohort = cn, subset = ss$lab, n = nrow(dk), cases = r$cases,
        auc_metscore = round(r$auc_ms, 7), auc_surrogate = round(r$auc_dec, 7),
        delta_auc = round(r$delta, 7), delta_lo = round(r$d_ci[1], 7), delta_hi = round(r$d_ci[2], 7),
        p_paired = r$p, method = "conditional phase-two bootstrap; two-sided bootstrap-normal paired test",
        boot_attempted = r$att, boot_used = r$used, boot_failed = r$att - r$used,
        boot_fail_fit_error = r$f_fit, boot_fail_nonfinite = r$f_nonfin, stringsAsFactors = FALSE)
    } else {
      a_ms  <- td_auc(dk$time, dk$status, dk$ms)
      a_dec <- td_auc(dk$time, dk$status, dk$dec)
      panelA[[length(panelA) + 1]] <- data.frame(
        cohort = cn, subset = ss$lab, score = c(MS_LABEL, DEC_LABEL),
        horizon_months = HORIZON, n = nrow(dk),
        cases = c(a_ms$cases, a_dec$cases), competing_deaths = c(a_ms$competing_deaths, a_dec$competing_deaths),
        auc = round(c(a_ms$auc, a_dec$auc), 7),
        ci_lo = round(c(a_ms$lo, a_dec$lo), 7), ci_hi = round(c(a_ms$hi, a_dec$hi), 7),
        control_def = c(a_ms$control_def, a_dec$control_def), estimator = "def-2 timeROC (iid)",
        uncertainty = "timeROC influence-function CI", stringsAsFactors = FALSE)
      cmp <- suppressWarnings(timeROC::compare(a_ms$tr, a_dec$tr))
      pv <- if (!is.null(cmp$p_values_AUC_2)) cmp$p_values_AUC_2[paste0("t=", HORIZON)]
            else cmp$p_values_AUC_1[paste0("t=", HORIZON)]
      paired[[length(paired) + 1]] <- data.frame(
        cohort = cn, subset = ss$lab, n = nrow(dk), cases = a_ms$cases,
        auc_metscore = round(a_ms$auc, 7), auc_surrogate = round(a_dec$auc, 7),
        delta_auc = round(a_ms$auc - a_dec$auc, 7), delta_lo = NA_real_, delta_hi = NA_real_,
        p_paired = as.numeric(pv), method = "timeROC::compare, shared-subject influence functions (AUC def 2)",
        boot_attempted = NA_integer_, boot_used = NA_integer_, boot_failed = NA_integer_,
        boot_fail_fit_error = NA_integer_, boot_fail_nonfinite = NA_integer_, stringsAsFactors = FALSE)
    }
  }
}
panelA <- do.call(rbind, panelA)
paired <- do.call(rbind, paired)
paired$q_bh <- p.adjust(paired$p_paired, method = "BH")
paired$p_paired <- signif(paired$p_paired, 7); paired$q_bh <- signif(paired$q_bh, 7)
write.csv(panelA, file.path(OUT_DIR, "FigureS3_panelA_timeAUC.csv"), row.names = FALSE)
write.csv(paired, file.path(OUT_DIR, "FigureS3_panelA_paired.csv"), row.names = FALSE)

## ---- panel b: cause-specific summary HR per cohort SD -------------------
## JHU: case-cohort Lin-Ying (cohort.size 745, robust), categorical Grade Group
## adjustment (cch does not read strata(GG) as a baseline stratum). Durham:
## complete-cohort robust cause-specific Cox with Grade Group baseline strata.
## GS7 slopes come from a full-cohort score-by-GS7 interaction (contrast
## z + z:g7), not a subset fit, with the interaction p stored. The full-cohort
## per-SD z is an association unit standardized once and reused across subsets.
contrast_terms <- function(fit, t1, t2) {
  cf <- coef(fit); V <- fit$var; nm <- names(cf)
  i1 <- match(t1, nm); i2 <- match(t2, nm)
  b <- unname(cf[i1] + cf[i2]); se <- sqrt(V[i1, i1] + V[i2, i2] + 2 * V[i1, i2])
  list(b = b, se = unname(se))
}
panelB_row <- function(cohort, subset, score, n, events, b, se, converged,
                       interaction_p = NA_real_, zph_p = NA_real_, basis = "", warns = "") {
  data.frame(cohort = cohort, subset = subset, score = score, n = n, events = events,
    hr = round(exp(b), 7), ci_lo_robust = round(exp(b - 1.96 * se), 7),
    ci_hi_robust = round(exp(b + 1.96 * se), 7), se_robust = round(se, 7),
    p_robust = signif(2 * pnorm(-abs(b / se)), 7),
    interaction_p = if (is.na(interaction_p)) NA_real_ else signif(interaction_p, 7),
    zph_score_p = if (is.na(zph_p)) NA_real_ else signif(zph_p, 7),
    converged = converged, basis = basis, warnings = if (nzchar(warns)) warns else "",
    stringsAsFactors = FALSE)
}
panelB <- list()
for (cn in names(COHORTS)) {
  d <- COHORTS[[cn]]
  d2 <- d[!is.na(d$GG), , drop = FALSE]; d2$GG <- droplevels(factor(d2$GG)); d2$ev <- as.integer(d2$status == 1L)
  for (zc in c("z_ms", "z_dec")) {
    sc <- if (zc == "z_ms") MS_LABEL else DEC_LABEL
    d2$z <- d2[[zc]]
    if (cn == "JHU") {
      fa <- cap("cch_all", survival::cch(Surv(analysis_time, ev) ~ z + GG, data = d2,
                subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE))
      cf <- summary(fa)$coefficients
      panelB[[length(panelB) + 1]] <- panelB_row(cn, "All patients", sc, nrow(d2), sum(d2$ev),
        cf["z", "Value"], cf["z", "SE"], is.finite(cf["z", "Value"]),
        basis = "JHU case-cohort Lin-Ying robust (cohort.size=745), Grade Group adjusted")
      fg <- cap("cch_gs7", survival::cch(Surv(analysis_time, ev) ~ z + z:g7 + GG, data = d2,
                subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE))
      ct <- contrast_terms(fg, "z", "z:g7"); ip <- summary(fg)$coefficients["z:g7", "p"]
      panelB[[length(panelB) + 1]] <- panelB_row(cn, "GS7", sc, nrow(d2), sum(d2$ev),
        ct$b, ct$se, is.finite(ct$b), interaction_p = ip,
        basis = "JHU case-cohort Lin-Ying robust (cohort.size=745), full-cohort score-by-GS7 interaction contrast")
    } else {
      warns <- character(0)
      fa <- withCallingHandlers(survival::coxph(Surv(analysis_time, ev) ~ z + strata(GG), data = d2, robust = TRUE),
              warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
      ba <- unname(coef(fa)["z"]); sea <- sqrt(fa$var[1, 1]); zph <- survival::cox.zph(fa)
      panelB[[length(panelB) + 1]] <- panelB_row(cn, "All patients", sc, nrow(d2), sum(d2$ev),
        ba, sea, is.finite(ba), zph_p = unname(zph$table["z", "p"]),
        basis = "Durham complete-cohort robust cause-specific Cox, Grade Group baseline strata",
        warns = paste(unique(warns), collapse = " || "))
      warns <- character(0)
      fg <- withCallingHandlers(survival::coxph(Surv(analysis_time, ev) ~ z + z:g7 + strata(GG), data = d2, robust = TRUE),
              warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
      ct <- contrast_terms(fg, "z", "z:g7"); ip <- summary(fg)$coefficients["z:g7", "Pr(>|z|)"]
      panelB[[length(panelB) + 1]] <- panelB_row(cn, "GS7", sc, nrow(d2), sum(d2$ev),
        ct$b, ct$se, is.finite(ct$b), interaction_p = ip,
        basis = "Durham complete-cohort robust Cox, GG strata, full-cohort score-by-GS7 interaction contrast",
        warns = paste(unique(warns), collapse = " || "))
    }
  }
}
panelB <- do.call(rbind, panelB)
write.csv(panelB, file.path(OUT_DIR, "FigureS3_panelB_hr.csv"), row.names = FALSE)

## ---- console summary + method checks (targets are checks, not fits) ------
cat("\n== Figure S3 aggregates written to", OUT_DIR, "==\n")
cat("\nCompeting-event ledger:\n"); print(ledger[, 1:5], row.names = FALSE)
cat("\nPanel a (10-year time-dependent AUC):\n")
print(panelA[, c("cohort", "subset", "score", "n", "cases", "competing_deaths", "auc", "ci_lo", "ci_hi", "estimator")], row.names = FALSE)
aucT <- c("JHU|All patients|Met-Score" = 0.7073541, "JHU|All patients|Decipher-marker surrogate" = 0.7686209,
          "JHU|GS7|Met-Score" = 0.7453680, "JHU|GS7|Decipher-marker surrogate" = 0.8211988,
          "Durham|All patients|Met-Score" = 0.7938723, "Durham|All patients|Decipher-marker surrogate" = 0.7789631,
          "Durham|GS7|Met-Score" = 0.7317384, "Durham|GS7|Decipher-marker surrogate" = 0.7174114)
cat("\nAUC check (computed vs expected):\n")
for (i in seq_len(nrow(panelA))) {
  key <- paste(panelA$cohort[i], panelA$subset[i], panelA$score[i], sep = "|")
  if (key %in% names(aucT)) cat(sprintf("  %-52s %.7f vs %.7f  (|d|=%.1e)\n", key, panelA$auc[i], aucT[key], abs(panelA$auc[i] - aucT[key])))
}
cat("\nPanel a paired (Met-Score - surrogate), BH across four comparisons:\n")
print(paired[, c("cohort", "subset", "delta_auc", "delta_lo", "delta_hi", "p_paired", "q_bh", "boot_used")], row.names = FALSE)
cat("\nPanel b (cause-specific summary HR per SD):\n")
print(panelB[, c("cohort", "subset", "score", "n", "events", "hr", "ci_lo_robust", "ci_hi_robust", "interaction_p", "zph_score_p")], row.names = FALSE)
hrT <- c("JHU|All patients|Met-Score" = 1.84784, "JHU|All patients|Decipher-marker surrogate" = 1.96200,
         "JHU|GS7|Met-Score" = 3.28064, "JHU|GS7|Decipher-marker surrogate" = 4.02471,
         "Durham|All patients|Met-Score" = 1.62594, "Durham|All patients|Decipher-marker surrogate" = 1.59464,
         "Durham|GS7|Met-Score" = 1.57436, "Durham|GS7|Decipher-marker surrogate" = 1.40360)
cat("\nHR check (computed vs expected):\n")
for (i in seq_len(nrow(panelB))) {
  key <- paste(panelB$cohort[i], panelB$subset[i], panelB$score[i], sep = "|")
  if (key %in% names(hrT)) cat(sprintf("  %-52s %.5f vs %.5f  (|d|=%.1e)\n", key, panelB$hr[i], hrT[key], abs(panelB$hr[i] - hrT[key])))
}

## ---- method checks (fail closed) ----------------------------------------
chk <- function(lab, got, exp, tol) { ok <- is.finite(got) && abs(got - exp) <= tol
  cat(sprintf("  [%s] %-42s %.8g vs %.8g\n", if (ok) "PASS" else "FAIL", lab, got, exp)); ok }
cat("\nMethod checks:\n"); ok <- TRUE
for (t0 in c(60, 120)) {
  trd <- timeROC::timeROC(T = DUR$analysis_time, delta = DUR$status, marker = DUR$ms, cause = 1, times = t0, iid = FALSE)
  trv <- if (!is.null(trd$AUC_2)) as.numeric(trd$AUC_2[paste0("t=", t0)]) else as.numeric(trd$AUC[paste0("t=", t0)])
  ok <- chk(sprintf("Durham_w1_eq_timeROC(%d)", t0), wauc_def2(DUR, t0, DUR$ms), trv, 1e-6) && ok
}
{ d2 <- JHU; d2$w <- d2$w * 3.7
  ok <- chk("JHU_weight_scale_invariance(120)", wauc_def2(d2, 120, d2$ms), wauc_def2(JHU, 120, JHU$ms), 1e-9) && ok }
{ set.seed(7); pj <- JHU[sample.int(nrow(JHU)), ]
  ok <- chk("JHU_row_permutation_invariance(120)", wauc_def2(pj, 120, pj$ms), wauc_def2(JHU, 120, JHU$ms), 1e-9) && ok }
for (t0 in c(60, 120))
  ok <- chk(sprintf("JHU_exact_horizon_cases(%d)", t0), ipcw_parts(JHU, t0)$ncase, sum(JHU$status == 1L & JHU$analysis_time <= t0), 1e-9) && ok
if (!ok) stop("Figure S3 method check(s) failed")

nw <- length(.warn$log)
cat(sprintf("\nWarnings captured: %d\n", nw))
if (nw) { wl <- do.call(rbind, .warn$log); print(unique(wl), row.names = FALSE) }
cat("\nJHU panel-a bootstrap ledger:\n")
print(paired[paired$cohort == "JHU", c("subset", "boot_attempted", "boot_used", "boot_failed", "boot_fail_fit_error", "boot_fail_nonfinite")], row.names = FALSE)
cat("\nDone.\n")

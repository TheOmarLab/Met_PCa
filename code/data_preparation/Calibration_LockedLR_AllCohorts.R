## Risk-stratification summary for the locked ridge-logistic Met-Score, JHU + Durham VA.
##
## The locked Met-Score is a binary classifier (cv.glmnet, family = "binomial") with no
## time horizon, so its output is a class probability, not an absolute 10-year risk. This
## script makes no probability-calibration claim. From one aligned competing-event table
## per cohort it reports, for the metastasis event and the locked HIGH/LOW threshold:
##   (a) cumulative/dynamic IPCW time-dependent AUC of the locked probability at 5 and 10y,
##   (b) 10-year metastasis cumulative incidence by HIGH/LOW group (Aalen-Johansen, death
##       without metastasis as the competing event; naive 1-S KM printed for reference),
##   (c) the cause-specific HIGH-vs-LOW hazard ratio,
##   (d) locked-cut cumulative/dynamic Se/Sp/PPV/NPV and ROC curves at 5 and 10y.
##
## JHU is a Ross two-phase case-cohort (source 745, random subcohort 265): the retained 239
## are 28 subcohort cases, 146 subcohort controls, and 65 outside-subcohort cases. All JHU
## estimators are phase-two weighted (subcohort controls 745/265, all cases 1). The high-vs-
## low association is a Lin-Ying case-cohort Cox (survival::cch); absolute-risk and operating-
## characteristic uncertainty is a fixed-seed conditional case-cohort bootstrap resampling
## within the three design strata (not phase-one/design-complete). Durham is a complete
## external cohort with all-one weights, ordinary Aalen-Johansen/Cox, and patient bootstrap.
##
## Competing-event table: status 1 = metastasis at its metastasis time; status 2 = death
## without prior metastasis at the death time when the death time is at or before the
## metastasis follow-up; otherwise status 0 = censored. analysis_time is the death time for
## status 2 and the metastasis follow-up otherwise. Horizons are evaluated immediately above
## 60/120 months so events at exactly the horizon are cases (score >= threshold cut).
##
## Output: outs/RiskStratification_LockedMetScore_AllCohorts.csv
##         outs/RiskStratification_AJ_curves_AllCohorts.csv
##         outs/FigureS6_ROC_curves.csv
##         outs/FigureS6_operating_characteristics.csv
##         outs/RiskStratification_method_checks.csv

suppressPackageStartupMessages({
  library(survival); library(timeROC)
})

out_dir <- "./outs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("./code/utils/locked_metscore.R")
LOCKED_THRESHOLD <- load_locked_metscore()$threshold
HORIZONS <- c(60, 120)
ALPHA    <- 265 / 745                      # source random-subcohort sampling fraction
B_BOOT   <- as.integer(Sys.getenv("METPCA_CALIBRATION_BOOTSTRAPS", "2000"))
SEED_JHU <- 20260814L
SEED_DUR <- 20260815L

# warning capture: warnings are logged with context, never silently dropped.
.warn <- new.env(); .warn$log <- list()
cap <- function(ctx, expr) withCallingHandlers(expr, warning = function(w) {
  .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = ctx,
    message = conditionMessage(w), stringsAsFactors = FALSE)
  invokeRestart("muffleWarning")
})

evtf <- function(e) factor(e, levels = c(0, 1, 2), labels = c("cens", "met", "death"))

# ---- one aligned competing-event table per cohort ------------------------
build_event_table <- function(score, mettime, metevent, death, deathtime, threshold,
                              cohort, expect_n, expect_mets, expect_cdeath, cchdef = NULL) {
  score <- as.numeric(score); mettime <- as.numeric(mettime); metevent <- as.integer(metevent)
  death     <- if (is.null(death))     rep(0L, length(score))       else as.integer(death)
  deathtime <- if (is.null(deathtime)) rep(NA_real_, length(score)) else as.numeric(deathtime)
  keep <- is.finite(score) & is.finite(mettime) & mettime > 0 & !is.na(metevent)
  score <- score[keep]; mettime <- mettime[keep]; metevent <- metevent[keep]
  death <- death[keep]; deathtime <- deathtime[keep]
  status <- ifelse(metevent == 1L, 1L,
             ifelse(death == 1L & !is.na(deathtime) & deathtime <= mettime, 2L, 0L))
  analysis_time <- ifelse(status == 2L, deathtime, mettime)
  group <- factor(ifelse(score >= threshold, "High", "Low"), levels = c("Low", "High"))
  d <- data.frame(score = score, group = group, analysis_time = as.numeric(analysis_time),
                  status = as.integer(status), stringsAsFactors = FALSE)
  if (!is.null(cchdef)) {
    cchv <- as.character(cchdef)[keep]
    d$cch   <- cchv
    d$w     <- ifelse(cchv == "Sub-cohort controls", 1 / ALPHA, 1)   # 745/265 controls, 1 cases
    d$insub <- as.integer(cchv %in% c("Sub-cohort cases", "Sub-cohort controls"))
  } else {
    d$cch <- "cohort"; d$w <- 1; d$insub <- 1L
  }
  d$id <- seq_len(nrow(d))
  stopifnot(all(is.finite(d$analysis_time)), all(d$analysis_time > 0),
            all(d$status %in% c(0L, 1L, 2L)), !any(is.na(d$group)))
  n <- nrow(d); mets <- sum(d$status == 1L); cdeath <- sum(d$status == 2L)
  if (!(n == expect_n && mets == expect_mets && cdeath == expect_cdeath))
    stop(sprintf("%s tally %d/%d/%d != expected %d/%d/%d", cohort, n, mets, cdeath,
                 expect_n, expect_mets, expect_cdeath))
  d$evt <- evtf(d$status)
  d
}

# ---- fast weighted Aalen-Johansen cause-1 CIF at t0 (counting-process) ----
wcif1 <- function(time, evt, w, t0) {
  o <- order(time); time <- time[o]; evt <- evt[o]; w <- w[o]
  Wrev <- rev(cumsum(rev(w)))
  et <- unique(time[evt != 0 & time <= t0]); if (!length(et)) return(0)
  S <- 1; cif <- 0
  for (u in et) {
    idx <- which(time == u); R <- Wrev[idx[1]]; if (R <= 0) next
    dNany <- sum(w[idx][evt[idx] != 0]); dN1 <- sum(w[idx][evt[idx] == 1])
    cif <- cif + S * dN1 / R; S <- S * (1 - dNany / R)
  }
  cif
}
cif_grp <- function(d, g, t0) {
  s <- d$group == g; if (!any(s)) return(NA_real_)
  wcif1(d$analysis_time[s], d$status[s], d$w[s], t0)
}

# ---- weighted cumulative/dynamic def-2 IPCW discrimination ----------------
# Cases: status 1 with analysis_time <= t0 (exact-horizon events included). Controls:
# event-free past t0, plus competing deaths by t0 (definition 2). IPCW uses the phase-two-
# weighted reverse Kaplan-Meier: cases and competing-death controls at G(T-), event-free
# controls at G(t0). At w = 1 this reduces to timeROC / SeSpPPVNPV.
ipcw_parts <- function(d, t0) {
  cens <- as.integer(d$status == 0L)
  km <- cap("censoringKM", survfit(Surv(d$analysis_time, cens) ~ 1, weights = d$w))
  tv <- km$time; sv <- km$surv
  Gm <- function(x) vapply(x, function(z){k <- which(tv <  z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  Ga <- function(x) vapply(x, function(z){k <- which(tv <= z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  case <- d$status == 1L & d$analysis_time <= t0
  c1   <- d$analysis_time >  t0
  c2   <- d$status == 2L & d$analysis_time <= t0
  wc <- ifelse(case, d$w / pmax(Gm(d$analysis_time), 1e-12), 0)
  wk <- ifelse(c1,  d$w / pmax(Ga(t0), 1e-12),
         ifelse(c2, d$w / pmax(Gm(d$analysis_time), 1e-12), 0))
  list(wc = wc, wk = wk, case = case, ncase = sum(case))
}
wauc_def2 <- function(d, t0, marker = d$score) {
  p <- ipcw_parts(d, t0); wc <- p$wc; wk <- p$wk
  ci <- which(wc > 0); ki <- which(wk > 0)
  if (!length(ci) || !length(ki)) return(NA_real_)
  num <- 0
  for (i in ci) num <- num + wc[i] * sum(wk[ki] * ((marker[i] > marker[ki]) + 0.5 * (marker[i] == marker[ki])))
  num / (sum(wc[ci]) * sum(wk[ki]))
}
# locked-cut cumulative/dynamic Se/Sp/PPV/NPV at threshold thr (score >= thr = positive)
oc_at <- function(d, t0, thr, pi_t0) {
  p <- ipcw_parts(d, t0); wc <- p$wc; wk <- p$wk
  pos <- d$score >= thr
  se <- sum(wc[pos]) / sum(wc)
  sp <- sum(wk[!pos]) / sum(wk)
  ppv <- (pi_t0 * se) / (pi_t0 * se + (1 - pi_t0) * (1 - sp))
  npv <- ((1 - pi_t0) * sp) / (pi_t0 * (1 - se) + (1 - pi_t0) * sp)
  c(Se = se, Sp = sp, PPV = ppv, NPV = npv)
}
# ROC curve (fpr = 1-Sp, tpr = Se) over a marker grid
roc_curve <- function(d, t0, grid) {
  p <- ipcw_parts(d, t0); wc <- p$wc; wk <- p$wk
  swc <- sum(wc); swk <- sum(wk)
  do.call(rbind, lapply(grid, function(c0) {
    pos <- d$score >= c0
    data.frame(threshold = c0, tpr = sum(wc[pos]) / swc, fpr = 1 - sum(wk[!pos]) / swk)
  }))
}

# ---- point-estimate weighted AJ CIF curves by group -----------------------
aj_curve_point <- function(d, grid = seq(0, 120, by = 1)) {
  sf <- cap("ajcurve", survfit(Surv(analysis_time, evt) ~ group, data = d, weights = w))
  s  <- summary(sf, times = grid, extend = TRUE)
  mi <- which(sf$states == "met")
  grp <- sub("^group=", "", as.character(s$strata))
  data.frame(RiskGroup = grp, time = s$time, cif = s$pstate[, mi],
             cif_lo = s$lower[, mi], cif_hi = s$upper[, mi], stringsAsFactors = FALSE)
}

run_cohort <- function(cohort, d, weighted, B, seed) {
  grid <- seq(0, 120, by = 1)
  thr  <- LOCKED_THRESHOLD
  pit  <- setNames(vapply(HORIZONS, function(t0) wcif1(d$analysis_time, d$status, d$w, t0), 0),
                   as.character(HORIZONS))
  # point estimates
  cif_grp_pt <- list()
  for (t0 in HORIZONS) cif_grp_pt[[as.character(t0)]] <-
    c(Low = cif_grp(d, "Low", t0), High = cif_grp(d, "High", t0))
  gap10 <- cif_grp_pt[["120"]]["High"] - cif_grp_pt[["120"]]["Low"]
  auc  <- setNames(vapply(HORIZONS, function(t0) wauc_def2(d, t0), 0), as.character(HORIZONS))
  oc   <- lapply(HORIZONS, function(t0) oc_at(d, t0, thr, pit[[as.character(t0)]]))
  names(oc) <- as.character(HORIZONS)
  mgrid <- sort(unique(c(-Inf, d$score, thr, Inf)))
  roc  <- lapply(HORIZONS, function(t0) { r <- roc_curve(d, t0, mgrid); r$Cohort <- cohort; r$horizon <- t0; r })
  names(roc) <- as.character(HORIZONS)
  curve <- aj_curve_point(d); curve$Cohort <- cohort
  # exact-horizon case-count ledger identity
  for (t0 in HORIZONS)
    stopifnot(ipcw_parts(d, t0)$ncase == sum(d$status == 1L & d$analysis_time <= t0))

  # bootstrap: JHU = conditional resample within design strata; Durham = patient resample
  strat <- if (weighted) split(seq_len(nrow(d)), d$cch) else list(all = seq_len(nrow(d)))
  set.seed(seed)
  bkeys <- c("auc60","auc120","gap120","Se60","Sp60","PPV60","NPV60","Se120","Sp120","PPV120","NPV120")
  acc <- setNames(vector("list", length(bkeys)), bkeys); for (k in bkeys) acc[[k]] <- numeric(0)
  band <- list("60" = list(Low=NULL,High=NULL), "120" = list(Low=NULL,High=NULL))  # curve band (grid CIFs)
  reasons <- c(absent_group = 0L, fit_error = 0L, nonfinite = 0L); used <- 0L
  for (b in seq_len(B)) {
    ix <- unlist(lapply(strat, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
    db <- d[ix, , drop = FALSE]
    gt <- table(db$group)
    if (!all(c("Low","High") %in% names(gt)) || any(gt[c("Low","High")] < 1)) {
      reasons[["absent_group"]] <- reasons[["absent_group"]] + 1L; next }
    r <- tryCatch({
      a60 <- wauc_def2(db, 60); a120 <- wauc_def2(db, 120)
      lo120 <- cif_grp(db, "Low", 120); hi120 <- cif_grp(db, "High", 120)
      p60 <- wcif1(db$analysis_time, db$status, db$w, 60); p120 <- wcif1(db$analysis_time, db$status, db$w, 120)
      o60 <- oc_at(db, 60, thr, p60); o120 <- oc_at(db, 120, thr, p120)
      cl <- vapply(grid, function(t) cif_grp(db, "Low", t), 0); ch <- vapply(grid, function(t) cif_grp(db, "High", t), 0)
      list(a60=a60,a120=a120,gap=hi120-lo120,o60=o60,o120=o120,cl=cl,ch=ch)
    }, error = function(e) NULL)
    if (is.null(r)) { reasons[["fit_error"]] <- reasons[["fit_error"]] + 1L; next }
    vals <- c(r$a60,r$a120,r$gap,r$o60,r$o120)
    if (any(!is.finite(vals))) { reasons[["nonfinite"]] <- reasons[["nonfinite"]] + 1L; next }
    used <- used + 1L
    acc$auc60 <- c(acc$auc60,r$a60); acc$auc120 <- c(acc$auc120,r$a120); acc$gap120 <- c(acc$gap120,r$gap)
    acc$Se60<-c(acc$Se60,r$o60["Se"]); acc$Sp60<-c(acc$Sp60,r$o60["Sp"]); acc$PPV60<-c(acc$PPV60,r$o60["PPV"]); acc$NPV60<-c(acc$NPV60,r$o60["NPV"])
    acc$Se120<-c(acc$Se120,r$o120["Se"]); acc$Sp120<-c(acc$Sp120,r$o120["Sp"]); acc$PPV120<-c(acc$PPV120,r$o120["PPV"]); acc$NPV120<-c(acc$NPV120,r$o120["NPV"])
    band[["120"]]$Low  <- rbind(band[["120"]]$Low,  r$cl); band[["120"]]$High <- rbind(band[["120"]]$High, r$ch)
  }
  q <- function(x) if (length(x) >= 50) unname(quantile(x, c(.025,.975))) else c(NA_real_, NA_real_)
  ci <- lapply(acc, q)
  list(cohort=cohort, weighted=weighted, pit=pit, cif_grp=cif_grp_pt, gap10=as.numeric(gap10),
       auc=auc, oc=oc, roc=roc, curve=curve, ci=ci, band=band,
       boot=list(B=B, used=used, failed=B-used, reasons=reasons),
       n=nrow(d), mets=sum(d$status==1L), cdeath=sum(d$status==2L))
}

# ---- high-vs-low association --------------------------------------------
hr_jhu_cch <- function(d) {
  fit <- cap("cch_high_low", survival::cch(Surv(analysis_time, status == 1L) ~ group, data = d,
             subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE))
  cf <- summary(fit)$coefficients            # single groupHigh row, unnamed
  b <- cf[1, "Value"]; se <- cf[1, "SE"]
  list(hr = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se), p = cf[1, "p"],
       basis = "case-cohort Lin-Ying robust (cohort.size=745)")
}
hr_cox <- function(d) {
  fit <- survival::coxph(Surv(analysis_time, as.integer(status == 1L)) ~ group, data = d)
  s <- summary(fit)
  list(hr = s$conf.int["groupHigh","exp(coef)"], lo = s$conf.int["groupHigh","lower .95"],
       hi = s$conf.int["groupHigh","upper .95"], p = s$coefficients["groupHigh","Pr(>|z|)"],
       basis = "complete-cohort cause-specific Cox")
}

# ==========================================================================
# JHU: weighted case-cohort
load("./outs/coxdata.rda")
djhu <- build_event_table(CoxData_jhu[["Met-Score prob"]], CoxData_jhu$met_time, CoxData_jhu$met,
                          CoxData_jhu$os, CoxData_jhu$os_time, LOCKED_THRESHOLD, "JHU",
                          239L, 93L, 6L, cchdef = CoxData_jhu[["post_rp_patients_cchdef"]])
stopifnot(sum(djhu$cch == "Sub-cohort cases") == 28L, sum(djhu$cch == "Sub-cohort controls") == 146L,
          sum(djhu$cch == "cases") == 65L)
RJ <- run_cohort("JHU", djhu, weighted = TRUE, B = B_BOOT, seed = SEED_JHU)
RJ$hr <- hr_jhu_cch(djhu)

# Durham VA: complete external cohort
load("./output/Durham/durham_metscore_batchcorrected.rda")
ddur <- build_event_table(clin_valid$MetScore_prob, clin_valid$surgmets, clin_valid$mets,
                          clin_valid$dead, clin_valid$limbo, LOCKED_THRESHOLD, "Durham VA",
                          555L, 40L, 167L)
RD <- run_cohort("Durham VA", ddur, weighted = FALSE, B = B_BOOT, seed = SEED_DUR)
RD$hr <- hr_cox(ddur)

# ==========================================================================
# method-equivalence / invariance checks (checks only, never fitting targets)
chk <- list()
addchk <- function(name, got, ref, tol) chk[[length(chk)+1]] <<-
  data.frame(check = name, value = got, reference = ref, pass = abs(got - ref) < tol, row.names = NULL)
# Durham all-one-weight AUC == timeROC def-2 (no exact-horizon events -> conventions coincide)
for (t0 in HORIZONS) {
  tr <- timeROC::timeROC(T = ddur$analysis_time, delta = ddur$status, marker = ddur$score,
                         cause = 1, times = t0, iid = FALSE)
  trauc <- if (!is.null(tr$AUC_2)) as.numeric(tr$AUC_2[paste0("t=",t0)]) else as.numeric(tr$AUC[paste0("t=",t0)])
  addchk(sprintf("Durham_w1_AUC_eq_timeROC(%d)", t0), wauc_def2(ddur, t0), trauc, 1e-6)
}
# Durham all-one-weight Se/Sp == timeROC::SeSpPPVNPV at the frozen cut
for (t0 in HORIZONS) {
  ss <- cap("SeSp", timeROC::SeSpPPVNPV(cutpoint = LOCKED_THRESHOLD, T = ddur$analysis_time,
            delta = ddur$status, marker = ddur$score, cause = 1, times = t0, weighting = "marginal"))
  se_ref <- as.numeric(ss$TP[paste0("t=",t0)]); sp_ref <- 1 - as.numeric(ss$FP_2[paste0("t=",t0)])
  o <- oc_at(ddur, t0, LOCKED_THRESHOLD, wcif1(ddur$analysis_time, ddur$status, ddur$w, t0))
  addchk(sprintf("Durham_w1_Se_eq_SeSpPPVNPV(%d)", t0), o["Se"], se_ref, 1e-6)
  addchk(sprintf("Durham_w1_Sp_eq_SeSpPPVNPV(%d)", t0), o["Sp"], sp_ref, 1e-6)
}
# no-censoring reduction: on rows with an event (no censoring), IPCW AUC == direct weighted AUC
{ de <- djhu[djhu$status != 0L, ]; t0 <- 120
  p <- ipcw_parts(de, t0); direct <- { wc <- ifelse(de$status==1L & de$analysis_time<=t0, de$w, 0)
    wk <- ifelse(de$analysis_time>t0 | (de$status==2L & de$analysis_time<=t0), de$w, 0)
    ci<-which(wc>0); ki<-which(wk>0); num<-0
    for (i in ci) num<-num+wc[i]*sum(wk[ki]*((de$score[i]>de$score[ki])+0.5*(de$score[i]==de$score[ki]))); num/(sum(wc[ci])*sum(wk[ki])) }
  addchk("no_censoring_reduces_to_direct(JHU,120)", wauc_def2(de, t0), direct, 1e-8) }
# phase-weight scale invariance: scaling all weights by a constant leaves weighted AUC unchanged
{ d2 <- djhu; d2$w <- d2$w * 3.7
  addchk("weight_scale_invariance(JHU,120)", wauc_def2(d2, 120), wauc_def2(djhu, 120), 1e-9) }
# row-permutation invariance
{ set.seed(1); pj <- djhu[sample.int(nrow(djhu)), ]
  addchk("row_permutation_invariance(JHU,120)", wauc_def2(pj, 120), wauc_def2(djhu, 120), 1e-9) }
# score orientation: reversing the marker gives AUC = 1 - AUC
addchk("score_orientation(JHU,120)", { d3 <- djhu; d3$score <- -d3$score; wauc_def2(d3, 120) },
       1 - wauc_def2(djhu, 120), 1e-9)
# locked-class identity: group partition equals score >= threshold
addchk("locked_class_identity(JHU)",
       mean(djhu$group == ifelse(djhu$score >= LOCKED_THRESHOLD, "High", "Low")), 1, 1e-12)
# exact-horizon case-count identity
addchk("exact_horizon_cases(JHU,60)",  ipcw_parts(djhu,60)$ncase,  sum(djhu$status==1L & djhu$analysis_time<=60), 1e-9)
addchk("exact_horizon_cases(JHU,120)", ipcw_parts(djhu,120)$ncase, sum(djhu$status==1L & djhu$analysis_time<=120), 1e-9)
# threshold-tie guard: no score exactly equals the cut (>= and > partitions coincide)
n_tie <- sum(djhu$score == LOCKED_THRESHOLD) + sum(ddur$score == LOCKED_THRESHOLD)
if (n_tie > 0) stop(sprintf("threshold ties present (%d): strict-cut equivalence not guaranteed", n_tie))
chk_df <- do.call(rbind, chk)
if (!all(chk_df$pass)) stop("method-equivalence check failed:\n",
                            paste(capture.output(print(chk_df[!chk_df$pass, ])), collapse = "\n"))

# ==========================================================================
# assemble outputs
summ_row <- function(R) {
  cg <- R$cif_grp[["120"]]
  do.call(rbind, lapply(c("Low","High"), function(g) data.frame(
    Cohort = R$cohort, RiskGroup = g,
    N = sum((if (R$cohort=="JHU") djhu else ddur)$group == g),
    Events = sum((if (R$cohort=="JHU") djhu else ddur)$group == g &
                 (if (R$cohort=="JHU") djhu else ddur)$status == 1L),
    AJ_CIF_10y = round(cg[[g]], 4),
    n_total = R$n, mets_total = R$mets, cdeath_total = R$cdeath,
    HR_high_vs_low = round(R$hr$hr, 4), HR_lo = round(R$hr$lo, 4), HR_hi = round(R$hr$hi, 4),
    Cox_p = signif(R$hr$p, 4),
    AJ_diff_10y = round(R$gap10, 4), AJ_diff_lo = round(R$ci$gap120[1], 4), AJ_diff_hi = round(R$ci$gap120[2], 4),
    AUC_5y = round(R$auc[["60"]], 4), AUC_5y_lo = round(R$ci$auc60[1], 4), AUC_5y_hi = round(R$ci$auc60[2], 4),
    AUC_10y = round(R$auc[["120"]], 4), AUC_10y_lo = round(R$ci$auc120[1], 4), AUC_10y_hi = round(R$ci$auc120[2], 4),
    boot_used = R$boot$used, boot_failed = R$boot$failed,
    basis = if (R$weighted) "JHU case-cohort: phase-two-weighted AJ/def-2 AUC, Lin-Ying cch HR, conditional case-cohort bootstrap"
            else "Durham complete cohort: ordinary AJ/Cox/def-2 AUC, patient bootstrap",
    stringsAsFactors = FALSE)))
}
tbl <- rbind(summ_row(RJ), summ_row(RD))
write.csv(tbl, file.path(out_dir, "RiskStratification_LockedMetScore_AllCohorts.csv"), row.names = FALSE)

# AJ CIF curves (JHU band = conditional bootstrap percentiles; Durham band = ordinary survfit)
curve_out <- function(R) {
  cu <- R$curve
  if (R$weighted) {
    band <- do.call(rbind, lapply(c("Low","High"), function(g) {
      M <- R$band[["120"]][[g]]
      lo <- apply(M, 2, function(x){x<-x[is.finite(x)]; if(length(x)>=50) quantile(x,.025) else NA_real_})
      hi <- apply(M, 2, function(x){x<-x[is.finite(x)]; if(length(x)>=50) quantile(x,.975) else NA_real_})
      data.frame(RiskGroup = g, time = seq(0,120,by=1), cif_lo = lo, cif_hi = hi, stringsAsFactors = FALSE)
    }))
    cu$cif_lo <- band$cif_lo[match(paste(cu$RiskGroup,cu$time), paste(band$RiskGroup,band$time))]
    cu$cif_hi <- band$cif_hi[match(paste(cu$RiskGroup,cu$time), paste(band$RiskGroup,band$time))]
  }
  d <- if (R$cohort=="JHU") djhu else ddur
  cu$n_at_risk <- mapply(function(g,t) sum(d$analysis_time[d$group==g] >= t), cu$RiskGroup, cu$time)
  cu[, c("Cohort","RiskGroup","time","cif","cif_lo","cif_hi","n_at_risk")]
}
curves <- rbind(curve_out(RJ), curve_out(RD))
curves[c("cif","cif_lo","cif_hi")] <- lapply(curves[c("cif","cif_lo","cif_hi")], function(x) round(x, 6))
write.csv(curves, file.path(out_dir, "RiskStratification_AJ_curves_AllCohorts.csv"), row.names = FALSE)

# Figure S6: ROC curves + operating characteristics
roc_all <- do.call(rbind, lapply(list(RJ, RD), function(R) do.call(rbind, lapply(HORIZONS, function(t0) {
  r <- R$roc[[as.character(t0)]]; r$is_frozen_cut <- r$threshold == LOCKED_THRESHOLD
  r[, c("Cohort","horizon","threshold","fpr","tpr","is_frozen_cut")] }))))
roc_all[c("threshold","fpr","tpr")] <- lapply(roc_all[c("threshold","fpr","tpr")], function(x) round(x, 6))
write.csv(roc_all, file.path(out_dir, "FigureS6_ROC_curves.csv"), row.names = FALSE)

oc_all <- do.call(rbind, lapply(list(RJ, RD), function(R) do.call(rbind, lapply(HORIZONS, function(t0) {
  h <- as.character(t0); o <- R$oc[[h]]
  do.call(rbind, lapply(c("Se","Sp","PPV","NPV"), function(m) data.frame(
    Cohort = R$cohort, horizon = t0, metric = m, estimate = round(o[[m]], 4),
    lo = round(R$ci[[paste0(m,h)]][1], 4), hi = round(R$ci[[paste0(m,h)]][2], 4),
    threshold = LOCKED_THRESHOLD,
    uncertainty = if (R$weighted) "conditional case-cohort bootstrap" else "patient bootstrap",
    stringsAsFactors = FALSE))) }))))
write.csv(oc_all, file.path(out_dir, "FigureS6_operating_characteristics.csv"), row.names = FALSE)

# method checks + bootstrap/warning ledger
warn_tab <- if (length(.warn$log)) do.call(rbind, .warn$log) else
  data.frame(context="none", message="no warnings captured", stringsAsFactors = FALSE)
wsum <- if (nrow(warn_tab)) as.data.frame(table(context = warn_tab$context)) else data.frame(context="none", Freq=0)
boot_led <- do.call(rbind, lapply(list(RJ, RD), function(R) data.frame(
  Cohort = R$cohort, B = R$boot$B, used = R$boot$used, failed = R$boot$failed,
  absent_group = R$boot$reasons[["absent_group"]], fit_error = R$boot$reasons[["fit_error"]],
  nonfinite = R$boot$reasons[["nonfinite"]], stringsAsFactors = FALSE)))
write.csv(chk_df, file.path(out_dir, "RiskStratification_method_checks.csv"), row.names = FALSE)

# ---- console report ------------------------------------------------------
cat(sprintf("\n=== Locked Met-Score risk stratification (JHU case-cohort corrected) ===\n"))
cat(sprintf("threshold = %.17f ; B = %d ; ties at cut = %d\n", LOCKED_THRESHOLD, B_BOOT, n_tie))
cat("JHU design strata: 28 subcohort cases / 146 subcohort controls / 65 outside cases ; control weight 745/265\n\n")
for (R in list(RJ, RD)) {
  cat(sprintf("%s: n=%d mets=%d competing-deaths=%d [%s]\n", R$cohort, R$n, R$mets, R$cdeath,
              if (R$weighted) "weighted case-cohort" else "complete cohort"))
  for (t0 in HORIZONS) cat(sprintf("  %d-mo AUC = %.4f (%.4f-%.4f)  cases<=h = %d\n", t0,
    R$auc[[as.character(t0)]], R$ci[[paste0("auc",t0)]][1], R$ci[[paste0("auc",t0)]][2],
    ipcw_parts(if(R$cohort=="JHU") djhu else ddur, t0)$ncase))
  cat(sprintf("  10y CIF High=%.4f Low=%.4f  gap=%.4f (%.4f-%.4f)  HR=%.4f (%.4f-%.4f) p=%s [%s]\n",
    R$cif_grp[["120"]][["High"]], R$cif_grp[["120"]][["Low"]], R$gap10, R$ci$gap120[1], R$ci$gap120[2],
    R$hr$hr, R$hr$lo, R$hr$hi, format(signif(R$hr$p,3)), R$hr$basis))
  cat(sprintf("  bootstrap used=%d/%d failed=%d (absent_group=%d fit_error=%d nonfinite=%d)\n\n",
    R$boot$used, R$boot$B, R$boot$failed, R$boot$reasons[["absent_group"]], R$boot$reasons[["fit_error"]], R$boot$reasons[["nonfinite"]]))
}
cat(sprintf("method-equivalence checks: %d/%d pass\n", sum(chk_df$pass), nrow(chk_df)))
cat("Written: RiskStratification_{LockedMetScore_AllCohorts,AJ_curves_AllCohorts,method_checks}.csv, FigureS6_{ROC_curves,operating_characteristics}.csv\n")
cat("=== DONE ===\n")

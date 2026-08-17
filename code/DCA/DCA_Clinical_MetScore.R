###############################################################################
# DCA_Clinical_MetScore.R
#
# JHU case-cohort development, frozen Durham external validation, and JHU-only
# CAPRA-S utility for the metastasis competing-risk endpoint at 60 and 120 months.
#
# JHU is a two-phase case-cohort (Ross et al., Eur Urol 2016): a random 35%
# subcohort (265/745) plus all out-of-subcohort metastatic cases. Absolute risk,
# calibration, and decision curves use inverse-probability (phase-two) weighting:
# subcohort noncases weight 745/265, all cases weight 1. The metastasis model is
# an inverse-probability-weighted Fine-Gray subdistribution model (survival::
# finegray + weighted coxph), which reproduces the published Ross Decipher
# absolute-risk anchor. A design-valid two-phase variance for the subdistribution
# model is not constructible from the sampled data alone (the phase-one frame is
# unavailable), so relative-effect inference is reported from a Lin-Ying
# case-cohort model and the weighted absolute-risk/DCA results are point estimates.
# Durham is a complete external cohort scored by the frozen JHU model with no
# weighting.
#
# Endpoint: metastasis competing-risk process. evt 1 = metastasis, 2 = death
# without metastasis inside the metastasis window, 0 = censored.
###############################################################################

suppressPackageStartupMessages({
  library(survival); library(cmprsk); library(prodlim); library(timeROC)
  library(riskRegression)
})
OUTD <- "./outs/DCA"; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

T5  <- 60L; T10 <- 120L; HORIZ <- c(T5, T10)
THR   <- seq(0.01, 0.20, by = 0.01)
IMPL  <- c(0.05, 0.10, 0.15, 0.20)
ALPHA <- 265 / 745                       # source random-subcohort sampling fraction
B     <- as.integer(Sys.getenv("METPCA_DCA_BOOTSTRAPS", "2000"))
set.seed(1L)

# warning capture (never suppressed); collected into a disposition table
.warn <- new.env(); .warn$log <- list()
capture_warn <- function(tag, expr) {
  withCallingHandlers(expr, warning = function(w) {
    .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = tag,
      message = conditionMessage(w), stringsAsFactors = FALSE)
    invokeRestart("muffleWarning")
  })
}

###############################################################################
# covariates
###############################################################################
gg_ord <- function(total_gs, primary) {
  g <- rep(NA_real_, length(total_gs))
  g[total_gs <= 6]                <- 1
  g[total_gs == 7 & primary == 3] <- 2
  g[total_gs == 7 & primary == 4] <- 3
  g[total_gs == 8]                <- 4
  g[total_gs >= 9]                <- 5
  g
}
pt_high <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[grepl("^T2", x)]      <- 0L
  out[grepl("^T3|^T4", x)]  <- 1L
  out
}

load("./outs/coxdata.rda")                                    # CoxData_jhu
load("./output/Durham/durham_metscore_batchcorrected.rda")    # clin_valid

build_jhu <- function(j) {
  met   <- as.integer(j$met); dead <- as.integer(j$os)
  mtime <- as.numeric(j$met_time); dtime <- as.numeric(j$os_time)
  evt   <- ifelse(met == 1, 1L,
             ifelse(dead == 1 & met == 0 & !is.na(dtime) & dtime <= mtime, 2L, 0L))
  etime <- ifelse(evt == 2L, dtime, mtime)
  cch   <- as.character(j[["post_rp_patients_cchdef"]])
  data.frame(
    time    = etime, evt = evt,
    GG_ord  = gg_ord(as.numeric(as.character(j[["Pathological GS"]])),
                     as.numeric(as.character(j$pathgs_p))),
    log2PSA = log2(as.numeric(j$preop_psa) + 1),
    pT_hi   = pt_high(j$pstage),
    score   = as.numeric(j[["Met-Score prob"]]),
    rf22    = as.numeric(j$rf22_scan),
    cch     = cch,
    w       = ifelse(cch == "Sub-cohort controls", 1 / ALPHA, 1),
    insub   = as.integer(cch %in% c("Sub-cohort cases", "Sub-cohort controls")),
    # CAPRA-S components (face value)
    psa     = as.numeric(j$preop_psa),
    gp      = as.integer(j$pathgs_p), gs = as.integer(j$pathgs_s),
    smpos   = as.integer(trimws(tolower(as.character(j[["re.PositiveMargin"]]))) == "yes"),
    ece     = as.integer(j$ece), svi = as.integer(j$svi), lni = as.integer(j$lni),
    stringsAsFactors = FALSE)
}
build_durham <- function(d) {
  met   <- as.integer(d$mets); dead <- as.integer(d$dead)
  mtime <- as.numeric(d$surgmets); dtime <- as.numeric(d$limbo)
  evt   <- ifelse(met == 1, 1L,
             ifelse(dead == 1 & met == 0 & !is.na(dtime) & dtime <= mtime, 2L, 0L))
  etime <- ifelse(evt == 2L, dtime, mtime)
  data.frame(
    time = etime, evt = evt,
    GG_ord = gg_ord(as.numeric(as.character(d$PathGleason)), as.numeric(as.character(d$pogl1))),
    log2PSA = log2(as.numeric(d$psapresurg) + 1),
    pT_hi = pt_high(d$stg),
    score = as.numeric(d$MetScore_prob),
    stringsAsFactors = FALSE)
}

jhu <- build_jhu(CoxData_jhu); dur <- build_durham(clin_valid)
jhu$id <- seq_len(nrow(jhu))

# JHU cohort anchors (case-cohort structure); the cohort is all 239 with 93 events
stopifnot(nrow(jhu) == 239L, sum(jhu$evt == 1L) == 93L)
stopifnot(sum(jhu$cch == "Sub-cohort cases") == 28L,
          sum(jhu$cch == "Sub-cohort controls") == 146L,
          sum(jhu$cch == "cases") == 65L)
stopifnot(sum(dur$evt == 1L) == 40L, nrow(dur) == 555L)

# design-weighted (phase-two) standardization of the continuous Met-Score, frozen
mean_JHU <- sum(jhu$w * jhu$score) / sum(jhu$w)
sd_JHU   <- sqrt(sum(jhu$w * (jhu$score - mean_JHU)^2) / sum(jhu$w))
jhu$ms_dev <- (jhu$score - mean_JHU) / sd_JHU
dur$ms_dev <- (dur$score - mean_JHU) / sd_JHU

# model-specific complete-case set (same covariates for both nested models); the
# JHU cohort remains 239/93, this is a per-model exclusion
cc_vars <- c("time", "evt", "GG_ord", "log2PSA", "pT_hi", "ms_dev")
jhu$cc  <- stats::complete.cases(jhu[, cc_vars]) & jhu$time > 0
dur$cc  <- stats::complete.cases(dur[, cc_vars]) & dur$time > 0
jhu_cc  <- jhu[jhu$cc, ]; dur_cc <- dur[dur$cc, ]
excl <- data.frame(
  variable = c("GG_ord", "pT_hi", "log2PSA", "ms_dev", "time_nonpositive"),
  n_missing_jhu = c(sum(is.na(jhu$GG_ord)), sum(is.na(jhu$pT_hi)),
                    sum(is.na(jhu$log2PSA)), sum(is.na(jhu$ms_dev)), sum(jhu$time <= 0)))

clin_terms <- c("GG_ord", "log2PSA", "pT_hi")
full_terms <- c("GG_ord", "log2PSA", "pT_hi", "ms_dev")

###############################################################################
# weighted Fine-Gray (IPW subdistribution) + absolute risk
###############################################################################
evtf <- function(e) factor(e, levels = c(0, 1, 2), labels = c("cens", "met", "death"))

fit_wfg <- function(df, terms) {
  d <- data.frame(time = df$time, evtf = evtf(df$evt), df[, terms, drop = FALSE],
                  w = df$w, id = df$id)
  fg <- capture_warn("finegray", survival::finegray(
          Surv(time, evtf) ~ ., data = d, etype = "met", weights = w, id = id))
  fit <- capture_warn("coxph_wfg", survival::coxph(
          stats::as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                                  paste(terms, collapse = "+"))),
          data = fg, weights = fgwt, cluster = id, robust = TRUE, ties = "breslow"))
  fit
}
# predicted subdistribution CIF at t0 for newdf, from a frozen weighted-FG fit
cif_wfg <- function(fit, newdf, terms, t0) {
  sf <- survival::survfit(fit, newdata = newdf[, terms, drop = FALSE])
  ti <- if (any(sf$time <= t0)) max(which(sf$time <= t0)) else 1L
  as.numeric(1 - sf$surv[ti, ])
}

M_clin <- fit_wfg(jhu_cc, clin_terms)
M_full <- fit_wfg(jhu_cc, full_terms)

for (t0 in HORIZ) {
  jhu_cc[[paste0("cif_clin_", t0)]] <- cif_wfg(M_clin, jhu_cc, clin_terms, t0)
  jhu_cc[[paste0("cif_full_", t0)]] <- cif_wfg(M_full, jhu_cc, full_terms, t0)
  dur_cc[[paste0("cif_clin_", t0)]] <- cif_wfg(M_clin, dur_cc, clin_terms, t0)
  dur_cc[[paste0("cif_full_", t0)]] <- cif_wfg(M_full, dur_cc, full_terms, t0)
}

###############################################################################
# weighted competing-risk estimators for net benefit
###############################################################################
# weighted Aalen-Johansen metastasis CIF at t0 within a subset (survfit reference)
aj_cif_w <- function(time, evt, w, keep, t0) {
  if (sum(keep) == 0) return(NA_real_)
  sf <- tryCatch(survfit(Surv(time[keep], evtf(evt[keep])) ~ 1, weights = w[keep]),
                 error = function(e) NULL)
  if (is.null(sf)) return(NA_real_)
  mi <- which(sf$states == "met"); if (!length(mi)) return(0)
  as.numeric(summary(sf, times = t0, extend = TRUE)$pstate[, mi])
}
# fast weighted Aalen-Johansen cause-1 CIF at t0 (counting-process form); equals
# aj_cif_w but avoids survfit so it is affordable inside the resampling loops
wcif1_fast <- function(time, evt, w, t0) {
  o <- order(time); time <- time[o]; evt <- evt[o]; w <- w[o]
  Wrev <- rev(cumsum(rev(w)))                 # risk-set weight = sum w[time >= u]
  et <- unique(time[evt != 0 & time <= t0])
  if (!length(et)) return(0)
  S <- 1; cif <- 0
  for (u in et) {
    idx <- which(time == u); R <- Wrev[idx[1]]
    if (R <= 0) next
    dNany <- sum(w[idx][evt[idx] != 0]); dN1 <- sum(w[idx][evt[idx] == 1])
    cif <- cif + S * dN1 / R; S <- S * (1 - dNany / R)
  }
  cif
}
cif_pos_fast <- function(time, evt, w, keep, t0) {
  if (!any(keep)) return(NA_real_)
  wcif1_fast(time[keep], evt[keep], w[keep], t0)
}
# competing-risk net benefit per 100, phase-two weighted, over THR
nb_curve_w <- function(pred, time, evt, w, thr = THR, t0) {
  Wtot <- sum(w)
  vapply(thr, function(pt) {
    pos <- pred >= pt; if (!any(pos)) return(0)
    pi_pos  <- sum(w[pos]) / Wtot
    cif_pos <- cif_pos_fast(time, evt, w, pos, t0)
    if (is.na(cif_pos)) return(NA_real_)
    100 * pi_pos * (cif_pos - (1 - cif_pos) * (pt / (1 - pt)))
  }, numeric(1))
}
nb_all_w <- function(time, evt, w, thr = THR, t0) {
  cif <- wcif1_fast(time, evt, w, t0)
  100 * (cif - (1 - cif) * (thr / (1 - thr)))
}

###############################################################################
# weighted performance metrics (reduce to Score / timeROC at w = 1)
###############################################################################
# explicit identifier-safe competing-risk IPCW Brier + IPA. Censoring is the
# phase-two-weighted reverse Kaplan-Meier; the Graf weighting evaluates events at
# G(T-) and subjects still at risk past t0 at G(t0). All vectors stay in the
# caller's row order (no residual re-sorting), and the null uses the identical
# target, censoring model, and weights before IPA. At w = 1 this reduces to
# riskRegression::Score up to the reverse-KM tie convention.
wbrier <- function(pred, time, evt, w, t0) {
  cens <- as.integer(evt == 0)
  km <- capture_warn("brier_censKM", survfit(Surv(time, cens) ~ 1, weights = w))
  tv <- km$time; sv <- km$surv
  Gm <- function(x) vapply(x, function(z){k <- which(tv <  z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  Ga <- function(x) vapply(x, function(z){k <- which(tv <= z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  ipcw <- ifelse(time <= t0 & evt != 0, 1 / pmax(Gm(time), 1e-12),
           ifelse(time > t0, 1 / pmax(Ga(t0), 1e-12), 0))
  Y  <- as.integer(time <= t0 & evt == 1)
  m0 <- wcif1_fast(time, evt, w, t0)               # weighted marginal null CIF
  bm <- sum(w * ipcw * (Y - pred)^2) / sum(w)
  bn <- sum(w * ipcw * (Y - m0)^2)  / sum(w)
  c(brier = bm, brier_null = bn, ipa = 1 - bm / bn)
}
# weighted competing-risk def-2 cumulative/dynamic AUC (0.5 for ties)
wauc_def2 <- function(time, evt, marker, w, t0) {
  cens <- as.integer(evt == 0)
  km <- survfit(Surv(time, cens) ~ 1, weights = w); tv <- km$time; sv <- km$surv
  Gminus <- function(x) vapply(x, function(z){k <- which(tv <  z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  Gat    <- function(x) vapply(x, function(z){k <- which(tv <= z); if (!length(k)) 1 else sv[max(k)]}, numeric(1))
  Gsub <- Gminus(time); Gt0 <- Gat(t0)
  case <- time <  t0 & evt == 1
  c1   <- time >  t0; c2 <- time < t0 & evt == 2
  wc <- ifelse(case, w / pmax(Gsub, 1e-12), 0)
  wk <- ifelse(c1, w / pmax(Gt0, 1e-12), ifelse(c2, w / pmax(Gsub, 1e-12), 0))
  ci <- which(wc > 0); ki <- which(wk > 0)
  if (!length(ci) || !length(ki)) return(NA_real_)
  num <- 0
  for (i in ci) num <- num + wc[i] * sum(wk[ki] * ((marker[i] > marker[ki]) + 0.5 * (marker[i] == marker[ki])))
  num / (sum(wc[ci]) * sum(wk[ki]))
}

cat(sprintf("JHU analytic n=%d events=%d (cohort 239/93; %d complete-case excluded)\n",
            nrow(jhu_cc), sum(jhu_cc$evt == 1L), 239L - nrow(jhu_cc)))
cat(sprintf("Durham analytic n=%d events=%d\n", nrow(dur_cc), sum(dur_cc$evt == 1L)))

# locked Met-Score threshold (unchanged, read from config)
.meta  <- read.csv("./config/metscore_locked_v1_metadata.csv", stringsAsFactors = FALSE)
MS_THR <- as.numeric(.meta$value[.meta$field == "threshold"])

###############################################################################
# CAPRA-S score (reconstructed face value), computed here so the Decipher anchor
# can use the exact Ross Table 3 Model 2 specification (CAPRA-S + Decipher)
###############################################################################
psa_pts <- function(x) ifelse(x<=6,0L, ifelse(x<=10,1L, ifelse(x<=20,2L,3L)))
gl_pts  <- function(gp, gs){ tot<-gp+gs; ifelse(tot<=6,0L, ifelse(gp==3 & gs==4,1L, ifelse(gp==4 & gs==3,2L, ifelse(tot>=8,3L, NA_integer_)))) }
jhu$capras_psa <- psa_pts(jhu$psa)
jhu$capras_gl  <- gl_pts(jhu$gp, jhu$gs)
jhu$capras_sm  <- 2L*jhu$smpos; jhu$capras_ece <- jhu$ece; jhu$capras_svi <- 2L*jhu$svi; jhu$capras_lni <- jhu$lni
jhu$capras <- jhu$capras_psa + jhu$capras_gl + jhu$capras_sm + jhu$capras_ece + jhu$capras_svi + jhu$capras_lni

###############################################################################
# method verification against Ross Decipher anchors; checks only, never fitting
# targets. The two weighted nonparametric CIFs check the population weighting. The
# Decipher association check reproduces the exact Ross Table 3 Model 2
# specification (CAPRA-S + Decipher per 0.1) with a robust Lin-Ying cause-specific
# case-cohort Cox on all model-complete JHU observations.
###############################################################################
lowD <- jhu$rf22 < 0.45; highD <- jhu$rf22 > 0.60
jhu_dec <- jhu[stats::complete.cases(jhu[, c("capras","rf22")]) & jhu$time > 0, ]
jhu_dec$rf10 <- jhu_dec$rf22 * 10
dec_fit <- tryCatch(survival::cch(Surv(time, evt == 1) ~ capras + rf10, data = jhu_dec,
             subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE),
           error = function(e) NULL)
dec_hr <- dec_lo <- dec_hi <- NA_real_
if (!is.null(dec_fit)) {
  dcf <- summary(dec_fit)$coefficients
  dec_hr <- exp(dcf["rf10","Value"])
  dec_lo <- exp(dcf["rf10","Value"] - 1.96 * dcf["rf10","SE"])
  dec_hi <- exp(dcf["rf10","Value"] + 1.96 * dcf["rf10","SE"])
}
anchor <- data.frame(
  quantity = c("weighted_10y_CIF_Decipher_low_nonparam","weighted_10y_CIF_Decipher_high_nonparam",
               "RossModel2_LinYing_robust_Decipher_HR_per0.1"),
  observed = c(wcif1_fast(jhu$time[lowD], jhu$evt[lowD], jhu$w[lowD], T10),
               wcif1_fast(jhu$time[highD], jhu$evt[highD], jhu$w[highD], T10), dec_hr),
  observed_lo = c(NA, NA, dec_lo), observed_hi = c(NA, NA, dec_hi),
  published = c(0.12, 0.47, 1.32), published_lo = c(NA, NA, 1.17), published_hi = c(NA, NA, 1.51),
  spec = c("weighted Aalen-Johansen CIF","weighted Aalen-Johansen CIF",
           "Ross Table 3 Model 2: CAPRA-S + Decipher/0.1, Lin-Ying robust cch"))

###############################################################################
# JHU internal: apparent weighted DCA + performance, design-aware optimism
###############################################################################
jperf <- list(); jnb <- list()
for (t0 in HORIZ) {
  cc <- jhu_cc[[paste0("cif_clin_", t0)]]; cf <- jhu_cc[[paste0("cif_full_", t0)]]
  jnb[[as.character(t0)]] <- list(
    clin = nb_curve_w(cc, jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, t0),
    full = nb_curve_w(cf, jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, t0),
    all  = nb_all_w(jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, t0))
  bc <- wbrier(cc, jhu_cc$time, jhu_cc$evt, jhu_cc$w, t0)
  bf <- wbrier(cf, jhu_cc$time, jhu_cc$evt, jhu_cc$w, t0)
  ac <- wauc_def2(jhu_cc$time, jhu_cc$evt, cc, jhu_cc$w, t0)
  af <- wauc_def2(jhu_cc$time, jhu_cc$evt, cf, jhu_cc$w, t0)
  jperf[[as.character(t0)]] <- list(brier_clin=bc, brier_full=bf, auc_clin=ac, auc_full=af)
}

# conditional stratified optimism correction: resample WITHIN the three observed
# sample groups (subcohort cases, subcohort controls, out-of-subcohort cases),
# refit, apparent-minus-test. This does not recreate the 265/745 phase-one
# sampling design; it is a point bias correction only, with no design-valid
# interval (the two-phase frame is unavailable).
strata_idx <- split(seq_len(nrow(jhu_cc)), jhu_cc$cch)
opt_nb <- list(); for (t0 in HORIZ) opt_nb[[as.character(t0)]] <-
  list(clin = matrix(NA_real_, B, length(THR)), full = matrix(NA_real_, B, length(THR)))
opt_m <- list(); for (t0 in HORIZ) opt_m[[as.character(t0)]] <-
  matrix(NA_real_, B, 2, dimnames = list(NULL, c("auc_clin","auc_full")))
# 10-year JHU clinical-impact within the same stratified resample (conditional
# internal; not phase-one design-valid): augmented-minus-clinical change per 100
jhu_imp_dsel <- matrix(NA_real_, B, length(IMPL)); jhu_imp_dmet <- matrix(NA_real_, B, length(IMPL))
nboot_ok <- 0L; nboot_fail <- 0L
for (b in seq_len(B)) {
  ix <- unlist(lapply(strata_idx, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
  cb <- jhu_cc[ix, ]; cb$id <- seq_len(nrow(cb))
  Mc <- tryCatch(fit_wfg(cb, clin_terms), error = function(e) NULL)
  Mf <- tryCatch(fit_wfg(cb, full_terms), error = function(e) NULL)
  if (is.null(Mc) || is.null(Mf)) { nboot_fail <- nboot_fail + 1L; next }
  ok <- TRUE
  for (t0 in HORIZ) {
    ccb <- tryCatch(cif_wfg(Mc, cb, clin_terms, t0), error=function(e) NULL)
    cfb <- tryCatch(cif_wfg(Mf, cb, full_terms, t0), error=function(e) NULL)
    cct <- tryCatch(cif_wfg(Mc, jhu_cc, clin_terms, t0), error=function(e) NULL)
    cft <- tryCatch(cif_wfg(Mf, jhu_cc, full_terms, t0), error=function(e) NULL)
    if (any(sapply(list(ccb,cfb,cct,cft), function(z) is.null(z) || any(!is.finite(z))))) { ok <- FALSE; break }
    app_c <- nb_curve_w(ccb, cb$time, cb$evt, cb$w, THR, t0); tst_c <- nb_curve_w(cct, jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, t0)
    app_f <- nb_curve_w(cfb, cb$time, cb$evt, cb$w, THR, t0); tst_f <- nb_curve_w(cft, jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, t0)
    opt_nb[[as.character(t0)]]$clin[b, ] <- app_c - tst_c
    opt_nb[[as.character(t0)]]$full[b, ] <- app_f - tst_f
    aa_c <- wauc_def2(cb$time, cb$evt, ccb, cb$w, t0); ta_c <- wauc_def2(jhu_cc$time, jhu_cc$evt, cct, jhu_cc$w, t0)
    aa_f <- wauc_def2(cb$time, cb$evt, cfb, cb$w, t0); ta_f <- wauc_def2(jhu_cc$time, jhu_cc$evt, cft, jhu_cc$w, t0)
    opt_m[[as.character(t0)]][b, ] <- c(aa_c - ta_c, aa_f - ta_f)
    if (t0 == T10) for (jj in seq_along(IMPL)) { pt <- IMPL[jj]
      sf <- cfb >= pt; sc <- ccb >= pt; Wt <- sum(cb$w)
      pif <- sum(cb$w[sf])/Wt; pic <- sum(cb$w[sc])/Wt
      cff <- if (any(sf)) cif_pos_fast(cb$time, cb$evt, cb$w, sf, T10) else 0
      cfc <- if (any(sc)) cif_pos_fast(cb$time, cb$evt, cb$w, sc, T10) else 0
      jhu_imp_dsel[b, jj] <- 100*(pif - pic); jhu_imp_dmet[b, jj] <- 100*(pif*cff - pic*cfc)
    }
  }
  if (ok) nboot_ok <- nboot_ok + 1L else nboot_fail <- nboot_fail + 1L
}
jnb_corr <- list()
for (t0 in HORIZ) {
  ch <- as.character(t0)
  jnb_corr[[ch]] <- list(
    clin = jnb[[ch]]$clin - colMeans(opt_nb[[ch]]$clin, na.rm = TRUE),
    full = jnb[[ch]]$full - colMeans(opt_nb[[ch]]$full, na.rm = TRUE))
  om <- colMeans(opt_m[[ch]], na.rm = TRUE)
  jperf[[ch]]$brier_clin_corr <- NA_real_   # apparent design-weighted Brier/IPA only
  jperf[[ch]]$brier_full_corr <- NA_real_
  jperf[[ch]]$auc_clin_corr   <- jperf[[ch]]$auc_clin - om["auc_clin"]
  jperf[[ch]]$auc_full_corr   <- jperf[[ch]]$auc_full - om["auc_full"]
}
cat(sprintf("JHU conditional stratified optimism correction: %d ok, %d failed (B=%d)\n", nboot_ok, nboot_fail, B))

###############################################################################
# Durham external (frozen JHU models, complete cohort, unweighted)
###############################################################################
w1 <- rep(1, nrow(dur_cc))
dnb <- list(); dperf <- list()
for (t0 in HORIZ) {
  cc <- dur_cc[[paste0("cif_clin_", t0)]]; cf <- dur_cc[[paste0("cif_full_", t0)]]
  dnb[[as.character(t0)]] <- list(
    clin = nb_curve_w(cc, dur_cc$time, dur_cc$evt, w1, THR, t0),
    full = nb_curve_w(cf, dur_cc$time, dur_cc$evt, w1, THR, t0),
    all  = nb_all_w(dur_cc$time, dur_cc$evt, w1, THR, t0))
  bc <- wbrier(cc, dur_cc$time, dur_cc$evt, w1, t0); bf <- wbrier(cf, dur_cc$time, dur_cc$evt, w1, t0)
  ac <- wauc_def2(dur_cc$time, dur_cc$evt, cc, w1, t0); af <- wauc_def2(dur_cc$time, dur_cc$evt, cf, w1, t0)
  dperf[[as.character(t0)]] <- list(brier_clin=bc, brier_full=bf, auc_clin=ac, auc_full=af)
}
# Durham paired patient-level percentile bootstrap for net benefit (pointwise) and
# for the 10-year clinical-impact quantities: augmented-minus-clinical change in
# patients selected and metastases captured per 100 at the fixed IMPL thresholds
nd <- nrow(dur_cc); dboot <- list()
for (t0 in HORIZ) dboot[[as.character(t0)]] <- list(clin=matrix(NA_real_,B,length(THR)),
                                                    full=matrix(NA_real_,B,length(THR)),
                                                    diff=matrix(NA_real_,B,length(THR)))
imp_dsel <- matrix(NA_real_, B, length(IMPL)); imp_dmet <- matrix(NA_real_, B, length(IMPL))
durF10 <- dur_cc[[paste0("cif_full_",T10)]]; durC10 <- dur_cc[[paste0("cif_clin_",T10)]]
for (b in seq_len(B)) {
  ix <- sample.int(nd, nd, replace = TRUE)
  for (t0 in HORIZ) {
    cvec <- nb_curve_w(dur_cc[[paste0("cif_clin_",t0)]][ix], dur_cc$time[ix], dur_cc$evt[ix], w1[ix], THR, t0)
    fvec <- nb_curve_w(dur_cc[[paste0("cif_full_",t0)]][ix], dur_cc$time[ix], dur_cc$evt[ix], w1[ix], THR, t0)
    dboot[[as.character(t0)]]$clin[b,] <- cvec
    dboot[[as.character(t0)]]$full[b,] <- fvec
    dboot[[as.character(t0)]]$diff[b,] <- fvec - cvec
  }
  tt <- dur_cc$time[ix]; ee <- dur_cc$evt[ix]; fp <- durF10[ix]; cp <- durC10[ix]; nn <- length(ix)
  for (jj in seq_along(IMPL)) { pt <- IMPL[jj]
    sf <- fp >= pt; sc <- cp >= pt; pif <- mean(sf); pic <- mean(sc)
    cff <- if (any(sf)) cif_pos_fast(tt, ee, rep(1,nn), sf, T10) else 0
    cfc <- if (any(sc)) cif_pos_fast(tt, ee, rep(1,nn), sc, T10) else 0
    imp_dsel[b, jj] <- 100*(pif - pic); imp_dmet[b, jj] <- 100*(pif*cff - pic*cfc)
  }
}
qcol <- function(m,p) apply(m, 2, function(v) quantile(v, p, na.rm=TRUE))

# Durham external calibration at 10y (frozen full model): O/E, cloglog CIL + slope.
# A bootstrap fit that hits a singular gradient is counted as a failed replicate
# and its coefficients are discarded; valid/failed counts are kept separately for
# the joint intercept/slope fit and the calibration-in-the-large fit. The flexible
# calibration band is a paired patient bootstrap on a fixed grid, not a loess
# standard-error band.
clog <- function(p) log(-log(1 - pmin(pmax(p, 1e-6), 1 - 1e-6))); cinv <- function(x) 1 - exp(-exp(x))
pv_jack <- function(time, evt, t0) as.numeric(prodlim::jackknife(
  prodlim::prodlim(prodlim::Hist(time, evt) ~ 1, data = data.frame(time, evt)), times = t0, cause = 1))
nls_fit <- function(form, start) {
  bad <- FALSE
  fit <- withCallingHandlers(
    tryCatch(nls(form, start = start, control = nls.control(maxiter = 200, warnOnly = TRUE)),
             error = function(e) NULL),
    warning = function(w) {
      .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = "calibration_nls",
        message = conditionMessage(w), stringsAsFactors = FALSE)
      if (grepl("singular|gradient|converg|step factor|minFactor", conditionMessage(w), ignore.case = TRUE)) bad <<- TRUE
      invokeRestart("muffleWarning")
    })
  if (is.null(fit) || bad) NULL else fit
}
obs_dur <- wcif1_fast(dur_cc$time, dur_cc$evt, w1, T10)
predD   <- dur_cc[[paste0("cif_full_", T10)]]
OE_dur  <- obs_dur / mean(predD)
pvD     <- pv_jack(dur_cc$time, dur_cc$evt, T10); lpD <- clog(predD)
fj  <- nls_fit(pvD ~ cinv(a + b * lpD), list(a = 0, b = 1))
fc0 <- nls_fit(pvD ~ cinv(a + lpD),     list(a = 0))
cal_int <- if (!is.null(fj))  coef(fj)["a"]  else NA_real_
cal_slp <- if (!is.null(fj))  coef(fj)["b"]  else NA_real_
cal_cil <- if (!is.null(fc0)) coef(fc0)["a"] else NA_real_
gp <- seq(min(predD), max(predD), length.out = 60)
cal_boot <- matrix(NA_real_, B, 3); flx_boot <- matrix(NA_real_, B, length(gp))
cal_joint_ok <- 0L; cal_joint_fail <- 0L; cal_cil_ok <- 0L; cal_cil_fail <- 0L
for (b in seq_len(B)) {
  ix <- sample.int(nd, nd, replace = TRUE)
  pvb <- tryCatch(pv_jack(dur_cc$time[ix], dur_cc$evt[ix], T10), error = function(e) NULL)
  if (is.null(pvb)) { cal_joint_fail <- cal_joint_fail + 1L; cal_cil_fail <- cal_cil_fail + 1L; next }
  lpb <- lpD[ix]
  fjb <- nls_fit(pvb ~ cinv(a + b * lpb), list(a = 0, b = 1))
  fcb <- nls_fit(pvb ~ cinv(a + lpb),     list(a = 0))
  if (!is.null(fjb)) { cal_boot[b, 1:2] <- coef(fjb); cal_joint_ok <- cal_joint_ok + 1L } else cal_joint_fail <- cal_joint_fail + 1L
  if (!is.null(fcb)) { cal_boot[b, 3]   <- coef(fcb); cal_cil_ok   <- cal_cil_ok   + 1L } else cal_cil_fail   <- cal_cil_fail   + 1L
  lob <- tryCatch(loess(pv ~ pr, data = data.frame(pv = pvb, pr = predD[ix]), span = 0.9, degree = 1), error = function(e) NULL)
  if (!is.null(lob)) flx_boot[b, ] <- suppressWarnings(predict(lob, newdata = data.frame(pr = gp)))
}
qv <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
lo_full <- loess(pv ~ pr, data = data.frame(pv = pvD, pr = predD), span = 0.9, degree = 1)
flexcal <- data.frame(pred = gp, obs = as.numeric(predict(lo_full, newdata = data.frame(pr = gp))),
                      lo = apply(flx_boot, 2, qv, 0.025), hi = apply(flx_boot, 2, qv, 0.975))
cal_conv <- data.frame(
  fit = c("joint_intercept_slope", "calibration_in_large"),
  point_available = c(!is.null(fj), !is.null(fc0)),
  boot_valid = c(cal_joint_ok, cal_cil_ok), boot_failed = c(cal_joint_fail, cal_cil_fail))

###############################################################################
# CAPRA-S (JHU only), reconstructed face value; component audit + rule DCA
###############################################################################
# CAPRA-S component scores are computed earlier (for the Decipher Model-2 anchor)
cap_cc <- jhu[!is.na(jhu$capras) & jhu$time > 0, ]
.n_below3 <- sum(jhu$capras < 3, na.rm = TRUE); .n_ok <- sum(!is.na(jhu$capras))
capras_audit <- data.frame(
  component = c("PSA(0-3)","Gleason(0-3,incl 4+3=2)","SM(0/2)","ECE(0/1)","SVI(0/2)","LNI(0/1)",
                "total","capras_below3_vs_RossTable4_low","n_nonmissing","n_missing_total"),
  summary = c(
    paste(names(table(jhu$capras_psa)), table(jhu$capras_psa), sep=":", collapse=" "),
    paste(names(table(jhu$capras_gl,  useNA="ifany")), table(jhu$capras_gl, useNA="ifany"), sep=":", collapse=" "),
    paste(names(table(jhu$capras_sm)), table(jhu$capras_sm), sep=":", collapse=" "),
    paste(names(table(jhu$capras_ece)),table(jhu$capras_ece),sep=":", collapse=" "),
    paste(names(table(jhu$capras_svi)),table(jhu$capras_svi),sep=":", collapse=" "),
    paste(names(table(jhu$capras_lni)),table(jhu$capras_lni),sep=":", collapse=" "),
    sprintf("median=%s range=%s", median(jhu$capras,na.rm=TRUE), paste(range(jhu$capras,na.rm=TRUE),collapse="-")),
    sprintf("%d (%.1f%%); Ross Table 4 CAPRA-S-low = 2.2%%", .n_below3, 100 * .n_below3 / .n_ok),
    as.character(.n_ok), as.character(sum(is.na(jhu$capras)))))
# rule strategies (weighted): CAPRA-S high (>=6) vs selective reflex
highpos <- cap_cc$capras >= 6
reflex  <- highpos | (cap_cc$capras >= 3 & cap_cc$capras <= 5 & cap_cc$score >= MS_THR)
rule_nb <- function(pos, t0){
  Wtot<-sum(cap_cc$w); pi_pos<-sum(cap_cc$w[pos])/Wtot; cifp<-cif_pos_fast(cap_cc$time,cap_cc$evt,cap_cc$w,pos,t0)
  100*pi_pos*(cifp - (1-cifp)*(THR/(1-THR)))
}
capras_rule <- do.call(rbind, lapply(HORIZ, function(t0){
  data.frame(horizon=t0, threshold=THR,
    nb_capras_high = rule_nb(highpos,t0), nb_reflex = rule_nb(reflex,t0),
    nb_treat_all = nb_all_w(cap_cc$time,cap_cc$evt,cap_cc$w,THR,t0))
}))
# reflex impact quantities (weighted)
Wtot<-sum(cap_cc$w); pct_tested <- 100*sum(cap_cc$w[cap_cc$capras>=3 & cap_cc$capras<=5])/Wtot
added <- 100*(sum(cap_cc$w[reflex])-sum(cap_cc$w[highpos]))/Wtot
capras_impact <- data.frame(
  quantity=c("pct_intermediate_requiring_MetScore_test","pct_added_by_reflex_vs_high",
             "n_capras_high","n_reflex_positive","MetScore_threshold"),
  value=c(pct_tested, added, sum(highpos), sum(reflex), MS_THR))

###############################################################################
# full-239 Lin-Ying case-cohort relative-association sensitivity (cause-specific)
###############################################################################
cch_row <- tryCatch({
  sc <- survival::cch(Surv(time, evt == 1) ~ GG_ord + log2PSA + pT_hi + ms_dev,
          data = jhu_cc, subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE)
  cf <- summary(sc)$coefficients
  data.frame(term = rownames(cf), coef = cf[,"Value"], robust_se = cf[,"SE"],
             HR = exp(cf[,"Value"]), lo = exp(cf[,"Value"]-1.96*cf[,"SE"]), hi = exp(cf[,"Value"]+1.96*cf[,"SE"]),
             robust_p = 2*pnorm(-abs(cf[,"Z"])),
             estimator = "LinYing_cause_specific_cch_cohort.size745_robust", row.names=NULL)
}, error = function(e) data.frame(term="UNAVAILABLE", coef=NA, robust_se=NA, HR=NA, lo=NA, hi=NA,
             robust_p=NA, estimator=paste("cch_failed:", conditionMessage(e))))
# weighted FG ms_dev coefficient: a noninferential working-model point HR retained
# only for the frozen absolute-risk prediction, not case-cohort inference
fg_ms <- as.numeric(coef(M_full)["ms_dev"])

cat("=== analysis complete; writing outputs ===\n")

###############################################################################
# frozen deployable artifact (no patient identifiers)
###############################################################################
base_cif <- function(fit, terms, t0) {
  z0 <- as.data.frame(as.list(setNames(rep(0, length(terms)), terms)))
  sf <- survfit(fit, newdata = z0); ti <- max(which(sf$time <= t0)); as.numeric(1 - sf$surv[ti])
}
frozen <- list(meta = list(horizons = HORIZ, alpha = ALPHA, clin_terms = clin_terms,
                 full_terms = full_terms, mean_JHU = mean_JHU, sd_JHU = sd_JHU,
                 note = "weighted Fine-Gray (IPW subdistribution); ms_dev on design-weighted JHU scale"),
               coef_clin = coef(M_clin), coef_full = coef(M_full),
               M_clin = M_clin, M_full = M_full)
saveRDS(frozen, file.path(OUTD, "frozen_combined_model.rds"))
fr <- readRDS(file.path(OUTD, "frozen_combined_model.rds"))
repro_err <- max(abs(cif_wfg(fr$M_full, dur_cc, full_terms, T10) - dur_cc[[paste0("cif_full_",T10)]]))

###############################################################################
# mathematical / reduction checks
###############################################################################
mc <- list()
addmc <- function(check, value, ref, tol) mc[[length(mc)+1L]] <<-
  data.frame(check=check, value=value, reference=ref, diff=value-ref, pass=abs(value-ref) <= tol)
# all-positive predictions reproduce treat-all; all-negative reproduce treat-none
nb_allpos <- nb_curve_w(rep(1, nrow(jhu_cc)), jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, T10)
nb_allneg <- nb_curve_w(rep(0, nrow(jhu_cc)), jhu_cc$time, jhu_cc$evt, jhu_cc$w, THR, T10)
addmc("all_positive_equals_treat_all_(max abs dev)", max(abs(nb_allpos - jnb[["120"]]$all)), 0, 1e-8)
addmc("all_negative_equals_treat_none_(max abs)", max(abs(nb_allneg)), 0, 1e-8)
addmc("fast_CIF_equals_survfit_AJ", wcif1_fast(jhu_cc$time,jhu_cc$evt,jhu_cc$w,T10),
      aj_cif_w(jhu_cc$time,jhu_cc$evt,jhu_cc$w,rep(TRUE,nrow(jhu_cc)),T10), 1e-8)
# w=1 Brier reduces to Score and w=1 AUC reproduces timeROC def-2, both at 1e-6
p10 <- dur_cc[[paste0("cif_full_",T10)]]
sc1 <- riskRegression::Score(list(m=p10), formula=Hist(time,evt)~1, data=dur_cc, times=T10, cause=1,
         metrics="brier", cens.model="km", conf.int=FALSE, null.model=FALSE)
addmc("w1_Brier_reduces_to_Score(Durham,120)", unname(wbrier(p10, dur_cc$time, dur_cc$evt, w1, T10)["brier"]),
      as.data.frame(sc1$Brier$score)$Brier[1], 1e-6)
tr1 <- timeROC::timeROC(T=dur_cc$time, delta=dur_cc$evt, marker=p10, cause=1, weighting="marginal", times=T10, iid=FALSE)
addmc("w1_AUC_reproduces_timeROC_def2(Durham,120)", wauc_def2(dur_cc$time, dur_cc$evt, p10, w1, T10),
      as.numeric(tail(tr1$AUC_2,1)), 1e-6)
addmc("anchor_weighted_CIF_Decipher_low_vs_Ross0.12", anchor$observed[1], 0.12, 0.03)
addmc("anchor_weighted_CIF_Decipher_high_vs_Ross0.47", anchor$observed[2], 0.47, 0.03)
addmc("frozen_reproduces_Durham_predictions", repro_err, 0, 1e-8)
# intercept-only weighted Fine-Gray reproduces the weighted Aalen-Johansen CIF:
# the baseline / estimating-equation equivalence gate on which absolute risk and
# the model-based DCA depend
fg0 <- survival::finegray(Surv(time, evtf) ~ ., etype = "met", weights = w, id = id,
        data = data.frame(time = jhu_cc$time, evtf = evtf(jhu_cc$evt), w = jhu_cc$w, id = jhu_cc$id))
M0  <- survival::coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, data = fg0, weights = fgwt, ties = "breslow")
bh0 <- survival::basehaz(M0, centered = FALSE)
fg0_cif <- function(t0){ H <- bh0$hazard[max(which(bh0$time <= t0))]; 1 - exp(-H) }
for (t0 in HORIZ)
  addmc(sprintf("interceptFG_reproduces_weightedAJ_CIF(%d)", t0), fg0_cif(t0),
        aj_cif_w(jhu_cc$time, jhu_cc$evt, jhu_cc$w, rep(TRUE, nrow(jhu_cc)), t0), 5e-3)
math_checks <- do.call(rbind, mc)

# estimator-equivalence + permutation-invariance: Brier/IPA/AUC/DCA are invariant
# to an arbitrary row permutation (the values are computed from the data, not
# hard-coded), and the all-one Brier reduces to Score
set.seed(20260813L); prm <- sample.int(nrow(jhu_cc)); pj <- jhu_cc[prm, ]; cf10col <- paste0("cif_full_", T10)
b_ref <- jperf[["120"]]$brier_full; a_ref <- jperf[["120"]]$auc_full; nb_ref <- jnb[["120"]]$full
b_prm <- wbrier(pj[[cf10col]], pj$time, pj$evt, pj$w, T10)
a_prm <- wauc_def2(pj$time, pj$evt, pj[[cf10col]], pj$w, T10)
nb_prm <- nb_curve_w(pj[[cf10col]], pj$time, pj$evt, pj$w, THR, T10)
perm_checks <- data.frame(
  check = c("Brier_permutation_invariant","IPA_permutation_invariant","AUC_permutation_invariant",
            "DCA_permutation_invariant_max_abs","allone_Brier_reduces_to_Score_120"),
  value = c(b_prm["brier"], b_prm["ipa"], a_prm, max(abs(nb_prm - nb_ref)),
            unname(wbrier(p10, dur_cc$time, dur_cc$evt, w1, T10)["brier"])),
  reference = c(b_ref["brier"], b_ref["ipa"], a_ref, 0, as.data.frame(sc1$Brier$score)$Brier[1]),
  tol = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-6), row.names = NULL)
perm_checks$diff <- perm_checks$value - perm_checks$reference
perm_checks$pass <- abs(perm_checks$diff) <= perm_checks$tol

###############################################################################
# write CSV outputs
###############################################################################
W <- function(df, name) write.csv(df, file.path(OUTD, name), row.names = FALSE)

W(data.frame(
  item = c("JHU_cohort_n","JHU_events","subcohort_cases","subcohort_controls","outside_cases",
           "JHU_analytic_n","JHU_analytic_events","complete_case_excluded",
           "excluded_GG_ord","excluded_pT_hi","sampling_fraction_265_745","control_weight_745_265",
           "design_weighted_score_mean","design_weighted_score_sd",
           "Durham_n","Durham_events"),
  value = c(239, 93, 28, 146, 65, nrow(jhu_cc), sum(jhu_cc$evt==1L), 239-nrow(jhu_cc),
            excl$n_missing_jhu[excl$variable=="GG_ord"], excl$n_missing_jhu[excl$variable=="pT_hi"],
            ALPHA, 1/ALPHA, mean_JHU, sd_JHU, nrow(dur_cc), sum(dur_cc$evt==1L))),
  "sample_event_ledger.csv")

# the weighted Fine-Gray coxph sandwich se/p are noninferential working-model
# quantities (not case-cohort inference); only coef/HR feed the frozen prediction
coef_tab <- function(fit, model){ s<-summary(fit)$coefficients
  data.frame(model=model, term=rownames(s), coef=s[,"coef"], HR=exp(s[,"coef"]),
             wfg_working_se_noninferential=s[,"robust se"], wfg_working_p_noninferential=s[,"Pr(>|z|)"], row.names=NULL) }
W(rbind(coef_tab(M_clin,"M_clin"), coef_tab(M_full,"M_full"),
        data.frame(model="standardization", term=c("mean_JHU","sd_JHU"),
                   coef=c(mean_JHU,sd_JHU), HR=NA, wfg_working_se_noninferential=NA, wfg_working_p_noninferential=NA),
        data.frame(model="baseline_CIF", term=c(paste0("clin_",HORIZ),paste0("full_",HORIZ)),
                   coef=c(sapply(HORIZ,function(t) base_cif(M_clin,clin_terms,t)),
                          sapply(HORIZ,function(t) base_cif(M_full,full_terms,t))),
                   HR=NA, wfg_working_se_noninferential=NA, wfg_working_p_noninferential=NA)),
  "model_coefficients_baselines.csv")

perf_rows <- do.call(rbind, lapply(HORIZ, function(t0){ ch<-as.character(t0); p<-jperf[[ch]]; d<-dperf[[ch]]
  rbind(
    data.frame(cohort="JHU", basis="apparent",  horizon=t0, model="clin", AUC=p$auc_clin, Brier=p$brier_clin["brier"], IPA=p$brier_clin["ipa"]),
    data.frame(cohort="JHU", basis="apparent",  horizon=t0, model="full", AUC=p$auc_full, Brier=p$brier_full["brier"], IPA=p$brier_full["ipa"]),
    data.frame(cohort="JHU", basis="optimism_corrected", horizon=t0, model="clin", AUC=p$auc_clin_corr, Brier=p$brier_clin_corr, IPA=NA),
    data.frame(cohort="JHU", basis="optimism_corrected", horizon=t0, model="full", AUC=p$auc_full_corr, Brier=p$brier_full_corr, IPA=NA),
    data.frame(cohort="Durham", basis="frozen_external", horizon=t0, model="clin", AUC=d$auc_clin, Brier=d$brier_clin["brier"], IPA=d$brier_clin["ipa"]),
    data.frame(cohort="Durham", basis="frozen_external", horizon=t0, model="full", AUC=d$auc_full, Brier=d$brier_full["brier"], IPA=d$brier_full["ipa"])) }))
W(perf_rows, "performance_5y10y.csv")

W(rbind(
  data.frame(metric=c("Durham_observed_AJ_CIF_120","Durham_mean_pred_120","OE_ratio","EO_factor",
                      "CIL_intercept","joint_intercept","joint_intercept_lo","joint_intercept_hi",
                      "joint_slope","joint_slope_lo","joint_slope_hi"),
             value=c(obs_dur, mean(predD), OE_dur, 1/OE_dur, cal_cil, cal_int, qv(cal_boot[,1],.025), qv(cal_boot[,1],.975),
                     cal_slp, qv(cal_boot[,2],.025), qv(cal_boot[,2],.975)), pred=NA, obs=NA, lo=NA, hi=NA),
  data.frame(metric="flexible_curve", value=NA, pred=flexcal$pred, obs=flexcal$obs, lo=flexcal$lo, hi=flexcal$hi)),
  "calibration_5y10y.csv")

dca_rows <- list()
for (t0 in HORIZ) { ch<-as.character(t0)
  dl<-qcol(dboot[[ch]]$clin,.025); dh<-qcol(dboot[[ch]]$clin,.975); fl<-qcol(dboot[[ch]]$full,.025); fh<-qcol(dboot[[ch]]$full,.975)
  dca_rows[[length(dca_rows)+1]] <- rbind(
    data.frame(cohort="JHU", basis="apparent", model="clin", horizon=t0, threshold=THR, nb=jnb[[ch]]$clin, lo=NA, hi=NA),
    data.frame(cohort="JHU", basis="apparent", model="full", horizon=t0, threshold=THR, nb=jnb[[ch]]$full, lo=NA, hi=NA),
    data.frame(cohort="JHU", basis="optimism_corrected", model="clin", horizon=t0, threshold=THR, nb=jnb_corr[[ch]]$clin, lo=NA, hi=NA),
    data.frame(cohort="JHU", basis="optimism_corrected", model="full", horizon=t0, threshold=THR, nb=jnb_corr[[ch]]$full, lo=NA, hi=NA),
    data.frame(cohort="JHU", basis="apparent", model="treat_all", horizon=t0, threshold=THR, nb=jnb[[ch]]$all, lo=NA, hi=NA),
    data.frame(cohort="JHU", basis="apparent", model="treat_none", horizon=t0, threshold=THR, nb=0, lo=NA, hi=NA),
    data.frame(cohort="Durham", basis="frozen", model="clin", horizon=t0, threshold=THR, nb=dnb[[ch]]$clin, lo=dl, hi=dh),
    data.frame(cohort="Durham", basis="frozen", model="full", horizon=t0, threshold=THR, nb=dnb[[ch]]$full, lo=fl, hi=fh),
    data.frame(cohort="Durham", basis="frozen", model="treat_all", horizon=t0, threshold=THR, nb=dnb[[ch]]$all, lo=NA, hi=NA),
    data.frame(cohort="Durham", basis="frozen", model="treat_none", horizon=t0, threshold=THR, nb=0, lo=NA, hi=NA)) }
W(do.call(rbind, dca_rows), "dca_curves.csv")

# JHU delta NB and full_exceeds_both use the conditional-stratified-optimism-
# corrected curves (jnb_corr), matching Figure S1 clinical-utility candidate; Durham uses frozen external
dnb_rows <- do.call(rbind, lapply(HORIZ, function(t0){ ch<-as.character(t0)
  jd <- jnb_corr[[ch]]$full - jnb_corr[[ch]]$clin
  ddl<-qcol(dboot[[ch]]$diff,.025); ddh<-qcol(dboot[[ch]]$diff,.975)
  rbind(
    data.frame(cohort="JHU", basis="conditional_stratified_optimism_corrected", horizon=t0, threshold=THR,
               delta_nb=jd, lo=NA, hi=NA, full_exceeds_both=jnb_corr[[ch]]$full > pmax(jnb[[ch]]$all,0)),
    data.frame(cohort="Durham", basis="frozen_external", horizon=t0, threshold=THR,
               delta_nb=dnb[[ch]]$full-dnb[[ch]]$clin, lo=ddl, hi=ddh,
               full_exceeds_both=dnb[[ch]]$full > pmax(dnb[[ch]]$all,0))) }))
W(dnb_rows, "delta_net_benefit.csv")

impl <- do.call(rbind, lapply(HORIZ, function(t0){ ch<-as.character(t0)
  do.call(rbind, lapply(c("JHU","Durham"), function(coh){
    if(coh=="JHU"){ tt<-jhu_cc$time; ee<-jhu_cc$evt; ww<-jhu_cc$w
      cc<-jhu_cc[[paste0("cif_clin_",t0)]]; cf<-jhu_cc[[paste0("cif_full_",t0)]]; boot<-NULL
    } else { tt<-dur_cc$time; ee<-dur_cc$evt; ww<-w1
      cc<-dur_cc[[paste0("cif_clin_",t0)]]; cf<-dur_cc[[paste0("cif_full_",t0)]]; boot<-dboot[[ch]]$diff }
    Wtot<-sum(ww); cif_all<-wcif1_fast(tt,ee,ww,t0)
    do.call(rbind, lapply(IMPL, function(pt){
      selc<-cf>=pt; selk<-cc>=pt
      pif<-sum(ww[selc])/Wtot; cff<-if(any(selc)) cif_pos_fast(tt,ee,ww,selc,t0) else NA
      pic<-sum(ww[selk])/Wtot; cfc<-if(any(selk)) cif_pos_fast(tt,ee,ww,selk,t0) else NA
      mets_f<-100*pif*cff; mets_c<-100*pic*cfc
      unnec_f<-100*pif*(1-cff); unnec_all<-100*(1-cif_all)
      nb_f<-100*pif*(cff-(1-cff)*pt/(1-pt)); nb_c<-100*pic*(cfc-(1-cfc)*pt/(1-pt))
      dnbv<-nb_f-nb_c
      ex_both <- nb_f > max(100*(cif_all-(1-cif_all)*pt/(1-pt)), 0)
      lo<-hi<-NA; if(!is.null(boot)){ j<-which.min(abs(THR-pt)); lo<-quantile(boot[,j],.025,na.rm=TRUE); hi<-quantile(boot[,j],.975,na.rm=TRUE) }
      tt_ratio <- if(!is.na(dnbv) && dnbv>0 && ex_both) (100*pif)/dnbv else NA
      data.frame(cohort=coh, basis=if(coh=="JHU") "internal_apparent_descriptive" else "frozen_external",
                 horizon=t0, threshold=pt, eligible=length(tt),
                 selected_per100_full=100*pif, selected_per100_clin=100*pic,
                 mets_captured_per100_full=mets_f, mets_captured_per100_clin=mets_c,
                 incremental_mets_per100=mets_f-mets_c, additional_evaluations_per100=100*(pif-pic),
                 unnecessary_evaluations_per100_full=unnec_f, net_reduction_unnecessary_per100=unnec_all-unnec_f,
                 delta_nb=dnbv, delta_nb_lo=lo, delta_nb_hi=hi, full_exceeds_both_defaults=ex_both,
                 tests_per_net_TP_equiv=tt_ratio) })) })) }))
W(impl, "implementation_table.csv")

W(capras_audit, "capras_component_audit.csv")
W(rbind(
  data.frame(type="rule_dca", capras_rule, quantity=NA, qval=NA),
  data.frame(type="impact", horizon=NA, threshold=NA, nb_capras_high=NA, nb_reflex=NA, nb_treat_all=NA,
             quantity=capras_impact$quantity, qval=capras_impact$value)),
  "capras_rule_utility.csv")

W(rbind(cch_row,
        data.frame(term="ms_dev_weighted_FG_subdistribution_noninferential", coef=fg_ms, robust_se=NA,
                   HR=exp(fg_ms), lo=NA, hi=NA, robust_p=NA, estimator="IPW_weighted_finegray_coxph_point_only")),
  "full239_association_sensitivity.csv")

W(data.frame(quantity=c("B","seed","JHU_optimism_ok","JHU_optimism_fail",
                        "Durham_paired_boot","design_variance","note"),
             value=c(B, 1, nboot_ok, nboot_fail, B, "two_phase_frame_unavailable",
                     "JHU DCA uncertainty limited to a conditional stratified optimism point correction; Lin-Ying cch gives design-valid relative inference")),
  "bootstrap_ledger.csv")

W(math_checks, "mathematical_checks.csv")
W(perm_checks, "estimator_equivalence_permutation_checks.csv")
W(cal_conv, "calibration_convergence_ledger.csv")
W(anchor, "decipher_anchor_check.csv")

# 10-year clinical impact / resource use, both cohorts (panel d source):
# augmented-minus-clinical change in patients selected and metastases captured per
# 100. JHU is case-cohort weighted with conditional-internal stratified-bootstrap
# intervals (not phase-one design-valid); Durham is frozen external with paired
# patient-bootstrap 95% intervals.
imp_durham <- do.call(rbind, lapply(seq_along(IMPL), function(jj){ pt<-IMPL[jj]
  sf<-durF10>=pt; sc<-durC10>=pt; pif<-mean(sf); pic<-mean(sc)
  cff<-if(any(sf)) cif_pos_fast(dur_cc$time,dur_cc$evt,w1,sf,T10) else 0
  cfc<-if(any(sc)) cif_pos_fast(dur_cc$time,dur_cc$evt,w1,sc,T10) else 0
  j<-which.min(abs(THR-pt))
  data.frame(cohort="Durham", basis="frozen external", horizon=T10, threshold=pt,
    selected_per100_full=100*pif, selected_per100_clin=100*pic, diff_selected=100*(pif-pic),
    diff_selected_lo=qv(imp_dsel[,jj],.025), diff_selected_hi=qv(imp_dsel[,jj],.975),
    mets_captured_per100_full=100*pif*cff, mets_captured_per100_clin=100*pic*cfc, diff_mets=100*(pif*cff-pic*cfc),
    diff_mets_lo=qv(imp_dmet[,jj],.025), diff_mets_hi=qv(imp_dmet[,jj],.975),
    delta_nb=dnb[["120"]]$full[j]-dnb[["120"]]$clin[j],
    delta_nb_lo=qv(dboot[["120"]]$diff[,j],.025), delta_nb_hi=qv(dboot[["120"]]$diff[,j],.975)) }))
imp_jhu <- do.call(rbind, lapply(seq_along(IMPL), function(jj){ pt<-IMPL[jj]
  cf<-jhu_cc[[paste0("cif_full_",T10)]]; cc<-jhu_cc[[paste0("cif_clin_",T10)]]; w<-jhu_cc$w; Wt<-sum(w)
  sf<-cf>=pt; sc<-cc>=pt; pif<-sum(w[sf])/Wt; pic<-sum(w[sc])/Wt
  cff<-if(any(sf)) cif_pos_fast(jhu_cc$time,jhu_cc$evt,w,sf,T10) else 0
  cfc<-if(any(sc)) cif_pos_fast(jhu_cc$time,jhu_cc$evt,w,sc,T10) else 0
  j<-which.min(abs(THR-pt))
  data.frame(cohort="JHU", basis="conditional internal", horizon=T10, threshold=pt,
    selected_per100_full=100*pif, selected_per100_clin=100*pic, diff_selected=100*(pif-pic),
    diff_selected_lo=qv(jhu_imp_dsel[,jj],.025), diff_selected_hi=qv(jhu_imp_dsel[,jj],.975),
    mets_captured_per100_full=100*pif*cff, mets_captured_per100_clin=100*pic*cfc, diff_mets=100*(pif*cff-pic*cfc),
    diff_mets_lo=qv(jhu_imp_dmet[,jj],.025), diff_mets_hi=qv(jhu_imp_dmet[,jj],.975),
    delta_nb=jnb_corr[["120"]]$full[j]-jnb_corr[["120"]]$clin[j], delta_nb_lo=NA, delta_nb_hi=NA) }))
imp10 <- rbind(imp_jhu, imp_durham)
W(imp10, "clinical_impact_10y.csv")

# frozen 10-year risk-band reclassification transitions (raw analytic counts) from
# the clinical vs clinical+Met-Score CIF predictions; bands <10%, 10-20%, >20%.
# Source of the parallel-flow panels c and d. No refit, rescale, or recalibration.
band3 <- function(x) cut(x, breaks = c(-Inf, 0.10, 0.20, Inf), labels = c("<10%","10-20%",">20%"))
tmat <- function(cc, cf) matrix(as.integer(table(band3(cc), band3(cf))), 3, 3,
          dimnames = list(clinical = c("<10%","10-20%",">20%"), full = c("<10%","10-20%",">20%")))
jhu_tm <- tmat(jhu_cc[[paste0("cif_clin_",T10)]], jhu_cc[[paste0("cif_full_",T10)]])
dur_tm <- tmat(dur_cc[[paste0("cif_clin_",T10)]], dur_cc[[paste0("cif_full_",T10)]])
tm_df <- function(m, coh) data.frame(cohort = coh,
           clinical_band = rep(rownames(m), 3), full_band = rep(colnames(m), each = 3), n = as.integer(m))
W(rbind(tm_df(jhu_tm, "JHU"), tm_df(dur_tm, "Durham")), "risk_band_transitions_10y.csv")

###############################################################################
# Figures are produced separately by code/DCA/figS1_clinical_utility_candidate.py,
# which reads the CSVs written above (dca_curves.csv, delta_net_benefit.csv,
# risk_band_transitions_10y.csv, capras_rule_utility.csv). This script writes
# only numerical results; it does not draw figures.
###############################################################################

# warning ledger written last so it captures warnings from the entire run
warn_df <- if (length(.warn$log)) do.call(rbind, .warn$log) else data.frame(context=character(), message=character())
warn_tab <- if (nrow(warn_df)) as.data.frame(table(warn_df$context, warn_df$message)) else
  data.frame(Var1="none", Var2="no warnings captured", Freq=0)
names(warn_tab) <- c("context","message","count")
W(warn_tab[warn_tab$count>0 | warn_tab$context=="none", ], "warning_disposition.csv")

cat(sprintf("wrote outputs to %s ; frozen reproduction err=%.2e ; math=%d/%d perm=%d/%d\n",
            OUTD, repro_err, sum(math_checks$pass), nrow(math_checks), sum(perm_checks$pass), nrow(perm_checks)))
cat("=== DONE ===\n")

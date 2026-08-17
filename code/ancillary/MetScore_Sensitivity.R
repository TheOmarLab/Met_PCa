############################################################################
# MetScore_Sensitivity.R -- Figure S6 aggregate producer (six robustness panels).
# Aggregate-only; writes identifier-free CSVs under outs/FigureS6/ and no figure.
# JHU is the Ross two-phase case-cohort: baseline/interaction models use
# survival::cch (Lin-Ying, cohort.size=745, robust); follow-up and salvage
# sensitivity use a Barlow start-stop case-cohort implementation (subcohort
# person-time at weight 745/265, outside cases at their failure time), validated
# against the Lin-Ying result. Durham is the complete external cohort with robust
# cause-specific Cox. Met-Score enters High-vs-Low (frozen class) in panels a/b and
# per cohort SD (association scale only) in panels c/d/e; the frozen probabilities,
# class, and threshold are never refit or rethresholded.
#   a follow-up robustness (full / 5y / 10y administrative truncation)
#   b 12-month salvage landmark (landmark, + salvage by 12 months)
#   c Met-Score x race interaction (Durham): White/Black adjusted Met-Score HR per SD
#   d multivariable covariate forest (GG, log2PSA, pT, Met-Score/SD)
#   e Met-Score HR under common-clinical vs + CCP gene-set adjustment
#   f alternative-cutoff sensitivity: the locked Youden threshold vs three
#     development-derived cutoffs (median / Sens90 / Spec90) applied unchanged
#     to JHU and Durham (frozen classifier; never refit or rethresholded). The
#     Met-Score vs CCP gene-set correlation is retained as a numerical output
#     (panelF_ccp_correlation.csv) but is no longer a manuscript panel.
############################################################################
rm(list = ls())
suppressWarnings(suppressMessages({ library(survival); library(readxl) }))

ROOT <- local({ e <- Sys.getenv("MET_PCA_ROOT", ""); if (nzchar(e) && dir.exists(e)) e else normalizePath(".") })
OUT <- file.path(ROOT, "outs", "FigureS6"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ALPHA  <- 265 / 745
B_BOOT <- as.integer(Sys.getenv("METPCA_SENSITIVITY_BOOTSTRAPS", "2000"))
SEED_F <- 20260814L
HORIZONS <- c(60L, 120L)

.warn <- new.env(); .warn$log <- list()
cap <- function(ctx, expr) withCallingHandlers(expr, warning = function(w) {
  .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = ctx, message = conditionMessage(w),
                                                     stringsAsFactors = FALSE)
  invokeRestart("muffleWarning") })

gg_from <- function(gs, pr) { g <- rep(NA_character_, length(gs))
  g[gs <= 6] <- "GG1"; g[gs == 7 & pr == 3] <- "GG2"; g[gs == 7 & pr == 4] <- "GG3"
  g[gs == 8] <- "GG4"; g[gs >= 9] <- "GG5"; g }
pt_collapse <- function(x) { x <- toupper(trimws(as.character(x))); o <- rep(NA_character_, length(x))
  o[grepl("^T2", x)] <- "T2"; o[grepl("^T3", x)] <- "T3"; o[grepl("^T4", x)] <- "T4"; o }
# cause-specific competing clock: metastasis = cause 1; death before metastasis is
# censored at its own time; analysis time is metastasis time otherwise.
clock <- function(met, mt, dead, dt) {
  st <- ifelse(met == 1L, 1L, ifelse(dead == 1L & met == 0L & is.finite(dt) & dt <= mt, 2L, 0L))
  list(status1 = as.integer(st == 1L), time = ifelse(st == 2L, dt, mt))
}

## ---- load cohorts --------------------------------------------------------
je <- new.env(); load(file.path(ROOT, "outs", "coxdata.rda"), envir = je); J0 <- get("CoxData_jhu", envir = je)
de <- new.env(); load(file.path(ROOT, "output", "Durham", "durham_metscore_batchcorrected.rda"), envir = de)
D0 <- get("clin_valid", envir = de)

jc <- clock(as.integer(J0$met), as.numeric(J0$met_time), as.integer(J0$os), as.numeric(J0$os_time))
cchv <- as.character(J0[["post_rp_patients_cchdef"]])
JHU <- data.frame(
  time = jc$time, ev = jc$status1,
  group = factor(ifelse(as.character(J0$MetScoreClass) == "High risk", "High", "Low"), levels = c("Low", "High")),
  ms_z = as.numeric(scale(as.numeric(J0[["Met-Score prob"]]))),
  GG = relevel(factor(gg_from(as.numeric(as.character(J0[["Pathological GS"]])), as.numeric(as.character(J0$pathgs_p)))), ref = "GG2"),
  log2PSA = log2(as.numeric(J0$preop_psa) + 1),
  pT = relevel(factor(pt_collapse(J0$pstage)), ref = "T2"),
  margin = as.integer(J0$sm), node = as.integer(J0$lni),
  insub = as.integer(cchv %in% c("Sub-cohort cases", "Sub-cohort controls")),
  w = ifelse(cchv == "Sub-cohort controls", 1 / ALPHA, 1),
  id = seq_len(nrow(J0)),
  race = as.character(J0$race.x), celfile = as.character(J0$celfile_name),
  rt_s = as.integer(J0$rt_s), adt_s = as.integer(J0$adt_s),
  rt_time = as.numeric(J0$rt_time), adt_time = as.numeric(J0$adt_time),
  met_time = as.numeric(J0$met_time),
  stringsAsFactors = FALSE)
JHU$race3 <- NA_character_
JHU$race3[JHU$race == "Caucasian"] <- "White"
JHU$race3[JHU$race == "African American"] <- "Black"
JHU$race3[JHU$race == "Other"] <- "Other"

dc <- clock(as.integer(D0$mets), as.numeric(D0$surgmets), as.integer(D0$dead), as.numeric(D0$limbo))
DUR <- data.frame(
  time = dc$time, ev = dc$status1,
  group = factor(ifelse(as.character(D0$MetScoreClass) == "High risk", "High", "Low"), levels = c("Low", "High")),
  ms_z = as.numeric(scale(as.numeric(D0$MetScore_prob))),
  GG = relevel(factor(gg_from(as.numeric(as.character(D0$PathGleason)), as.numeric(as.character(D0$pogl1)))), ref = "GG2"),
  log2PSA = log2(as.numeric(D0$psapresurg) + 1),
  pT = relevel(factor(pt_collapse(D0$stg)), ref = "T2"),
  margin = as.integer(D0$positivesurgicalmargins),
  # Durham node codebook: 0=negative, 1=positive, 2=unknown, 3=not done -> NA
  node = ifelse(as.numeric(D0$lymphnodeinvolvement) == 1, 1L,
                ifelse(as.numeric(D0$lymphnodeinvolvement) == 0, 0L, NA_integer_)),
  sample_id = as.character(D0$sample_id),
  salvagerx = as.character(D0$salvagerx), surgsalvagerx = as.numeric(D0$surgsalvagerx),
  stringsAsFactors = FALSE)

## ---- CCP gene-set score from expression (documented list + coverage) ------
CCP <- c("FOXM1","ASPM","TK1","PRC1","CDC20","BUB1B","PBK","DTL","CDKN3","RRM2","ASF1B",
         "CEP55","CDC2","DLGAP5","C18orf24","RAD51","KIF11","BIRC5","RAD54L","CENPM",
         "KIAA0101","KIF20A","PTTG1","CDCA8","NUSAP1","PLK1","CDCA3","ORC6L","CENPF","TOP2A","MCM10")
ALIAS <- list(CDC2 = c("CDC2","CDK1"), KIAA0101 = c("KIAA0101","PCLAF"),
              ORC6L = c("ORC6L","ORC6"), C18orf24 = c("C18orf24","RSL24D1","RSRC1"))
resolve_ccp <- function(rn) {
  used <- character(0); canon <- character(0)
  for (g in CCP) { cand <- if (!is.null(ALIAS[[g]])) ALIAS[[g]] else g
    hit <- cand[cand %in% rn][1]
    if (!is.na(hit)) { used <- c(used, hit); canon <- c(canon, g) } }
  list(used = used, canon = canon) }
ccp_score <- function(expr) {  # per-gene z across samples, mean, then cohort z
  z <- t(apply(expr, 1, function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }))
  as.numeric(scale(colMeans(z, na.rm = TRUE))) }

load(file.path(ROOT, "data", "Dataset7.rda"))            # Dataset7$expr (gene x sample, JHU)
jx <- Dataset7$expr
jr <- resolve_ccp(rownames(jx)); jccp <- ccp_score(jx[jr$used, , drop = FALSE]); names(jccp) <- colnames(jx)
JHU$ccp <- as.numeric(scale(jccp[JHU$celfile]))

DUR_XLSX <- file.path(ROOT, "data", "Durham_cohort_and_GRID_cohort", "Durham_cohort_011526.xlsx")
dz <- suppressMessages(read_excel(DUR_XLSX, sheet = "eset_gene_filtered"))
dsym <- dz$Symbol; dm <- as.matrix(dz[, !(colnames(dz) %in% c("ENTREZID", "Symbol"))]); rownames(dm) <- dsym
dr <- resolve_ccp(rownames(dm)); dccp <- ccp_score(dm[dr$used, , drop = FALSE]); names(dccp) <- colnames(dm)
DUR$ccp <- as.numeric(scale(dccp[DUR$sample_id]))

ccp_cov <- data.frame(
  cohort = c(rep("JHU", length(CCP)), rep("Durham", length(CCP))),
  canonical_gene = c(CCP, CCP),
  matched_symbol = c(ifelse(CCP %in% jr$canon, jr$used[match(CCP, jr$canon)], NA_character_),
                     ifelse(CCP %in% dr$canon, dr$used[match(CCP, dr$canon)], NA_character_)),
  matched = c(CCP %in% jr$canon, CCP %in% dr$canon), stringsAsFactors = FALSE)
write.csv(ccp_cov, file.path(OUT, "panelE_ccp_gene_coverage.csv"), row.names = FALSE)

## ---- design + event ledger (fail closed) ---------------------------------
stopifnot(nrow(JHU) == 239L, sum(JHU$ev) == 93L,
          sum(cchv == "Sub-cohort cases") == 28L, sum(cchv == "Sub-cohort controls") == 146L, sum(cchv == "cases") == 65L,
          nrow(DUR) == 555L, sum(DUR$ev) == 40L)
ledger <- data.frame(
  cohort = c("JHU", "Durham"), n = c(nrow(JHU), nrow(DUR)),
  metastasis_events = c(sum(JHU$ev), sum(DUR$ev)),
  subcohort_cases = c(28L, NA), subcohort_controls = c(146L, NA), outside_cases = c(65L, NA),
  design = c("Ross two-phase case-cohort (cohort.size=745, alpha=265/745)", "complete external cohort"),
  ccp_coverage = c(sprintf("%d/31", length(jr$used)), sprintf("%d/31", length(dr$used))),
  stringsAsFactors = FALSE)
write.csv(ledger, file.path(OUT, "sample_event_design_ledger.csv"), row.names = FALSE)

## ---- estimators ----------------------------------------------------------
mvchecks <- list()
addchk <- function(name, got, exp, tol) {
  mvchecks[[length(mvchecks) + 1L]] <<- data.frame(check = name, value = signif(got, 7),
    reference = signif(exp, 7), abs_diff = signif(abs(got - exp), 3),
    pass = isTRUE(abs(got - exp) <= tol), stringsAsFactors = FALSE) }

# JHU Lin-Ying cch high-vs-low (reference for panel a validation)
cch_hilo <- function(d) {
  f <- cap("cch_hilo", survival::cch(Surv(time, ev) ~ group, data = d, subcoh = ~ insub, id = ~ id,
                                     cohort.size = 745L, method = "LinYing", robust = TRUE))
  cf <- summary(f)$coefficients; v <- unname(cf[1, "Value"]); se <- unname(cf[1, "SE"])
  c(hr = exp(v), b = v, se = se, p = unname(cf[1, "p"])) }

# Barlow start-stop case-cohort high-vs-low, administrative truncation at T (Inf = full)
barlow_hilo <- function(d, T = Inf) {
  eps <- 0.5; rows <- list()
  for (i in seq_len(nrow(d))) {
    ti <- d$time[i]; ev <- d$ev[i]; ins <- d$insub[i]; ce <- ev == 1L & ti <= T; tt <- min(ti, T)
    if (ins == 1L) {
      if (ce) { if (tt - eps > 0) rows[[length(rows) + 1L]] <- data.frame(id = d$id[i], s = 0, e = tt - eps, z = 0L, group = d$group[i], w = 1 / ALPHA)
                rows[[length(rows) + 1L]] <- data.frame(id = d$id[i], s = max(tt - eps, 0), e = tt, z = 1L, group = d$group[i], w = 1)
      } else rows[[length(rows) + 1L]] <- data.frame(id = d$id[i], s = 0, e = tt, z = 0L, group = d$group[i], w = 1 / ALPHA)
    } else if (ev == 1L && ti <= T) {
      rows[[length(rows) + 1L]] <- data.frame(id = d$id[i], s = max(ti - eps, 0), e = ti, z = 1L, group = d$group[i], w = 1)
    }
  }
  ss <- do.call(rbind, rows)
  f <- cap("barlow", survival::coxph(Surv(s, e, z) ~ group, data = ss, weights = w, cluster = id,
                                     robust = TRUE, ties = "breslow", timefix = FALSE))
  cf <- summary(f)$coefficients; b <- cf["groupHigh", "coef"]; se <- cf["groupHigh", "robust se"]
  n_ev <- sum(ss$z); n_at <- length(unique(ss$id))
  c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = cf["groupHigh", "Pr(>|z|)"], events = n_ev, n = n_at) }

# Durham complete-cohort robust Cox high-vs-low, administrative truncation at T
dur_hilo <- function(d, T = Inf) {
  ev <- as.integer(d$ev == 1L & d$time <= T); tt <- pmin(d$time, T)
  f <- cap("dur_cox", survival::coxph(Surv(tt, ev) ~ group, data = d, robust = TRUE))
  cf <- summary(f)$coefficients; b <- cf["groupHigh", "coef"]; se <- cf["groupHigh", "robust se"]
  c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = cf["groupHigh", "Pr(>|z|)"],
    events = sum(ev), n = nrow(d)) }

## ---- panel a: follow-up robustness --------------------------------------
Ja <- JHU[is.finite(JHU$time) & JHU$time > 0, ]
Da <- DUR[is.finite(DUR$time) & DUR$time > 0, ]
lin <- cch_hilo(Ja); pren <- {
  f <- cap("cch_prentice", survival::cch(Surv(time, ev) ~ group, data = Ja, subcoh = ~ insub, id = ~ id,
                                         cohort.size = 745L, method = "Prentice")); exp(coef(f)[1]) }
b_full <- barlow_hilo(Ja, Inf)
addchk("panelA_startstop_full_vs_LinYing", unname(b_full["hr"]), unname(lin["hr"]), 0.6)
addchk("panelA_startstop_full_vs_Prentice", unname(b_full["hr"]), unname(pren), 0.2)
SCEN <- c(Inf, 60, 120); SLAB <- c("Full follow-up", "Administrative truncation at 5 years", "Administrative truncation at 10 years")
panelA <- list()
for (k in seq_along(SCEN)) { T <- SCEN[k]
  bj <- barlow_hilo(Ja, T); dj <- dur_hilo(Da, T)
  panelA[[length(panelA) + 1L]] <- data.frame(cohort = "JHU", scenario = SLAB[k], horizon_months = ifelse(is.infinite(T), NA, T),
    n = bj["n"], events = bj["events"], hr = round(bj["hr"], 5), ci_lo = round(bj["lo"], 5), ci_hi = round(bj["hi"], 5),
    p = signif(bj["p"], 5), estimator = "Barlow start-stop case-cohort (admin-censored)", stringsAsFactors = FALSE)
  panelA[[length(panelA) + 1L]] <- data.frame(cohort = "Durham", scenario = SLAB[k], horizon_months = ifelse(is.infinite(T), NA, T),
    n = dj["n"], events = dj["events"], hr = round(dj["hr"], 5), ci_lo = round(dj["lo"], 5), ci_hi = round(dj["hi"], 5),
    p = signif(dj["p"], 5), estimator = "complete-cohort robust Cox (admin-censored)", stringsAsFactors = FALSE) }
panelA <- do.call(rbind, panelA)
write.csv(panelA, file.path(OUT, "panelA_followup_robustness.csv"), row.names = FALSE)

## ---- panel b: 12-month salvage landmark ---------------------------------
LM <- 12
# JHU salvage-by-12mo from RT/ADT timestamps; unavailable if timing cannot support it
Ja$salv_time <- suppressWarnings(pmin(ifelse(Ja$rt_s %in% 1, Ja$rt_time, Inf),
                                      ifelse(Ja$adt_s %in% 1, Ja$adt_time, Inf), na.rm = TRUE))
Ja$salv_by_lm <- as.integer(is.finite(Ja$salv_time) & Ja$salv_time <= LM)
jlm <- Ja[Ja$time > LM, ]; jlm$time <- jlm$time - LM        # re-origin at landmark
# feasibility: salvage-by-12mo must precede metastasis for exposed metastatic cases
jlm_bad <- sum(jlm$salv_by_lm == 1L & jlm$ev == 1L & jlm$salv_time >= (jlm$met_time))
jlm_exposed_events <- sum(jlm$salv_by_lm == 1L & jlm$ev == 1L)
jhu_salv_ok <- jlm_exposed_events >= 5L && jlm_bad == 0L
panelB <- list()
bjl <- barlow_hilo(jlm, Inf)
panelB[[1]] <- data.frame(cohort = "JHU", model = "Landmark 12 months", n = bjl["n"], events = bjl["events"],
  hr = round(bjl["hr"], 5), ci_lo = round(bjl["lo"], 5), ci_hi = round(bjl["hi"], 5), p = signif(bjl["p"], 5),
  estimator = "Barlow start-stop case-cohort, 12-month landmark", available = TRUE, stringsAsFactors = FALSE)
if (jhu_salv_ok) {
  # start-stop high-vs-low with salvage_by_lm as an additional covariate
  eps <- 0.5; rows <- list()
  for (i in seq_len(nrow(jlm))) { ti <- jlm$time[i]; ev <- jlm$ev[i]; ins <- jlm$insub[i]; sv <- jlm$salv_by_lm[i]
    if (ins == 1L) { if (ev == 1L) { if (ti - eps > 0) rows[[length(rows)+1L]] <- data.frame(id=jlm$id[i],s=0,e=ti-eps,z=0L,group=jlm$group[i],salv=sv,w=1/ALPHA)
        rows[[length(rows)+1L]] <- data.frame(id=jlm$id[i],s=max(ti-eps,0),e=ti,z=1L,group=jlm$group[i],salv=sv,w=1)
      } else rows[[length(rows)+1L]] <- data.frame(id=jlm$id[i],s=0,e=ti,z=0L,group=jlm$group[i],salv=sv,w=1/ALPHA)
    } else if (ev == 1L) rows[[length(rows)+1L]] <- data.frame(id=jlm$id[i],s=max(ti-eps,0),e=ti,z=1L,group=jlm$group[i],salv=sv,w=1) }
  ssl <- do.call(rbind, rows)
  fjs <- cap("barlow_salv", survival::coxph(Surv(s,e,z) ~ group + salv, data = ssl, weights = w, cluster = id, robust = TRUE, ties = "breslow", timefix = FALSE))
  cf <- summary(fjs)$coefficients; b <- cf["groupHigh","coef"]; se <- cf["groupHigh","robust se"]
  panelB[[2]] <- data.frame(cohort = "JHU", model = "+ salvage by 12 months", n = length(unique(ssl$id)), events = sum(ssl$z),
    hr = round(exp(b),5), ci_lo = round(exp(b-1.96*se),5), ci_hi = round(exp(b+1.96*se),5), p = signif(cf["groupHigh","Pr(>|z|)"],5),
    estimator = "Barlow start-stop case-cohort + salvage-by-12mo", available = TRUE, stringsAsFactors = FALSE)
} else {
  panelB[[2]] <- data.frame(cohort = "JHU", model = "+ salvage by 12 months", n = NA, events = NA,
    hr = NA, ci_lo = NA, ci_hi = NA, p = NA,
    estimator = "unavailable: JHU salvage timestamps do not validly establish salvage-by-12mo", available = FALSE, stringsAsFactors = FALSE)
}
# Durham salvage-by-12mo from surgsalvagerx timestamp
DUR$salv_by_lm <- as.integer(DUR$salvagerx != "0" & is.finite(DUR$surgsalvagerx) & DUR$surgsalvagerx <= LM)
dlm <- DUR[DUR$time > LM, ]; dlm$time <- dlm$time - LM
d_base <- survival::coxph(Surv(time, ev) ~ group, data = dlm, robust = TRUE)
cf <- summary(d_base)$coefficients; b <- cf["groupHigh","coef"]; se <- cf["groupHigh","robust se"]
panelB[[3]] <- data.frame(cohort = "Durham", model = "Landmark 12 months", n = nrow(dlm), events = sum(dlm$ev),
  hr = round(exp(b),5), ci_lo = round(exp(b-1.96*se),5), ci_hi = round(exp(b+1.96*se),5), p = signif(cf["groupHigh","Pr(>|z|)"],5),
  estimator = "complete-cohort robust Cox, 12-month landmark", available = TRUE, stringsAsFactors = FALSE)
d_salv <- cap("dur_salv", survival::coxph(Surv(time, ev) ~ group + salv_by_lm, data = dlm, robust = TRUE))
cf <- summary(d_salv)$coefficients; b <- cf["groupHigh","coef"]; se <- cf["groupHigh","robust se"]
panelB[[4]] <- data.frame(cohort = "Durham", model = "+ salvage by 12 months", n = nrow(dlm), events = sum(dlm$ev),
  hr = round(exp(b),5), ci_lo = round(exp(b-1.96*se),5), ci_hi = round(exp(b+1.96*se),5), p = signif(cf["groupHigh","Pr(>|z|)"],5),
  estimator = "complete-cohort robust Cox + salvage-by-12mo", available = TRUE, stringsAsFactors = FALSE)
panelB <- do.call(rbind, panelB)
write.csv(panelB, file.path(OUT, "panelB_salvage_landmark.csv"), row.names = FALSE)

## ---- panel c: Met-Score x race interaction (Durham) ---------------------
# Race codes per Durham_cohort_clinical_data_022526.xlsx, sheet clin_codebook,
# row 6: 1=Black, 3=White, 10=Asian/Pacific Islander, 13=American Indian/Alaska
# Native. Formal cause-specific inference is limited to Black and White (all 40
# metastases); race-specific Met-Score HR/SD are model contrasts from one robust
# Cox model. Frozen probabilities/classes are used as-is.
RACE_MAP <- c("1"="Black","3"="White","10"="Asian/Pacific Islander","13"="American Indian/Alaska Native")
ETH_MAP  <- c("1"="Hispanic/Latino","2"="Not Hispanic/Latino")
gg_ord <- function(gs, pr) { o <- rep(NA_real_, length(gs)); o[gs<=6]<-1; o[gs==7&pr==3]<-2
  o[gs==7&pr==4]<-3; o[gs==8]<-4; o[gs>=9]<-5; o }
Rc <- data.frame(time = DUR$time, ev = DUR$ev, ms_z = DUR$ms_z, log2PSA = DUR$log2PSA,
  gg = gg_ord(as.numeric(as.character(D0$PathGleason)), as.numeric(as.character(D0$pogl1))),
  pT34 = ifelse(DUR$pT %in% c("T3","T4"), 1L, ifelse(DUR$pT == "T2", 0L, NA_integer_)),
  race = unname(RACE_MAP[as.character(D0$race)]), eth = unname(ETH_MAP[as.character(D0$ethnicity)]),
  stringsAsFactors = FALSE)
# race / ethnicity ledger (formal inference only where events exist)
raceL <- do.call(rbind, lapply(c("Black","White","Asian/Pacific Islander","American Indian/Alaska Native"), function(r) {
  s <- Rc[which(Rc$race == r), ]; data.frame(group_type="race", group=r, n=nrow(s), events=sum(s$ev),
    inference=ifelse(sum(s$ev)>0, "formal", "counts only (no events)"), stringsAsFactors=FALSE) }))
ethL <- do.call(rbind, lapply(c("Hispanic/Latino","Not Hispanic/Latino"), function(x) { s <- Rc[which(Rc$eth == x), ]
  data.frame(group_type="ethnicity", group=x, n=nrow(s), events=sum(s$ev), inference="counts only", stringsAsFactors=FALSE) }))
write.csv(rbind(raceL, ethL), file.path(OUT, "panelC_race_ledger.csv"), row.names = FALSE)
# one robust cause-specific Cox on Black + White with a Met-Score x race interaction
Rbw <- Rc[Rc$race %in% c("Black","White") & is.finite(Rc$time) & Rc$time > 0 &
          is.finite(Rc$ms_z) & is.finite(Rc$gg) & is.finite(Rc$log2PSA) & !is.na(Rc$pT34), ]
Rbw$race <- relevel(factor(Rbw$race), ref = "White")
fc <- cap("race_dur", survival::coxph(Surv(time, ev) ~ ms_z * race + gg + log2PSA + pT34, data = Rbw, robust = TRUE))
bC <- coef(fc); VC <- vcov(fc); im <- which(names(bC)=="ms_z"); ix <- grep("ms_z:raceBlack", names(bC))
white_b <- unname(bC[im]); white_se <- unname(sqrt(VC[im,im]))
black_b <- unname(bC[im] + bC[ix]); black_se <- unname(sqrt(VC[im,im] + VC[ix,ix] + 2*VC[im,ix]))
int_b <- unname(bC[ix]); int_se <- unname(sqrt(VC[ix,ix]))
ph_global <- unname(cap("race_zph", survival::cox.zph(fc))$table["GLOBAL","p"])
converged <- all(is.finite(bC)) && all(is.finite(diag(VC)))
mkrow <- function(rc, lab, b, se, n, ev) data.frame(cohort="Durham", race=rc, label=lab, n=n, events=ev,
  hr=round(exp(b),5), ci_lo=round(exp(b-1.96*se),5), ci_hi=round(exp(b+1.96*se),5),
  interaction_ratio=round(exp(int_b),5), interaction_lo=round(exp(int_b-1.96*int_se),5),
  interaction_hi=round(exp(int_b+1.96*int_se),5), interaction_p=signif(2*pnorm(-abs(int_b/int_se)),5),
  global_ph_p=signif(ph_global,5), converged=converged,
  estimator="robust cause-specific Cox, Met-Score x race interaction (White ref.)", stringsAsFactors=FALSE)
panelC <- rbind(mkrow("White","White (ref.)", white_b, white_se, sum(Rbw$race=="White"), sum(Rbw$ev[Rbw$race=="White"])),
                mkrow("Black","Black",        black_b, black_se, sum(Rbw$race=="Black"), sum(Rbw$ev[Rbw$race=="Black"])))
write.csv(panelC, file.path(OUT, "panelC_race_interaction.csv"), row.names = FALSE)

## ---- panel d: multivariable covariate forest ----------------------------
TERMS <- c("GGGG3","GGGG4","GGGG5","log2PSA","pTT3","pTT4","ms_z")
TLAB  <- c(GGGG3="GG3 (4+3)", GGGG4="GG4", GGGG5="GG5", log2PSA="log2(PSA+1)", pTT3="pT3", pTT4="pT4",
           ms_z="Met-Score (per SD)", node="Node +")
mv_rows <- function(fit, robust_cch) {
  cf <- summary(fit)$coefficients; nmv <- rownames(cf)
  bb <- if (robust_cch) cf[, "Value"] else cf[, "coef"]
  se <- if (robust_cch) cf[, "SE"] else cf[, "robust se"]
  out <- list()
  for (k in seq_along(nmv)) {
    t <- nmv[k]; b <- unname(bb[k]); s <- unname(se[k])
    lab <- if (t %in% names(TLAB)) unname(TLAB[[t]]) else t
    out[[length(out)+1L]] <- data.frame(term = t, label = lab,
      hr = round(exp(b),5), ci_lo = round(exp(b-1.96*s),5), ci_hi = round(exp(b+1.96*s),5),
      se = round(s,5), p = signif(2*pnorm(-abs(b/s)),5),
      estimable = is.finite(b) & is.finite(s) & abs(b) < 10 & s < 5, stringsAsFactors = FALSE) }
  do.call(rbind, out) }
# JHU common model (cch)
Jd <- JHU[stats::complete.cases(JHU[, c("time","ev","GG","log2PSA","pT","ms_z")]) & JHU$time > 0, ]
fjd <- cap("mv_jhu", survival::cch(Surv(time, ev) ~ GG + log2PSA + pT + ms_z, data = Jd, subcoh = ~ insub, id = ~ id,
                                   cohort.size = 745L, method = "LinYing", robust = TRUE))
rj <- mv_rows(fjd, TRUE); rj$cohort <- "JHU"; rj$model <- "common"; rj$display <- TRUE
# JHU + node sensitivity: the node-adjusted model exports both the Met-Score
# coefficient after node adjustment (displayed) and the node-status coefficient
# (exported, not displayed as the Met-Score result).
Jd2 <- JHU[stats::complete.cases(JHU[, c("time","ev","GG","log2PSA","pT","ms_z","node")]) & JHU$time > 0, ]
fjd2 <- tryCatch(cap("mv_jhu_node", survival::cch(Surv(time, ev) ~ GG + log2PSA + pT + ms_z + node, data = Jd2,
             subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE)), error = function(e) NULL)
rj_mn <- if (!is.null(fjd2)) { r <- mv_rows(fjd2, TRUE)
  msa <- r[r$term == "ms_z", ]; msa$label <- "Met-Score (+ node adjustment)"; msa$cohort <- "JHU"; msa$model <- "JHU sensitivity (+node)"; msa$display <- TRUE
  ndr <- r[r$term == "node", ]; ndr$label <- "Node status +"; ndr$cohort <- "JHU"; ndr$model <- "JHU sensitivity (+node)"; ndr$display <- FALSE
  rbind(msa, ndr) } else NULL
# Durham common model (robust Cox)
Dd <- DUR[stats::complete.cases(DUR[, c("time","ev","GG","log2PSA","pT","ms_z")]) & DUR$time > 0, ]
fdd <- cap("mv_dur", survival::coxph(Surv(time, ev) ~ GG + log2PSA + pT + ms_z, data = Dd, robust = TRUE))
rd <- mv_rows(fdd, FALSE); rd$cohort <- "Durham"; rd$model <- "common"; rd$display <- TRUE
panelD <- rbind(rj, rj_mn, rd)
panelD$n <- ifelse(panelD$cohort == "JHU" & panelD$model == "common", nrow(Jd),
             ifelse(panelD$cohort == "JHU", nrow(Jd2), nrow(Dd)))
panelD$events <- ifelse(panelD$cohort == "JHU" & panelD$model == "common", sum(Jd$ev),
                 ifelse(panelD$cohort == "JHU", sum(Jd2$ev), sum(Dd$ev)))
panelD <- panelD[, c("cohort","model","term","label","hr","ci_lo","ci_hi","se","p","estimable","display","n","events")]
write.csv(panelD, file.path(OUT, "panelD_multivariable.csv"), row.names = FALSE)

## ---- canonical Main Table 1 (multivariable MFS) from the same result --------
# Manuscript-facing view of the accepted multivariable MFS models, built from the
# in-memory panelD (no refit). Same HR/CI/p/n/events; public-facing labels only.
mt1_cmp <- c(GGGG3 = "GG3 (4+3) vs GG2 (3+4)", GGGG4 = "GG4 vs GG2", GGGG5 = "GG5 vs GG2",
             GGGG1 = "GG1 / GS≤6 vs GG2", log2PSA = "log2(PSA+1), per unit",
             pTT3 = "pT3 vs pT2", pTT4 = "pT4 vs pT2", ms_z = "Met-Score (cohort z; per 1 SD)")
mt1_var <- c(GGGG3 = "Grade Group", GGGG4 = "Grade Group", GGGG5 = "Grade Group", GGGG1 = "Grade Group",
             log2PSA = "log2(PSA+1)", pTT3 = "Pathological pT", pTT4 = "Pathological pT",
             ms_z = "Met-Score", node = "Node status")
mt1_est <- c(JHU = "Lin-Ying case-cohort (cohort.size=745, robust)", Durham = "robust cause-specific Cox")
mt1_adj_common <- "Grade Group + log2(PSA+1) + pathological pT + Met-Score (cohort z per SD); no surgical margin"
mt1_adj_node   <- "common covariates (Grade Group + log2(PSA+1) + pathological pT + Met-Score) retained, plus node status; no surgical margin"
is_node <- panelD$model == "JHU sensitivity (+node)"
mt1 <- data.frame(
  Cohort = ifelse(panelD$cohort == "JHU", "JHU", "Durham VA"),
  Model = ifelse(panelD$model == "common", "Common multivariable MFS", "JHU node-adjusted sensitivity"),
  Endpoint = "MFS",
  Variable = ifelse(is_node & panelD$term == "node", "Node status",
             ifelse(is_node & panelD$term == "ms_z", "Met-Score", unname(mt1_var[panelD$term]))),
  Comparison = ifelse(is_node & panelD$term == "node", "Node-positive vs node-negative",
               ifelse(is_node & panelD$term == "ms_z", "Met-Score after node adjustment", unname(mt1_cmp[panelD$term]))),
  HR    = ifelse(panelD$estimable, round(panelD$hr, 4), NA_real_),
  CI_lo = ifelse(panelD$estimable, round(panelD$ci_lo, 4), NA_real_),
  CI_hi = ifelse(panelD$estimable, round(panelD$ci_hi, 4), NA_real_),
  p     = ifelse(panelD$estimable, panelD$p, NA_real_),
  n = panelD$n, events = panelD$events,
  Status = ifelse(panelD$estimable, "estimated", "non-estimable: separation"),
  Estimator = unname(mt1_est[panelD$cohort]),
  Time_origin = "radical prostatectomy; metastasis event, death before metastasis censored at death time",
  Adjustment_set = ifelse(is_node, mt1_adj_node, mt1_adj_common),
  stringsAsFactors = FALSE)
mt1 <- mt1[, c("Cohort","Model","Endpoint","Variable","Comparison","HR","CI_lo","CI_hi","p",
               "n","events","Status","Estimator","Time_origin","Adjustment_set")]
write.csv(mt1, file.path(dirname(OUT), "MainTable1_multivariable_MFS.csv"), row.names = FALSE, na = "")
cat(sprintf("Wrote outs/MainTable1_multivariable_MFS.csv (%d rows)\n", nrow(mt1)))

## ---- panel e: Met-Score HR under common-clinical vs + CCP ----------------
mse <- function(dsub, robust_cch, add_ccp) {
  f <- if (!add_ccp) Surv(time, ev) ~ ms_z + GG + log2PSA + pT else Surv(time, ev) ~ ms_z + GG + log2PSA + pT + ccp
  if (robust_cch) { fit <- cap("e_cch", survival::cch(f, data = dsub, subcoh = ~ insub, id = ~ id, cohort.size = 745L, method = "LinYing", robust = TRUE))
    cf <- summary(fit)$coefficients; c(hr = exp(cf["ms_z","Value"]), lo = exp(cf["ms_z","Value"] - 1.96*cf["ms_z","SE"]), hi = exp(cf["ms_z","Value"] + 1.96*cf["ms_z","SE"]), p = cf["ms_z","p"]) }
  else { fit <- cap("e_cox", survival::coxph(f, data = dsub, robust = TRUE))
    cf <- summary(fit)$coefficients; c(hr = exp(cf["ms_z","coef"]), lo = exp(cf["ms_z","coef"] - 1.96*cf["ms_z","robust se"]), hi = exp(cf["ms_z","coef"] + 1.96*cf["ms_z","robust se"]), p = cf["ms_z","Pr(>|z|)"]) }
}
Je <- JHU[stats::complete.cases(JHU[, c("time","ev","ms_z","GG","log2PSA","pT","ccp")]) & JHU$time > 0, ]
Dce <- DUR[stats::complete.cases(DUR[, c("time","ev","ms_z","GG","log2PSA","pT","ccp")]) & DUR$time > 0, ]
panelE <- list()
for (cn in c("JHU","Durham")) { dsub <- if (cn == "JHU") Je else Dce; rc <- cn == "JHU"
  for (adj in c(FALSE, TRUE)) { r <- mse(dsub, rc, adj)
    panelE[[length(panelE)+1L]] <- data.frame(cohort = cn, adjustment = ifelse(adj, "+ CCP gene-set score", "Common clinical"),
      n = nrow(dsub), events = sum(dsub$ev), hr = round(r["hr"],5), ci_lo = round(r["lo"],5), ci_hi = round(r["hi"],5),
      p = signif(r["p"],5), estimator = ifelse(rc, "Lin-Ying case-cohort cch", "complete-cohort robust Cox"), stringsAsFactors = FALSE) } }
panelE <- do.call(rbind, panelE)
write.csv(panelE, file.path(OUT, "panelE_ccp_adjustment.csv"), row.names = FALSE)

## ---- Met-Score vs CCP gene-set score correlation (numerical output only) --
# Retained per revision as a numerical output; no longer a manuscript panel.
# JHU: design-weighted rank correlation + conditional design-stratified bootstrap.
# Durham: complete-cohort Spearman + patient bootstrap. B fixed-seed replicates.
wt_spearman <- function(x, y, w) { rx <- rank(x); ry <- rank(y)
  mx <- sum(w * rx) / sum(w); my <- sum(w * ry) / sum(w)
  cxy <- sum(w * (rx - mx) * (ry - my)); cxx <- sum(w * (rx - mx)^2); cyy <- sum(w * (ry - my)^2)
  cxy / sqrt(cxx * cyy) }
Jf <- JHU[is.finite(JHU$ms_z) & is.finite(JHU$ccp), ]
Df <- DUR[is.finite(DUR$ms_z) & is.finite(DUR$ccp), ]
rho_j <- wt_spearman(Jf$ms_z, Jf$ccp, Jf$w)
rho_d <- suppressWarnings(cor(Df$ms_z, Df$ccp, method = "spearman"))
set.seed(SEED_F)
strat <- split(seq_len(nrow(Jf)), as.character(J0[["post_rp_patients_cchdef"]])[match(Jf$id, JHU$id)])
bj <- numeric(0); jf_ok <- 0L; jf_fail <- 0L
for (b in seq_len(B_BOOT)) { ix <- unlist(lapply(strat, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
  v <- tryCatch(wt_spearman(Jf$ms_z[ix], Jf$ccp[ix], Jf$w[ix]), error = function(e) NA_real_)
  if (is.finite(v)) { bj <- c(bj, v); jf_ok <- jf_ok + 1L } else jf_fail <- jf_fail + 1L }
set.seed(SEED_F)
bd <- numeric(0); df_ok <- 0L; df_fail <- 0L; nd <- nrow(Df)
for (b in seq_len(B_BOOT)) { ix <- sample.int(nd, nd, replace = TRUE)
  v <- tryCatch(suppressWarnings(cor(Df$ms_z[ix], Df$ccp[ix], method = "spearman")), error = function(e) NA_real_)
  if (is.finite(v)) { bd <- c(bd, v); df_ok <- df_ok + 1L } else df_fail <- df_fail + 1L }
qj <- quantile(bj, c(.025, .975)); qd <- quantile(bd, c(.025, .975))
panelF <- data.frame(
  cohort = c("JHU", "Durham"), comparator = "CCP gene-set score",
  n = c(nrow(Jf), nrow(Df)), rho = round(c(rho_j, rho_d), 5),
  ci_lo = round(c(qj[1], qd[1]), 5), ci_hi = round(c(qj[2], qd[2]), 5),
  estimator = c("design-weighted Spearman", "Spearman"),
  resampling = c("conditional design-stratified bootstrap", "patient bootstrap"),
  boot_seed = SEED_F, boot_used = c(jf_ok, df_ok), boot_failed = c(jf_fail, df_fail), stringsAsFactors = FALSE)
write.csv(panelF, file.path(OUT, "panelF_ccp_correlation.csv"), row.names = FALSE)

## ---- panel f: alternative-cutoff sensitivity ----------------------------
# Four cutoffs on the probability scale: the locked Youden threshold (primary,
# from metadata) and three cutoffs derived from the locked classifier's
# development predictions only (never validation outcomes). All four are applied
# UNCHANGED to JHU (case-cohort, phase-two weights, conditional design bootstrap)
# and Durham (complete cohort, patient bootstrap) at the same competing-risk
# event, 5/10-year horizons, and censoring as the accepted Figure 2 analysis.
source(file.path(ROOT, "code", "utils", "locked_metscore.R"))
LM_MODEL   <- load_locked_metscore(config_dir = file.path(ROOT, "config"))
THR_YOUDEN <- LM_MODEL$threshold
me <- new.env(); load(file.path(ROOT, "outs", "MetastasisData_JHUOut.rda"), envir = me)
devX <- t(as.matrix(me$trainMat)[LM_MODEL$feature_names, , drop = FALSE])   # samples x features
devp <- plogis(as.numeric(LM_MODEL$intercept + devX %*% LM_MODEL$beta[LM_MODEL$feature_names]))
devy <- as.integer(as.character(me$trainGroup) == "Mets")
dev_se <- function(c0) sum(devp >= c0 & devy == 1L) / sum(devy == 1L)
dev_sp <- function(c0) sum(devp <  c0 & devy == 0L) / sum(devy == 0L)
# deterministic tie handling: candidate cutoffs are the sorted unique development
# probabilities, positive if prob >= cutoff; Se/Sp are step functions of the cutoff.
cand <- sort(unique(devp)); se_c <- vapply(cand, dev_se, 0); sp_c <- vapply(cand, dev_sp, 0)
cut_med    <- as.numeric(median(devp))          # median development probability
cut_sens90 <- max(cand[se_c >= 0.90])           # largest dev cutoff with Se>=0.90 (max Sp s.t. Se>=0.90)
cut_spec90 <- min(cand[sp_c >= 0.90])           # smallest dev cutoff with Sp>=0.90 (max Se s.t. Sp>=0.90)
you_resub  <- cand[which.max(se_c + sp_c - 1)]  # resubstitution Youden (validation only)
CUTS <- data.frame(
  cutoff_name  = c("Youden (locked)", "Median", "Sens90", "Spec90"),
  cutoff_value = c(THR_YOUDEN, cut_med, cut_sens90, cut_spec90),
  derivation   = c("locked classifier metadata (primary cutoff)",
                   "median development probability",
                   "max specificity s.t. development sensitivity >= 0.90 (largest qualifying dev-probability cutoff)",
                   "max sensitivity s.t. development specificity >= 0.90 (smallest qualifying dev-probability cutoff)"),
  is_primary   = c(TRUE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE)
CUTS$dev_sensitivity <- round(vapply(CUTS$cutoff_value, dev_se, 0), 5)
CUTS$dev_specificity <- round(vapply(CUTS$cutoff_value, dev_sp, 0), 5)
CUTS$dev_n <- length(devy); CUTS$dev_events <- sum(devy); CUTS$dev_prevalence <- round(mean(devy), 5)
CUTS$cutoff_value_out <- round(CUTS$cutoff_value, 6)
write.csv(CUTS[, c("cutoff_name","cutoff_value_out","derivation","is_primary",
                   "dev_n","dev_events","dev_prevalence","dev_sensitivity","dev_specificity")],
          file.path(OUT, "panelF_cutoff_definitions.csv"), row.names = FALSE)

# competing-event table: locked probability marker, cause-1 metastasis, def-2 controls
bet <- function(score, mt, met, dead, dt, cchdef, n_exp, m_exp, cd_exp) {
  score <- as.numeric(score); mt <- as.numeric(mt); met <- as.integer(met)
  dead <- if (is.null(dead)) rep(0L, length(score)) else as.integer(dead)
  dt   <- if (is.null(dt)) rep(NA_real_, length(score)) else as.numeric(dt)
  keep <- is.finite(score) & is.finite(mt) & mt > 0 & !is.na(met)
  score<-score[keep]; mt<-mt[keep]; met<-met[keep]; dead<-dead[keep]; dt<-dt[keep]
  status <- ifelse(met==1L, 1L, ifelse(dead==1L & !is.na(dt) & dt<=mt, 2L, 0L))
  d <- data.frame(score=score, atime=as.numeric(ifelse(status==2L, dt, mt)),
                  status=as.integer(status), stringsAsFactors=FALSE)
  if (!is.null(cchdef)) { cv <- as.character(cchdef)[keep]
    d$cch <- cv; d$w <- ifelse(cv=="Sub-cohort controls", 1/ALPHA, 1)
    d$insub <- as.integer(cv %in% c("Sub-cohort cases","Sub-cohort controls"))
  } else { d$cch <- "cohort"; d$w <- 1; d$insub <- 1L }
  d$id <- seq_len(nrow(d))
  stopifnot(nrow(d)==n_exp, sum(d$status==1L)==m_exp, sum(d$status==2L)==cd_exp,
            all(is.finite(d$atime)), all(d$atime > 0))
  d }
wcif1 <- function(time, evt, w, t0) { o <- order(time); time<-time[o]; evt<-evt[o]; w<-w[o]
  Wrev <- rev(cumsum(rev(w))); et <- unique(time[evt != 0 & time <= t0]); if (!length(et)) return(0)
  S <- 1; cif <- 0
  for (u in et) { idx <- which(time==u); R <- Wrev[idx[1]]; if (R <= 0) next
    cif <- cif + S * sum(w[idx][evt[idx]==1]) / R; S <- S * (1 - sum(w[idx][evt[idx]!=0]) / R) }
  cif }
ipcw_parts <- function(d, t0) {
  cens <- as.integer(d$status == 0L)
  km <- cap("f_censKM", survival::survfit(Surv(d$atime, cens) ~ 1, weights = d$w)); tv <- km$time; sv <- km$surv
  Gm <- function(x) vapply(x, function(z){k<-which(tv< z); if(!length(k)) 1 else sv[max(k)]}, numeric(1))
  Ga <- function(x) vapply(x, function(z){k<-which(tv<=z); if(!length(k)) 1 else sv[max(k)]}, numeric(1))
  case <- d$status==1L & d$atime<=t0; c1 <- d$atime>t0; c2 <- d$status==2L & d$atime<=t0
  wc <- ifelse(case, d$w / pmax(Gm(d$atime), 1e-12), 0)
  wk <- ifelse(c1, d$w / pmax(Ga(t0), 1e-12), ifelse(c2, d$w / pmax(Gm(d$atime), 1e-12), 0))
  list(wc = wc, wk = wk) }
oc_at2 <- function(d, t0, pos, pit) { p <- ipcw_parts(d, t0); wc<-p$wc; wk<-p$wk
  se <- sum(wc[pos])/sum(wc); sp <- sum(wk[!pos])/sum(wk)
  c(Se=se, Sp=sp, PPV=(pit*se)/(pit*se+(1-pit)*(1-sp)), NPV=((1-pit)*sp)/(pit*(1-se)+(1-pit)*sp)) }
wauc_bin <- function(d, t0, pos) { p <- ipcw_parts(d, t0); wc<-p$wc; wk<-p$wk
  ci<-which(wc>0); ki<-which(wk>0); if(!length(ci)||!length(ki)) return(NA_real_)
  m <- as.numeric(pos); num <- 0
  for (i in ci) num <- num + wc[i]*sum(wk[ki]*((m[i]>m[ki]) + 0.5*(m[i]==m[ki]))); num/(sum(wc[ci])*sum(wk[ki])) }
hr_cch_pos <- function(d, pos) { dd <- d; dd$grp <- factor(ifelse(pos,"High","Low"), levels=c("Low","High"))
  f <- cap("f_cch", survival::cch(Surv(atime, status==1L) ~ grp, data=dd, subcoh=~insub, id=~id,
                                  cohort.size=745L, method="LinYing", robust=TRUE))
  cf <- summary(f)$coefficients; b<-unname(cf[1,"Value"]); se<-unname(cf[1,"SE"])
  c(hr=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=unname(cf[1,"p"])) }
hr_cox_pos <- function(d, pos) { dd <- d; dd$grp <- factor(ifelse(pos,"High","Low"), levels=c("Low","High"))
  f <- cap("f_cox", survival::coxph(Surv(atime, as.integer(status==1L)) ~ grp, data=dd, robust=TRUE))
  cf <- summary(f)$coefficients; b<-cf["grpHigh","coef"]; se<-cf["grpHigh","robust se"]
  c(hr=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=cf["grpHigh","Pr(>|z|)"]) }

run_cuts <- function(d, weighted, seed, cohort) {
  K <- nrow(CUTS); cutvals <- CUTS$cutoff_value; cutnm <- CUTS$cutoff_name
  keys <- expand.grid(hi = seq_along(HORIZONS), k = seq_len(K)); nk <- nrow(keys)
  pit <- setNames(vapply(HORIZONS, function(t0) wcif1(d$atime, d$status, d$w, t0), 0), as.character(HORIZONS))
  pt <- vector("list", nk)
  for (j in seq_len(nk)) { t0 <- HORIZONS[keys$hi[j]]; k <- keys$k[j]; pos <- d$score >= cutvals[k]
    oc <- oc_at2(d, t0, pos, pit[[as.character(t0)]]); au <- wauc_bin(d, t0, pos)
    hrf <- if (weighted) sum(d$w[pos])/sum(d$w) else mean(pos)
    hr <- tryCatch(if (weighted) hr_cch_pos(d, pos) else hr_cox_pos(d, pos),
                   error = function(e) c(hr=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_))
    pt[[j]] <- list(t0=t0, k=k, oc=oc, auc=au, hrf=hrf, hr=hr,
                    n=nrow(d), events=sum(d$status==1L & d$atime<=t0)) }
  strat <- if (weighted) split(seq_len(nrow(d)), d$cch) else list(all = seq_len(nrow(d)))
  set.seed(seed)
  Se<-Sp<-PPV<-NPV<-AUC<- lapply(seq_len(nk), function(i) numeric(0))
  used <- 0L; f_absent <- 0L; f_err <- 0L; f_nonfin <- 0L
  for (bi in seq_len(B_BOOT)) {
    ix <- unlist(lapply(strat, function(g) g[sample.int(length(g), length(g), replace=TRUE)]))
    db <- d[ix, , drop=FALSE]
    r <- tryCatch({
      pitb <- setNames(vapply(HORIZONS, function(t0) wcif1(db$atime, db$status, db$w, t0), 0), as.character(HORIZONS))
      absent <- FALSE; vv <- vector("list", nk)
      for (j in seq_len(nk)) { t0 <- HORIZONS[keys$hi[j]]; k <- keys$k[j]; pos <- db$score >= cutvals[k]
        if (all(pos) || !any(pos)) { absent <- TRUE; break }
        vv[[j]] <- c(oc_at2(db, t0, pos, pitb[[as.character(t0)]]), AUC = wauc_bin(db, t0, pos)) }
      if (absent) NULL else vv
    }, error = function(e) "err")
    if (identical(r, "err")) { f_err <- f_err + 1L; next }
    if (is.null(r)) { f_absent <- f_absent + 1L; next }
    if (!all(vapply(r, function(v) all(is.finite(v)), logical(1)))) { f_nonfin <- f_nonfin + 1L; next }
    used <- used + 1L
    for (j in seq_len(nk)) { v <- r[[j]]
      Se[[j]]<-c(Se[[j]],v["Se"]); Sp[[j]]<-c(Sp[[j]],v["Sp"]); PPV[[j]]<-c(PPV[[j]],v["PPV"])
      NPV[[j]]<-c(NPV[[j]],v["NPV"]); AUC[[j]]<-c(AUC[[j]],v["AUC"]) } }
  qq <- function(x) if (length(x) >= 50) unname(quantile(x, c(.025, .975))) else c(NA_real_, NA_real_)
  rows <- list()
  for (j in seq_len(nk)) { P <- pt[[j]]; k <- keys$k[j]
    seCI<-qq(Se[[j]]); spCI<-qq(Sp[[j]]); ppvCI<-qq(PPV[[j]]); npvCI<-qq(NPV[[j]]); aucCI<-qq(AUC[[j]])
    rows[[length(rows)+1L]] <- data.frame(
      cohort = cohort, horizon_months = P$t0, cutoff_name = cutnm[k], cutoff_value = round(cutvals[k], 6),
      derivation = CUTS$derivation[k],
      high_risk_fraction = round(P$hrf, 5),
      high_risk_fraction_basis = if (weighted) "phase-two-weighted" else "observed",
      n = P$n, events = P$events,
      Se = round(P$oc["Se"], 5),  Se_lo = round(seCI[1], 5),  Se_hi = round(seCI[2], 5),
      Sp = round(P$oc["Sp"], 5),  Sp_lo = round(spCI[1], 5),  Sp_hi = round(spCI[2], 5),
      PPV = round(P$oc["PPV"], 5), PPV_lo = round(ppvCI[1], 5), PPV_hi = round(ppvCI[2], 5),
      NPV = round(P$oc["NPV"], 5), NPV_lo = round(npvCI[1], 5), NPV_hi = round(npvCI[2], 5),
      AUC_binary = round(P$auc, 5), AUC_lo = round(aucCI[1], 5), AUC_hi = round(aucCI[2], 5),
      HR = round(P$hr["hr"], 5), HR_lo = round(P$hr["lo"], 5), HR_hi = round(P$hr["hi"], 5), HR_p = signif(P$hr["p"], 5),
      HR_estimator = if (weighted) "Lin-Ying case-cohort cch (analytic robust CI)" else "complete-cohort robust Cox (analytic CI)",
      uncertainty = if (weighted) "conditional design-stratified bootstrap" else "patient bootstrap",
      boot_seed = seed, boot_used = used, boot_failed = B_BOOT - used, stringsAsFactors = FALSE) }
  do.call(rbind, rows) }

djhu_f <- bet(J0[["Met-Score prob"]], J0$met_time, J0$met, J0$os, J0$os_time,
              J0[["post_rp_patients_cchdef"]], 239L, 93L, 6L)
ddur_f <- bet(D0$MetScore_prob, D0$surgmets, D0$mets, D0$dead, D0$limbo, NULL, 555L, 40L, 167L)
panelF_thr <- rbind(run_cuts(djhu_f, TRUE, SEED_F, "JHU"), run_cuts(ddur_f, FALSE, SEED_F + 1L, "Durham"))
write.csv(panelF_thr, file.path(OUT, "panelF_threshold_strategies.csv"), row.names = FALSE)
# cutoff-analysis validation checks
addchk("panelF_dev_prevalence", mean(devy), 0.306, 0.001)
addchk("panelF_youden_resub_vs_locked", you_resub, THR_YOUDEN, 0.05)   # near-agreement (locked = between-sample midpoint)
{ posY <- ddur_f$score >= THR_YOUDEN; pit120 <- wcif1(ddur_f$atime, ddur_f$status, ddur_f$w, 120)
  ocY <- oc_at2(ddur_f, 120, posY, pit120)
  addchk("panelF_binaryAUC_eq_meanSeSp(Durham,120,Youden)", unname(wauc_bin(ddur_f, 120, posY)),
         unname((ocY["Se"] + ocY["Sp"]) / 2), 1e-8) }

## ---- validation + warning ledgers ---------------------------------------
mvchecks[[length(mvchecks)+1L]] <- data.frame(check = "panelA_LinYing_reference_HR", value = signif(unname(lin["hr"]),7),
  reference = signif(unname(lin["hr"]),7), abs_diff = 0, pass = TRUE, stringsAsFactors = FALSE)
mvtab <- do.call(rbind, mvchecks)
write.csv(mvtab, file.path(OUT, "method_validation_checks.csv"), row.names = FALSE)
wl <- if (length(.warn$log)) do.call(rbind, .warn$log) else data.frame(context = character(0), message = character(0))
write.csv(unique(wl), file.path(OUT, "warning_convergence_ledger.csv"), row.names = FALSE)

## ---- console summary -----------------------------------------------------
cat("\n== Figure S6 aggregates ->", OUT, "==\n")
cat(sprintf("panel a JHU start-stop full HR %.4f (LinYing %.4f, Prentice %.4f)\n", b_full["hr"], lin["hr"], pren))
cat("panel a rows:\n"); print(panelA[, c("cohort","scenario","n","events","hr","ci_lo","ci_hi")], row.names = FALSE)
cat("\npanel b:\n"); print(panelB[, c("cohort","model","n","events","hr","ci_lo","ci_hi","available")], row.names = FALSE)
cat(sprintf("\npanel c Durham Met-Score x race interaction ratio=%.4f (%.4f-%.4f) p=%.4g; global PH p=%.4g\n",
  panelC$interaction_ratio[1], panelC$interaction_lo[1], panelC$interaction_hi[1], panelC$interaction_p[1], panelC$global_ph_p[1]))
print(panelC[, c("race","label","n","events","hr","ci_lo","ci_hi")], row.names = FALSE)
cat("\npanel d:\n"); print(panelD[, c("cohort","term","hr","ci_lo","ci_hi","estimable")], row.names = FALSE)
cat("\npanel e:\n"); print(panelE[, c("cohort","adjustment","hr","ci_lo","ci_hi")], row.names = FALSE)
cat("\npanel f cutoffs (development-derived):\n")
print(CUTS[, c("cutoff_name","cutoff_value_out","dev_sensitivity","dev_specificity","is_primary")], row.names = FALSE)
cat("panel f external validation at 10y (120mo):\n")
print(panelF_thr[panelF_thr$horizon_months == 120L, c("cohort","cutoff_name","high_risk_fraction","Se","Sp","AUC_binary","HR")], row.names = FALSE)
cat(sprintf("resubstitution Youden = %.4f (locked = %.4f); dev prevalence = %.3f\n", you_resub, THR_YOUDEN, mean(devy)))
cat("(CCP correlation retained as numerical output panelF_ccp_correlation.csv)\n")
cat(sprintf("\nCCP coverage JHU %d/31, Durham %d/31; warnings captured %d\n", length(jr$used), length(dr$used), nrow(unique(wl))))
cat("Done.\n")

############################################################################
# Figure S5 panels c/d producer: exploratory signature-robustness modeling.
#
# This script does NOT redevelop, refit, retune, or alter the frozen 41-feature
# Met-Score. The five modeling strategies compared here are exploratory. Panel c
# is a nested internal-external cross-validation over the six discovery cohorts
# (development-strategy robustness). Panel d is an external discrimination
# sensitivity analysis in JHU and Durham where M1 is the actual frozen Met-Score
# and the four alternatives are training-fitted. It writes identifier-free
# aggregates to outs/FigureS5/ and never uses JHU/Durham outcomes or feature
# availability for gene selection, tuning, scaling, or imputation.
############################################################################
suppressWarnings(suppressMessages({
  library(glmnet); library(limma); library(timeROC); library(survival); library(readxl); library(pROC)
}))
if (!("--figure-s5-only" %in% commandArgs(trailingOnly = TRUE)))
  message("Met_PCa_LASSO_vs_Ridge.R produces Figure S5 panels c/d; expected flag --figure-s5-only.")

ROOT <- local({ v <- Sys.getenv("MET_PCA_ROOT", ""); if (nzchar(v) && dir.exists(v)) v else normalizePath(".") })
setwd(ROOT)
OUT <- file.path(ROOT, "outs", "FigureS5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source("code/utils/locked_metscore.R"); MODEL <- load_locked_metscore(); MG <- MODEL$feature_names  # 41 deployed
HORIZON <- 120L
LGRID <- exp(seq(log(10), log(1e-4), length.out = 60))       # broad, decreasing
WGRID <- exp(seq(log(1e4), log(1e-8), length.out = 100))     # one-time widen grid

## ---- helpers -------------------------------------------------------------
qn <- function(m) normalizeBetweenArrays(m, method = "quantile")
bridge41 <- function(raw41, ref41) {  # per-gene rank -> training-reference quantile
  out <- matrix(NA_real_, nrow(raw41), ncol(raw41), dimnames = dimnames(raw41))
  for (g in rownames(raw41)) {
    rk <- rank(raw41[g, ], ties.method = "average") / (ncol(raw41) + 1)
    out[g, ] <- quantile(ref41[g, ], probs = rk)
  }
  out
}
hedges_g_pooled <- function(mat, y, coh) {  # inverse-variance fixed pool of per-cohort Hedges g
  G <- nrow(mat); num <- rep(0, G); den <- rep(0, G)
  for (ct in unique(coh)) {
    x1 <- mat[, coh == ct & y == 1, drop = FALSE]; x0 <- mat[, coh == ct & y == 0, drop = FALSE]
    n1 <- ncol(x1); n0 <- ncol(x0); if (n1 < 2 || n0 < 2) next
    sp <- sqrt(((n1 - 1) * apply(x1, 1, var) + (n0 - 1) * apply(x0, 1, var)) / (n1 + n0 - 2))
    d <- (rowMeans(x1) - rowMeans(x0)) / sp
    g <- (1 - 3 / (4 * (n1 + n0) - 9)) * d
    vg <- (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0)); w <- 1 / vg
    ok <- is.finite(g) & is.finite(w); num[ok] <- num[ok] + w[ok] * g[ok]; den[ok] <- den[ok] + w[ok]
  }
  p <- ifelse(den > 0, num / den, 0); names(p) <- rownames(mat); p
}
safe_auc <- function(y, p) tryCatch(as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))),
                                     error = function(e) NA_real_)
auc_delong <- function(y, p) {
  r <- pROC::roc(y, p, quiet = TRUE, direction = "<"); ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  c(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3])
}

## ---- discovery data ------------------------------------------------------
e <- new.env(); load("outs/MetastasisData_JHUOut.rda", envir = e)
trainMat <- get("trainMat", e); trainGroup <- get("trainGroup", e)
testMat <- get("testMat", e)                                   # JHU expression, training gene space
y_all <- as.integer(trainGroup == "Mets")
map <- read.csv("outs/train_sample_to_gse.csv", stringsAsFactors = FALSE)
cohort <- map$gse[match(colnames(trainMat), map$sample_id)]
stopifnot(!any(is.na(cohort)))
EXP <- list(GSE116918 = c(226, 22), GSE41408 = c(39, 9), GSE46691 = c(333, 212),
            GSE51066 = c(34, 51), GSE55935 = c(36, 8), GSE70769 = c(26, 4))
for (ct in names(EXP)) stopifnot(sum(cohort == ct & y_all == 0) == EXP[[ct]][1],
                                 sum(cohort == ct & y_all == 1) == EXP[[ct]][2])
stopifnot(ncol(trainMat) == 1000L, sum(y_all == 1) == 306L, sum(y_all == 0) == 694L)
stopifnot(all(MG %in% rownames(trainMat)))
COHORTS <- names(EXP)

## QN reference cache keyed by the sorted training-cohort set; only the 41
## deployed rows are kept, so the full transcriptome is never held in memory.
.ref_cache <- new.env()
get_ref <- function(cohs) {
  key <- paste(sort(cohs), collapse = "_")
  if (!is.null(.ref_cache[[key]])) return(.ref_cache[[key]])
  idx <- which(cohort %in% cohs)
  full <- qn(trainMat[, idx]); qr41 <- full[MG, , drop = FALSE]; rm(full)
  ranked <- names(sort(abs(hedges_g_pooled(qr41, y_all[idx], cohort[idx])), decreasing = TRUE))
  val <- list(qr41 = qr41, ranked = ranked, idx = idx)
  .ref_cache[[key]] <- val; val
}

MODELS <- list(
  list(key = "full_ridge",  label = "Full-panel ridge (IECV)", alpha = 0,   feat = "all"),
  list(key = "top10_ridge", label = "Top-10 ridge",            alpha = 0,   feat = "top10"),
  list(key = "top20_ridge", label = "Top-20 ridge",            alpha = 0,   feat = "top20"),
  list(key = "lasso",       label = "LASSO (alpha=1)",         alpha = 1,   feat = "all"),
  list(key = "enet",        label = "Elastic net (alpha=0.5)", alpha = 0.5, feat = "all"))
pick_feats <- function(spec, ranked)
  if (spec$feat == "all") MG else if (spec$feat == "top10") ranked[1:10] else ranked[1:20]

## leave-one-cohort-out tuning over the cohorts in `cohs`; grid decreasing so
## the smallest index is the strongest penalty (deterministic tie rule).
loco_tune <- function(spec, cohs, grid) {
  amat <- matrix(NA_real_, length(grid), length(cohs), dimnames = list(NULL, cohs))
  audit <- list()
  for (ci in seq_along(cohs)) {
    ref <- get_ref(setdiff(cohs, cohs[ci])); feats <- pick_feats(spec, ref$ranked)
    te <- which(cohort == cohs[ci])
    br <- bridge41(trainMat[MG, te, drop = FALSE], ref$qr41)
    w <- character(0)
    fit <- withCallingHandlers(
      glmnet(t(ref$qr41[feats, , drop = FALSE]), y_all[ref$idx], family = "binomial",
             alpha = spec$alpha, lambda = grid, standardize = FALSE),
      warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
    pp <- predict(fit, newx = t(br[feats, , drop = FALSE]), s = grid, type = "response")
    col <- vapply(seq_along(grid), function(li) safe_auc(y_all[te], pp[, li]), numeric(1))
    # fail closed: no silent NA averaging
    if (any(fit$jerr != 0))
      stop(sprintf("loco_tune %s / inner %s: glmnet jerr=%s", spec$label, cohs[ci], paste(fit$jerr, collapse = ",")))
    if (!all(is.finite(col)))
      stop(sprintf("loco_tune %s / inner %s: non-finite held-out AUC at %d/%d grid points",
                   spec$label, cohs[ci], sum(!is.finite(col)), length(col)))
    amat[, ci] <- col
    audit[[cohs[ci]]] <- list(jerr = as.integer(any(fit$jerr != 0)),
                              warnings = if (length(w)) paste(unique(w), collapse = " || ") else "")
  }
  m <- rowMeans(amat); sel <- min(which(m >= max(m) - 1e-12))
  list(mean_auc = m, sel = sel, lambda = grid[sel], amat = amat,
       boundary = (sel == 1L || sel == length(grid)), grid = grid, audit = audit)
}
# the actual tuning grid is passed in so a widened grid is honored by the fit
fit_final <- function(spec, cohs, lambda, grid) {
  ref <- get_ref(cohs); feats <- pick_feats(spec, ref$ranked)
  fit <- glmnet(t(ref$qr41[feats, , drop = FALSE]), y_all[ref$idx], family = "binomial",
                alpha = spec$alpha, lambda = grid, standardize = FALSE)
  co <- as.numeric(coef(fit, s = lambda))[-1]
  list(fit = fit, feats = feats, coef = setNames(co, feats), ref = ref)
}

## ---- PANEL C: nested internal-external cross-validation ------------------
cat("Panel c: nested internal-external CV over 6 discovery cohorts ...\n")
c_outer <- list(); c_inner <- list(); c_tune <- list(); c_models <- list()
for (om in COHORTS) {
  tr_cohs <- setdiff(COHORTS, om); te <- which(cohort == om)
  for (spec in MODELS) {
    tu <- loco_tune(spec, tr_cohs, LGRID); widened <- FALSE
    if (tu$boundary) { tu <- loco_tune(spec, tr_cohs, WGRID); widened <- TRUE }
    ff <- fit_final(spec, tr_cohs, tu$lambda, tu$grid)
    br <- bridge41(trainMat[MG, te, drop = FALSE], ff$ref$qr41)
    warns <- character(0)
    p <- withCallingHandlers(
      as.numeric(predict(ff$fit, newx = t(br[ff$feats, , drop = FALSE]), s = tu$lambda, type = "response")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    a <- auc_delong(y_all[te], p)
    active <- sum(abs(ff$coef) > 0)
    conv <- fit_conv <- all(ff$fit$jerr == 0) && is.finite(a["auc"])
    c_outer[[length(c_outer) + 1]] <- data.frame(
      fold_omitted = om, model = spec$label, n = length(te),
      cases = sum(y_all[te] == 1), controls = sum(y_all[te] == 0), n_genes = length(ff$feats),
      selected_lambda = signif(tu$lambda, 8), active_features = active, boundary_widened = widened,
      converged = conv, auc = round(a["auc"], 7), auc_lo = round(a["lo"], 7), auc_hi = round(a["hi"], 7),
      warnings = if (length(warns)) paste(unique(warns), collapse = " || ") else "", stringsAsFactors = FALSE)
    inner_at_sel <- tu$amat[tu$sel, ]
    for (ic in names(inner_at_sel))
      c_inner[[length(c_inner) + 1]] <- data.frame(fold_omitted = om, model = spec$label,
        inner_cohort = ic, selected_lambda = signif(tu$lambda, 8),
        inner_auc = round(inner_at_sel[[ic]], 7),
        inner_jerr = tu$audit[[ic]]$jerr, inner_warnings = tu$audit[[ic]]$warnings,
        stringsAsFactors = FALSE)
    c_tune[[length(c_tune) + 1]] <- data.frame(fold_omitted = om, model = spec$label,
      lambda = signif(tu$grid, 8), mean_inner_auc = round(tu$mean_auc, 7),
      selected = tu$grid == tu$lambda, stringsAsFactors = FALSE)
    c_models[[length(c_models) + 1]] <- data.frame(fold_omitted = om, model = spec$label,
      n_genes = length(ff$feats), active_features = active,
      selected_genes = paste(ff$feats, collapse = ";"),
      coefficients = paste(sprintf("%s:%.6g", names(ff$coef), ff$coef), collapse = ";"),
      stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, c_outer),  file.path(OUT, "FigureS5_panelC_outer_auc.csv"), row.names = FALSE)
write.csv(do.call(rbind, c_inner),  file.path(OUT, "FigureS5_panelC_inner_auc.csv"), row.names = FALSE)
write.csv(do.call(rbind, c_tune),   file.path(OUT, "FigureS5_panelC_tuning.csv"), row.names = FALSE)
write.csv(do.call(rbind, c_models), file.path(OUT, "FigureS5_panelC_models.csv"), row.names = FALSE)

## ---- PANEL D: external discrimination sensitivity ------------------------
cat("Panel d: external time-dependent AUC in JHU and Durham ...\n")
ALTS <- MODELS[sapply(MODELS, function(m) m$key != "full_ridge")]   # 4 alternatives (M1 = frozen)
ref6 <- get_ref(COHORTS)                                            # all-1000 QN reference (41 rows)
# final alternative models: tune via 6-cohort LOCO, fit on all 1000
final <- list(); d_final_rows <- list()
for (spec in ALTS) {
  tu <- loco_tune(spec, COHORTS, LGRID); widened <- FALSE
  if (tu$boundary) { tu <- loco_tune(spec, COHORTS, WGRID); widened <- TRUE }
  warns <- character(0)
  ff <- withCallingHandlers(fit_final(spec, COHORTS, tu$lambda, tu$grid),
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  final[[spec$key]] <- list(spec = spec, lambda = tu$lambda, ff = ff)
  d_final_rows[[length(d_final_rows) + 1]] <- data.frame(
    model = spec$label, n_train = 1000L, selected_lambda = signif(tu$lambda, 8),
    boundary_widened = widened, n_genes = length(ff$feats), active_features = sum(abs(ff$coef) > 0),
    converged = all(ff$fit$jerr == 0), selected_genes = paste(ff$feats, collapse = ";"),
    warnings = if (length(warns)) paste(unique(warns), collapse = " || ") else "", stringsAsFactors = FALSE)
}

# external bridged 41-gene matrices (all-training QN reference), JHU then Durham
jhu_ids <- sub(".*?(JHU[0-9]+).*", "\\1", colnames(testMat))
jhu_br <- bridge41(testMat[MG, , drop = FALSE], ref6$qr41); colnames(jhu_br) <- jhu_ids
jhu_X <- t(jhu_br)                                                  # samples x 41

dur_raw <- as.data.frame(read_excel("data/Durham_cohort_and_GRID_cohort/Durham_cohort_011526.xlsx",
                                     sheet = "eset_gene_filtered"))
drn <- dur_raw$Symbol; dm <- apply(dur_raw[, -c(1, 2)], 2, as.numeric); rownames(dm) <- drn
dv <- new.env(); load("output/Durham/durham_metscore_batchcorrected.rda", envir = dv); clin <- get("clin_valid", dv)
cx <- new.env(); load("outs/coxdata.rda", envir = cx); coxj <- get("CoxData_jhu", cx)
dur_ids <- as.character(clin$sample_id)
present <- MG[MG %in% rownames(dm)]; missing_d <- setdiff(MG, present)
stopifnot(identical(missing_d, "KCTD14"))                          # only KCTD14 may be absent
# fail-closed integrity checks before external scoring
stopifnot(!anyDuplicated(jhu_ids), !anyDuplicated(dur_ids))         # unique cohort IDs
stopifnot(setequal(jhu_ids, as.character(coxj$sample_id)))         # JHU marker/event ID sets align
stopifnot(all(dur_ids %in% colnames(dm)))                          # Durham clinical IDs present in expression
.dup_req <- MG[vapply(MG, function(g) sum(rownames(dm) == g) > 1L, logical(1))]
stopifnot(length(.dup_req) == 0L)                                  # no duplicate required gene rows
stopifnot(all(is.finite(dm[present, dur_ids, drop = FALSE])), all(is.finite(testMat[MG, , drop = FALSE])))
dur_br <- matrix(NA_real_, length(MG), length(dur_ids), dimnames = list(MG, dur_ids))
dur_br[present, ] <- bridge41(dm[present, dur_ids, drop = FALSE], ref6$qr41[present, , drop = FALSE])
dur_br["KCTD14", ] <- mean(ref6$qr41["KCTD14", ])                  # training-reference mean imputation
dur_X <- t(dur_br)
stopifnot(all(is.finite(jhu_X)), all(is.finite(dur_X)))            # finite bridged matrices

# frozen-score round trip: reconstruct M1 from the bridged matrices, verify vs saved
p_jhu_frozen_saved <- setNames(as.numeric(coxj[["Met-Score prob"]]), as.character(coxj$sample_id))[jhu_ids]
p_dur_frozen_saved <- setNames(as.numeric(clin$MetScore_prob), as.character(clin$sample_id))[dur_ids]
rt_jhu <- locked_metscore_score(jhu_X, MODEL)$prob
rt_dur <- locked_metscore_score(dur_X, MODEL)$prob
rt_j <- max(abs(rt_jhu - p_jhu_frozen_saved)); rt_d <- max(abs(rt_dur - p_dur_frozen_saved))
cm_j <- sum((rt_jhu >= MODEL$threshold) != (p_jhu_frozen_saved >= MODEL$threshold))
cm_d <- sum((rt_dur >= MODEL$threshold) != (p_dur_frozen_saved >= MODEL$threshold))
cat(sprintf("Frozen round trip: JHU max|diff|=%.3e (class mismatch %d); Durham max|diff|=%.3e (class mismatch %d)\n",
            rt_j, cm_j, rt_d, cm_d))
if (max(rt_j, rt_d) > 1e-12 || (cm_j + cm_d) > 0)
  stop("Frozen Met-Score round trip exceeded tolerance; bridged matrices do not reproduce the saved scores.")

# canonical competing-event clocks (must equal the signature-robustness event ledger)
build_ce <- function(met, tm, dth, td, nmet, ndth) {
  st <- ifelse(met == 1L, 1L, ifelse(dth == 1L & is.finite(td) & td <= tm, 2L, 0L))
  tt <- ifelse(st == 2L, td, tm)
  stopifnot(sum(st == 1L) == nmet, sum(st == 2L) == ndth, all(is.finite(tt)), all(tt > 0))
  list(status = as.integer(st), time = as.numeric(tt))
}
jce <- build_ce(as.integer(coxj$met), as.numeric(coxj$met_time), as.integer(coxj$os), as.numeric(coxj$os_time), 93L, 6L)
dce <- build_ce(as.integer(clin$mets), as.numeric(clin$surgmets), as.integer(clin$dead), as.numeric(clin$limbo), 40L, 167L)
names(jce$status) <- names(jce$time) <- as.character(coxj$sample_id)
names(dce$status) <- names(dce$time) <- as.character(clin$sample_id)
led <- data.frame(cohort = c("JHU", "Durham"), n = c(239L, 555L),
                  metastases = c(sum(jce$status == 1), sum(dce$status == 1)),
                  competing_deaths = c(sum(jce$status == 2), sum(dce$status == 2)),
                  censored = c(sum(jce$status == 0), sum(dce$status == 0)), stringsAsFactors = FALSE)
write.csv(led, file.path(OUT, "FigureS5_panelD_event_ledger.csv"), row.names = FALSE)

# marker set per cohort: M1 frozen (saved) + 4 alternatives (predicted)
marker_jhu <- list("Frozen Met-Score" = p_jhu_frozen_saved)
marker_dur <- list("Frozen Met-Score" = p_dur_frozen_saved)
for (spec in ALTS) {
  fj <- final[[spec$key]]
  marker_jhu[[spec$label]] <- setNames(as.numeric(predict(fj$ff$fit, newx = jhu_X[, fj$ff$feats, drop = FALSE],
                                                           s = fj$lambda, type = "response")), jhu_ids)
  marker_dur[[spec$label]] <- setNames(as.numeric(predict(fj$ff$fit, newx = dur_X[, fj$ff$feats, drop = FALSE],
                                                           s = fj$lambda, type = "response")), dur_ids)
}
# External discrimination uses the accepted Figure-2 estimator: competing-risk,
# phase-two-weighted, exact-10-year def-2 IPCW AUC. JHU carries the case-cohort
# phase-two weights (subcohort controls 745/265, cases 1) and a conditional
# design-stratified bootstrap; Durham is the complete cohort (w=1) with a patient
# bootstrap. All five predictions are evaluated together on each shared resample so
# comparisons with the frozen model are paired.
ALPHA_S5D <- 265 / 745
B_S5D <- as.integer(Sys.getenv("S5D_B", "2000"))
SEED_S5D <- c(JHU = 20260814L, Durham = 20260815L)
ipcw_parts <- function(status, atime, w, t0) {
  cens <- as.integer(status == 0L)
  km <- survival::survfit(Surv(atime, cens) ~ 1, weights = w); tv <- km$time; sv <- km$surv
  Gm <- function(x) vapply(x, function(z){k<-which(tv< z); if(!length(k)) 1 else sv[max(k)]}, numeric(1))
  Ga <- function(x) vapply(x, function(z){k<-which(tv<=z); if(!length(k)) 1 else sv[max(k)]}, numeric(1))
  case <- status==1L & atime<=t0; c1 <- atime>t0; c2 <- status==2L & atime<=t0
  wc <- ifelse(case, w/pmax(Gm(atime),1e-12), 0)
  wk <- ifelse(c1, w/pmax(Ga(t0),1e-12), ifelse(c2, w/pmax(Gm(atime),1e-12), 0))
  list(wc = wc, wk = wk) }
wauc_from <- function(wc, wk, m) {   # weighted def-2 AUC from precomputed IPCW weights
  ci <- which(wc>0); ki <- which(wk>0); if(!length(ci)||!length(ki)) return(NA_real_)
  num <- 0; for (i in ci) num <- num + wc[i]*sum(wk[ki]*((m[i]>m[ki]) + 0.5*(m[i]==m[ki]))); num/(sum(wc[ci])*sum(wk[ki])) }

design <- list(
  JHU = list(status=jce$status, atime=jce$time, ids=names(jce$status),
             w=ifelse(as.character(coxj[["post_rp_patients_cchdef"]])=="Sub-cohort controls", 1/ALPHA_S5D, 1),
             strat=as.character(coxj[["post_rp_patients_cchdef"]]),
             basis="JHU case-cohort: phase-two-weighted def-2 IPCW AUC; conditional design-stratified bootstrap"),
  Durham = list(status=dce$status, atime=dce$time, ids=names(dce$status),
             w=rep(1, length(dce$status)), strat=rep("cohort", length(dce$status)),
             basis="Durham complete cohort: def-2 IPCW AUC; patient bootstrap"))
d_auc <- list(); d_pair <- list(); boot_led <- list()
for (co in c("JHU", "Durham")) {
  D <- design[[co]]; mk <- if (co == "JHU") marker_jhu else marker_dur; ids <- D$ids; mnames <- names(mk)
  M <- sapply(mnames, function(nm) as.numeric(mk[[nm]][ids])); rownames(M) <- ids   # subjects x models
  pp <- ipcw_parts(D$status, D$atime, D$w, HORIZON)
  auc_pt <- setNames(vapply(mnames, function(nm) wauc_from(pp$wc, pp$wk, M[, nm]), numeric(1)), mnames)
  idx_by <- split(seq_along(ids), D$strat); set.seed(SEED_S5D[[co]])
  bA <- matrix(NA_real_, B_S5D, length(mnames), dimnames = list(NULL, mnames))
  succ <- 0L; f_case <- 0L; f_ctrl <- 0L; f_nf <- 0L
  for (b in seq_len(B_S5D)) {
    ix <- unlist(lapply(idx_by, function(g) g[sample.int(length(g), length(g), replace = TRUE)]))
    pb <- ipcw_parts(D$status[ix], D$atime[ix], D$w[ix], HORIZON)
    if (!any(pb$wc > 0)) { f_case <- f_case + 1L; next }
    if (!any(pb$wk > 0)) { f_ctrl <- f_ctrl + 1L; next }
    row <- vapply(mnames, function(nm) wauc_from(pb$wc, pb$wk, M[ix, nm]), numeric(1))
    if (any(!is.finite(row))) { f_nf <- f_nf + 1L; next }
    succ <- succ + 1L; bA[succ, ] <- row }
  bA <- bA[seq_len(succ), , drop = FALSE]
  qci <- function(x){ x <- x[is.finite(x)]; if (length(x) >= 50) unname(quantile(x, c(.025, .975))) else c(NA_real_, NA_real_) }
  for (nm in mnames) {
    ci <- qci(bA[, nm])
    d_auc[[length(d_auc) + 1]] <- data.frame(cohort = co, model = nm, is_frozen = (nm == "Frozen Met-Score"),
      frozen_or_trained = ifelse(nm == "Frozen Met-Score", "frozen locked classifier", "trained in development data only"),
      n = length(ids), cases_by_120 = sum(D$status == 1L & D$atime <= HORIZON),
      competing_deaths_by_120 = sum(D$status == 2L & D$atime <= HORIZON), beyond_120 = sum(D$atime > HORIZON),
      censored_before_120 = sum(D$status == 0L & D$atime <= HORIZON), horizon_months = HORIZON,
      auc = round(auc_pt[[nm]], 7), ci_lo = round(ci[1], 7), ci_hi = round(ci[2], 7),
      boot_attempts = B_S5D, boot_success = succ, boot_fail = B_S5D - succ,
      control_def = "definition 2 (competing death as control)", estimator = D$basis, stringsAsFactors = FALSE) }
  for (spec in ALTS) {
    nm <- spec$label; dd <- bA[, nm] - bA[, "Frozen Met-Score"]; dd <- dd[is.finite(dd)]
    dci <- qci(dd)
    pv <- if (length(dd) >= 50) min(1, 2 * min((sum(dd <= 0) + 1) / (length(dd) + 1), (sum(dd >= 0) + 1) / (length(dd) + 1))) else NA_real_
    d_pair[[length(d_pair) + 1]] <- data.frame(cohort = co, model = nm,
      delta_auc = round(auc_pt[[nm]] - auc_pt[["Frozen Met-Score"]], 7),
      delta_lo = round(dci[1], 7), delta_hi = round(dci[2], 7), p_paired = pv, boot_success = length(dd),
      method = paste0("paired ", if (co == "JHU") "conditional design-stratified" else "patient", " bootstrap vs frozen (def-2 IPCW AUC)"),
      stringsAsFactors = FALSE) }
  boot_led[[length(boot_led) + 1]] <- data.frame(cohort = co, B = B_S5D, success = succ,
    fail_no_cases = f_case, fail_no_controls = f_ctrl, fail_nonfinite = f_nf, seed = SEED_S5D[[co]], stringsAsFactors = FALSE) }
d_pair <- do.call(rbind, d_pair); d_pair$q_bh <- signif(p.adjust(d_pair$p_paired, method = "BH"), 7); d_pair$p_paired <- signif(d_pair$p_paired, 7)
write.csv(do.call(rbind, d_auc), file.path(OUT, "FigureS5_panelD_timeAUC.csv"), row.names = FALSE)
write.csv(d_pair, file.path(OUT, "FigureS5_panelD_paired.csv"), row.names = FALSE)
write.csv(do.call(rbind, boot_led), file.path(OUT, "FigureS5_panelD_bootstrap_ledger.csv"), row.names = FALSE)
fm <- do.call(rbind, d_final_rows)
fm$roundtrip_jhu_maxabs <- signif(rt_j, 3); fm$roundtrip_durham_maxabs <- signif(rt_d, 3)
write.csv(fm, file.path(OUT, "FigureS5_panelD_final_models.csv"), row.names = FALSE)

cat("\nPanel c outer AUCs:\n"); print(do.call(rbind, c_outer)[, c("fold_omitted","model","auc","selected_lambda","active_features","boundary_widened")], row.names = FALSE)
cat("\nPanel d time-dependent AUCs:\n"); print(do.call(rbind, d_auc)[, c("cohort","model","auc","ci_lo","ci_hi","cases_by_120","competing_deaths_by_120","beyond_120","censored_before_120")], row.names = FALSE)
cat("\nPanel d paired (alternative - frozen):\n"); print(d_pair, row.names = FALSE)
cat("\nFigure S5 panels c/d written to", OUT, "\n")

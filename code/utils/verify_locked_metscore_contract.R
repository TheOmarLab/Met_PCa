## Verifier for the frozen Met-Score locked contract.
##
## Usage:
##   Rscript --vanilla code/utils/verify_locked_metscore_contract.R \
##     --output-dir <dir> --r2-deposit-dir <dir>
##
## Checks: live artifact SHA-256 vs the pinned value and metadata; tracked config
## vs the legacy artifact and the R2 deposit; 45-gene vs 41-deployed identity
## (exact row-wise); helper round-trip vs legacy predict() on the complete 1000-
## record training matrix and the complete 239-record JHU matrix; helper
## failure-mode tests; config-tamper invariants (temporary copies); successful-
## behaviour tests. Durham is recorded UNAVAILABLE (no provenance-bound post-bridge
## input exists). Result is PARTIAL (exit 2) when required checks pass but Durham
## has not executed; FAIL (exit 1) if any required check fails.

options(warn = 1)
library(glmnet)   # load + predict() the legacy cv.glmnet object
library(limma)    # normalizeBetweenArrays for the training reference

## ---- args ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 1L && i < length(args)) return(args[i + 1L])
  NULL
}
output_dir     <- get_arg("--output-dir")
r2_deposit_dir <- get_arg("--r2-deposit-dir")
if (is.null(output_dir))     stop("verify: --output-dir <dir> is required.")
if (is.null(r2_deposit_dir)) stop("verify: --r2-deposit-dir <dir> is required.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(r2_deposit_dir))
  stop(sprintf("verify: --r2-deposit-dir does not exist: %s", r2_deposit_dir))

## ---- paths ---------------------------------------------------------------
REPO     <- getwd()
ARTIFACT <- file.path(REPO, "outs/metscore_logReg_genes_.rda")
FSG_RDA  <- file.path(REPO, "outs/filtersiggenes_MetaScore.rda")
JHU_RDA  <- file.path(REPO, "outs/MetastasisData_JHUOut.rda")
HELPER   <- file.path(REPO, "code/utils/locked_metscore.R")
DEP_GENE <- file.path(r2_deposit_dir, "metscore_45gene_list.csv")
DEP_COEF <- file.path(r2_deposit_dir, "metscore_coefficients.csv")
for (f in c(DEP_GENE, DEP_COEF))
  if (!file.exists(f)) stop(sprintf("verify: required deposit file missing: %s", f))
source(HELPER)

AUTH_SHA <- "713e5d9260527de3c9e5c25e7a70d5d29c0b2bfaeebe0f277a9db0f5e080c574"
PRED_TOL <- 1e-12   # numeric tolerance for prediction / legacy-double comparison only

## ---- accumulators --------------------------------------------------------
CV <- data.frame(check=character(), detail=character(), observed=character(),
                 expected=character(), status=character(), stringsAsFactors=FALSE)
add <- function(check, detail, observed, expected, status)
  CV <<- rbind(CV, data.frame(check=check, detail=detail,
    observed=as.character(observed), expected=as.character(expected),
    status=status, stringsAsFactors=FALSE))
FM  <- data.frame(test=character(), expectation=character(), errored=logical(),
                  message=character(), status=character(), stringsAsFactors=FALSE)
INV <- FM[0, ]
BEH <- data.frame(test=character(), expectation=character(),
                  observed=character(), status=character(), stringsAsFactors=FALSE)
.hardfail <- function(test, expectation, expr) {
  errored <- FALSE; msg <- ""
  tryCatch(force(expr), error=function(ee){errored<<-TRUE; msg<<-conditionMessage(ee)})
  cat(sprintf("  [hardfail] %-30s errored=%-5s -> %s\n", test, errored,
              if (errored) "PASS" else "FAIL"))
  data.frame(test=test, expectation=expectation, errored=errored, message=msg,
             status=if (errored) "PASS" else "FAIL", stringsAsFactors=FALSE)
}
addfm  <- function(test, expectation, expr) FM  <<- rbind(FM,  .hardfail(test, expectation, expr))
addinv <- function(test, expr)             INV <<- rbind(INV, .hardfail(test, "tampered config must be rejected", expr))
addbeh <- function(test, expectation, ok, observed) {
  BEH <<- rbind(BEH, data.frame(test=test, expectation=expectation,
    observed=as.character(observed), status=if (isTRUE(ok)) "PASS" else "FAIL",
    stringsAsFactors=FALSE))
  cat(sprintf("  [behaviour] %-32s -> %s\n", test, if (isTRUE(ok)) "PASS" else "FAIL"))
}

## ---- SHA-256 via an installed utility ------------------------------------
sha256_of <- function(path) {
  tool <- NA_character_; targ <- character(0)
  if (nzchar(Sys.which("shasum")))        { tool <- "shasum";    targ <- c("-a","256", path) }
  else if (nzchar(Sys.which("sha256sum"))){ tool <- "sha256sum"; targ <- path }
  else stop("verify: no SHA-256 utility (shasum/sha256sum) found on PATH.")
  res <- system2(tool, targ, stdout = TRUE, stderr = TRUE)
  st  <- attr(res, "status")
  if (!is.null(st) && st != 0L)
    stop(sprintf("verify: %s failed (status %s): %s", tool, st, paste(res, collapse=" ")))
  hit <- grep("^[0-9a-f]{64}", res, value = TRUE)
  if (length(hit) != 1L)
    stop(sprintf("verify: %s did not return one SHA-256: %s", tool, paste(res, collapse=" ")))
  sub("^([0-9a-f]{64}).*", "\\1", hit[1])
}

## ---- 0. load config through the helper -----------------------------------
model <- load_locked_metscore()
feats <- model$feature_names
add("config_load", "load_locked_metscore() feature count", length(feats), "41",
    if (length(feats) == 41L) "PASS" else "FAIL")

## ---- 1. live SHA-256 vs pinned + metadata --------------------------------
live_sha <- sha256_of(ARTIFACT)
add("sha256", "legacy artifact live SHA-256 computed", live_sha, "<64 hex>",
    if (grepl("^[0-9a-f]{64}$", live_sha)) "PASS" else "FAIL")
add("sha256", "live SHA-256 == pinned SHA", live_sha, AUTH_SHA,
    if (identical(live_sha, AUTH_SHA)) "PASS" else "FAIL")
add("sha256", "live SHA-256 == metadata source_artifact_sha256",
    model$source_artifact_sha256, live_sha,
    if (identical(live_sha, model$source_artifact_sha256)) "PASS" else "FAIL")

## ---- 2. config vs legacy artifact ----------------------------------------
e <- new.env(); load(ARTIFACT, envir=e)
m_legacy <- e$logReg_model; lam_leg <- e$LAMBDA_LOCKED; thr_leg <- e$LOCKED_THRESHOLD
ms_leg   <- as.character(e$Meta_Score)
co_leg   <- as.matrix(coef(m_legacy, s=lam_leg))
int_leg  <- co_leg[1,1]; gen_leg <- rownames(co_leg)[-1]
beta_leg <- co_leg[-1,1]; names(beta_leg) <- gen_leg
add("config_vs_legacy", "deployed gene order identical to coef(logReg_model)",
    identical(feats, gen_leg), "TRUE", if (identical(feats, gen_leg)) "PASS" else "FAIL")
add("config_vs_legacy", "deployed gene set identical to artifact Meta_Score",
    identical(feats, ms_leg), "TRUE", if (identical(feats, ms_leg)) "PASS" else "FAIL")
dbeta <- max(abs(model$beta[gen_leg] - beta_leg))
add("config_vs_legacy", "max abs coefficient diff vs legacy (numeric)",
    sprintf("%.3e", dbeta), "<= 1e-12", if (dbeta <= PRED_TOL) "PASS" else "FAIL")
dint <- abs(model$intercept - int_leg)
add("config_vs_legacy", "abs intercept diff vs legacy", sprintf("%.3e", dint),
    "<= 1e-12", if (dint <= PRED_TOL) "PASS" else "FAIL")
dlam <- abs(model$lambda - lam_leg)
add("config_vs_legacy", "abs lambda diff vs legacy", sprintf("%.3e", dlam),
    "<= 1e-12", if (dlam <= PRED_TOL) "PASS" else "FAIL")
dthr <- abs(model$threshold - thr_leg)
add("config_vs_legacy", "abs threshold diff vs legacy", sprintf("%.3e", dthr),
    "<= 1e-12", if (dthr <= PRED_TOL) "PASS" else "FAIL")

## ---- 3. config vs R2 deposit ---------------------------------------------
dep_coef <- read.csv(DEP_COEF, stringsAsFactors=FALSE)
dep_int  <- dep_coef$beta[dep_coef$gene == "(Intercept)"]
dep_only <- dep_coef[dep_coef$gene != "(Intercept)", ]
idx      <- match(feats, dep_only$gene)
depb_ok  <- !any(is.na(idx))
add("config_vs_deposit", "all 41 deployed genes present in deposit coefficients",
    depb_ok, "TRUE", if (depb_ok) "PASS" else "FAIL")
if (depb_ok) {
  ddb <- max(abs(model$beta[feats] - dep_only$beta[idx]))
  add("config_vs_deposit", "max abs coefficient diff vs deposit",
      sprintf("%.3e", ddb), "<= 1e-12", if (ddb <= PRED_TOL) "PASS" else "FAIL")
}
ddi <- abs(model$intercept - dep_int)
add("config_vs_deposit", "abs intercept diff vs deposit", sprintf("%.3e", ddi),
    "<= 1e-12", if (ddi <= PRED_TOL) "PASS" else "FAIL")
cvl <- data.frame(gene=feats, config_coef=model$beta[feats], legacy_coef=beta_leg[feats],
  deposit_coef=if (depb_ok) dep_only$beta[idx] else NA_real_,
  abs_diff_legacy=abs(model$beta[feats]-beta_leg[feats]),
  abs_diff_deposit=if (depb_ok) abs(model$beta[feats]-dep_only$beta[idx]) else NA_real_,
  frozen_order_ok=identical(feats, gen_leg), stringsAsFactors=FALSE)
write.csv(cvl, file.path(output_dir, "config_vs_legacy_coefficients.csv"), row.names=FALSE)

## ---- 4. 45-gene vs 41-deployed identity (exact row-wise) -----------------
fe <- new.env(); load(FSG_RDA, envir=fe); FSG <- fe[[ls(fe)[1]]]
sig_cfg <- read.csv(file.path(model$config_dir, "metscore_signature_v1.csv"), colClasses="character")
dep_gene <- read.csv(DEP_GENE, colClasses="character")
add("sig_45_41", "signature config has 45 rows", nrow(sig_cfg), "45",
    if (nrow(sig_cfg) == 45L) "PASS" else "FAIL")
cols3 <- c("gene","direction","used_in_locked_model")
rowwise_ok <- all(cols3 %in% names(sig_cfg)) && all(cols3 %in% names(dep_gene)) &&
              identical(sig_cfg[, cols3], dep_gene[, cols3])
add("sig_45_41", "signature rows identical to deposit (gene,direction,used; exact order)",
    rowwise_ok, "TRUE", if (rowwise_ok) "PASS" else "FAIL")
add("sig_45_41", "signature genes == repo Filter_SignatureGenes (45, as set)",
    setequal(sig_cfg$gene, FSG), "TRUE", if (setequal(sig_cfg$gene, FSG)) "PASS" else "FAIL")
used_true  <- sig_cfg$gene[sig_cfg$used_in_locked_model == "TRUE"]
used_false <- sig_cfg$gene[sig_cfg$used_in_locked_model == "FALSE"]
add("sig_45_41", "used==TRUE set == 41 deployed genes",
    setequal(used_true, feats), "TRUE", if (setequal(used_true, feats)) "PASS" else "FAIL")
add("sig_45_41", "used==FALSE set == 4 excluded genes",
    paste(sort(used_false), collapse=","), "ARL6IP1,KIAA1210,SEM1,TMEM121B",
    if (setequal(used_false, c("ARL6IP1","SEM1","KIAA1210","TMEM121B"))) "PASS" else "FAIL")

## ---- 5. round-trip on complete authorized matrices -----------------------
qn_to_train <- function(test_mat, train_ref) {
  common <- intersect(rownames(test_mat), rownames(train_ref)); out <- test_mat
  for (g in common) {
    rk <- rank(test_mat[g, ], ties.method="average") / (ncol(test_mat) + 1)
    out[g, ] <- quantile(train_ref[g, ], probs=rk, names=FALSE, type=7)
  }
  out
}
je <- new.env(); load(JHU_RDA, envir=je)
usedTrainMat <- normalizeBetweenArrays(je$trainMat, method="quantile")
expr_jhu     <- qn_to_train(je$testMat, usedTrainMat)
X_train_full <- t(usedTrainMat)   # 1000 x genes, complete
X_jhu_full   <- t(expr_jhu)       # 239  x genes, complete
roundtrip <- function(cohort, Xfull) {
  Xf <- Xfull[, feats, drop=FALSE]
  legacy <- as.numeric(predict(m_legacy, newx=Xf, s=lam_leg, type="response"))
  manual <- as.numeric(locked_metscore_score(Xfull, model)$prob)
  d <- abs(manual - legacy); maxd <- max(d); meand <- mean(d)
  mis <- sum((manual >= model$threshold) != (legacy >= thr_leg))
  status <- if (maxd <= PRED_TOL && mis == 0L) "PASS" else "FAIL"
  add("roundtrip", sprintf("%s n", cohort), nrow(Xfull), "-", status)
  add("roundtrip", sprintf("%s max abs prob diff", cohort), sprintf("%.3e", maxd), "<= 1e-12", status)
  add("roundtrip", sprintf("%s mean abs prob diff", cohort), sprintf("%.3e", meand), "-", status)
  add("roundtrip", sprintf("%s threshold-class mismatches", cohort), mis, "0", status)
  cat(sprintf("  [roundtrip %s] n=%d max=%.3e mean=%.3e mismatch=%d -> %s\n",
              cohort, nrow(Xfull), maxd, meand, mis, status))
  status
}
st_train <- roundtrip("TRAINING", X_train_full)
st_jhu   <- roundtrip("JHU",      X_jhu_full)

## ---- 6. helper failure-mode tests (real 1000-record matrix) --------------
tmp_empty <- tempfile("lm_empty_"); dir.create(tmp_empty)
addfm("missing_config", "load errors when config dir empty",
      load_locked_metscore(config_dir = tmp_empty))
addfm("missing_required_gene", "score errors when a locked gene column absent",
      locked_metscore_score(X_train_full[, colnames(X_train_full) != feats[1], drop=FALSE], model))
addfm("duplicated_required_gene", "score errors when a locked gene duplicated",
      locked_metscore_score(cbind(X_train_full, X_train_full[, feats[1], drop=FALSE]), model))
Xnf <- X_train_full; Xnf[1, feats[1]] <- NA_real_
addfm("nonfinite_required_value", "score errors on non-finite required value",
      locked_metscore_score(Xnf, model))

## ---- 7. config-tamper invariants (single-token edits on temp copies) -----
CFG <- model$config_dir
COEF <- "metscore_locked_v1_coefficients.csv"; META <- "metscore_locked_v1_metadata.csv"
copy_cfg <- function() {
  d <- tempfile("lmcfg_"); dir.create(d)
  file.copy(file.path(CFG, COEF), file.path(d, COEF))
  file.copy(file.path(CFG, META), file.path(d, META))
  d
}
coef_tok <- function(d, dataline) {  # value token of a data line (1 = first gene)
  ln <- readLines(file.path(d, COEF)); strsplit(ln[dataline + 1L], ",", fixed=TRUE)[[1]][2]
}
set_coef <- function(d, dataline, gene=NULL, value=NULL) {
  p <- file.path(d, COEF); ln <- readLines(p); t <- strsplit(ln[dataline + 1L], ",", fixed=TRUE)[[1]]
  if (!is.null(gene))  t[1] <- gene
  if (!is.null(value)) t[2] <- value
  ln[dataline + 1L] <- paste(t, collapse=","); writeLines(ln, p)
}
swap_coef <- function(d, a, b) {
  p <- file.path(d, COEF); ln <- readLines(p); i <- a + 1L; j <- b + 1L
  tmp <- ln[i]; ln[i] <- ln[j]; ln[j] <- tmp; writeLines(ln, p)
}
meta_val <- function(d, field) {
  ln <- readLines(file.path(d, META)); i <- grep(sprintf('^"%s",', field), ln)
  sub(sprintf('^"%s","(.*)"$', field), "\\1", ln[i])
}
set_meta <- function(d, field, value) {
  p <- file.path(d, META); ln <- readLines(p); i <- grep(sprintf('^"%s",', field), ln)
  ln[i] <- sprintf('"%s","%s"', field, value); writeLines(ln, p)
}
dup_meta <- function(d, field) {
  p <- file.path(d, META); ln <- readLines(p); i <- grep(sprintf('^"%s",', field), ln)
  writeLines(append(ln, ln[i]), p)
}
plus1e14 <- function(tok) sprintf("%.17g", as.numeric(tok) + 1e-14)
inv <- function(name, mutate) { d <- copy_cfg(); mutate(d); addinv(name, load_locked_metscore(config_dir = d)) }

inv("changed_gene_identity",  function(d) set_coef(d, 1, gene = "NOTAGENE"))
inv("changed_gene_order",     function(d) swap_coef(d, 1, 2))
inv("changed_coefficient",    function(d) set_coef(d, 1, value = sprintf("%.17g", as.numeric(coef_tok(d,1)) + 0.01)))
inv("coefficient_plus_1e14",  function(d) set_coef(d, 1, value = plus1e14(coef_tok(d, 1))))
inv("equivalent_numeric_reformat", function(d) set_coef(d, 1, value = paste0(coef_tok(d, 1), "0")))
inv("wrong_version",          function(d) set_meta(d, "version", "metscore_locked_v2"))
inv("intercept_plus_1e14",    function(d) set_meta(d, "intercept", plus1e14(meta_val(d, "intercept"))))
inv("lambda_plus_1e14",       function(d) set_meta(d, "lambda", plus1e14(meta_val(d, "lambda"))))
inv("threshold_plus_1e14",    function(d) set_meta(d, "threshold", plus1e14(meta_val(d, "threshold"))))
inv("lambda_le_zero",         function(d) set_meta(d, "lambda", "-1"))
inv("changed_lambda",         function(d) set_meta(d, "lambda", "0.5"))
inv("threshold_out_of_range", function(d) set_meta(d, "threshold", "1.5"))
inv("changed_threshold",      function(d) set_meta(d, "threshold", "0.5"))
inv("wrong_bio_sig_count",    function(d) set_meta(d, "biological_signature_count", "44"))
inv("malformed_sha",          function(d) set_meta(d, "source_artifact_sha256", "not-a-sha"))
inv("wrong_sha",              function(d) set_meta(d, "source_artifact_sha256", paste(rep("a",64), collapse="")))
inv("duplicated_metadata_field", function(d) dup_meta(d, "version"))

## ---- 8. successful-behaviour tests (complete training matrix) ------------
p_ref  <- locked_metscore_score(X_train_full[, feats, drop=FALSE], model)$prob
p_shuf <- locked_metscore_score(X_train_full[, rev(feats), drop=FALSE], model)$prob
addbeh("reordered_columns_identical", "reordered required columns -> identical probs",
       isTRUE(all.equal(unname(p_ref), unname(p_shuf), tolerance=0)),
       sprintf("max abs diff %.3e", max(abs(p_ref - p_shuf))))
p_full <- locked_metscore_score(X_train_full, model)$prob
addbeh("extra_columns_ignored", "unrelated columns ignored -> identical probs",
       isTRUE(all.equal(unname(p_ref), unname(p_full), tolerance=0)),
       sprintf("max abs diff %.3e", max(abs(p_ref - p_full))))
addbeh("rownames_preserved", "sample names preserved on output",
       identical(names(p_full), rownames(X_train_full)),
       sprintf("names identical = %s", identical(names(p_full), rownames(X_train_full))))
ext_ok <- FALSE; ext_msg <- ""; old_wd <- getwd()
tryCatch({
  ext <- tempfile("lm_ext_"); dir.create(ext); setwd(ext)
  ev <- new.env(); source(HELPER, local = ev)     # absolute path, sets ofile
  m2 <- ev$load_locked_metscore()                  # no config_dir supplied
  ext_ok  <- identical(m2$feature_names, feats)
  ext_msg <- sprintf("loaded %d genes from %s", length(m2$feature_names), m2$config_dir)
}, error=function(ee){ ext_msg <<- conditionMessage(ee) }, finally = setwd(old_wd))
addbeh("external_cwd_absolute_source", "config resolves when sourced by absolute path outside repo",
       ext_ok, ext_msg)
writeLines(c(sprintf("external-cwd absolute-path source load: %s", if (ext_ok) "PASS" else "FAIL"),
             sprintf("detail: %s", ext_msg)),
           file.path(output_dir, "..", "logs", "path_resolution.log"))

## ---- 9. Durham: recorded UNAVAILABLE -------------------------------------
durham_detail <- paste0(
  "No provenance-bound Durham post-bridge model-input artifact exists. The saved ",
  "object output/Durham/durham_metscore_batchcorrected.rda holds clin_valid ",
  "(sample-level score), not the 41-gene model-input matrix; matrix shape alone ",
  "cannot establish Durham identity or preprocessing, so no generic input is ",
  "accepted and the matrix is not reconstructed or persisted here.")
add("roundtrip", "DURHAM", "UNAVAILABLE", durham_detail, "UNAVAILABLE")
cat("  [roundtrip DURHAM] UNAVAILABLE\n")

## ---- write results + status ----------------------------------------------
write.csv(CV,  file.path(output_dir, "contract_verification.csv"), row.names=FALSE)
write.csv(FM,  file.path(output_dir, "failure_mode_tests.csv"),    row.names=FALSE)
write.csv(INV, file.path(output_dir, "invariant_tests.csv"),       row.names=FALSE)
write.csv(BEH, file.path(output_dir, "behavior_tests.csv"),        row.names=FALSE)
tracked <- c("code/utils/locked_metscore.R","code/utils/verify_locked_metscore_contract.R",
             "code/survival_analysis/Met_PCa_Survival.R",
             "code/validation/Durham_MetScore_Validation_BatchCorrected.R",
             "code/data_preparation/Calibration_LockedLR_AllCohorts.R",
             "config/metscore_locked_v1_coefficients.csv","config/metscore_locked_v1_metadata.csv",
             "config/metscore_signature_v1.csv")
fh <- vapply(tracked, function(p){ fp<-file.path(REPO,p); if (file.exists(fp)) sha256_of(fp) else "MISSING" }, character(1))
writeLines(c(sprintf("%s  %s", fh, tracked),
             sprintf("%s  %s", live_sha, "outs/metscore_logReg_genes_.rda (legacy artifact, live)")),
           file.path(output_dir, "file_sha256.txt"))

req_fail <- sum(CV$status[!(CV$check=="roundtrip" & grepl("^DURHAM", CV$detail))] == "FAIL") +
            sum(FM$status == "FAIL") + sum(INV$status == "FAIL") + sum(BEH$status == "FAIL")
cat("\n================ SUMMARY ================\n")
cat(sprintf("contract: %d PASS / %d FAIL / %d UNAVAILABLE\n",
    sum(CV$status=="PASS"), sum(CV$status=="FAIL"), sum(CV$status=="UNAVAILABLE")))
cat(sprintf("failure-mode %d/%d | invariants %d/%d | behaviour %d/%d\n",
    sum(FM$status=="PASS"), nrow(FM), sum(INV$status=="PASS"), nrow(INV),
    sum(BEH$status=="PASS"), nrow(BEH)))
cat(sprintf("round-trip: TRAINING=%s JHU=%s DURHAM=UNAVAILABLE\n", st_train, st_jhu))
if (req_fail > 0L) { cat("VERIFY_RESULT: FAIL\n"); quit(status = 1L, save = "no") }
cat("VERIFY_RESULT: PARTIAL\n")
quit(status = 2L, save = "no")

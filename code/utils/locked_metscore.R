## Loader and scorer for the frozen Met-Score binary classifier.
## Reads the tracked config in config/, checks it against the canonical values
## pinned below, and scores an already-preprocessed matrix as
##   p = plogis(intercept + X %*% beta)
## with the fixed decision threshold. It performs no fitting, feature selection,
## scaling, imputation, or recalibration.
##
## API:
##   load_locked_metscore(config_dir = NULL) -> list(feature_names, beta,
##     intercept, lambda, threshold, version, deployed_feature_count,
##     biological_signature_count, source_artifact_sha256, config_dir)
##   locked_metscore_score(x, model) -> list(prob, class, threshold, feature_names)

## ---- canonical values (must match the tracked config exactly) ------------
## 41 deployed genes, in model order.
.LM_GENES <- c("TMSB10","ENSA","ASPN","YWHAZ","HES6","STC2","ASNS","HAVCR2","F5",
               "RFTN1","SOX4","PTPN9","ALDH1A1","MRPL11","GABRD","RC3H2","CST2",
               "CXCR4","FOXH1","KIF7","BARD1","CADPS","RNF19A","CAMK2N1","GPR37",
               "KCTD14","AZGP1","PART1","CHRNA2","DPT","EDN3","LTF","SIDT1","CBLL1",
               "PTN","CCK","UFM1","CPA3","CDC42EP5","AKAP7","KLF4")
## coefficient strings, compared as exact text (order matches .LM_GENES).
.LM_COEF <- c("0.0585442650366939","0.0303321200252018","0.0353307063171115",
              "0.0397202587325218","0.0484974690199927","0.0407681758498952",
              "0.0359744747413487","0.0224268657383136","0.0241858483224475",
              "0.0444216408194979","0.0170105746186685","0.0468832359340348",
              "0.0259414068737554","0.0223652699444301","0.0521568100111632",
              "0.035243445128123","0.0418095946610025","0.0464842750709513",
              "0.0491269321447875","0.0244343282085666","0.0220146925245359",
              "0.0416804820494348","0.0180300719990522","0.067194186487649",
              "0.0262324623411712","-0.0495967151775115","-0.0869310323092144",
              "-0.0921387926175466","-0.0440260708605766","-0.0476894404300216",
              "-0.026563178655553","-0.063835182193674","-0.0333937932588824",
              "-0.0504418792189937","-0.0345302897373036","-0.0330979128065671",
              "-0.0439450141672567","-0.0437304373986636","-0.0600743006693073",
              "-0.0221268095259599","-0.046557750286702")
## metadata values, compared as exact text.
.LM_META <- c(
  version                    = "metscore_locked_v1",
  intercept                  = "-0.89165217196777258",
  lambda                     = "0.63286147959398531",
  threshold                  = "0.27035556168582586",
  deployed_feature_count     = "41",
  biological_signature_count = "45",
  source_artifact_sha256     = "713e5d9260527de3c9e5c25e7a70d5d29c0b2bfaeebe0f277a9db0f5e080c574")

## accepts optional sign, decimal, and exponent; rejects anything else.
.lm_num_ok <- function(s)
  grepl("^[+-]?(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))([eE][+-]?[0-9]+)?$", s)

## ---- config resolution ---------------------------------------------------
## Directory of this file, captured when it is sourced (source() sets ofile),
## so config resolves even when sourced by absolute path from another cwd.
.LM_SELF_DIR <- local({
  d <- NA_character_
  fr <- sys.frames()
  for (i in seq_along(fr)) {
    of <- fr[[i]]$ofile
    if (!is.null(of) && nzchar(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.na(d)) {
    ca <- commandArgs(FALSE); m <- grep("^--file=", ca, value = TRUE)
    if (length(m)) d <- dirname(normalizePath(sub("^--file=", "", m[1]), mustWork = FALSE))
  }
  d
})

.locked_config_dir <- function(config_dir = NULL) {
  if (!is.null(config_dir)) {
    if (!dir.exists(config_dir))
      stop(sprintf("locked_metscore: config_dir does not exist: %s", config_dir))
    return(normalizePath(config_dir))
  }
  cand <- character(0)
  if (!is.na(.LM_SELF_DIR))
    cand <- c(cand, file.path(dirname(dirname(.LM_SELF_DIR)), "config"))
  cand <- c(cand, file.path(getwd(), "config"), "config")
  for (cd in cand) if (dir.exists(cd)) return(normalizePath(cd))
  stop(sprintf("locked_metscore: could not locate config/ (searched: %s).",
               paste(cand, collapse = " ; ")))
}

## ---- load + validate -----------------------------------------------------
load_locked_metscore <- function(config_dir = NULL) {
  cd <- .locked_config_dir(config_dir)
  coef_path <- file.path(cd, "metscore_locked_v1_coefficients.csv")
  meta_path <- file.path(cd, "metscore_locked_v1_metadata.csv")
  if (!file.exists(coef_path))
    stop(sprintf("locked_metscore: missing coefficients config: %s", coef_path))
  if (!file.exists(meta_path))
    stop(sprintf("locked_metscore: missing metadata config: %s", meta_path))

  ## coefficients: exact gene and coefficient text, then numeric conversion
  lines <- readLines(coef_path)
  if (length(lines) < 2L) stop("locked_metscore: coefficients config is empty.")
  if (!identical(strsplit(lines[1], ",", fixed = TRUE)[[1]], c("gene", "coefficient")))
    stop("locked_metscore: coefficients header must be 'gene,coefficient'.")
  body  <- lines[-1][nzchar(trimws(lines[-1]))]
  parts <- strsplit(body, ",", fixed = TRUE)
  if (any(lengths(parts) != 2L)) stop("locked_metscore: malformed coefficient row(s).")
  genes <- vapply(parts, `[`, character(1), 1L)
  cstr  <- vapply(parts, `[`, character(1), 2L)
  if (!identical(genes, .LM_GENES))
    stop("locked_metscore: gene identity or order does not match the locked configuration.")
  if (!identical(cstr, .LM_COEF))
    stop("locked_metscore: coefficient values do not match the locked configuration.")
  if (!all(.lm_num_ok(cstr)))
    stop("locked_metscore: non-numeric coefficient text in config.")
  beta <- as.numeric(cstr); names(beta) <- genes

  ## metadata: exact text for each pinned field, then range checks
  meta <- read.csv(meta_path, colClasses = "character")
  if (!all(c("field", "value") %in% names(meta)))
    stop("locked_metscore: metadata config must have columns 'field','value'.")
  if (any(duplicated(meta$field)))
    stop("locked_metscore: duplicated metadata field(s).")
  getf <- function(k) {
    v <- meta$value[meta$field == k]
    if (length(v) != 1L)
      stop(sprintf("locked_metscore: metadata field '%s' missing or duplicated.", k))
    v
  }
  for (k in names(.LM_META))
    if (!identical(getf(k), unname(.LM_META[k])))
      stop(sprintf("locked_metscore: metadata '%s' does not match the locked configuration.", k))
  int_s <- getf("intercept"); lam_s <- getf("lambda"); thr_s <- getf("threshold")
  if (!all(.lm_num_ok(c(int_s, lam_s, thr_s))))
    stop("locked_metscore: non-numeric intercept/lambda/threshold in metadata.")
  if (!grepl("^[0-9]+$", getf("deployed_feature_count")) ||
      !grepl("^[0-9]+$", getf("biological_signature_count")))
    stop("locked_metscore: non-integer count field in metadata.")
  intercept <- as.numeric(int_s); lambda <- as.numeric(lam_s); threshold <- as.numeric(thr_s)
  if (!(lambda > 0))               stop("locked_metscore: lambda must be > 0.")
  if (threshold < 0 || threshold > 1) stop("locked_metscore: threshold must be within [0,1].")

  list(feature_names = genes, beta = beta, intercept = intercept, lambda = lambda,
       threshold = threshold, version = getf("version"),
       deployed_feature_count = as.integer(getf("deployed_feature_count")),
       biological_signature_count = as.integer(getf("biological_signature_count")),
       source_artifact_sha256 = getf("source_artifact_sha256"), config_dir = cd)
}

## ---- score ---------------------------------------------------------------
## x: samples in rows, genes in columns (glmnet newx orientation). Extra numeric
## columns are ignored; each locked gene must be present exactly once. Column
## names are required; missing or duplicated required genes and non-finite
## required values are errors.
locked_metscore_score <- function(x, model = load_locked_metscore()) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("locked_metscore: 'x' must be a matrix or data.frame (samples x genes).")
  x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "double"
  cn <- colnames(x)
  if (is.null(cn) || any(is.na(cn)) || any(!nzchar(cn)))
    stop("locked_metscore: input matrix must have non-empty column (gene) names.")
  feats <- model$feature_names
  miss <- setdiff(feats, cn)
  if (length(miss) > 0)
    stop(sprintf("locked_metscore: %d required gene(s) absent from input: %s",
                 length(miss), paste(miss, collapse = ", ")))
  dup <- feats[vapply(feats, function(g) sum(cn == g) > 1L, logical(1))]
  if (length(dup) > 0)
    stop(sprintf("locked_metscore: required gene(s) present more than once: %s",
                 paste(unique(dup), collapse = ", ")))
  X <- x[, feats, drop = FALSE]
  if (any(!is.finite(X)))
    stop("locked_metscore: non-finite value(s) in required model-input columns.")
  lp   <- as.numeric(X %*% model$beta) + model$intercept
  prob <- as.numeric(stats::plogis(lp))
  cls  <- ifelse(prob >= model$threshold, "High", "Low")
  nm <- rownames(x)
  if (!is.null(nm)) { names(prob) <- nm; names(cls) <- nm }
  list(prob = prob, class = cls, threshold = model$threshold, feature_names = feats)
}

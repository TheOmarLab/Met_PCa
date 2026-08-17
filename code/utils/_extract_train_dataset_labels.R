# ============================================================
# Extract GSM -> GSE mapping from the per-Dataset rda files,
# write to outs/train_sample_to_gse.csv. Read by
# code/signature_discovery/Met_PCa_LASSO_vs_Ridge.R to supply the
# discovery-cohort labels used in Figure S5 c/d.
#
# Each Dataset is an R list with $pheno (rownames = GSM IDs) and
# $expr (colnames = GSM IDs). We harvest both for redundancy.
# ============================================================

suppressPackageStartupMessages({
  library(here)
})

# Locate project root robustly
ROOT <- tryCatch(here::here(), error = function(e) ".")
if (!file.exists(file.path(ROOT, "data", "ProstateData.rda"))) {
  # Fallback to the working directory if here() did not resolve the root
  candidates <- c(
    getwd(),
    "."
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "data", "ProstateData.rda"))) {
      ROOT <- cand
      break
    }
  }
}
cat("ROOT =", ROOT, "\n")

# ----- DATASET → GSE mapping (from prostate_data_collection.R line 444) -----
ds2gse <- list(
  Dataset1  = "GSE116918",
  Dataset2  = "GSE55935",
  Dataset3  = "GSE51066",
  Dataset4  = "GSE46691",
  Dataset5  = "GSE41408",
  Dataset10 = "GSE70769"
)

# ----- Load ProstateData.rda (Datasets 2/3/4/5) ------------------------------
load(file.path(ROOT, "data", "ProstateData.rda"))
# Load Dataset1 + Dataset10 from their own rda files
load(file.path(ROOT, "data", "Dataset1.rda"))
load(file.path(ROOT, "data", "Dataset10.rda"))

extract_ids <- function(ds) {
  ids <- character(0)
  if (!is.null(ds$pheno)) {
    ids <- unique(c(ids, rownames(ds$pheno)))
  }
  if (!is.null(ds$expr)) {
    ids <- unique(c(ids, colnames(ds$expr)))
  }
  ids
}

rows <- list()
for (ds_name in names(ds2gse)) {
  obj <- get(ds_name)
  ids <- extract_ids(obj)
  cat(sprintf("%-10s = %-10s n_ids=%d\n", ds_name, ds2gse[[ds_name]], length(ids)))
  if (length(ids) > 0) {
    rows[[ds_name]] <- data.frame(
      sample_id = ids,
      gse       = ds2gse[[ds_name]],
      stringsAsFactors = FALSE
    )
  }
}

map_df <- do.call(rbind, rows)
out_path <- file.path(ROOT, "outs", "train_sample_to_gse.csv")
write.csv(map_df, out_path, row.names = FALSE)
cat(sprintf("\nWrote %d sample→GSE rows → %s\n", nrow(map_df), out_path))

# Quick distribution print
tab <- table(map_df$gse)
cat("\nDistribution of mapped sample IDs by GSE:\n")
print(tab)

#!/usr/bin/env bash
# setup.sh — one-time environment setup for the Met-Score project.
# Run from the project root: bash setup.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

# ---------------------------------------------------------------------------
# 1. Version checks
# ---------------------------------------------------------------------------
echo ""
echo "=== Checking versions ==="

R_VERSION=$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "NOT FOUND")
echo "R:      $R_VERSION  (required: 4.4 or 4.5)"
if [[ "$R_VERSION" == "NOT FOUND" ]]; then
  echo "  ERROR: R not found. Install from https://cran.r-project.org"
fi

PY_VERSION=$(python3 --version 2>/dev/null || echo "NOT FOUND")
echo "Python: $PY_VERSION  (required: 3.9.6, the verified analysis environment)"

# ---------------------------------------------------------------------------
# 2. Python conda environment
# ---------------------------------------------------------------------------
echo ""
echo "=== Setting up Python conda environment ==="

if ! command -v conda &>/dev/null; then
  echo "  ERROR: conda not found. Install Miniconda from https://docs.conda.io/en/latest/miniconda.html"
  echo "  Then re-run this script."
  exit 1
fi

if conda env list | grep -q "^metpca "; then
  echo "  Environment 'metpca' already exists — updating..."
  conda env update -n metpca -f "$PROJECT_ROOT/environment.yml" --prune
else
  echo "  Creating environment 'metpca'..."
  conda env create -f "$PROJECT_ROOT/environment.yml"
fi

echo "  Python environment ready."
echo "  Activate with: conda activate metpca"

# ---------------------------------------------------------------------------
# 3. R dependencies
# ---------------------------------------------------------------------------
echo ""
echo "=== Installing R dependencies ==="

Rscript - <<'EOF'
pkgs_cran <- c(
  # Data wrangling
  "data.table", "dplyr", "tidyverse", "reshape", "readr", "purrr",
  "forcats", "broom", "scales", "xtable", "readxl", "here", "matrixStats",
  # Layout
  "gridExtra",
  # Statistics / modelling
  "caret", "glmnet", "mltools", "metafor", "sampling",
  # Survival
  "survival", "survminer", "pROC", "precrec", "timeROC",
  "survRM2", "cmprsk", "rms", "effsize", "splines",
  # Visualisation
  "ggplot2", "ggpubr", "pheatmap", "RColorBrewer",
  "ragg", "patchwork", "plotROC", "cowplot",
  # Decision curve / risk / concordance
  "riskRegression", "pec", "prodlim",
  # Single-cell
  "Seurat", "Matrix",
  # Mixed models (AXL script)
  "lme4", "lmerTest",
  # Functional enrichment (functional_enrichment_MetScore.R)
  "msigdbr", "digest",
  # Excel output (Met_Score_Gene_Consistency.R)
  "openxlsx"
)

pkgs_bioc <- c(
  "GEOquery", "Biobase", "limma", "edgeR", "genefilter",
  "fgsea", "DESeq2",
  "ComplexHeatmap", "circlize",
  "GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db", "AnnotationDbi",
  "survcomp"
)

missing_cran <- pkgs_cran[!pkgs_cran %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- pkgs_bioc[!pkgs_bioc %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# MetaIntegrator is archived from CRAN — install its dependencies first, then the package itself
if (!requireNamespace("MetaIntegrator", quietly = TRUE)) {
  message("Installing MetaIntegrator dependencies...")

  # Step A: Bioconductor deps (includes COCONUT's own deps: sva, BiocParallel)
  BiocManager::install(
    c("multtest", "GEOmetadb", "sva", "BiocParallel", "limma"),
    ask = FALSE, update = FALSE
  )

  # Step B: CRAN deps for the MetaIntegrator chain (manhattanly before MetaIntegrator)
  meta_cran_deps <- c("pander", "Rmisc", "pracma", "Metrics", "manhattanly", "DT", "HGNChelper", "gtools")
  missing_meta <- meta_cran_deps[!meta_cran_deps %in% rownames(installed.packages())]
  if (length(missing_meta) > 0)
    install.packages(missing_meta, repos = "https://cloud.r-project.org")

  # Step C: RMySQL — needs system MySQL headers, best-effort only
  if (!requireNamespace("RMySQL", quietly = TRUE)) {
    tryCatch(
      install.packages("RMySQL", repos = "https://cloud.r-project.org"),
      error = function(e) message("  Skipping RMySQL (MySQL headers not found — not required for analysis)")
    )
  }

  # Step D: COCONUT — archived (latest: 1.0.2); needs sva/BiocParallel/limma (done in Step A)
  if (!requireNamespace("COCONUT", quietly = TRUE)) {
    message("Installing COCONUT from CRAN archive...")
    install.packages(
      "https://cran.r-project.org/src/contrib/Archive/COCONUT/COCONUT_1.0.2.tar.gz",
      repos = NULL, type = "source"
    )
    if (!requireNamespace("COCONUT", quietly = TRUE))
      stop("COCONUT failed to install. Check that sva and BiocParallel are available.")
  }

  # Step D2: manhattanly — archived (latest: 0.3.0); in MetaIntegrator Suggests
  if (!requireNamespace("manhattanly", quietly = TRUE)) {
    message("Installing manhattanly from CRAN archive...")
    install.packages(
      "https://cran.r-project.org/src/contrib/Archive/manhattanly/manhattanly_0.3.0.tar.gz",
      repos = NULL, type = "source"
    )
  }

  message("Installing MetaIntegrator from CRAN archive...")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/MetaIntegrator/MetaIntegrator_2.1.3.tar.gz",
    repos = NULL, type = "source"
  )
}

# ESTIMATE — tumour purity for Purity_ESTIMATE.R (Figure S1a); hosted on R-Forge, not CRAN
if (!requireNamespace("estimate", quietly = TRUE)) {
  message("Installing estimate from R-Forge...")
  install.packages("estimate", repos = "http://R-Forge.R-project.org")
  if (!requireNamespace("estimate", quietly = TRUE))
    message("  WARNING: estimate failed to install (only needed for Purity_ESTIMATE.R / Figure S1a)")
}

message("R packages ready.")
EOF

# ---------------------------------------------------------------------------
# 4. Create output subdirectories that scripts expect
# ---------------------------------------------------------------------------
# Scripts locate the project root from their own file path (or the MET_PCA_ROOT
# environment variable); setup does not modify any tracked source.
echo ""
echo "=== Creating output directories ==="
mkdir -p "$PROJECT_ROOT/outs/Decipher"
mkdir -p "$PROJECT_ROOT/outs/Figure3"
mkdir -p "$PROJECT_ROOT/outs/FigureS3"
mkdir -p "$PROJECT_ROOT/outs/FigureS5"
mkdir -p "$PROJECT_ROOT/outs/FigureS6"
mkdir -p "$PROJECT_ROOT/outs/comparison"
mkdir -p "$PROJECT_ROOT/outs/kfoury"
mkdir -p "$PROJECT_ROOT/figures/scRNAseq/kfoury"
mkdir -p "$PROJECT_ROOT/figures/Durham"
mkdir -p "$PROJECT_ROOT/figures/survival"
mkdir -p "$PROJECT_ROOT/figures/time_dependent"
mkdir -p "$PROJECT_ROOT/figures/DCA"
mkdir -p "$PROJECT_ROOT/outs/functional_enrichment"
mkdir -p "$PROJECT_ROOT/figures/functional_enrichment"
echo "Output directories created."

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup complete ==="
echo ""
echo "Next steps:"
echo "  1. Open Met_PCa.Rproj in RStudio (or run Rscript from the project root)."
echo "  2. Follow the execution order in README.md, starting with:"
echo "       Phase 1:  code/data_preparation/prostate_data_collection.R"
echo "                 code/data_preparation/process_GSE268308_GSE268309.R"
echo "       Phase 2:  code/signature_discovery/MET_PCa_MetaIntegrator.R"
echo "       Phase 3:  code/survival_analysis/Met_PCa_Survival.R   *** central script ***"
echo "       Phase 4:  code/validation/Durham_MetScore_Validation_BatchCorrected.R"
echo "       Phase 5:  python3 code/decipher/decipher_canonical_stage1_durham_full.py"
echo "                 python3 code/decipher/decipher_canonical.py"
echo "       Phase 6:  python3 code/figures/Figure2_unified_validation.py"
echo "                 code/ancillary/MetScore_Sensitivity.R                  # Figure S6 aggregates"
echo "                 python3 code/figures/FigureS6_sensitivity.py            # Figure S6 renderer"
echo "                 code/survival_analysis/Met_PCa_Survival_Multivariate.R"
echo "                 python3 code/figures/Figure3_GS7.py"
echo "                 Rscript code/ancillary/Met_PCa_Survival_DECIPHER.R   # Figure S3 aggregates"
echo "                 python3 code/figures/FigureS3_Decipher.py            # Figure S3 renderer"

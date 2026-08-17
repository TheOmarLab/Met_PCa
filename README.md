# MET_PCa: Met-Score for Prostate Cancer Metastasis Risk

This repository includes the analysis scripts used to produce the results reported in the manuscript entitled **“A Conserved Metastatic Competence Signature from Primary Prostate Tumors”**.

## Overview

The code covers the discovery, implementation, and validation of the locked Met-Score classifier, and generation of the manuscript figures and supplementary tables.

## Repository structure

- `code/` — analysis and figure-generation scripts.
- `config/` — locked Met-Score and Decipher-surrogate configuration files.
- `setup.sh` and `environment.yml` — software setup.

## Setup

```
bash setup.sh
conda activate metpca
```

Python dependencies are listed in `environment.yml`, the R dependencies are
installed by `setup.sh`, and analyses should be run from the repository root.

## Main analysis scripts

- `code/signature_discovery/MET_PCa_MetaIntegrator.R` — signature discovery and Figure 1.
- `code/survival_analysis/Met_PCa_Survival.R` — frozen scoring and the JHU analysis object.
- `code/validation/Durham_MetScore_Validation_BatchCorrected.R` — Durham validation.
- `code/data_preparation/Calibration_LockedLR_AllCohorts.R` and
  `code/figures/Figure2_unified_validation.py` — Figure 2.
- `code/survival_analysis/Met_PCa_Survival_Multivariate.R` and
  `code/figures/Figure3_GS7.py` — Figure 3 and multivariable results.
- `code/figures/Figure4_Biopsy_Validation.R` — Figure 4.
- `code/figures/FigureS1_purity_clinical_utility.py` through
  `code/figures/FigureS6_sensitivity.py` — supplementary figures.
- `code/supplementary/generate_supplementary_tables.py` — supplementary tables.

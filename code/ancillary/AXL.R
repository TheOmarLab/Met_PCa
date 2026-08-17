#######################################################################
## Bulk RNA-seq: GSE281461 – AXL inhibition in bone PDX
## Goal: Does AXL blockade suppress Met-Score?
##
## Score definition: MetScore here is the MetaIntegrator calculateScore over the
## 45-gene Met-Score signature (filter$posGeneNames / filter$negGeneNames), a
## directional module z-score (positive-gene z-scores minus negative-gene
## z-scores). It is not the frozen 41-feature ridge-logistic Met-Score classifier.
#######################################################################

rm(list = ls())

## Packages
library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(MetaIntegrator)
library(edgeR)

## Load Met-Score objects (original signature)
load("./outs/PP_filter_MetaScore.rda")          ## -> filter
load("./outs/filtersiggenes_MetaScore.rda")     ## -> Filter_SignatureGenes

## Rebuild Positive / Negative gene lists from filter
PositiveGenes <- filter$posGeneNames
NegativeGenes <- filter$negGeneNames

#######################################################################
## 1. Download / load GSE281461
#######################################################################

## If you already saved expression + pheno as .rda, you can load those here
## and skip getGEO. Otherwise:

gse281_list <- GEOquery::getGEO("GSE281461", GSEMatrix = TRUE)
length(gse281_list)

## Use the first (or appropriate) ExpressionSet – adjust index if needed
gse281 <- gse281_list[[1]]
pheno281_raw <- Biobase::pData(gse281)

#######
# get expression data
# This downloads *all* supplementary files for the series
getGEOSuppFiles("GSE281461", makeDirectory = FALSE)
list.files(pattern = "GSE281461")
counts281 <- read.delim(
  gzfile("GSE281461_gene_counts.xls.gz"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(counts281)

gene_symbol_col <- 31L
sample_cols     <- 2:30

expr281_counts <- as.matrix(counts281[, sample_cols])
rownames(expr281_counts) <- counts281[[gene_symbol_col]]

# Drop rows without a symbol or duplicated symbols (keeps first occurrence)
valid <- !is.na(rownames(expr281_counts)) & rownames(expr281_counts) != ""
expr281_counts <- expr281_counts[valid, , drop = FALSE]
expr281_counts <- expr281_counts[!duplicated(rownames(expr281_counts)), , drop = FALSE]

dim(expr281_counts)
head(rownames(expr281_counts))
head(colnames(expr281_counts))

#######################################################################
## 2. Build clean phenotype table
##    *** YOU MUST ADAPT THE grepl() PATTERNS TO MATCH GSE281461 META ***
#######################################################################

pheno281 <- pheno281_raw %>%
  mutate(
    SampleID       = geo_accession,   # GSM...
    SampleTitle    = title,           # "147_AXL 1", etc.
    PDX_line       = `cell line:ch1`,
    Treatment_raw  = `treatment:ch1`,
    Treatment_group = case_when(
      grepl("batiraxcept", Treatment_raw, ignore.case = TRUE) ~ "AXL",
      grepl("vehicle",     Treatment_raw, ignore.case = TRUE) ~ "Vehicle",
      TRUE ~ NA_character_
    ),
    PDX_family = case_when(
      grepl("LuCaP 147CR", PDX_line, ignore.case = TRUE) ~ "LuCaP_147CR",
      grepl("LuCaP 147",   PDX_line, ignore.case = TRUE) ~ "LuCaP_147",
      grepl("LuCaP 35CR",  PDX_line, ignore.case = TRUE) ~ "LuCaP_35CR",
      grepl("LuCaP 35",    PDX_line, ignore.case = TRUE) ~ "LuCaP_35",
      TRUE ~ PDX_line
    )
  )

# Keep only AXL vs Vehicle for now
pheno281_axl <- pheno281 %>%
  filter(Treatment_group %in% c("AXL", "Vehicle"))

## Keep only AXL vs vehicle
pheno281_axl <- pheno281 %>%
  filter(!is.na(Treatment_group))

table(pheno281_axl$PDX_family, pheno281_axl$Treatment_group)


# Align expression matrix to phenotype 
# Samples present in both expression and pheno
common_titles <- intersect(colnames(expr281_counts), pheno281_axl$SampleTitle)
length(common_titles)
common_titles

# Subset expression
expr281_counts_sub <- expr281_counts[, common_titles, drop = FALSE]

# Subset and reorder pheno
pheno281_axl_sub <- pheno281_axl %>%
  filter(SampleTitle %in% common_titles) %>%
  arrange(match(SampleTitle, common_titles))

# Reorder expression to match pheno row order
expr281_counts_sub <- expr281_counts_sub[, pheno281_axl_sub$SampleTitle, drop = FALSE]

stopifnot(all(colnames(expr281_counts_sub) == pheno281_axl_sub$SampleTitle))

#######################################################################
## 3. Normalize expression
#######################################################################

# Filter very lowly expressed genes (optional but sensible)
keep <- rowSums(expr281_counts_sub > 1) >= 2
expr281_counts_filt <- expr281_counts_sub[keep, ]

# Library-size–normalized logCPM
expr281_logcpm <- cpm(expr281_counts_filt, log = TRUE, prior.count = 1)

dim(expr281_logcpm)
range(expr281_logcpm)


## Keep only genes present in the Met-Score signature
sig_genes <- intersect(Filter_SignatureGenes, rownames(expr281_logcpm))
length(sig_genes)

expr281_sig <- expr281_logcpm[sig_genes, , drop = FALSE]
dim(expr281_sig)


## Match
colnames(expr281_counts_sub)
rownames(pheno281_axl_sub)

## Build a pheno object that matches expr columns exactly
pheno281_axl_matched <- pheno281_axl_sub[
  match(colnames(expr281_sig), pheno281_axl_sub$SampleTitle),
]

# Sanity check: no NAs in match
stopifnot(!any(is.na(rownames(pheno281_axl_matched))))

# Now set rownames to match expr column names
rownames(pheno281_axl_matched) <- colnames(expr281_sig)

# Final consistency check
stopifnot(all(colnames(expr281_sig) == rownames(pheno281_axl_matched)))

#######################################################################
## 4. Compute Met-Score on GSE281461 bone PDXs
#######################################################################

Dataset_GSE281461 <- list(
  expr          = expr281_sig,
  pheno         = pheno281_axl_matched,
  keys          = rownames(expr281_sig),
  formattedName = "GSE281461_AXL_vs_Vehicle"
)


MetScore_vec <- calculateScore(
  filterObject  = filter,
  datasetObject = Dataset_GSE281461
)

pheno281_axl_matched$MetScore <- MetScore_vec

## Save for downstream / survival etc.
save(pheno281_axl_matched, expr281_sig, file = "./outs/GSE281461_bonePDX_MetScore.rda")



#######################################################################
## Statistical tests: does Met-Score drop with AXL inhibition?
## Pre-specified direction: AXL inhibition is expected to LOWER Met-Score
## (Vehicle > AXL). The tests below are therefore one-sided by design: the
## Wilcoxon uses alternative = "greater" (Vehicle > AXL) and the mixed-model
## contrast uses a one-sided z-test on the treatment coefficient (beta < 0).
#######################################################################

pheno_clean <- pheno281_axl_matched %>%
  mutate(
    Treatment_clean = case_when(
      grepl("vehicle", Treatment_raw,  ignore.case = TRUE) ~ "Vehicle",
      grepl("batiraxcept", Treatment_raw, ignore.case = TRUE) &
        !grepl("doc", Treatment_raw, ignore.case = TRUE)     ~ "AXL",   # monotherapy only
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment_clean))

## Drop PDX families that don’t have both Vehicle and AXL
fam_tab <- pheno_clean %>%
  count(PDX_family, Treatment_clean) %>%
  pivot_wider(names_from = Treatment_clean, values_from = n, values_fill = 0)

keep_fams <- fam_tab %>%
  filter(Vehicle > 0, AXL > 0) %>%
  pull(PDX_family)

pheno_clean <- pheno_clean %>%
  filter(PDX_family %in% keep_fams) %>%
  mutate(
    Treatment_clean = factor(Treatment_clean, levels = c("Vehicle", "AXL"))
  )

table(pheno_clean$PDX_family, pheno_clean$Treatment_clean)

#######################################################################
## Figure: Met-Score in Control vs AXL inhibitor
#######################################################################

cols_trt <- c("Vehicle" = "#2b2eb5", "AXL" = "#ed6905")

p_plot <- ggplot(pheno_clean,
       aes(x = Treatment_clean, y = MetScore, fill = Treatment_clean)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.9, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1.6, alpha = 0.7) +
  facet_wrap(~ PDX_family, nrow = 1) +
  labs(x = NULL, y = "Met-Score") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 16, face = "bold"),
    axis.text       = element_text(size = 14, face = "bold"),
    strip.text      = element_text(size = 14, face = "bold")
  )

ggsave("./figures/Fig_AXL_GSE281461_MetScore_bonePDX.pdf",
       plot = p_plot, width = 6, height = 4.6, useDingbats = FALSE)


by(pheno_clean, pheno_clean$PDX_family, function(df) {
  wilcox.test(
    MetScore ~ Treatment_clean,
    data = df,
    exact = FALSE,
    alternative = "greater"   # Vehicle > AXL if you expect AXL ↓ Met-Score
  )
})


summ_by_fam <- pheno_clean %>%
  dplyr::group_by(PDX_family, Treatment_clean) %>%
  dplyr::summarise(
    mean_MetScore = mean(MetScore),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = Treatment_clean,
    values_from = mean_MetScore
  ) %>%
  dplyr::mutate(
    diff_AXL_vs_Vehicle = AXL - Vehicle
  )

summ_by_fam


ggplot(summ_by_fam, aes(x = PDX_family, y = diff_AXL_vs_Vehicle)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(aes(xend = PDX_family, y = 0, yend = diff_AXL_vs_Vehicle)) +
  geom_point(size = 3) +
  labs(y = "Δ Met-Score (AXL − Vehicle)", x = NULL) +
  theme_classic(base_size = 14)


library(lme4)
library(lmerTest)

m1 <- lmer(MetScore ~ Treatment_clean + (1 | PDX_family),
           data = pheno_clean)

summary(m1)
confint(m1, parm = "Treatment_cleanAXL", method = "Wald")


beta     <- fixef(m1)["Treatment_cleanAXL"]
se_beta  <- sqrt(vcov(m1)["Treatment_cleanAXL", "Treatment_cleanAXL"])
z        <- beta / se_beta
# One-sided by the pre-specified direction (AXL lowers Met-Score, so beta < 0):
# pnorm(z) is the lower-tail p-value for the AXL-vs-Vehicle coefficient.
p_one_sided <- pnorm(z)

beta
p_one_sided



#######################################################################
## Supplementary figure: AXL inhibition lowers Met-Score in LuCaP_147CR
#######################################################################

## 1) Subset to LuCaP_147CR only
pheno_147cr <- pheno_clean %>%
  dplyr::filter(PDX_family == "LuCaP_147CR")

## 2) One-sided Wilcoxon test (Vehicle > AXL, expecting AXL ↓ Met-Score)
wil_147cr <- wilcox.test(
  MetScore ~ Treatment_clean,
  data        = pheno_147cr,
  exact       = FALSE,
  alternative = "greater"
)

wil_147cr

p_label <- paste0("Wilcoxon one-sided P = ",
                  formatC(wil_147cr$p.value, format = "e", digits = 2))

## 3) Colors as before
cols_trt <- c("Vehicle" = "#2b2eb5", "AXL" = "#ed6905")

## 4) Decide y-position for annotation
y_max <- max(pheno_147cr$MetScore, na.rm = TRUE)

p_147cr <- ggplot(pheno_147cr,
                  aes(x = Treatment_clean,
                      y = MetScore,
                      fill = Treatment_clean)) +
  # optional violin for shape
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.9, color = NA) +
  geom_boxplot(width = 0.20,
               outlier.shape = NA,
               alpha = 0.9,
               color = "black") +
  geom_jitter(width = 0.08,
              size  = 2.0,
              alpha = 0.8,
              color = "black") +
  scale_fill_manual(values = cols_trt) +
  labs(x = NULL, y = "Met-Score",
       title = "LuCaP 147 CR bone PDX: Vehicle vs AXL inhibitor") +
  annotate("text",
           x = 1.8,
           y = y_max + 0.3,
           label = p_label,
           size = 3.8,
           fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12, face = "bold"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin     = margin(5.5, 15, 5.5, 5.5)
  )

p_147cr

ggsave("./figures/SuppFig_AXL_LuCaP147CR_MetScore.pdf",
       plot = p_147cr,
       width = 6, height = 6, useDingbats = FALSE)

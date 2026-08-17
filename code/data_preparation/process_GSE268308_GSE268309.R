############################################################################
## Clean work space
rm(list = ls())

## Load  packages
library(GEOquery)
library(data.table)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

## robust GEO downloads (libcurl; longer timeout for large series/GPL files)
options(download.file.method = "libcurl", timeout = 600)

######################################################################################################
# get the esets from GEO
######################################################################################################

# GEOquery's bundled httr2 downloader is unreliable here; fetch each series
# matrix with libcurl and parse the local file. Only sample phenotype is needed
# (getGPL = FALSE); gene annotation is done downstream via org.Hs.eg.db.
gse308_sm <- tempfile(fileext = ".txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268308/matrix/GSE268308_series_matrix.txt.gz",
              gse308_sm, method = "libcurl", mode = "wb", quiet = TRUE)
GSE268308 <- getGEO(filename = gse308_sm, getGPL = FALSE)


gse309_sm <- tempfile(fileext = ".txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268309/matrix/GSE268309_series_matrix.txt.gz",
              gse309_sm, method = "libcurl", mode = "wb", quiet = TRUE)
GSE268309 <- getGEO(filename = gse309_sm, getGPL = FALSE)


######################################################################################################
# work on expression data
######################################################################################################

####################################
# GSE268308
###################################

# get the raw expression for GSE268308
counts_GSE268308 <- fread("./data/GSE268308_FullTable_RawCounts_Discovery_dataset.txt.gz")


# Gene symbol annotation
ensembl_ids_counts_GSE268308 <- counts_GSE268308$V1
count_mat_GSE268308 <- as.matrix(counts_GSE268308[, -1, with = FALSE])
rownames(count_mat_GSE268308) <- ensembl_ids_counts_GSE268308

# Clean Ensembl IDs
ensembl_clean_GSE268308 <- sub("\\.\\d+$", "", rownames(count_mat_GSE268308))

symbols_GSE268308 <- mapIds(org.Hs.eg.db,
                      keys = ensembl_clean_GSE268308,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")


# Keep only rows with a mapped symbol
valid_idx_GSE268308 <- !is.na(symbols_GSE268308)
counts_GSE268308_sym <- count_mat_GSE268308[valid_idx_GSE268308, ]
symbols_GSE268308_nonNA <- symbols_GSE268308[valid_idx_GSE268308]

# Set rownames to symbols
rownames(counts_GSE268308_sym) <- symbols_GSE268308_nonNA

# Collapse duplicated symbols by SUMMING counts
counts_GSE268308_sym <- rowsum(counts_GSE268308_sym,
                         group = rownames(counts_GSE268308_sym),
                         reorder = FALSE)

# mapped gene-symbol universe before expression filtering (coverage provenance)
mapped_symbols_GSE268308 <- rownames(counts_GSE268308_sym)

# normalization
dge_GSE268308 <- DGEList(counts = counts_GSE268308_sym)

keep_GSE268308 <- filterByExpr(dge_GSE268308)
dge_GSE268308 <- dge_GSE268308[keep_GSE268308, ]

dge_GSE268308 <- calcNormFactors(dge_GSE268308, method = "TMM")

expr_GSE268308 <- cpm(dge_GSE268308, log = TRUE, prior.count = 1)

range(expr_GSE268308)

####################################
# GSE268309
###################################

# get the raw expression for GSE268309.
# NOTE: the raw-counts file has a malformed (off-by-one) header: its first row holds
# the 32 sample names (S_01..S_32) but the data rows carry 33 columns (a leading
# Ensembl-gene-id column + 32 sample columns). A naive fread() therefore shifts every
# sample label by one and orphans the 32nd sample under the auto-name "V33", which
# later breaks the pheno match. We read the header and data separately so column 1 is
# the gene id and the 32 sample columns keep their correct names.
gse309_hdr <- strsplit(
  readLines("./data/GSE268309_FullTable_RawCounts_Validation_dataset.txt.gz", n = 1),
  "\t")[[1]]
gse309_sample_names <- gse309_hdr[gse309_hdr != ""]           # 32 sample names
counts_GSE268309 <- fread("./data/GSE268309_FullTable_RawCounts_Validation_dataset.txt.gz",
                          header = FALSE, skip = 1)
# one gene-id column + exactly 32 nonempty, unique sample columns
stopifnot(length(gse309_sample_names) == 32L,
          all(nzchar(gse309_sample_names)),
          !anyDuplicated(gse309_sample_names),
          ncol(counts_GSE268309) == length(gse309_sample_names) + 1L)


# Gene symbol annotation
ensembl_ids_counts_GSE268309 <- counts_GSE268309[[1]]
count_mat_GSE268309 <- as.matrix(counts_GSE268309[, -1, with = FALSE])
rownames(count_mat_GSE268309) <- ensembl_ids_counts_GSE268309
colnames(count_mat_GSE268309) <- gse309_sample_names

# Clean Ensembl IDs
ensembl_clean_GSE268309 <- sub("\\.\\d+$", "", rownames(count_mat_GSE268309))

symbols_GSE268309 <- mapIds(org.Hs.eg.db,
                            keys = ensembl_clean_GSE268309,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")


# Keep only rows with a mapped symbol
valid_idx_GSE268309 <- !is.na(symbols_GSE268309)
counts_GSE268309_sym <- count_mat_GSE268309[valid_idx_GSE268309, ]
symbols_GSE268309_nonNA <- symbols_GSE268309[valid_idx_GSE268309]

# Set rownames to symbols
rownames(counts_GSE268309_sym) <- symbols_GSE268309_nonNA

# Collapse duplicated symbols by SUMMING counts
counts_GSE268309_sym <- rowsum(counts_GSE268309_sym,
                               group = rownames(counts_GSE268309_sym),
                               reorder = FALSE)

# mapped gene-symbol universe before expression filtering (coverage provenance)
mapped_symbols_GSE268309 <- rownames(counts_GSE268309_sym)

# normalization
dge_GSE268309 <- DGEList(counts = counts_GSE268309_sym)

keep_GSE268309 <- filterByExpr(dge_GSE268309)
dge_GSE268309 <- dge_GSE268309[keep_GSE268309, ]

dge_GSE268309 <- calcNormFactors(dge_GSE268309, method = "TMM")

expr_GSE268309 <- cpm(dge_GSE268309, log = TRUE, prior.count = 1)

range(expr_GSE268309)



######################################################################################################
# work on Pheno data
######################################################################################################

############################
# process pheno for GSE268308
############################

pheno_GSE268308 <- pData(GSE268308)

idx_GSE268308 <- match(colnames(expr_GSE268308), pheno_GSE268308$description)
stopifnot(!any(is.na(idx_GSE268308)))
# reorder pheno to match counts columns
pheno_GSE268308 <- pheno_GSE268308[idx_GSE268308, ]

# enforce rownames(pheno) == colnames(counts)
rownames(pheno_GSE268308) <- colnames(expr_GSE268308)

# check alignment
stopifnot(all(rownames(pheno_GSE268308) == colnames(expr_GSE268308)))


table(pheno_GSE268308$`cell line:ch1`)

pheno_GSE268308$Metastasis <- ifelse( pheno_GSE268308$`cell line:ch1` == "metastatic hormone naive prostate cancer", "Mets", "No_Mets")

pheno_GSE268308$Metastasis <- factor(pheno_GSE268308$Metastasis,
                                     levels = c("No_Mets", "Mets"))


table(pheno_GSE268308$Metastasis)

group_GSE268308 <- pheno_GSE268308$Metastasis
# structural check: n=78, 47 localized (No_Mets), 31 mHNPC (Mets)
stopifnot(ncol(expr_GSE268308) == 78L,
          sum(group_GSE268308 == "No_Mets") == 47L,
          sum(group_GSE268308 == "Mets") == 31L)


############################
# process pheno for GSE268309
############################

pheno_GSE268309 <- pData(GSE268309)

idx_GSE268309 <- match(colnames(expr_GSE268309), pheno_GSE268309$description)
stopifnot(!any(is.na(idx_GSE268309)))
# reorder pheno to match counts columns
pheno_GSE268309 <- pheno_GSE268309[idx_GSE268309, ]

# enforce rownames(pheno) == colnames(counts)
rownames(pheno_GSE268309) <- colnames(expr_GSE268309)

# check alignment
stopifnot(all(rownames(pheno_GSE268309) == colnames(expr_GSE268309)))


table(pheno_GSE268309$`disease:ch1`)

pheno_GSE268309$Metastasis <- ifelse( pheno_GSE268309$`disease:ch1` == "metastatic hormone naive prostate cancer", "Mets", "No_Mets")

pheno_GSE268309$Metastasis <- factor(pheno_GSE268309$Metastasis,
                                     levels = c("No_Mets", "Mets"))


table(pheno_GSE268309$Metastasis)

group_GSE268309 <- pheno_GSE268309$Metastasis
# structural check: n=32, 17 localized (No_Mets), 15 mHNPC (Mets)
stopifnot(ncol(expr_GSE268309) == 32L,
          sum(group_GSE268309 == "No_Mets") == 17L,
          sum(group_GSE268309 == "Mets") == 15L)


######################################################################################################
# save: expression matrices, mapped-symbol universe, minimal phenotype, named group
######################################################################################################
# name group vectors by sample and build minimal (Metastasis-only) phenotype
names(group_GSE268308) <- colnames(expr_GSE268308)
names(group_GSE268309) <- colnames(expr_GSE268309)
pheno_GSE268308 <- data.frame(Metastasis = group_GSE268308,
                              row.names = colnames(expr_GSE268308))
pheno_GSE268309 <- data.frame(Metastasis = group_GSE268309,
                              row.names = colnames(expr_GSE268309))

# fail-closed: nonempty, unique names identical across expr / phenotype / group
for (cn in list(colnames(expr_GSE268308), colnames(expr_GSE268309)))
  stopifnot(length(cn) > 0L, all(nzchar(cn)), !anyDuplicated(cn))
stopifnot(identical(colnames(expr_GSE268308), rownames(pheno_GSE268308)),
          identical(colnames(expr_GSE268308), names(group_GSE268308)),
          identical(colnames(expr_GSE268309), rownames(pheno_GSE268309)),
          identical(colnames(expr_GSE268309), names(group_GSE268309)))

save(expr_GSE268308, mapped_symbols_GSE268308, pheno_GSE268308, group_GSE268308,
     expr_GSE268309, mapped_symbols_GSE268309, pheno_GSE268309, group_GSE268309,
     file = "./data/GSE268308_GSE268309.rda")




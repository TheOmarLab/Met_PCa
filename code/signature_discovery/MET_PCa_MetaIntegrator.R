############################################################################
#code for the filtered gene signature - Meta Score
## Clean work space
rm(list = ls())

## Load  packages
library(MetaIntegrator)
library(GEOquery)
library(caret)
library(genefilter)
library(mltools)
library(pROC)
library(data.table)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ragg)
library(scales)
library(metafor)

###########################################################################
#set path

## Load training data sets 
load("./data/ProstateData.rda")
load("./data/Dataset1.rda")
load("./data/Dataset10.rda")

## Load Testing data set
# This is JHU
load("./data/Dataset7.rda")


##########################################################################

## Getting the phenotype data for each data set
#Training
pheno1 <- Dataset1$pheno
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno
pheno5 <- Dataset5$pheno
pheno10 <- Dataset10$pheno
#JHU for testing
pheno7 <- Dataset7$pheno

#########################################################################

##Gene expression matrix for each data set
#Training
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr
expr5 <- Dataset5$expr
expr10 <- Dataset10$expr
#JHU for testing
expr7 <- Dataset7$expr


###################################################################
## Checking if the expression data are normalized and log2 transformed
# boxplot(expr3[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)

###################################################################

## Annotate expression
#log 2 

## Expr1
# head(rownames(expr1))
# rownames(expr1) <- Dataset1$keys
# expr1 <- expr1[!is.na(rownames(expr1)), ]
# 
# rownames(expr1) <- gsub("\\,.+", "", rownames(expr1))
# rownames(expr1) <- gsub("\\-.+", "", rownames(expr1))
# 
# expr1 <- aggregate(expr1[,], list(Gene = rownames(expr1)), FUN = mean)
# rownames(expr1) <- expr1$Gene
# expr1$Gene <- NULL
# expr1 <- as.matrix(expr1)
# dim(expr1)
# X1 <- expr1
# ffun <- filterfun(pOverA(p = 0.5, A = 50))
# filt1 <- genefilter(2^X1, ffun)
# expr1 <- expr1[filt1, ]
# dim(expr1)

#expr1 <- t(scale(t(expr1), center = T, scale = T))
#####################
## expr2
#head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)

X2 <- expr2
ffun <- filterfun(cv(a = 0.1, b = 10))
filt2 <- genefilter(2^X2, ffun)
table(filt2)
expr2 <- expr2[filt2, ]
dim(expr2)

#expr2 <- t(scale(t(expr2), center = T, scale = T))

#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)

X3 <- expr3
#ffun <- filterfun(cv(a = 0.1, b = 1))
filt3 <- genefilter(2^X3, ffun)
table(filt3)
expr3 <- expr3[filt3, ]

#expr3 <- t(scale(t(expr3), center = T, scale = T))

# #######################
dim(expr4)
X4 <- expr4
filt4 <- genefilter(2^X4, ffun)
table(filt4)
expr4 <- expr4[filt4, ]
# 
#expr4 <- t(scale(t(expr4), center = T, scale = T))

######################
## expr5
head(rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)
# 
X5 <- expr5
filt5 <- genefilter(2^X5, ffun)
table(filt5)
expr5 <- expr5[filt5, ]

#expr5 <- t(scale(t(expr5), center = T, scale = T))

#####################
## Expr7
# Processed
#range(expr7)
#expr7 <- expr7+2

# dim(expr7)
# X7 <- expr7
# filt7 <- genefilter(2^X7, ffun)
# expr7 <- expr7[filt7, ]

expr7 <- expr7+2

#expr7 <- t(scale(t(expr7), center = T, scale = T))


#######################

########################
## Expr10 (Processed)
# rownames(expr10) <- Dataset10$keys
# summary(is.na(rownames(expr10)))
# 
# expr10 <- expr10[!is.na(rownames(expr10)), ]
dim(expr10)
sel = which(apply(expr10, 1, function(x) all(is.finite(x)) ))
expr10 <- expr10[sel,]

dim(expr10)
X10 <- expr10
filt10 <- genefilter(2^X10, ffun)
table(filt10)
expr10 <- expr10[filt10, ]

#expr10 <- t(scale(t(expr10), center = T, scale = T))

####################################################################
##Adding Mets/No_Mets to eaach data set
pheno1$Metastasis <- pheno1$`met event (1=yes, 0=no):ch1`
pheno1$Metastasis[pheno1$Metastasis == 0] <- "No_Mets"
pheno1$Metastasis[pheno1$Metastasis == 1] <- "Mets"
pheno1$Metastasis <- as.factor(pheno1$Metastasis)
table(pheno1$Metastasis)
all(rownames(pheno1) == colnames(expr1))

Dataset1$pheno <- pheno1
Dataset1$expr <- expr1
Dataset1$keys <- rownames(expr1)


## Modify pheno2
# Remove cell lines
pheno2 <- pheno2[-c(47:54), ]
pheno2$Metastasis <- pheno2$`lymph node metastasis status:ch1`
pheno2$Metastasis[pheno2$Metastasis == 0] <- "No_Mets"
pheno2$Metastasis[pheno2$Metastasis == 1] <- "Mets"
pheno2 <- pheno2[!(pheno2$Metastasis == "NA"), ]

pheno2$Metastasis <- as.factor(pheno2$Metastasis)
table(pheno2$Metastasis)
# Modify expr2
expr2 <- expr2[,colnames(expr2) %in% rownames(pheno2)]
all(rownames(pheno2) == colnames(expr2))
## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset2$expr <- expr2
Dataset2$pheno <- pheno2

####### 
## Modify pheno4
pheno4$Metastasis <- pheno4$`metastatic event:ch1`
pheno4$Metastasis[pheno4$Metastasis == 0] <- "No_Mets"
pheno4$Metastasis[pheno4$Metastasis == 1] <- "Mets"
pheno4 <- pheno4[!(pheno4$Metastasis == "NA"), ]

pheno4$Metastasis <- as.factor(pheno4$Metastasis)
table(pheno4$Metastasis)
levels(pheno4$Metastasis)

## Modify sample names to match sample names of pheno4
head(colnames(expr4))
colnames(expr4) <- gsub(".CEL", "", colnames(expr4))
colnames(expr4) <- gsub(".+\\.", "", colnames(expr4))
expr4 <- expr4[,order(colnames(expr4))]

rownames(pheno4) <- pheno4$description
head(rownames(pheno4))
rownames(pheno4) <- gsub(".+\\.", "", rownames(pheno4))
pheno4 <- pheno4[order(rownames(pheno4)), ]

all(rownames(pheno4) == colnames(expr4))

##Replace the expression and phenotype data in the dataset with the new modified versions
Dataset4$expr <- expr4
Dataset4$pheno <- pheno4

######
# Modify pheno3
pheno3$Metastasis <- pheno3$`metastatic event:ch1`
pheno3$Metastasis[pheno3$Metastasis == 0] <- "No_Mets"
pheno3$Metastasis[pheno3$Metastasis == 1] <- "Mets"

pheno3$Metastasis <- as.factor(pheno3$Metastasis)
table(pheno3$Metastasis)
levels(pheno3$Metastasis)

all(rownames(pheno3) == colnames(expr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset3$pheno <- pheno3
Dataset3$expr <- expr3

#####
## Modify Pheno5
table(pheno5$`development of metastasis:ch1`)

pheno5$Metastasis <- pheno5$`development of metastasis:ch1`

pheno5$Metastasis[pheno5$Metastasis == "no"] <- "No_Mets"
pheno5$Metastasis[pheno5$Metastasis == "yes"] <- "Mets"

pheno5$Metastasis <- as.factor(pheno5$Metastasis)
table(pheno5$Metastasis)
all(rownames(pheno5) == colnames(expr5))

Dataset5$expr <- expr5
Dataset5$pheno <- pheno5


########
# Pheno7
# Processed
pheno7$Metastasis <- as.factor(pheno7$Metastasis) 

Dataset7$expr <- expr7
Dataset7$pheno <- pheno7
Dataset7$keys <- rownames(expr7)
#########


###########
## Pheno10
# Processed
#NOTE: we lose a good amount of samples due to lack of Mets information
pheno10$Metastasis <- pheno10$`pathology stage:ch1`
pheno10$Metastasis <- gsub("T.+M", "M", pheno10$Metastasis)
pheno10$Metastasis <- gsub("N.+", "", pheno10$Metastasis)

table(pheno10$Metastasis)
pheno10 <- pheno10[!(pheno10$Metastasis == "Mx"), ]
pheno10 <- pheno10[!(pheno10$Metastasis == "pMx"), ]
pheno10 <- pheno10[!(pheno10$Metastasis == "U"), ]
pheno10$Metastasis <- as.factor(pheno10$Metastasis)
levels(pheno10$Metastasis) <- c("No_Mets", "No_Mets", "Mets", "Mets")
table(pheno10$Metastasis)

expr10 <- expr10[, colnames(expr10) %in% rownames(pheno10)]
all(rownames(pheno10) == colnames(expr10))

Dataset10$expr <- expr10
Dataset10$pheno <- pheno10
Dataset10$keys <- rownames(expr10)

############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘met’ or 0 for ‘no met’)
Dataset1 <- classFunction(Dataset1, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset2 <- classFunction(Dataset2, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset3 <- classFunction(Dataset3, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset4 <- classFunction(Dataset4, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset5 <- classFunction(Dataset5, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset10 <- classFunction(Dataset10, column = "Metastasis", diseaseTerms = c("Mets"))
#JHU
Dataset7 <- classFunction(Dataset7, column = "Metastasis", diseaseTerms = c("Mets"))

#########################################################################
## Metaanalysis

## Creating the metaobject
AllDataSets <- list(Dataset1, Dataset2, Dataset3, Dataset4, Dataset5, Dataset10)
names(AllDataSets) <- c(Dataset1$formattedName, Dataset2$formattedName,
                        Dataset3$formattedName, Dataset4$formattedName,
                        Dataset5$formattedName, Dataset10$formattedName)


Prostate_meta <- list()
Prostate_meta$originalData <- AllDataSets

## Replace keys within each data set
Prostate_meta$originalData$GSE116918$keys <- rownames(expr1)
Prostate_meta$originalData$GSE55935$keys <- rownames(expr2)
Prostate_meta$originalData$GSE51066$keys <- rownames(expr3)
Prostate_meta$originalData$GSE46691$keys <- rownames(expr4)
Prostate_meta$originalData$GSE41408$keys <- rownames(expr5)
Prostate_meta$originalData$GSE70769$keys <- rownames(expr10)


## Check the meta object before the metaanalysis
#checkDataObject(Prostate_meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis

Prostate_meta <- geneSymbolCorrection(Prostate_meta)

## Run the meta analysis algorithm
Prostate_metaanalysis <- runMetaAnalysis(Prostate_meta, 
                                         runLeaveOneOutAnalysis = T, maxCores = 3)

## Filter out significant genes from the metaanalysis results 
#this will be the gene signature that separates Mets from No_Mets
Prostate_metaanalysis <- filterGenes(Prostate_metaanalysis, isLeaveOneOut = F, 
                                     effectSizeThresh = 0.2, FDRThresh = 0.05, 
                                     numberStudiesThresh = 4)
## Assigning a name to the filter
filter <- Prostate_metaanalysis$filterResult$FDR0.05_es0.2_nStudies4_looaFALSE_hetero0
filter

#list of upregulated and downregulated genes in the signature
PositiveGenes <- filter$posGeneNames
NegativeGenes <- filter$negGeneNames
#save(PositiveGenes, NegativeGenes, file = "./Objs/PP_SigGenes_Pre.rda")

## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, 
                                         getMostRecentFilter(Prostate_metaanalysis))

## Save the filter gene signature
Filter_SignatureGenes <- c(PositiveGenes, NegativeGenes)
save(Filter_SignatureGenes, file = "./outs/filtersiggenes_MetaScore.rda")

save(filter, file = "./outs/PP_filter_MetaScore.rda")

## Save the genome-wide meta-analysis object (per-gene pooled effect sizes/FDR).
## Consumed by Met_Score_Gene_Consistency.R and the functional-enrichment GSEA ranking.
meta_analysis_results <- Prostate_metaanalysis
save(meta_analysis_results, file = "./outs/meta_analysis_results.rda")

## Save as a table
#write.table(filter_summary$pos, file = "./Objs/Meta/Positive_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
#write.table(filter_summary$neg, file = "./Objs/Meta/Negative_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")

write.table(PositiveGenes, file="./outs/MetScore_pos_genes.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(NegativeGenes, file="./outs/MetScore_neg_genes.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#############################################################
## Filter the gene signature for more accuracy and AUC
#############################################################
# Using forward search 
# New_filter <- forwardSearch(metaObject = Prostate_metaanalysis, filterObject = filter)
# 
# ## Replace the old filter with the new smaller one
# Prostate_metaanalysis$filterResults$FDR0.05_es0.2_nStudies4_looaFALSE_hetero0$posGeneNames <- New_filter$posGeneNames
# Prostate_metaanalysis$filterResults$FDR0.05_es0.2_nStudies4_looaFALSE_hetero0$negGeneNames <- New_filter$negGeneNames
# 
# New_filter <- Prostate_metaanalysis$filterResults[[1]]
# New_filter_summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, getMostRecentFilter(Prostate_metaanalysis))
# 
# ## Save the tables of positive and negative genes
# #write.table(New_filter_summary$pos, file = "./outs/NewFilter_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)
# #write.table(New_filter_summary$neg, file = "./outs/NewFilter_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)
# 
# ## Gene names
# PositiveGenes2 <- New_filter$posGeneNames
# PositiveGenes2
# NegativeGenes2 <- New_filter$negGeneNames
# NegativeGenes2
# 
# Filter_SignatureGenes2 <- c(PositiveGenes2, NegativeGenes2)
# save(Filter_SignatureGenes2, file = "./outs/filtersiggenes_MetaScore2.rda")
# 
# save(New_filter, file = "./outs/PP_filter_MetaScore2.rda")

##########################################################################################################################

## Create a summary ROC curve of the training sets
set.seed(333)
png(filename = "./figures/Fig1C_training_cohort_ROC.png",
    width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = Prostate_metaanalysis, filterObject = filter)
dev.off()

##########################################################################################################################

###############################################################################
# Figure 1 — Met-Score gene expression heatmap
###############################################################################

sig <- c(filter$posGeneNames, filter$negGeneNames)

es_mat <- Prostate_metaanalysis$metaAnalysis$datasetEffectSizes  # genes x studies
es_sig <- es_mat[sig, , drop = FALSE]

# order genes: positives then negatives
es_sig <- es_sig[c(filter$posGeneNames, filter$negGeneNames), , drop = FALSE]

row_split <- factor(
  c(rep("Up in Mets", length(filter$posGeneNames)),
    rep("Down in Mets", length(filter$negGeneNames))),
  levels = c("Up in Mets", "Down in Mets")
)

# symmetric color scale around 0
mx <- max(abs(es_sig), na.rm = TRUE)
col_fun <- colorRamp2(c(-mx, 0, mx), c("#2166ac", "white", "#b2182b"))

ht_es <- Heatmap(
  es_sig,
  name = "Hedges' g",
  col = col_fun,
  na_col = "grey95",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  row_split = row_split,
  row_gap = unit(3, "mm"),
  
  # =====================
  # Font harmonization
  # =====================
  row_title_gp     = gpar(fontsize = 12, fontface = "bold"),
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 10),
  column_names_rot = 45,
  column_title     = "Training cohorts",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

pdf("./figures/Fig1A_effectSize_heatmap.pdf",
    width = 5.5,
    height = 0.16 * nrow(es_sig) + 2)  # scale height to gene count
draw(ht_es, heatmap_legend_side = "right")
dev.off()




##########################################
# Meta effect-size lollipop (pooled ES per gene)
meta_tab <- as.data.frame(Prostate_metaanalysis$metaAnalysis$pooledResults)
meta_tab$gene <- rownames(meta_tab)

sig <- c(filter$posGeneNames, filter$negGeneNames)

meta_sig <- meta_tab %>%
  filter(gene %in% sig) %>%
  mutate(direction = ifelse(gene %in% filter$posGeneNames, "Up in Mets", "Down in Mets")) %>%
  arrange(factor(direction, levels=c("Up in Mets","Down in Mets")),
          desc(effectSize))

meta_sig$gene <- factor(meta_sig$gene, levels = rev(meta_sig$gene))

# lock gene order for plotting (top->bottom after flip)
meta_sig$gene <- factor(meta_sig$gene, levels = rev(meta_sig$gene))

cols_dir <- c("Up in Mets"   = "#d7301f",
              "Down in Mets" = "#2b8cbe")

p_lolli <- ggplot(meta_sig, aes(x = gene, y = effectSize, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey65") +
  geom_pointrange(
    aes(ymin = effectSize - 1.96 * effectSizeStandardError,
        ymax = effectSize + 1.96 * effectSizeStandardError),
    linewidth = 0.45, fatten = 2.2
  ) +
  coord_flip(clip = "off") +
  scale_color_manual(values = cols_dir) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 9),
    
    legend.position = c(0.72, 0.12),             # Inside bottom-right
    legend.direction = "vertical",
    legend.box.background = element_rect(fill="white", color=NA),
    legend.background = element_rect(fill="white", color=NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    
    plot.margin  = margin(5.5, 5.5, 5.5, 5.5)     # Remove huge right margin
  ) + 
  
  labs(x = NULL, y = "Pooled effect size (Hedges' g)",
       title = "Met-Score pooled effect sizes")

p_lolli

# save with height scaling to number of genes -> no overlap
ggsave("./figures/Fig1B_metaEffectSizes_lollipop.pdf",
       plot = p_lolli,
       width = 4,
       height = 6,  
       limitsize = FALSE)





###############################################################################
# Internal cross-study validation (LOSO: leave-one-study-out)
###############################################################################
# Output:
#   - ./outs/LOSO_internal_validation_AUC.csv
#   - ./figures/LOSO/LOSO_AUC_forest.pdf/.tiff
#   - ./figures/LOSO/LOSO_ROC_<heldout>.pdf  (optional per-cohort ROC)


################################################################
# pooled AUC across held-out cohorts
################################################################
dir.create("./figures/LOSO", showWarnings = FALSE, recursive = TRUE)

# ---- helper: get most recent filter robustly
extract_filter_obj <- function(ma_obj) {
  # Newer MetaIntegrator versions store in $filterResult (a named list)
  if (!is.null(ma_obj$filterResult)) {
    fr <- ma_obj$filterResult
    if (is.list(fr) && length(fr) >= 1) {
      # usually one element; take first
      f1 <- fr[[1]]
      if (is.list(f1) && !is.null(f1$posGeneNames) && !is.null(f1$negGeneNames)) return(f1)
    }
  }
  
  # Some versions store in $filterResults
  if (!is.null(ma_obj$filterResults)) {
    fr <- ma_obj$filterResults
    if (is.list(fr) && length(fr) >= 1) {
      f1 <- fr[[1]]
      if (is.list(f1) && !is.null(f1$posGeneNames) && !is.null(f1$negGeneNames)) return(f1)
    }
  }
  
  stop("Could not extract filter object (posGeneNames/negGeneNames) from MetaIntegrator output.")
}

extract_filter_obj <- function(ma_obj,
                               filter_name = "FDR0.05_es0.2_nStudies4_looaFALSE_hetero0") {
  
  fr <- ma_obj$filterResult
  if (is.null(fr)) fr <- ma_obj$filterResults
  if (is.null(fr)) stop("No filterResult/filterResults found.")
  
  # If the intended filter exists, use it
  if (!is.null(fr[[filter_name]])) return(fr[[filter_name]])
  
  # Otherwise: find any filter that matches your thresholds (fallback)
  # (useful if MetaIntegrator changes naming slightly)
  nm <- names(fr)
  hit <- grep("FDR0\\.05.*es0\\.2.*nStudies4", nm, value = TRUE)
  if (length(hit) >= 1) return(fr[[hit[1]]])
  
  # Last resort: take first element but warn
  warning("Could not find expected filter name; falling back to first filter element.")
  return(fr[[1]])
}

# ---- helper: compute AUC + DeLong CI
auc_ci <- function(y, score) {
  df <- data.frame(y=y, score=score) %>% filter(!is.na(y), !is.na(score))
  # MetaIntegrator dataset$class is usually 0/1 numeric
  df$y <- as.numeric(df$y)
  
  if (length(unique(df$y)) < 2) {
    return(list(auc=NA, ci_low=NA, ci_high=NA, roc=NULL))
  }
  
  r <- pROC::roc(response=df$y, predictor=df$score, levels=c(0,1), direction="<", quiet=TRUE)
  list(
    auc = as.numeric(pROC::auc(r)),
    ci_low = as.numeric(pROC::ci.auc(r))[1],
    ci_high = as.numeric(pROC::ci.auc(r))[3],
    roc = r
  )
}


# ---- MAIN LOSO
run_LOSO <- function(ds_all,
                     effectSizeThresh = 0.2,
                     FDRThresh = 0.05,
                     numberStudiesThresh = 4,
                     maxCores = 3,
                     save_rocs = TRUE) {
  
  # Ensure datasets named
  if (is.null(names(ds_all)) || any(names(ds_all) == "")) {
    names(ds_all) <- vapply(ds_all, function(d) d$formattedName, character(1))
  }
  
  res <- lapply(names(ds_all), function(heldout_name) {
    
    message("LOSO held-out cohort: ", heldout_name)
    
    train_list <- ds_all[names(ds_all) != heldout_name]
    test_ds    <- ds_all[[heldout_name]]

    # MetaIntegrator requires keys aligned to expression row names in every split.
    train_list <- lapply(train_list, function(d) {
      d$keys <- rownames(d$expr)
      d
    })
    test_ds$keys <- rownames(test_ds$expr)

    # meta-object
    meta_obj <- list(originalData = train_list)

    # Harmonize gene symbols across the five training cohorts.
    meta_obj <- geneSymbolCorrection(meta_obj)
    
    # run meta-analysis (no internal LOO here; LOSO is external)
    ma <- runMetaAnalysis(meta_obj, runLeaveOneOutAnalysis = FALSE, maxCores = maxCores)
    
    # filter genes with your thresholds
    ma <- filterGenes(ma,
                      isLeaveOneOut = FALSE,
                      effectSizeThresh = effectSizeThresh,
                      FDRThresh = FDRThresh,
                      numberStudiesThresh = min(numberStudiesThresh, length(train_list)))
    
    filt <- extract_filter_obj(ma)
    
    sig <- c(filt$posGeneNames, filt$negGeneNames)
    if (length(sig) == 0) {
      warning("No genes selected for held-out: ", heldout_name)
      return(data.frame(
        heldout = heldout_name,
        n = length(test_ds$class),
        n_event = sum(test_ds$class == 1, na.rm=TRUE),
        n_pos = 0, n_neg = 0,
        auc = NA, ci_low = NA, ci_high = NA
      ))
    }
    
    # score held-out cohort using the filter learned on the OTHER 5 cohorts
    score <- calculateScore(filterObject = filt, datasetObject = test_ds)
    a <- auc_ci(test_ds$class, score)
    
    # optional: save ROC plot for each held-out cohort
    if (save_rocs && !is.null(a$roc)) {
      roc_df <- data.frame(
        FPR = 1 - a$roc$specificities,
        TPR = a$roc$sensitivities
      )
      
      p <- ggplot(roc_df, aes(FPR, TPR)) +
        geom_line(linewidth=1.2, color="#2b2eb5") +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
        coord_equal() +
        annotate("text", x=0.55, y=0.1,
                 label=paste0("AUC = ", sprintf("%.2f", a$auc),
                              "\n95% CI ", sprintf("%.2f", a$ci_low), "–", sprintf("%.2f", a$ci_high)),
                 size=5, hjust=0) +
        labs(
          title = paste0("LOSO ROC (held-out: ", heldout_name, ")"),
          x = "False Positive Rate (1–Specificity)",
          y = "True Positive Rate (Sensitivity)"
        ) +
        theme_classic(base_size=14) +
        theme(
          plot.title = element_text(face="bold", hjust=0.5),
          axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold")
        )
      
      ggsave(file.path("./figures/LOSO", paste0("LOSO_ROC_", heldout_name, ".pdf")),
             p, width=4.8, height=4.8, useDingbats=FALSE)
    }
    
    data.frame(
      heldout = heldout_name,
      n = length(test_ds$class),
      n_event = sum(test_ds$class == 1, na.rm=TRUE),
      n_pos = length(filt$posGeneNames),
      n_neg = length(filt$negGeneNames),
      auc = a$auc,
      ci_low = a$ci_low,
      ci_high = a$ci_high
    )
  })
  
  bind_rows(res)
}

#################
# ---- Run it
#################
# --- define ds_all from existing objects ---
ds_all <- AllDataSets

# ensure names exist
if (is.null(names(ds_all)) || any(names(ds_all) == "")) {
  names(ds_all) <- vapply(ds_all, function(d) d$formattedName, character(1))
}

# sanity check
names(ds_all)

sapply(ds_all, function(d) length(d$class))

set.seed(333)
loso_df <- run_LOSO(ds_all,
                    effectSizeThresh = 0.2,
                    FDRThresh = 0.05,
                    numberStudiesThresh = 4,
                    maxCores = 3,
                    save_rocs = TRUE)

write.csv(loso_df, "./outs/LOSO_internal_validation_AUC.csv", row.names=FALSE)
print(loso_df)

# Forest-style plot
p_loso <- ggplot(loso_df, aes(x=reorder(heldout, auc), y=auc)) +
  geom_point(size=3, color="#2b2eb5") +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=0.15, color="#2b2eb5") +
  coord_flip() +
  labs(x="Held-out cohort", y="AUC (DeLong 95% CI)",
       title="Internal cross-study validation (LOSO)") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.title=element_text(face="bold"),
        axis.text=element_text(face="bold"))

ggsave("./figures/LOSO/LOSO_AUC_forest.pdf", p_loso,
       width=6.5, height=4.5, useDingbats=FALSE)
ggsave("./figures/LOSO/LOSO_AUC_forest.tiff", p_loso,
       width=6.5, height=4.5, dpi=450, compression="lzw")


##########################################################################################################################
### TESTING
##########################################################################################################################

###################################################################################################
# process jhu
###################################################################################################

#Rename testing data set
jhu <- Dataset7

# Acess the phenotype and expression matrix of the testing data set
pheno_jhu <- pheno7
expr_jhu <- expr7

## Modify testing data set to match the testing
rownames(expr_jhu) <- jhu$keys
summary(is.na(rownames(expr_jhu)))
expr_jhu <- expr_jhu[!is.na(rownames(expr_jhu)), ]
dim(expr_jhu)


## Modify val_pheno
pheno_jhu$Metastasis <- pheno_jhu$met
pheno_jhu$Metastasis[pheno_jhu$Metastasis == 0] <- "No_Mets"
pheno_jhu$Metastasis[pheno_jhu$Metastasis == 1] <- "Mets"
pheno_jhu$Metastasis <- as.factor(pheno_jhu$Metastasis)
table(pheno_jhu$Metastasis)
all(rownames(pheno_jhu) == colnames(expr_jhu))

jhu$pheno <- pheno_jhu
jhu$expr <- expr_jhu
jhu$keys <- rownames(expr_jhu)

## Label the samples
jhu <- classFunction(jhu, column = "Metastasis", diseaseTerms = c("Mets"))

# ROC curve
png(filename = "./figures/ROC_Test_JHU.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = jhu, filterObject = filter, title = "JHU Nat. History")
dev.off()


#### Calculate a signature score (Z score) and add it to the phenotype table
pheno_jhu$score <- calculateScore(filterObject = filter, datasetObject = jhu)

## save val_pheno1 with z-score for futher survival analysis
save(pheno_jhu, PositiveGenes, NegativeGenes, file = "./outs/jhu_pheno_filter_MetaScorer_Zscore.rda")


############## 
# plot ROC
##############

df_jhu <- pheno_jhu %>% 
  dplyr::select(Metastasis, score) %>%
  filter(!is.na(Metastasis), !is.na(score))

# ROC
roc_obj <- pROC::roc(
  response  = df_jhu$Metastasis,
  predictor = df_jhu$score,
  levels    = c("No_Mets","Mets"),
  direction = "<"
)

auc_val <- as.numeric(pROC::auc(roc_obj))
ci_auc  <- as.numeric(pROC::ci.auc(roc_obj))

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

fs <- 14  

p_jhu_roc <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth=1.2, color="#2b2eb5") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  coord_equal() +
  annotate("text", x=0.55, y=0.1,
           label=paste0(
             "AUC = ", sprintf("%.2f", auc_val), "\n",
             "95% CI ", sprintf("%.2f", ci_auc[1]), "–", sprintf("%.2f", ci_auc[3])
           ),
           size=5, hjust=0) +
  labs(
    title = "",
    x = "False Positive Rate (1–Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_classic(base_size = fs) +
  theme(
    axis.title  = element_text(face="bold", size=fs),
    axis.text   = element_text(face="bold", size=fs),
    plot.title  = element_text(face="bold", hjust=0.5, size=fs)
  )

ggsave("./figures/JHU_binary_ROC_diagnostic.pdf",
       plot=p_jhu_roc, width=4.8, height=4.8, useDingbats=FALSE)
ggsave("./figures/JHU_binary_ROC_diagnostic.tiff",
       plot=p_jhu_roc, width=4.8, height=4.8, dpi=400, compression="lzw")


###############
# Met-Score Violin Plot
################
cols2 <- c("No_Mets"="#2b2eb5", "Mets"="#ed6905")
group_col <- "Metastasis"
score_col <- "score"
df_jhu[[group_col]] <- as.character(df_jhu[[group_col]])

df_jhu[[group_col]] <- dplyr::recode(
  df_jhu[[group_col]],
  "No_Mets" = "No Metastasis",
  "no_mets" = "No Metastasis",
  "Localized" = "No Metastasis",
  "Mets" = "Metastasis",
  "mets" = "Metastasis",
  "mHNPC" = "Metastasis"
)

# enforce left → right order in plot
df_jhu[[group_col]] <- factor(df_jhu[[group_col]],
                              levels = c("No Metastasis", "Metastasis"))

p_jhu_violin <- ggplot(df_jhu, aes_string(x=group_col, y=score_col, fill=group_col)) +
  geom_violin(trim=FALSE, alpha=0.6, width=0.9, color=NA) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.9) +
  geom_jitter(width=0.08, size=1.6, alpha=0.7) +
  stat_compare_means(
    method="wilcox.test",
    label="p.format",
    size=5,
    label.y = max(df_jhu[[score_col]]) + 0.35
  ) +
  scale_fill_manual(values = c("No Metastasis"="#2b2eb5",
                               "Metastasis"="#ed6905")) +
  labs(x=NULL, y="Met-Score (z)", title="") +
  theme_classic(base_size=16) +
  theme(
    legend.position="none",
    axis.title.y = element_text(size=16, face="bold"),
    axis.text   = element_text(size=16, face="bold")
  )

ggsave("./figures/JHU_MetScore_distribution_diagnostic.pdf",
       plot=p_jhu_violin, width=4.6, height=4.6, useDingbats=FALSE)
ggsave("./figures/JHU_MetScore_distribution_diagnostic.tiff",
       plot=p_jhu_violin, width=4.6, height=4.6, dpi=400, compression="lzw")


##########################
# Calibration plot for JHU cohort (binary diagnostic)
##########################
# Build logistic regression using Met-Score only
calib_df <- df_jhu %>%
  mutate(
    Metastasis_bin = ifelse(Metastasis == "Mets", 1, 0)
  )

table(calib_df$Metastasis_bin, calib_df$score > median(calib_df$score))

calib_mod <- glm(Metastasis_bin ~ score, data = calib_df, family = "binomial")

# Predicted probabilities
calib_df$pred_prob <- predict(calib_mod, type = "response")

# Bin into deciles of predicted risk
calib_df$bin <- ntile(calib_df$pred_prob, 10)

calib_summary <- calib_df %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(
    mean_pred = mean(pred_prob),
    mean_obs  = mean(Metastasis_bin),
    n = dplyr::n()
  )

# Calibration plot
p_calib <- ggplot(calib_summary, aes(x = mean_pred, y = mean_obs)) +
  geom_point(size = 3, color = "#2b2eb5") +
  geom_line(linewidth = 1, color = "#2b2eb5") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", linewidth = 1, color = "gray50") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Predicted probability of metastasis",
    y = "Observed metastasis rate",
    title = ""
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 16, face = "bold")
  )

# Save
ggsave("./figures/JHU_binary_calibration_diagnostic.pdf",
       plot = p_calib, width = 4.6, height = 4.6, useDingbats = FALSE)
ggsave("./figures/JHU_binary_calibration_diagnostic.tiff",
       plot = p_calib, width = 4.6, height = 4.6, dpi = 400, compression = "lzw")



########################
# observed effect sizes per gene in JHU
#########################
missing <- setdiff(sig, rownames(expr_jhu))
missing

sig_jhu <- intersect(sig, rownames(expr_jhu))

df_jhu <- pheno_jhu %>%
  dplyr::select(Metastasis, score = score)

# Convert to No_Mets / Mets coding for consistency
grp <- pheno_jhu$Metastasis

# For each gene, compute effect size
compute_es <- function(gene) {
  x <- as.numeric(expr_jhu[gene, ])
  g <- grp
  g0 <- x[g == "No_Mets"]
  g1 <- x[g == "Mets"]
  
  # Hedges' g
  n0 <- length(g0); n1 <- length(g1)
  s0 <- sd(g0);     s1 <- sd(g1)
  sp <- sqrt(((n0-1)*s0^2 + (n1-1)*s1^2)/(n0+n1-2))
  g_raw <- (mean(g1) - mean(g0)) / sp
  J <- 1 - 3/(4*(n0+n1)-9)
  g_hedges <- g_raw * J
  
  # Standard error approximation
  se <- sqrt( (n0+n1)/(n0*n1) + (g_hedges^2)/(2*(n0+n1)) )
  
  return(data.frame(
    gene = gene,
    effectSize = g_hedges,
    effectSizeSE = se
  ))
}

# Apply to all signature genes
es_list <- lapply(sig_jhu, compute_es)
es_tab  <- do.call(rbind, es_list)

# Add direction
es_tab$direction <- ifelse(es_tab$gene %in% filter$posGeneNames,
                           "Up in Mets", "Down in Mets")

# Order genes: Up first, then Down
es_tab <- es_tab %>%
  arrange(factor(direction, levels = c("Up in Mets", "Down in Mets")),
          desc(effectSize))

# Correct factor order for plotting
es_tab$gene <- factor(es_tab$gene, levels = rev(es_tab$gene))

cols_dir <- c("Up in Mets"   = "#d7301f",
              "Down in Mets" = "#2b8cbe")

p_jhu_lolli <- ggplot(es_tab,
                      aes(x = gene, y = effectSize, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.5, color = "grey65") +
  geom_pointrange(
    aes(ymin = effectSize - 1.96 * effectSizeSE,
        ymax = effectSize + 1.96 * effectSizeSE),
    linewidth = 0.45, fatten = 2.2
  ) +
  coord_flip(clip = "off") +
  scale_color_manual(values = cols_dir) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 9),
    legend.title = element_blank(),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Observed effect size (Hedges' g)",
    title = "JHU Cohort: Observed Effect Sizes per Gene"
  )

ggsave("./figures/JHU_effectSizes_lollipop.pdf",
       plot = p_jhu_lolli,
       width = 4, height = 6, useDingbats = FALSE)

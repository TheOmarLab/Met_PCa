############################################################################
# Builds the canonical scored JHU survival object (outs/coxdata.rda).
# Applies the frozen Met-Score classifier after training-reference quantile
# mapping; performs no classifier fitting and no outcome analysis.
############################################################################
rm(list = ls())
suppressPackageStartupMessages(library(limma))   # normalizeBetweenArrays (reference mapping)

# ---- authorized inputs ---------------------------------------------------
load("./data/Dataset7.rda")                        # JHU phenotype
pheno_jhu <- Dataset7$pheno
load("./outs/MetastasisData_JHUOut.rda")           # trainMat/testMat + trainGroup/testGroup

# ---- training-reference quantile mapping ---------------------------------
# Quantile-normalize the training matrix, then map each JHU gene onto the
# training quantile distribution so every cohort is scored on one per-gene
# scale. This is reference mapping, not model fitting.
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
qn_to_train <- function(test_mat, train_ref) {
  common <- intersect(rownames(test_mat), rownames(train_ref))
  out <- test_mat
  for (g in common) {
    train_vals <- train_ref[g, ]
    test_ranks <- rank(test_mat[g, ], ties.method = "average") / (ncol(test_mat) + 1)
    out[g, ]   <- quantile(train_vals, probs = test_ranks, names = FALSE, type = 7)
  }
  out
}
expr_jhu <- qn_to_train(testMat, usedTrainMat)

# ---- frozen Met-Score contract -------------------------------------------
# The loader validates the tracked contract (41 genes, frozen order, coefficient
# text, metadata SHA, threshold); scoring applies it without any refit.
source("code/utils/locked_metscore.R")
.locked_model    <- load_locked_metscore()
Meta_Score       <- .locked_model$feature_names    # 41 deployed genes, frozen order
LOCKED_THRESHOLD <- .locked_model$threshold
cat(sprintf("Loaded frozen Met-Score contract (version=%s, threshold=%.6f, %d genes)\n",
            .locked_model$version, LOCKED_THRESHOLD, length(Meta_Score)))

# ---- score the complete training and JHU matrices ------------------------
x_train           <- t(usedTrainMat)
Train_prob_logReg <- as.numeric(locked_metscore_score(x_train, .locked_model)$prob)
x_jhu             <- t(expr_jhu)
jhu_prob_logReg   <- as.numeric(locked_metscore_score(x_jhu, .locked_model)$prob)

# ---- apply the frozen threshold ------------------------------------------
jhu_predClasses_logReg <- ifelse(jhu_prob_logReg >= LOCKED_THRESHOLD, "1", "0")
jhu_predClasses_logReg <- factor(jhu_predClasses_logReg, levels = c("0", "1"))

# ---- build the canonical CoxData_jhu object ------------------------------
pheno_jhu$pathgs <- factor(pheno_jhu$pathgs, levels = c("7", "8", "9"))
CoxData_jhu <- data.frame(pheno_jhu)
CoxData_jhu$jhu_prob_logReg <- jhu_prob_logReg
CoxData_jhu$jhu_predClasses_logReg <- jhu_predClasses_logReg
colnames(CoxData_jhu)[colnames(CoxData_jhu) == "jhu_prob_logReg"] <- "Met-Score prob"
CoxData_jhu$MetScoreClass <- factor(
  CoxData_jhu$jhu_predClasses_logReg,
  levels = c("0", "1"),
  labels = c("Low risk", "High risk"))
colnames(CoxData_jhu)[colnames(CoxData_jhu) == "pathgs"] <- "Pathological GS"

# ---- invariants ----------------------------------------------------------
stopifnot(
  # row identity: JHU phenotype rows align 1:1 and in order with scored columns
  identical(rownames(pheno_jhu), colnames(expr_jhu)),
  # sample + event ledger
  nrow(CoxData_jhu) == 239L,
  sum(as.integer(as.character(CoxData_jhu$met)) == 1L, na.rm = TRUE) == 93L,
  # feature order: frozen 41-gene contract, all present in the scored matrices
  identical(Meta_Score, .locked_model$feature_names),
  all(Meta_Score %in% rownames(expr_jhu)),
  all(Meta_Score %in% rownames(usedTrainMat)),
  # score: valid probabilities from the frozen contract (no refit)
  length(jhu_prob_logReg) == nrow(CoxData_jhu),
  all(is.finite(jhu_prob_logReg)), all(jhu_prob_logReg >= 0 & jhu_prob_logReg <= 1),
  all(is.finite(Train_prob_logReg)),
  # class: frozen-threshold labels and counts
  identical(levels(CoxData_jhu$MetScoreClass), c("Low risk", "High risk")),
  identical(as.integer(table(CoxData_jhu$MetScoreClass)[c("Low risk", "High risk")]), c(101L, 138L))
)

save(CoxData_jhu, file = "./outs/coxdata.rda")
cat(sprintf("Wrote outs/coxdata.rda  (CoxData_jhu: n=%d, cols=%d, mets=%d, MetScoreClass Low/High=%d/%d)\n",
            nrow(CoxData_jhu), ncol(CoxData_jhu),
            sum(as.integer(as.character(CoxData_jhu$met)) == 1L, na.rm = TRUE),
            sum(CoxData_jhu$MetScoreClass == "Low risk"),
            sum(CoxData_jhu$MetScoreClass == "High risk")))

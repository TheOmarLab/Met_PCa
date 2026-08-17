###############################################################################
# functional_enrichment_MetScore.R
#
# Functional context for the Met-Score metastasis signature. Two distinct,
# complementary analyses run from the genome-wide discovery meta object.
#
# (A) ORA of the frozen 45-gene signature against the 50 MSigDB Hallmarks.
#     Direct annotation of the selected signature: which Hallmarks the POS
#     (up in metastasis) and NEG (down) genes fall into. Right-tail
#     hypergeometric against the assayed discovery genes annotated in Hallmark,
#     computed for every Hallmark so zero-overlap terms are reported too.
#     This is annotation of the signature, not a positive enrichment claim.
#
# (B) Genome-wide GSEA: rank every assayed gene by its signed pooled meta
#     statistic (Hedges' g / SE) and test the 50 Hallmarks with fgsea. This
#     describes the broader metastasis-association program the signature is
#     drawn from; it is program-level context, not enrichment or validation of
#     the classifier.
#
# input : outs/meta_analysis_results.rda  (genome-wide pooled effect sizes)
#         outs/PP_filter_MetaScore.rda     (POS/NEG signature membership)
#         config/metscore_signature_v1.csv (frozen 45-gene/direction contract)
# run   : Rscript code/gsea/functional_enrichment_MetScore.R   (from project root)
###############################################################################

# warning capture: warnings are logged with context, never silently dropped.
.warn <- new.env(); .warn$log <- list()
cap <- function(ctx, expr) withCallingHandlers(expr, warning = function(w) {
  .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = ctx,
    message = conditionMessage(w), stringsAsFactors = FALSE)
  invokeRestart("muffleWarning")
})

cap("package_load", suppressPackageStartupMessages({
  library(fgsea); library(msigdbr); library(dplyr)
  library(ggplot2); library(data.table); library(digest)
}))

META_RDA <- "./outs/meta_analysis_results.rda"
FILT_RDA <- "./outs/PP_filter_MetaScore.rda"
SIG_CSV  <- "./config/metscore_signature_v1.csv"
OUT_DIR  <- "./outs/functional_enrichment"
FIG_DIR  <- "./figures/functional_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
SEED <- 42L
set.seed(SEED)

# ---- 1) genome-wide signed meta rank -----------------------------------------
cap("load_inputs", { load(META_RDA); load(FILT_RDA) })
pr <- meta_analysis_results$metaAnalysis$pooledResults
gw <- data.frame(Gene = toupper(rownames(pr)),
                 g  = pr$effectSize,
                 se = pr$effectSizeStandardError,
                 stringsAsFactors = FALSE)
gw <- gw[!is.na(gw$Gene) & is.finite(gw$g) & is.finite(gw$se) & gw$se > 0, ]
n_dup_symbol <- sum(duplicated(gw$Gene))
gw <- gw[!duplicated(gw$Gene), ]                 # unique uppercase symbols
gw$z <- gw$g / gw$se                             # signed pooled meta statistic
assayed <- gw$Gene
n_ties <- sum(duplicated(gw$z))                  # tied ranks (reported)

# ---- 2) frozen signature contract, fail closed -------------------------------
sig_cfg <- read.csv(SIG_CSV, stringsAsFactors = FALSE)
sig_cfg$gene <- toupper(sig_cfg$gene)
pos_cfg <- sig_cfg$gene[grepl("^POS", sig_cfg$direction)]
neg_cfg <- sig_cfg$gene[grepl("^NEG", sig_cfg$direction)]
pos_flt <- toupper(filter$posGeneNames)
neg_flt <- toupper(filter$negGeneNames)
stopifnot(
  nrow(sig_cfg) == 45L,
  length(unique(sig_cfg$gene)) == 45L,
  length(pos_cfg) == 27L, length(neg_cfg) == 18L,
  setequal(pos_cfg, pos_flt), setequal(neg_cfg, neg_flt),
  all(c(pos_cfg, neg_cfg) %in% assayed)
)
pos <- pos_cfg; neg <- neg_cfg; sig <- c(pos, neg)

# ---- 3) Hallmark collection + annotated background ---------------------------
mdf <- cap("msigdbr", msigdbr(species = "Homo sapiens", collection = "H"))
t2g <- unique(data.frame(gs = mdf$gs_name, gene = toupper(mdf$gene_symbol),
                         stringsAsFactors = FALSE))
db_version <- as.character(mdf$db_version[1])
terms <- sort(unique(t2g$gs))                    # 50 Hallmarks, fixed order
hall  <- split(t2g$gene, t2g$gs)
universe <- intersect(assayed, unique(t2g$gene)) # assayed genes annotated in Hallmark
N <- length(universe)
# sorted TERM2GENE snapshot hash (order-independent provenance)
t2g_sorted <- t2g[order(t2g$gs, t2g$gene), ]
t2g_sha <- digest(paste(t2g_sorted$gs, t2g_sorted$gene, sep = "\t", collapse = "\n"),
                  algo = "sha256", serialize = FALSE)

# ---- 4) ORA: right-tail hypergeometric, every Hallmark x every family --------
# effective background is the annotated universe; genes with no Hallmark
# annotation cannot land in any term and are outside the hypergeometric frame.
ora_family <- function(gset, fam) {
  gu <- intersect(gset, universe); nn <- length(gu)
  rbindlist(lapply(terms, function(tm) {
    tg <- intersect(hall[[tm]], universe); M <- length(tg)
    ov <- intersect(gu, tg); k <- length(ov)
    p  <- if (k == 0L) 1 else phyper(k - 1L, M, N - M, nn, lower.tail = FALSE)
    data.table(family = fam, hallmark = tm, term_size = M, set_size = nn,
               overlap = k, overlap_genes = paste(sort(ov), collapse = ";"),
               raw_p = p)
  }))
}
ora <- rbindlist(lapply(list(c("Signature"), c("POS"), c("NEG")),
                        function(f) ora_family(switch(f, Signature = sig, POS = pos, NEG = neg), f)))
ora[, q_family := p.adjust(raw_p, "BH"), by = family]     # BH within each 50-test family
ora[, q_global := p.adjust(raw_p, "BH")]                  # BH across all 150 tests
ora[, family := factor(family, levels = c("Signature", "POS", "NEG"))]
setorder(ora, family, raw_p, hallmark)
fwrite(ora, file.path(OUT_DIR, "MetScore_Hallmark_ORA.csv"))

# ---- 5) 45-gene Hallmark membership (complete; genes with none kept as NA) ----
g2t <- split(t2g$gs, t2g$gene)
memb <- rbindlist(lapply(seq_len(nrow(sig_cfg)), function(i) {
  gsym <- sig_cfg$gene[i]; hits <- sort(unique(g2t[[gsym]]))
  dirn <- if (grepl("^POS", sig_cfg$direction[i])) "POS" else "NEG"
  if (length(hits) == 0L) data.table(gene = gsym, direction = dirn, hallmark = NA_character_)
  else data.table(gene = gsym, direction = dirn, hallmark = hits)
}))
memb[, in_universe := gene %in% universe]
fwrite(memb, file.path(OUT_DIR, "MetScore_Hallmark_membership.csv"))

# ---- 6) genome-wide GSEA -----------------------------------------------------
ranks <- sort(setNames(gw$z, gw$Gene), decreasing = TRUE)
set.seed(SEED)
gsea <- cap("fgsea", fgseaMultilevel(pathways = hall, stats = ranks,
                                     minSize = 10, maxSize = 500, eps = 0))
gsea[, padj := p.adjust(pval, "BH")]              # BH across all tested Hallmarks
setorder(gsea, padj, pval)
gsea_out <- copy(gsea)
gsea_out[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
setcolorder(gsea_out, c("pathway", "NES", "ES", "pval", "padj", "size", "leadingEdge"))
fwrite(gsea_out, file.path(OUT_DIR, "MetScore_genomewide_Hallmark_GSEA.csv"))

# ---- 7) Figure S4: signed NES lollipop, FDR-supported Hallmarks only ---------
tidy <- function(x) {
  s <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x)))
  gsub("\\bv1\\b", "V1", gsub("\\bv2\\b", "V2",
    gsub("\\b(e2f|g2m|emt|uv|tgf|myc|dna|mtorc1)\\b", "\\U\\1", s, perl = TRUE)))
}
fig <- gsea[padj < 0.05]
fig[, label := tidy(pathway)]
fig[, dirn := ifelse(NES > 0, "Higher (NES > 0)", "Lower (NES < 0)")]
fig[, dirn := factor(dirn, levels = c("Higher (NES > 0)", "Lower (NES < 0)"))]
fig[, neglogq := -log10(padj)]
setorder(fig, NES)
fig[, label := factor(label, levels = label)]
fwrite(fig[, .(pathway, label, NES, ES, pval, padj, neglogq, size, direction = dirn)],
       file.path(OUT_DIR, "FigureS4_plot_data.csv"))

col_dir <- c("Higher (NES > 0)" = "#C0392B", "Lower (NES < 0)" = "#2E86C1")
shp_dir <- c("Higher (NES > 0)" = 16, "Lower (NES < 0)" = 17)
p_fig <- ggplot(fig, aes(x = NES, y = label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES, yend = label, colour = dirn), linewidth = 0.7) +
  geom_point(aes(colour = dirn, shape = dirn, size = neglogq)) +
  scale_colour_manual(values = col_dir, name = "Enrichment direction") +
  scale_shape_manual(values = shp_dir, name = "Enrichment direction") +
  scale_size_continuous(name = expression(-log[10]~"FDR"), range = c(3, 8)) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 11),
        legend.title = element_text(size = 10), legend.text = element_text(size = 9.5),
        legend.key.height = grid::unit(13, "pt"),
        plot.margin = margin(6, 12, 6, 6)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3.4)),
         shape = guide_legend(order = 1), size = guide_legend(order = 2))
ggsave(file.path(FIG_DIR, "FigureS4_functional_enrichment.pdf"), p_fig,
       width = 7.6, height = 5.6, useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "FigureS4_functional_enrichment.tiff"), p_fig,
       width = 7.6, height = 5.6, dpi = 500, compression = "lzw")
ggsave(file.path(FIG_DIR, "FigureS4_functional_enrichment.png"), p_fig,
       width = 7.6, height = 5.6, dpi = 300)

# ---- 8) manifest: freezes, versions, counts, warnings ------------------------
warn_tab <- if (length(.warn$log)) do.call(rbind, .warn$log) else
  data.frame(context = "none", message = "no warnings captured", stringsAsFactors = FALSE)
warn_str <- paste(sprintf("[%s] %s", warn_tab$context, warn_tab$message), collapse = " | ")
pkg_v <- function(p) as.character(packageVersion(p))
manifest <- data.frame(
  key = c("meta_object_sha256", "signature_config_sha256", "term2gene_sorted_sha256",
          "R_version", "fgsea_version", "msigdbr_version", "msigdbr_db_version",
          "dplyr_version", "ggplot2_version", "data.table_version", "digest_version",
          "seed", "assayed_universe", "annotated_universe", "duplicate_symbols_dropped",
          "tied_ranks", "signature_genes", "POS_genes", "NEG_genes",
          "signature_mapped", "POS_mapped", "NEG_mapped",
          "ora_rows", "ora_positive_overlap_rows", "gsea_tested", "gsea_sig_q05",
          "warnings"),
  value = c(digest(file = META_RDA, algo = "sha256"),
            digest(file = SIG_CSV, algo = "sha256"), t2g_sha,
            R.version.string, pkg_v("fgsea"), pkg_v("msigdbr"), db_version,
            pkg_v("dplyr"), pkg_v("ggplot2"), pkg_v("data.table"), pkg_v("digest"),
            SEED, length(assayed), N, n_dup_symbol, n_ties,
            length(sig), length(pos), length(neg),
            length(intersect(sig, universe)), length(intersect(pos, universe)),
            length(intersect(neg, universe)),
            nrow(ora), sum(ora$overlap > 0), nrow(gsea), sum(gsea$padj < 0.05),
            warn_str),
  stringsAsFactors = FALSE)
fwrite(manifest, file.path(OUT_DIR, "functional_enrichment_manifest.csv"))

cat(sprintf("ORA rows=%d (positive-overlap=%d) | GSEA tested=%d (q<0.05=%d) | universe=%d/%d | ties=%d\n",
            nrow(ora), sum(ora$overlap > 0), nrow(gsea), sum(gsea$padj < 0.05), N, length(assayed), n_ties))
cat("done. tables in outs/functional_enrichment/, Figure S4 in figures/functional_enrichment/\n")

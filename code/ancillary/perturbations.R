############################################################
## Met-Score perturbation validation panel (RNA-seq)
## GSE288991, GSE285692, GSE296237, GSE287409
## - DESeq2 differential expression -> ranked stats
## - fgseaMultilevel on MetScore_POS / MetScore_NEG
## - Summary dot plot (rows = contrasts, cols = gene sets)
##
## EXPLORATORY mechanistic analysis. It tests directional enrichment of the
## 45-gene Met-Score signature (27 positive / 18 negative genes) after
## perturbation, and does not use or refit the frozen 41-feature ridge-logistic
## Met-Score classifier.
############################################################
rm(list = ls())

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)  
  library(matrixStats)
  library(DESeq2)
  library(limma)
  library(Matrix)
  library(fgsea)
  library(ggplot2)
})

# -----------------------------
# 0) Met-Score gene sets
# -----------------------------
met_pos <- c("TMSB10","ENSA","ASPN","YWHAZ","HES6","STC2","ASNS","HAVCR2","ARL6IP1","F5",
             "RFTN1","SOX4","PTPN9","ALDH1A1","MRPL11","GABRD","RC3H2","CST2","CXCR4",
             "SEM1","FOXH1","KIF7","BARD1","CADPS","RNF19A","CAMK2N1","GPR37")

met_neg <- c("KCTD14","AZGP1","PART1","CHRNA2","DPT","EDN3","KIAA1210","LTF","SIDT1","CBLL1",
             "PTN","CCK","UFM1","CPA3","CDC42EP5","TMEM121B","AKAP7","KLF4")

all_syms <- unique(toupper(mapIds(org.Hs.eg.db,
                                  keys=keys(org.Hs.eg.db, keytype="SYMBOL"),
                                  keytype="SYMBOL",
                                  column="SYMBOL",
                                  multiVals="first")))

met_pos_u <- unique(toupper(met_pos))
met_neg_u <- unique(toupper(met_neg))

message("MetScore_POS present in org.Hs.eg.db: ",
        sum(met_pos_u %in% all_syms), "/", length(met_pos_u))
message("MetScore_NEG present in org.Hs.eg.db: ",
        sum(met_neg_u %in% all_syms), "/", length(met_neg_u))

pathways <- list(
  MetScore_POS = met_pos_u,
  MetScore_NEG = met_neg_u
)


# -----------------------------
# 1) Helpers
# -----------------------------
out_dir <- "./outs/perturbation_met_score"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

build_coldata_from_counts_cols <- function(counts_cols) {
  s <- counts_cols
  s_low <- tolower(s)
  
  # Case A: DHT_NC1, CSS_SI2
  media_a <- ifelse(grepl("^dht_", s_low), "DHT",
                    ifelse(grepl("^css_", s_low), "CSS", NA))
  cond_a  <- ifelse(grepl("_nc\\d+$", s_low), "siNC",
                    ifelse(grepl("_si\\d+$", s_low), "siKMT2D", NA))
  
  # Case B: LNCaP_DN1, LNCaP_DS2, LNCaP_CN1, LNCaP_CS2
  #   D = DHT, C = CSS
  #   N = siNC, S = siKMT2D
  media_b <- ifelse(grepl("^lncap_d[ns]\\d+$", s_low), "DHT",
                    ifelse(grepl("^lncap_c[ns]\\d+$", s_low), "CSS", NA))
  cond_b  <- ifelse(grepl("^lncap_[dc]n\\d+$", s_low), "siNC",
                    ifelse(grepl("^lncap_[dc]s\\d+$", s_low), "siKMT2D", NA))
  
  media <- ifelse(!is.na(media_a), media_a, media_b)
  condition <- ifelse(!is.na(cond_a), cond_a, cond_b)
  
  out <- data.frame(
    media = factor(media),
    condition = factor(condition),
    row.names = counts_cols,
    stringsAsFactors = FALSE
  )
  out
}

rename_Sxx_to_GSM <- function(counts_obj, coldata, verbose = TRUE) {
  gsm_ids <- rownames(coldata)
  if (is.null(gsm_ids) || length(gsm_ids) == 0) stop("coldata rownames (GSM IDs) are missing.")
  
  # Helper: rename on a character vector of column names
  rename_vec <- function(cn) {
    sample_cols <- grep("^S_\\d+$", cn, value = TRUE)
    if (length(sample_cols) == 0) sample_cols <- grep("^S_\\d+", cn, value = TRUE)
    
    if (length(sample_cols) == 0) {
      if (verbose) message("No S_XX-like sample columns detected. Skipping renaming.")
      return(cn)
    }
    if (length(sample_cols) != length(gsm_ids)) {
      stop("Cannot order-map S_XX columns to GSMs: n(S_XX)=", length(sample_cols),
           " vs n(GSMs)=", length(gsm_ids), ". Confirm order or use GSM<->run mapping.")
    }
    
    if (verbose) {
      message("Renaming ", length(sample_cols), " columns (S_XX -> GSM IDs) using SeriesMatrix order.")
      message("Example: ", sample_cols[1], " -> ", gsm_ids[1])
    }
    
    idx <- match(sample_cols, cn)
    cn[idx] <- gsm_ids
    cn
  }
  
  # data.table / data.frame
  if (is.data.table(counts_obj) || is.data.frame(counts_obj)) {
    setnames(counts_obj, names(counts_obj), rename_vec(names(counts_obj)))
    return(counts_obj)
  }
  
  # matrix
  if (is.matrix(counts_obj)) {
    colnames(counts_obj) <- rename_vec(colnames(counts_obj))
    return(counts_obj)
  }
  
  stop("rename_Sxx_to_GSM: unsupported counts type: ", paste(class(counts_obj), collapse = ", "))
}

is_interval_counts <- function(counts_dt) {
  req <- c("Chr", "Start", "End")
  all(req %in% colnames(counts_dt))
}

read_counts_dt <- function(path) {
  message("Reading counts as data.table: ", path)
  dt <- data.table::fread(path, data.table = TRUE)
  return(dt)
}


intervals_to_gene_counts <- function(counts_dt, sample_ids_gsm,
                                     txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                     mapping = c("gene_body", "promoter"),
                                     promoter_up = 2000, promoter_down = 2000,
                                     keep_chr = paste0("chr", c(1:22, "X", "Y", "M")),
                                     verbose = TRUE) {
  mapping <- match.arg(mapping)
  
  stopifnot(all(c("Chr", "Start", "End") %in% colnames(counts_dt)))
  stopifnot(all(sample_ids_gsm %in% colnames(counts_dt)))
  
  setDT(counts_dt)
  
  # keep canonical chromosomes only (avoids alt contigs)
  counts_dt <- counts_dt[Chr %in% keep_chr]
  if (nrow(counts_dt) == 0) stop("After filtering to canonical chromosomes, no rows remain.")
  
  # GRanges for intervals
  gr <- GRanges(
    seqnames = counts_dt$Chr,
    ranges   = IRanges(start = as.integer(counts_dt$Start),
                       end   = as.integer(counts_dt$End))
  )
  
  # reference ranges
  if (mapping == "gene_body") {
    ref_gr <- genes(txdb)
  } else {
    ref_gr <- promoters(genes(txdb), upstream = promoter_up, downstream = promoter_down)
  }
  ref_gr <- keepStandardChromosomes(ref_gr, pruning.mode = "coarse")
  
  hits <- findOverlaps(gr, ref_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    stop("No overlaps between intervals and reference (", mapping, "). Check genome build / coordinates.")
  }
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  entrez <- names(ref_gr)[s]
  sym <- mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                column = "SYMBOL", multiVals = "first")
  sym <- toupper(as.character(sym))
  
  keep <- !is.na(sym) & sym != ""
  q <- q[keep]
  sym <- sym[keep]
  if (length(q) == 0) stop("Overlaps found but no mappable SYMBOLs after org.Hs.eg.db mapping.")
  
  # counts matrix for GSM samples (interval x sample)
  X <- as.matrix(counts_dt[, ..sample_ids_gsm])
  storage.mode(X) <- "numeric"
  
  # Build sparse incidence matrix intervals->genes and aggregate counts
  gene_levels <- unique(sym)
  gene_idx <- match(sym, gene_levels)
  
  G <- sparseMatrix(i = q, j = gene_idx, x = 1,
                    dims = c(nrow(X), length(gene_levels)))
  
  # gene x sample
  X_gene <- t(G) %*% Matrix(X, sparse = TRUE)
  rownames(X_gene) <- gene_levels
  colnames(X_gene) <- colnames(X)
  
  if (verbose) {
    message("Interval->", mapping, " aggregation: ",
            nrow(X), " intervals -> ", nrow(X_gene), " genes (", length(sample_ids_gsm), " samples).")
  }
  return(X_gene)
}


# Robust rank vector from DESeq2 result (recommended): signed -log10(padj)
make_rank <- function(res_tbl, gene_col="gene", lfc_col="log2FoldChange", padj_col="padj") {
  df <- res_tbl %>%
    dplyr::select(all_of(c(gene_col, lfc_col, padj_col))) %>%
    dplyr::filter(!is.na(.data[[gene_col]]), !is.na(.data[[lfc_col]]), !is.na(.data[[padj_col]])) %>%
    dplyr::mutate(
      gene = toupper(.data[[gene_col]]),
      padj = pmax(.data[[padj_col]], 1e-300),
      lfc  = .data[[lfc_col]],
      score = sign(lfc) * (-log10(padj))
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(score = score[which.max(abs(score))], .groups="drop")
  
  # break ties slightly (fgsea warns about ties)
  set.seed(7)
  df$score <- df$score + rnorm(nrow(df), 0, 1e-8)
  
  ranks <- df$score
  names(ranks) <- df$gene
  ranks <- ranks[!is.na(names(ranks))]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  ranks
}

make_rank_stat <- function(res_tbl, gene_col="gene", stat_col="stat") {
  df <- res_tbl %>%
    dplyr::select(all_of(c(gene_col, stat_col))) %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[stat_col]])) %>%
    mutate(
      gene = toupper(.data[[gene_col]]),
      score = .data[[stat_col]]
    ) %>%
    group_by(gene) %>%
    summarise(score = score[which.max(abs(score))], .groups="drop")
  
  set.seed(7)
  df$score <- df$score + rnorm(nrow(df), 0, 1e-8)
  
  ranks <- df$score
  names(ranks) <- df$gene
  ranks <- ranks[!duplicated(names(ranks))]
  sort(ranks, decreasing = TRUE)
}

make_rank_from_limma <- function(tt, gene_col = "gene") {
  stopifnot(all(c(gene_col, "logFC", "P.Value") %in% colnames(tt)))
  r <- with(tt, sign(logFC) * (-log10(pmax(P.Value, 1e-300))))
  names(r) <- toupper(tt[[gene_col]])
  # collapse duplicate genes by max abs score
  r <- tapply(r, names(r), function(x) x[which.max(abs(x))])
  r <- as.numeric(r); names(r) <- names(tapply(sign(tt$logFC), toupper(tt[[gene_col]]), length)) # reset names safely
  # The line above is clunky; better:
  r <- tapply(with(tt, sign(logFC) * (-log10(pmax(P.Value, 1e-300)))),
              toupper(tt[[gene_col]]),
              function(x) x[which.max(abs(x))])
  r <- sort(r, decreasing = TRUE)
  return(r)
}


run_fgsea_two_sets <- function(ranks, pathways, minSize=5, maxSize=5000) {
  # prefilter pathways to those with enough overlap
  pw2 <- lapply(pathways, function(gs) intersect(toupper(gs), names(ranks)))
  pw2 <- pw2[lengths(pw2) >= minSize]
  if (length(pw2) == 0) return(tibble())
  
  fg <- fgseaMultilevel(pathways=pw2, stats=ranks, minSize=minSize, maxSize=maxSize) %>%
    as.data.frame() %>% as_tibble() %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=";")))
  fg
}

# Reads a GEO *supplementary* counts matrix in common formats:
# - TSV/CSV: first column is gene, remaining are samples
# - Some counts are gzipped
read_counts_any <- function(path) {
  message("Reading counts: ", path)
  dt <- fread(path)
  # Heuristic: first col = gene id/symbol
  gene_col <- names(dt)[1]
  genes <- dt[[gene_col]]
  mat <- as.matrix(dt[, -1, with=FALSE])
  rownames(mat) <- genes
  mode(mat) <- "numeric"
  mat
}

# Download SeriesMatrix + supp files
get_geo_series_and_supp <- function(gse, destdir) {
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  gset <- getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
  if (length(gset) > 1) {
    # pick first platform by default; refine if needed
    gset <- gset[[1]]
  } else {
    gset <- gset[[1]]
  }
  supp <- getGEOSuppFiles(GEO = gse, baseDir = destdir, makeDirectory = FALSE)
  list(gset=gset, supp=supp)
}

# Run DESeq2 for a specified contrast
run_deseq2_contrast <- function(counts, coldata, design_formula, contrast_vec,
                                min_total_count=1,
                                met_pos=NULL, met_neg=NULL) {
  # ----met-score filter audit
  if (!is.null(met_pos) || !is.null(met_neg)) {
    met_all <- unique(toupper(c(met_pos, met_neg)))
    rn_up <- toupper(rownames(counts))
    met_in <- intersect(rn_up, met_all)
    if (length(met_in) > 0) {
      rs <- rowSums(counts)
      keep_vec <- rs >= min_total_count
      keep_map <- setNames(keep_vec, rn_up)
      filtered <- met_in[!keep_map[met_in]]
      message("MetScore genes filtered by rowSums<", min_total_count, ": ",
              ifelse(length(filtered)==0, "(none)", paste(filtered, collapse=", ")))
    }
  }
  
  keep <- rowSums(counts) >= min_total_count
  counts <- counts[keep, , drop=FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = coldata,
    design    = design_formula
  )
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = contrast_vec, independentFiltering = FALSE)
  res_tbl <- as.data.frame(res) %>% tibble::rownames_to_column("gene")
  list(dds=dds, res_tbl=res_tbl)
}

write_fgsea_csv <- function(fg, path) {
  fg2 <- as.data.frame(fg)
  if ("leadingEdge" %in% colnames(fg2)) {
    fg2$leadingEdge <- vapply(fg2$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  }
  write.csv(fg2, path, row.names = FALSE)
}

# Read one GSM FPKM file (2 columns: gene + value)
read_gsm_fpkm <- function(path) {
  dt <- data.table::fread(path, data.table = TRUE)
  
  # drop comment rows if present
  if (nrow(dt) > 0 && is.character(dt[[1]]) && any(grepl("^#", dt[[1]]))) {
    dt <- dt[!grepl("^#", dt[[1]]), ]
  }
  
  if (ncol(dt) < 2) stop("FPKM file has <2 columns: ", path)
  
  id_raw <- as.character(dt[[1]])
  val    <- suppressWarnings(as.numeric(dt[[2]]))
  
  # --- Map to SYMBOL ---
  id_up <- toupper(id_raw)
  
  # If already looks like SYMBOL (letters/numbers, not ENSG), keep
  looks_ensg <- grepl("^ENSG", id_up)
  
  sym <- id_up
  if (any(looks_ensg)) {
    ens <- sub("\\..*$", "", id_up[looks_ensg])  # drop version suffix
    mapped <- mapIds(
      org.Hs.eg.db,
      keys = ens,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    sym[looks_ensg] <- toupper(as.character(mapped))
  }
  
  out <- data.frame(gene = sym, value = val, stringsAsFactors = FALSE)
  out <- out[!is.na(out$gene) & out$gene != "" & is.finite(out$value), ]
  
  out
}

# Build expression matrix (genes x samples) from GSM*_FPKM files
build_expr_matrix_from_gsm_files <- function(files) {
  stopifnot(length(files) >= 4)
  
  sample_ids <- sub("^(GSM\\d+).*", "\\1", basename(files))
  sample_ids <- make.unique(sample_ids)
  
  lst <- vector("list", length(files))
  for (i in seq_along(files)) {
    x <- read_gsm_fpkm(files[i])
    x$sample <- sample_ids[i]
    lst[[i]] <- x
  }
  all <- dplyr::bind_rows(lst)
  
  wide <- tidyr::pivot_wider(all, names_from = sample, values_from = value)
  wide <- as.data.frame(wide)
  rownames(wide) <- wide$gene
  wide$gene <- NULL
  
  mat <- as.matrix(wide)
  storage.mode(mat) <- "numeric"
  
  # collapse duplicate genes by mean (FPKM)
  if (any(duplicated(rownames(mat)))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) / as.vector(table(rownames(mat)))
  }
  
  mat
}

# limma DE on log2(FPKM+1), returns ranks (t-stat) + fgsea
run_limma_contrast_fgsea <- function(expr_mat, coldata, factor_name, A, B, pathways) {
  stopifnot(all(colnames(expr_mat) == rownames(coldata)))
  
  y <- log2(expr_mat + 1)
  v <- matrixStats::rowVars(y)
  y <- y[v > 0, , drop = FALSE]
  
  grp <- droplevels(coldata[[factor_name]])
  grp <- relevel(grp, ref = B)
  design <- model.matrix(~ grp)
  
  fit <- lmFit(y, design)
  fit <- eBayes(fit)
  
  # coefficient corresponds to grpA vs grpB
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- toupper(rownames(tt))
  
  ranks <- tt$t
  names(ranks) <- tt$gene
  
  # drop NA / duplicate names
  keep <- !is.na(names(ranks)) & names(ranks) != ""
  ranks <- ranks[keep]
  ranks <- ranks[!duplicated(names(ranks))]
  set.seed(7)
  ranks <- ranks + rnorm(length(ranks), 0, 1e-8)
  ranks <- sort(ranks, decreasing = TRUE)
  fg <- run_fgsea_two_sets(ranks, pathways)
  
  list(tt = tt, ranks = ranks, fg = fg)
}

# -----------------------------
# 2) Dataset-specific configuration
#    (Edit ONLY this block when you add more datasets/contrasts)
# -----------------------------

# For each GSE: specify where counts file is (pattern match) and define contrasts using sample titles/characteristics.

configs <- list(
  
  # GSE287409: counts file is explicitly "GSE287409_AC-71_counts.txt.gz"
  # Design: shTCF19 vs Scramble in PC3 orthotopic tumors (6 per group). 
  GSE287409 = list(
    counts_pattern = "counts.*\\.txt(\\.gz)?$",
    build_groups = function(pdat) {
      geno <- as.character(pdat[["genotype:ch1"]])
      geno <- factor(geno)
      
      # sanity check
      print(table(geno, useNA="ifany"))
      
      data.frame(condition = geno, row.names = rownames(pdat))
    },
    contrasts = list(
      list(name="TCF19_KD_vs_Ctrl", factor="condition", A="shTCF19", B="shScramble")
      
    )
  ),
  
  # GSE285692: siKMT2D vs siNC in LNCaP and C4, DHT/CSS contexts. 
  GSE285692 = list(
    counts_pattern = "count|counts|matrix|txt|tsv|csv",  # flexible; refine after download listing
    build_groups = function(pdat) {
      data.frame(dummy = rep(1, nrow(pdat)), row.names = rownames(pdat))
    },
    contrasts = list(
      list(name="KMT2D_KD_vs_Ctrl_DHT", subset="media=='DHT'", factor="condition", A="siKMT2D", B="siNC"),
      list(name="KMT2D_KD_vs_Ctrl_CSS", subset="media=='CSS'", factor="condition", A="siKMT2D", B="siNC")
    )
  ),
  
  # GSE296237: shHOXB13 +/- CCS1477 (+ some p300/CBP KD arms). 
  GSE296237 = list(
    counts_pattern = "RAW\\.tar$|counts|featureCounts|\\.txt(\\.gz)?$",
    build_groups = function(pdat) {
      tit <- tolower(pdat$title)
      hox <- ifelse(str_detect(tit, "shhoxb13"), "shHOXB13",
                    ifelse(str_detect(tit, "shctrl"), "shCtrl", NA))
      drug <- dplyr::case_when(
        str_detect(tit, "ccs\\s*1477|ccs1477") ~ "CCS1477",
        str_detect(tit, "dmso") ~ "DMSO",
        TRUE ~ NA_character_
      )
      # optional: p300/CBP perturbations; many titles include "+shp300" etc
      p300 <- ifelse(str_detect(tit, "shp300"), "shp300", "none")
      cbp  <- ifelse(str_detect(tit, "shcbp"),  "shCBP",  "none")
      data.frame(
        hox = factor(hox),
        drug = factor(drug),
        p300 = factor(p300),
        cbp  = factor(cbp),
        row.names = rownames(pdat)
      )
    },
    contrasts = list(
      # Primary: HOXB13 KD effect under DMSO
      list(name="HOXB13_KD_vs_Ctrl_DMSO", subset="drug=='DMSO' & p300=='none' & cbp=='none'",
           factor="hox", A="shHOXB13", B="shCtrl"),
      # Primary: CCS1477 effect within shHOXB13
      list(name="CCS1477_vs_DMSO_in_shHOXB13", subset="hox=='shHOXB13' & p300=='none' & cbp=='none'",
           factor="drug", A="CCS1477", B="DMSO"),
      # Optional: CCS1477 effect within shCtrl
      list(name="CCS1477_vs_DMSO_in_shCtrl", subset="hox=='shCtrl' & p300=='none' & cbp=='none'",
           factor="drug", A="CCS1477", B="DMSO")
    )
  ),
  
  # GSE288991: EMD KO vs controls, RNA-seq triplicates, counts + RSEM available. 
  GSE288991 = list(
    # Prefer the STAR/RSEM gene-level files; avoid the KC-... COUNTS file that mixes many designs.
    counts_pattern = "STAR_RSEM_All_EMD_.*\\.csv\\.gz$",
    
    # Not used when we build coldata from columns, but keep it harmless
    build_groups = function(pdat) {
      data.frame(dummy = rep(1, nrow(pdat)), row.names = rownames(pdat))
    },
    
    # Parse sample column names like: DU145CO-1, DU145shEMD-2, 22RV1CO_1, C42BshEMD-3, etc.
    build_coldata_from_counts_cols = function(sample_cols_all) {
      
      s <- as.character(sample_cols_all)
      s_low <- tolower(s)
      
      # condition parser that handles:
      #  - ctrl1_RV1 / EMD2_RV1  (22RV1 file)
      #  - Con1 / EMD2           (C42B file)
      #  - DU145CO-1 / DU145shEMD-2 (if ever present)
      cond <- dplyr::case_when(
        grepl("shem d|shem d", gsub("_", "", s_low)) ~ "shEMD",          # ultra-defensive
        grepl("shem d", s_low) ~ "shEMD",
        grepl("shem d", gsub("-", "", s_low)) ~ "shEMD",
        grepl("\\bemd\\d+\\b", s_low) ~ "shEMD",                        # EMD1 / EMD2
        grepl("^emd\\d+$", s_low) ~ "shEMD",
        grepl("^emd\\d+_", s_low) ~ "shEMD",
        grepl("ctrl\\d", s_low) ~ "CO",                                 # ctrl1_RV1
        grepl("^con\\d+$", s_low) ~ "CO",                               # Con1
        grepl("\\bco\\b", s_low) ~ "CO",                                # ...CO...
        grepl("scramble|control", s_low) ~ "CO",
        TRUE ~ NA_character_
      )
      
      # cell line parser:
      #  - If RV1 appears, call it 22RV1 (matches your file)
      #  - Otherwise, for C42B file, force C42B
      #  - If DU145 appears, call it DU145
      cellline <- dplyr::case_when(
        grepl("rv1|22rv1", s_low) ~ "22RV1",
        grepl("du145", s_low) ~ "DU145",
        TRUE ~ "C42B"   # <- important: C42B file has Con1/EMD1 with no cellline token
      )
      
      out <- data.frame(
        cellline  = factor(cellline),
        condition = factor(cond, levels = c("CO","shEMD")),
        row.names = s,
        stringsAsFactors = FALSE
      )
      
      # keep only successfully parsed samples
      out <- out[!is.na(out$condition) & !is.na(out$cellline), , drop=FALSE]
      
      return(out)
    },
    
    # contrasts, stratified by cellline
    contrasts = list(
      list(name="EMD_KD_vs_Ctrl", factor="condition", A="shEMD", B="CO")
    )
  )
)

audit_met_overlap <- function(counts_mat, res_tbl, ranks, met_pos, met_neg) {
  met_pos <- toupper(met_pos); met_neg <- toupper(met_neg)
  
  in_counts_pos <- intersect(rownames(counts_mat), met_pos)
  in_counts_neg <- intersect(rownames(counts_mat), met_neg)
  
  # after DESeq2 results table exists
  res_genes <- toupper(res_tbl$gene)
  in_res_pos <- intersect(res_genes, met_pos)
  in_res_neg <- intersect(res_genes, met_neg)
  
  # after ranks built
  rank_genes <- names(ranks)
  in_rank_pos <- intersect(rank_genes, met_pos)
  in_rank_neg <- intersect(rank_genes, met_neg)
  
  cat("\n===== MetScore overlap audit =====\n")
  cat("POS in counts:", length(in_counts_pos), " | ", paste(in_counts_pos, collapse=", "), "\n")
  cat("NEG in counts:", length(in_counts_neg), " | ", paste(in_counts_neg, collapse=", "), "\n\n")
  
  cat("POS in DESeq2 res_tbl:", length(in_res_pos), " | ", paste(in_res_pos, collapse=", "), "\n")
  cat("NEG in DESeq2 res_tbl:", length(in_res_neg), " | ", paste(in_res_neg, collapse=", "), "\n\n")
  
  cat("POS in ranks:", length(in_rank_pos), " | ", paste(in_rank_pos, collapse=", "), "\n")
  cat("NEG in ranks:", length(in_rank_neg), " | ", paste(in_rank_neg, collapse=", "), "\n\n")
  
  # diagnose which ones got lost due to filtering / NA padj
  lost_pos <- setdiff(in_counts_pos, in_rank_pos)
  lost_neg <- setdiff(in_counts_neg, in_rank_neg)
  
  if (length(lost_pos) > 0) {
    cat("POS lost between counts -> ranks:", paste(lost_pos, collapse=", "), "\n")
    sub <- res_tbl %>% mutate(GENE=toupper(gene)) %>% filter(GENE %in% lost_pos) %>%
      dplyr::select(GENE, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
    print(sub)
    if (nrow(sub) == 0) {
      cat("  (not in DESeq2 results; likely filtered out before DESeq2 due to min_total_count)\n")
    }
  }
  if (length(lost_neg) > 0) {
    cat("NEG lost between counts -> ranks:", paste(lost_neg, collapse=", "), "\n")
    sub <- res_tbl %>% mutate(GENE=toupper(gene)) %>% filter(GENE %in% lost_neg) %>%
      dplyr::select(GENE, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
    print(sub)
    if (nrow(sub) == 0) {
      cat("  (not in DESeq2 results; likely filtered out before DESeq2 due to min_total_count)\n")
    }
  }
}

# -----------------------------
# 3) Main runner for a GSE
# -----------------------------
analyze_gse <- function(gse, cfg, pathways, out_dir) {
  
  gse_dir <- file.path(out_dir, gse)
  dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)
  
  geo <- get_geo_series_and_supp(gse, destdir = gse_dir)
  gset <- geo$gset
  pdat <- pData(gset)
  
  # Build colData (groups)
  coldata <- cfg$build_groups(pdat)
  
  # Identify a likely counts file from supplemental downloads
  supp_files <- rownames(geo$supp)
  message(gse, " supp files:\n", paste(supp_files, collapse="\n"))
  
  # Downloaded files live under gse_dir/GSE*/filename; GEOquery returns paths in geo$supp
  pick <- supp_files[str_detect(tolower(supp_files), tolower(cfg$counts_pattern))]
  if (length(pick) == 0) stop("No supplemental file matched counts_pattern for ", gse)
  
  # -----------------------------------------
  # LOOP over all matched counts files (pick)
  # -----------------------------------------
  fg_all_files  <- list()
  de_all_files  <- list()
  
  for (counts_path in pick) {
    
    # Robust: if path returned by GEOquery isn't directly readable, fall back to basename under gse_dir
    if (!file.exists(counts_path)) {
      counts_path2 <- file.path(gse_dir, basename(counts_path))
      if (file.exists(counts_path2)) counts_path <- counts_path2
    }
    if (!file.exists(counts_path)) {
      stop("Counts file does not exist after fallback: ", counts_path)
    }
    
    # Identify file label (used to namespace outputs)
    counts_file <- basename(counts_path)
    label_suffix <- tools::file_path_sans_ext(counts_file)
    label_suffix <- sub("\\.csv$", "", label_suffix, ignore.case = TRUE)
    label_suffix <- sub("\\.txt$", "", label_suffix, ignore.case = TRUE)
    label_suffix <- sub("\\.tsv$", "", label_suffix, ignore.case = TRUE)
    
    message("\n--- Processing counts file: ", counts_file, " ---")
    
    gse_subdir <- file.path(gse_dir, label_suffix)
    dir.create(gse_subdir, showWarnings = FALSE, recursive = TRUE)
    
    # Untar if needed (some GSEs use RAW.tar)
    if (str_detect(tolower(counts_path), "\\.tar$")) {
      message("Untarring: ", counts_path)
      untar(counts_path, exdir = file.path(gse_subdir, "untar"))
      
      all_files <- list.files(file.path(gse_subdir, "untar"),
                              recursive = TRUE, full.names = TRUE)
      
      # Prefer GSM*_FPKM files
      fpkm_files <- all_files[str_detect(basename(all_files), "^GSM\\d+_.*FPKM.*\\.(txt|tsv|csv)(\\.gz)?$")]
      if (length(fpkm_files) >= 4) {
        message("Detected per-sample FPKM files (n=", length(fpkm_files), "). Building expression matrix.")
        
        expr_mat <- build_expr_matrix_from_gsm_files(fpkm_files)
        
        # Build coldata from GEO metadata by GSM ID
        gsm_ids <- colnames(expr_mat)
        pdat_gsm <- pdat[match(gsm_ids, rownames(pdat)), , drop = FALSE]
        if (any(is.na(rownames(pdat_gsm)))) stop("Some GSMs not found in pData().")
        
        # Use YOUR existing build_groups() to create condition/drug/etc.
        coldata_use <- cfg$build_groups(pdat_gsm)
        # IMPORTANT: rownames must be GSM IDs matching expr_mat columns
        stopifnot(all(rownames(coldata_use) == colnames(expr_mat)))
        
        # Now run your configured contrasts using LIMMA instead of DESeq2
        fg_this <- list()
        de_this <- list()
        
        coldata_use <- coldata_use[!is.na(coldata_use$drug), , drop = FALSE]
        expr_mat <- expr_mat[, rownames(coldata_use), drop = FALSE]
        
        for (con in cfg$contrasts) {
          sub_idx <- rep(TRUE, nrow(coldata_use))
          if (!is.null(con$subset)) {
            sub_idx <- with(coldata_use, eval(parse(text = con$subset)))
            sub_idx <- as.logical(sub_idx)
          }
          
          cd_sub <- coldata_use[sub_idx, , drop = FALSE]
          keep_s <- rownames(cd_sub)
          if (length(keep_s) < 4) {
            message("[SKIP] ", con$name, " too few samples after subsetting.")
            next
          }
          
          f <- con$factor
          if (!(con$A %in% levels(cd_sub[[f]]) && con$B %in% levels(cd_sub[[f]]))) {
            message("[SKIP] ", con$name, " missing needed levels: ", con$A, " / ", con$B)
            next
          }
          
          expr_sub <- expr_mat[, keep_s, drop = FALSE]
          cd_sub <- cd_sub[keep_s, , drop = FALSE]
          
          message("\n[", gse, " | FPKM-matrix] Contrast: ", con$name)
          print(table(cd_sub[[f]], useNA = "ifany"))
          
          lim <- run_limma_contrast_fgsea(expr_sub, cd_sub, f, con$A, con$B, pathways)
          
          fg <- lim$fg %>% mutate(gse = gse, contrast = con$name, counts_file = "FPKM_from_TAR")
          tt <- lim$tt %>% mutate(gse = gse, contrast = con$name, counts_file = "FPKM_from_TAR")
          
          fg_this[[con$name]] <- fg
          de_this[[con$name]] <- tt
          
          fwrite(tt, file.path(gse_subdir, paste0(con$name, "_limma_FPKM_results.csv")))
          fwrite(fg, file.path(gse_subdir, paste0(con$name, "_fgsea_MetScore.csv")))
        }
        
        fg_file_all <- bind_rows(fg_this)
        de_file_all <- bind_rows(de_this)
        
        fwrite(fg_file_all, file.path(gse_subdir, "ALL_CONTRASTS_fgsea_MetScore.csv"))
        fwrite(de_file_all, file.path(gse_subdir, "ALL_CONTRASTS_limma_results.csv"))
        
        fg_all_files[[label_suffix]] <- fg_file_all
        de_all_files[[label_suffix]] <- de_file_all
        
        next  # IMPORTANT: skip the rest of the non-tar parsing for this counts_path
      }
      
      # If no FPKM files, fall back to your old “single matrix file” logic:
      # (keep your existing candidate-selection code here if you still want it)
      stop("TAR extracted, but no GSM*_FPKM files found to build a multi-sample matrix.")
    }
    
    # Read counts table
    counts_dt <- read_counts_dt(counts_path)
    
    # If S_XX columns exist, try to rename S_XX -> GSM order (works for GSE287409)
    counts_dt <- tryCatch(
      rename_Sxx_to_GSM(counts_dt, coldata, verbose = TRUE),
      error = function(e) {
        message("rename_Sxx_to_GSM failed (ok for some datasets): ", conditionMessage(e))
        counts_dt
      }
    )
    
    # Decide sample columns + coldata alignment
    common <- intersect(colnames(counts_dt), rownames(coldata))
    
    if (length(common) < 4) {
      
      ann_cols <- c(
        "gene_id","gene","gene_name","gene_symbol","symbol","V1",
        "Chr","Start","End","Strand","Length",
        "gene_chr","gene_start","gene_end","gene_strand","gene_length",
        "gene_biotype","gene_description","tf_family"
      )
      
      sample_cols_all <- setdiff(colnames(counts_dt), ann_cols)
      
      # Dataset-specific parser if provided (e.g. GSE288991)
      if (!is.null(cfg$build_coldata_from_counts_cols)) {
        message("No GSM overlap; using cfg$build_coldata_from_counts_cols on counts columns (n=",
                length(sample_cols_all), ").")
        
        message("Sample-like columns passed to parser:\n", paste(sample_cols_all, collapse=", "))
        coldata2 <- cfg$build_coldata_from_counts_cols(sample_cols_all)
      } else {
        # Default parser (works for KMT2D DHT/CSS files in GSE285692)
        sample_cols <- sample_cols_all[
          grepl("^(dht|css)_(nc|si)\\d+$", tolower(sample_cols_all)) |
            grepl("^lncap_[dc][ns]\\d+$", tolower(sample_cols_all))
        ]
        dropped <- setdiff(sample_cols_all, sample_cols)
        if (length(dropped) > 0) message("Dropped non-sample cols: ", paste(dropped, collapse=", "))
        message("No GSM overlap; building coldata from counts column names (n=", length(sample_cols), ").")
        coldata2 <- build_coldata_from_counts_cols(sample_cols)
      }
      
      if (nrow(coldata2) < 4) {
        stop("Too few samples in counts after parsing sample columns for: ", counts_file,
             "\nParsed coldata rows = ", nrow(coldata2),
             "\nCounts columns head: ", paste(head(sample_cols_all), collapse=", "))
      }
      
      coldata_use <- coldata2
      common <- rownames(coldata_use)
      
    } else {
      # We have GSM overlap; use the GEO-derived coldata
      coldata_use <- coldata
    }
    
    # Build gene-level count matrix
    counts_dt <- as.data.table(counts_dt)
    
    if (is_interval_counts(counts_dt)) {
      
      message("Detected interval-count format (Chr/Start/End). Aggregating to gene-level counts.")
      counts_mat <- intervals_to_gene_counts(
        counts_dt,
        sample_ids_gsm = common,
        mapping = "gene_body",
        promoter_up = 2000,
        promoter_down = 2000,
        verbose = TRUE
      )
      counts_mat <- as.matrix(counts_mat)
      
    } else {
      
      message("Detected gene-level count table. Building gene x sample matrix.")
      
      # choose best gene identifier column
      if ("gene_name" %in% names(counts_dt)) {
        gene_vec <- counts_dt[["gene_name"]]
      } else if ("gene_symbol" %in% names(counts_dt)) {
        gene_vec <- counts_dt[["gene_symbol"]]
      } else if (any(tolower(names(counts_dt)) == "symbol")) {
        sym_col <- names(counts_dt)[tolower(names(counts_dt)) == "symbol"][1]
        gene_vec <- counts_dt[[sym_col]]
      } else {
        gene_id_col <- if ("gene_id" %in% names(counts_dt)) "gene_id" else names(counts_dt)[1]
        ens <- as.character(counts_dt[[gene_id_col]])
        ens <- sub("\\..*$", "", ens)  # drop Ensembl version suffix
        gene_vec <- mapIds(org.Hs.eg.db, keys = ens,
                           keytype = "ENSEMBL", column = "SYMBOL",
                           multiVals = "first")
      }
      
      counts_mat <- as.matrix(counts_dt[, ..common])
      rownames(counts_mat) <- toupper(as.character(gene_vec))
      
      keep_rows <- !is.na(rownames(counts_mat)) & rownames(counts_mat) != ""
      counts_mat <- counts_mat[keep_rows, , drop = FALSE]
      storage.mode(counts_mat) <- "numeric"
      
      if (any(duplicated(rownames(counts_mat)))) {
        counts_mat <- rowsum(counts_mat, group = rownames(counts_mat), reorder = FALSE)
      }
    }
    
    # Align samples
    counts_mat <- counts_mat[, rownames(coldata_use), drop = FALSE]
    stopifnot(all(colnames(counts_mat) == rownames(coldata_use)))
    
    message("Counts dim (genes x samples): ", paste(dim(counts_mat), collapse = " x "))
    message("Coldata dim: ", paste(dim(coldata_use), collapse = " x "))
    
    # -----------------------------------------
    # Run contrasts for THIS counts file
    # -----------------------------------------
    fg_this  <- list()
    de_this  <- list()
    
    for (con in cfg$contrasts) {
      
      sub_idx <- rep(TRUE, nrow(coldata_use))
      if (!is.null(con$subset)) {
        sub_idx <- with(coldata_use, eval(parse(text = con$subset)))
        sub_idx <- as.logical(sub_idx)
      }
      
      cd_sub <- coldata_use[sub_idx, , drop = FALSE]
      keep_s <- intersect(rownames(cd_sub), colnames(counts_mat))
      cd_sub <- cd_sub[keep_s, , drop = FALSE]
      ct_sub <- counts_mat[, keep_s, drop = FALSE]
      
      f <- con$factor
      message("\n[", gse, " | ", label_suffix, "] Contrast: ", con$name)
      message("Levels for ", f, ": ", paste(levels(cd_sub[[f]]), collapse=", "))
      print(table(cd_sub[[f]], useNA="ifany"))
      
      if (nrow(cd_sub) < 4) {
        message("[SKIP] too few samples after subsetting.")
        next
      }
      if (!(con$A %in% levels(cd_sub[[f]]) && con$B %in% levels(cd_sub[[f]]))) {
        message("[SKIP] missing needed levels: ", con$A, " / ", con$B)
        next
      }
      
      cd_sub[[f]] <- droplevels(cd_sub[[f]])
      cd_sub[[f]] <- relevel(cd_sub[[f]], ref = con$B)
      message("Subset samples: ", paste(rownames(cd_sub), collapse=", "))
      
      design <- as.formula(paste0("~ ", f))
      de <- run_deseq2_contrast(
        counts = ct_sub,
        coldata = cd_sub,
        design_formula = design,
        contrast_vec = c(f, con$A, con$B),
        min_total_count = 1,
        met_pos = pathways$MetScore_POS,
        met_neg = pathways$MetScore_NEG
      )
      
      res_tbl <- de$res_tbl %>%
        mutate(gse = gse, contrast = con$name, counts_file = label_suffix)
      
      ranks <- make_rank_stat(res_tbl)
      audit_met_overlap(ct_sub, res_tbl, ranks, pathways$MetScore_POS, pathways$MetScore_NEG)
      
      fg <- run_fgsea_two_sets(ranks, pathways) %>%
        mutate(gse = gse, contrast = con$name, counts_file = label_suffix)
      
      fg_this[[con$name]] <- fg
      de_this[[con$name]] <- res_tbl
      
      # Save per-contrast outputs under subdir
      fwrite(res_tbl, file.path(gse_subdir, paste0(con$name, "_DESeq2_results.csv")))
      fwrite(fg,      file.path(gse_subdir, paste0(con$name, "_fgsea_MetScore.csv")))
    }
    
    fg_file_all <- bind_rows(fg_this)
    de_file_all <- bind_rows(de_this)
    
    fwrite(fg_file_all, file.path(gse_subdir, "ALL_CONTRASTS_fgsea_MetScore.csv"))
    fwrite(de_file_all, file.path(gse_subdir, "ALL_CONTRASTS_DESeq2_results.csv"))
    
    fg_all_files[[label_suffix]] <- fg_file_all
    de_all_files[[label_suffix]] <- de_file_all
  }
  
  # After processing all files:
  fg_all <- bind_rows(fg_all_files)
  de_all <- bind_rows(de_all_files)
  
  fwrite(fg_all, file.path(gse_dir, "ALL_FILES_fgsea_MetScore.csv"))
  fwrite(de_all, file.path(gse_dir, "ALL_FILES_DESeq2_results.csv"))
  
  return(list(fg = fg_all, de = de_all))
}

# -----------------------------
# 4) Run all GSEs
# -----------------------------

gse <- "GSE287409"
tmp <- analyze_gse(gse, configs[[gse]], pathways, out_dir)
table(tmp$de$contrast)
tmp$fg

gse <- "GSE285692"
tmp <- analyze_gse(gse, configs[[gse]], pathways, out_dir)
table(tmp$de$contrast)
tmp$fg

gse <- "GSE288991"
tmp <- analyze_gse(gse, configs[[gse]], pathways, out_dir)
table(tmp$de$contrast)
tmp$fg

gse <- "GSE296237"
tmp <- analyze_gse(gse, configs[[gse]], pathways, out_dir)
table(tmp$de$contrast)
tmp$fg


all_runs <- list()
for (gse in names(configs)) {
  message("\n====================\nRunning ", gse, "\n====================")
  all_runs[[gse]] <- analyze_gse(gse, configs[[gse]], pathways, out_dir)
}

fg_panel <- bind_rows(lapply(all_runs, function(x) x$fg)) %>%
  dplyr::mutate(
    counts_file = dplyr::coalesce(counts_file, "NA"),
    contrast_id = paste0(gse, " | ", counts_file, " | ", contrast),
    neglog10fdr = -log10(padj + 1e-300)
  )

fwrite(fg_panel, file.path(out_dir, "MetScore_perturbation_panel_fgsea.csv"))

# -----------------------------
# 5) Panel plot: NES (color) and -log10(FDR) (size)
# -----------------------------
# Put MetScore_NEG on the left and MetScore_POS on the right (explicit)
fg_panel2 <- fg_panel %>%
  mutate(
    pathway = factor(pathway, levels = c("MetScore_NEG", "MetScore_POS")),
    neglog10fdr = -log10(pmax(padj, 1e-300)),
    
    # Short row label that reads like a figure panel
    row_label = case_when(
      gse == "GSE287409" ~ "TCF19 KD (PC3 orthotopic)",
      gse == "GSE285692" & str_detect(counts_file, "C4")    & str_detect(contrast, "DHT") ~ "KMT2D KD (C4, DHT)",
      gse == "GSE285692" & str_detect(counts_file, "C4")    & str_detect(contrast, "CSS") ~ "KMT2D KD (C4, CSS)",
      gse == "GSE285692" & str_detect(counts_file, "LNCaP") & str_detect(contrast, "DHT") ~ "KMT2D KD (LNCaP, DHT)",
      gse == "GSE285692" & str_detect(counts_file, "LNCaP") & str_detect(contrast, "CSS") ~ "KMT2D KD (LNCaP, CSS)",
      gse == "GSE288991" & str_detect(counts_file, "22RV1") ~ "EMD KD (22Rv1)",
      gse == "GSE288991" & str_detect(counts_file, "C42B")  ~ "EMD KD (C42B)",
      gse == "GSE296237" & str_detect(contrast, "HOXB13_KD") ~ "HOXB13 KD (DMSO)",
      gse == "GSE296237" & str_detect(contrast, "in_shHOXB13") ~ "CCS1477 vs DMSO (shHOXB13)",
      gse == "GSE296237" & str_detect(contrast, "in_shCtrl")   ~ "CCS1477 vs DMSO (shCtrl)",
      TRUE ~ paste0(gse, " | ", contrast)
    )
    
    # Optional: append data-type tag for transparency (keeps it short)
    #row_label = case_when(
    #  gse == "GSE296237" ~ paste0(row_label, " [FPKM/limma]"),
    #  TRUE ~ paste0(row_label, " [counts/DESeq2]")
    #)
  )

# Order rows (group by perturbation family)
row_order <- c(
  "TCF19 KD (PC3 orthotopic)",
  "KMT2D KD (C4, DHT)",
  "KMT2D KD (C4, CSS)",
  "KMT2D KD (LNCaP, DHT)",
  "KMT2D KD (LNCaP, CSS)",
  "EMD KD (22Rv1)",
  "EMD KD (C42B)",
  "HOXB13 KD (DMSO)",
  "CCS1477 vs DMSO (shHOXB13)",
  "CCS1477 vs DMSO (shCtrl)"
)

fg_panel2 <- fg_panel2 %>%
  mutate(row_label = factor(row_label, levels = rev(row_order)))

# 2) Better NES scaling: symmetric around 0
nes_lim <- max(abs(fg_panel2$NES), na.rm = TRUE)

# 3) Plot (wide panel area, legends compact, readable labels)
p <- ggplot(fg_panel2, aes(x = pathway, y = row_label)) +
  geom_point(aes(size = neglog10fdr, color = NES), alpha = 0.95) +
  scale_x_discrete(
    labels = c(
      MetScore_NEG = "Met-Score NEG",
      MetScore_POS = "Met-Score POS"
    )
  ) +
  scale_color_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    limits = c(-nes_lim, nes_lim),
    name = "NES"
  ) +
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(2.5, 10),
    breaks = c(0, 1, 2, 3, 4, 5),
    labels = c("0", "1", "2", "3", "4", "5")
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Met-Score enrichment across metastasis-relevant perturbations",
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
    
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11),
    
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    
    panel.grid = element_blank(),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    
    plot.margin = margin(10, 20, 10, 10)
  )

# 4) Export with sane dimensions
h <- max(5.5, 0.45 * nlevels(fg_panel2$row_label))

ggsave(file.path(out_dir, "MetScore_perturbation_panel_dotplot.pdf"),
       plot = p, width = 11, height = h, useDingbats = FALSE)

ggsave(file.path(out_dir, "MetScore_perturbation_panel_dotplot.tiff"),
       plot = p, width = 11, height = h, dpi = 600, compression = "lzw")

ggsave(file.path(out_dir, "MetScore_perturbation_panel_dotplot.png"),
       plot = p, width = 11, height = h, dpi = 600)

p



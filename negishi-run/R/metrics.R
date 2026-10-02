# metrics.R: kit-only measurements on the saved episode outputs (not learner code).
# Runs in the OOD image after 05, 05b, and 05-biomart. Reads the learner directory $W and
# writes $KIT_RECORDS/metrics/*.tsv for py/summarize.py. Missing inputs are skipped.
suppressPackageStartupMessages({ library(DESeq2) })
W <- Sys.getenv("W"); rec <- Sys.getenv("KIT_RECORDS")
out <- file.path(rec, "metrics"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
wt <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)
anchors <- c("Cdkn1a", "Mdm2", "Bax", "Bbc3", "Pmaip1", "Gadd45a")
lfc_cut <- log2(1.5)
tracks <- list(
  genome = list(tab = "results/deseq2/DESeq2_results_joined.tsv", dds = "results/deseq2/dds_featurecounts.rds"),
  transcript = list(tab = "results/deseq2_kallisto/DESeq2_kallisto_results.tsv", dds = "results/deseq2_kallisto/dds_kallisto.rds"))
tabs <- list()
anc <- list(); pca <- list(); lib <- list(); top <- list()
for (tr in names(tracks)) {
  f <- file.path(W, tracks[[tr]]$tab)
  if (!file.exists(f)) { message("metrics: missing ", f); next }
  start <- as.numeric(Sys.getenv("KIT_RUN_START", "0"))
  if (as.numeric(file.mtime(f)) < start) { message("metrics: ignoring ", f, " (older than this run, not produced by it)"); next }
  x <- read.delim(f, check.names = FALSE)
  tabs[[tr]] <- x
  a <- x[match(anchors, x$external_gene_name), c("external_gene_name", "baseMean", "log2FoldChange", "padj")]
  a$external_gene_name <- anchors
  a$call <- ifelse(is.na(a$padj), "NA", ifelse(a$padj < 0.05 & a$log2FoldChange > 0, "UP",
                   ifelse(a$padj < 0.05, "DOWN", "ns")))
  a$passes_episode_cutoffs <- !is.na(a$padj) & a$padj <= 0.05 & abs(a$log2FoldChange) >= lfc_cut
  anc[[tr]] <- cbind(track = tr, a)
  sig <- x[!is.na(x$padj) & x$padj <= 0.05 & abs(x$log2FoldChange) >= lfc_cut, ]
  up <- head(sig[order(-sig$log2FoldChange), c("external_gene_name", "log2FoldChange", "padj")], 10)
  dn <- head(sig[order(sig$log2FoldChange), c("external_gene_name", "log2FoldChange", "padj")], 10)
  top[[tr]] <- rbind(cbind(track = tr, direction = "up", up), cbind(track = tr, direction = "down", dn))
  d <- file.path(W, tracks[[tr]]$dds)
  if (file.exists(d)) {
    dds <- readRDS(d)
    vsd <- vst(dds, blind = TRUE)
    rv <- matrixStats::rowVars(assay(vsd))
    sel <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
    p <- prcomp(t(assay(vsd)[sel, ]))
    pv <- p$sdev^2 / sum(p$sdev^2)
    cond <- as.character(dds$condition)
    pc1 <- p$x[, 1]
    sep <- (all(pc1[cond == "WT_IR"] > 0) && all(pc1[cond == "WT_mock"] < 0)) ||
           (all(pc1[cond == "WT_IR"] < 0) && all(pc1[cond == "WT_mock"] > 0))
    pca[[tr]] <- data.frame(track = tr, PC1_pct = round(100 * pv[1], 1), PC2_pct = round(100 * pv[2], 1),
                            PC1_separates_condition = sep)
    cs <- colSums(counts(dds))
    sf <- if (is.null(sizeFactors(dds))) colMeans(normalizationFactors(dds)) else sizeFactors(dds)
    lib[[tr]] <- data.frame(track = tr, sample = names(cs), counts = cs, size_factor = round(sf, 3))
  }
}
if (length(anc)) wt(do.call(rbind, anc), "anchors.tsv")
if (length(top)) wt(do.call(rbind, top), "top_genes.tsv")
if (length(pca)) wt(do.call(rbind, pca), "pca.tsv")
if (length(lib)) wt(do.call(rbind, lib), "libsize.tsv")
if (length(tabs) == 2) {
  g <- tabs$genome; t <- tabs$transcript
  sg <- g$ensembl_gene_id_version[!is.na(g$padj) & g$padj <= 0.05 & abs(g$log2FoldChange) >= lfc_cut]
  st <- t$ensembl_gene_id_version[!is.na(t$padj) & t$padj <= 0.05 & abs(t$log2FoldChange) >= lfc_cut]
  tg <- head(g$ensembl_gene_id_version[order(g$padj)], 50); tt <- head(t$ensembl_gene_id_version[order(t$padj)], 50)
  sh <- intersect(g$ensembl_gene_id_version, t$ensembl_gene_id_version)
  r <- cor(g$log2FoldChange[match(sh, g$ensembl_gene_id_version)], t$log2FoldChange[match(sh, t$ensembl_gene_id_version)],
           method = "spearman", use = "complete.obs")
  wt(data.frame(genes_tested_genome = nrow(g), genes_tested_transcript = nrow(t), shared_tested = length(sh),
                sig_genome = length(sg), sig_transcript = length(st), sig_both = length(intersect(sg, st)),
                jaccard = round(length(intersect(sg, st)) / length(union(sg, st)), 3),
                top50_overlap = length(intersect(tg, tt)), lfc_spearman_shared = round(r, 3)),
     "track_comparison.tsv")
}
# biomaRt spoiler vs staged annotation
fb <- file.path(W, "data/mart_biomart.tsv"); fs <- file.path(W, "data/mart.tsv")
if (file.exists(fb) && file.exists(fs)) {
  b <- read.delim(fb); s <- read.delim(fs)
  cts <- read.delim(file.path(W, "results/counts/gene_counts_clean.txt"), row.names = 1)
  ids <- rownames(cts)
  common <- intersect(b$ensembl_gene_id_version, s$ensembl_gene_id_version)
  wt(data.frame(columns_identical = identical(names(b), names(s)), rows_biomart = nrow(b), rows_staged = nrow(s),
                count_ids = length(ids), count_ids_in_biomart = sum(ids %in% b$ensembl_gene_id_version),
                count_ids_in_staged = sum(ids %in% s$ensembl_gene_id_version),
                symbol_agreement = round(mean(b$external_gene_name[match(common, b$ensembl_gene_id_version)] ==
                                              s$external_gene_name[match(common, s$ensembl_gene_id_version)], na.rm = TRUE), 4)),
     "biomart_vs_staged.tsv")
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
if (length(tabs) == 0) {
  cat("metrics: no DE tables yet (05 and 05b did not finish); failing so the step is not marked done\n")
  quit(status = 1, save = "no")
}
cat("metrics: wrote", paste(list.files(out), collapse = ", "), "\n")

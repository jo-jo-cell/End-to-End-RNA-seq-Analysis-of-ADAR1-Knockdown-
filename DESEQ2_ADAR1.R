# ==============================================================================
# RNA-seq differential expression: siADAR1 vs Scramble (DESeq2 + clusterProfiler)
# Positive log2FC = higher in siADAR1 (reference level = Scramble)
# ==============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
# Needed packages: ashr (lfcShrink), msigdbr (Hallmark GSEA, optional)

outdir <- "results"
dir.create(outdir, showWarnings = FALSE)
out <- function(f) file.path(outdir, f)

save_plot <- function(file, plot, w = 8, h = 7) {
  pdf(out(file), width = w, height = h)
  print(plot)
  dev.off()
}

# ---- 1. Load data ------------------------------------------------------------
counts   <- read.delim("counts_matrix.tsv", row.names = 1, check.names = FALSE)
metadata <- read.csv("experiment_design.csv", stringsAsFactors = FALSE)

# ---- 2. Sanitise gene identifiers -------------------------------------------
raw_symbols <- sub("_[0-9]+$", "", rownames(counts))
message(sprintf("Stripped suffixes from %d / %d genes",
                sum(raw_symbols != rownames(counts)), nrow(counts)))

# Keep the highest-count row per symbol, and log what was collapsed so
# multi-copy gene families can be inspected
ord  <- order(raw_symbols, -rowSums(counts))
keep <- ord[!duplicated(raw_symbols[ord])]
dups <- sort(unique(raw_symbols[duplicated(raw_symbols)]))
message(sprintf("Collapsed %d duplicated rows (%d symbols)",
                nrow(counts) - length(keep), length(dups)))
write.csv(data.frame(symbol = dups), out("collapsed_duplicate_symbols.csv"),
          row.names = FALSE)

counts <- counts[keep, , drop = FALSE]
rownames(counts) <- raw_symbols[keep]

# ---- 3. Wrangle metadata (Scramble = reference) -----------------------------
rownames(metadata) <- metadata$Run
metadata <- metadata[, "treatment", drop = FALSE]
metadata$treatment <- factor(metadata$treatment, levels = c("Scramble", "siADAR1"))

# ---- 4. Align counts and metadata -------------------------------------------
stopifnot(all(rownames(metadata) %in% colnames(counts)))
counts <- counts[, rownames(metadata), drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(counts)))
stopifnot(!anyDuplicated(rownames(counts)))

counts_all <- counts   # unfiltered copy for sex / knockdown checks

# ---- 5. Build DESeqDataSet --------------------------------------------------
# If sex, donor or batch differs between groups, add it here (e.g. ~ sex + treatment)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                              design = ~ treatment)
min_samples <- min(table(metadata$treatment))
dds <- dds[rowSums(counts(dds) >= 10) >= min_samples, ]

# ---- 6. Run DESeq -----------------------------------------------------------
dds <- DESeq(dds)
res     <- results(dds, contrast = c("treatment", "siADAR1", "Scramble"), alpha = 0.05)
res_shr <- lfcShrink(dds, contrast = c("treatment", "siADAR1", "Scramble"),
                     type = "ashr")
summary(res)

# ---- 7. Sample-level diagnostics (run these before interpreting anything) ----
cpm <- t(t(counts_all) / colSums(counts_all)) * 1e6
chk_genes <- c("ADAR", "XIST", "DDX3Y", "UTY", "KDM5D", "PCDH11Y", "UGT2B17")
chk <- cpm[intersect(chk_genes, rownames(cpm)), , drop = FALSE]
message("CPM of knockdown / sex / CNV-prone genes:")
print(round(chk, 1))
write.csv(round(chk, 2), out("diagnostic_gene_cpm.csv"))

# ---- 8. QC plots ------------------------------------------------------------
pdf(out("01_MA_plot.pdf"), width = 7, height = 6)
plotMA(res_shr, ylim = c(-3, 3), main = "MA plot (shrunken LFC)")
dev.off()

pdf(out("02_pvalue_histogram.pdf"), width = 7, height = 6)
hist(res$pvalue, breaks = 50, main = "p-value distribution",
     xlab = "p-value", col = "grey")
dev.off()

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "treatment", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = treatment, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw()
save_plot("03_PCA_plot.pdf", p_pca, 7, 6)

# ---- 9. Results table (built once) ------------------------------------------
res_df <- as.data.frame(res)
res_df$symbol     <- rownames(res_df)
res_df$lfc_shrunk <- res_shr[rownames(res_df), "log2FoldChange"]
res_df$padj_plot  <- pmax(ifelse(is.na(res_df$padj), 1, res_df$padj), 1e-300)

# Saved now, before any rows are dropped for GSEA
write.csv(res_df[order(res_df$padj), ], out("DE_all_results.csv"), row.names = FALSE)

sig           <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
up_siADAR1    <- subset(sig, log2FoldChange > 0)   # higher in siADAR1
up_Scramble   <- subset(sig, log2FoldChange < 0)   # higher in Scramble
message(sprintf("Significant: %d total | %d up in siADAR1 | %d up in Scramble",
                nrow(sig), nrow(up_siADAR1), nrow(up_Scramble)))

write.csv(sig[order(sig$padj), ],         out("DE_significant.csv"),            row.names = FALSE)
write.csv(up_siADAR1[order(up_siADAR1$padj), ],   out("DE_up_in_siADAR1.csv"),  row.names = FALSE)
write.csv(up_Scramble[order(up_Scramble$padj), ], out("DE_up_in_Scramble.csv"), row.names = FALSE)

# Interferon / dsRNA-sensing check (expected ADAR1 phenotype)
isg <- c("ISG15", "IFIT1", "IFIT3", "OAS1", "OAS2", "MX1", "IFI44L",
         "DDX58", "IFIH1", "STAT1")
isg_tab <- res_df[res_df$symbol %in% isg, c("symbol", "log2FoldChange", "padj")]
print(isg_tab)
write.csv(isg_tab, out("ISG_check.csv"), row.names = FALSE)

# ---- 10. Volcano plot -------------------------------------------------------
top_lab <- head(sig$symbol[order(sig$padj)], 20)
p_volc <- EnhancedVolcano(res_df,
                          lab      = res_df$symbol,
                          x        = "lfc_shrunk",
                          y        = "padj_plot",
                          selectLab = top_lab,
                          pCutoff  = 1e-50,
                          FCcutoff = 2.5,
                          pointSize = 2,
                          ylab     = "-Log10 adjusted P",
                          title    = "siADAR1 vs Scramble",
                          subtitle = "positive = higher in siADAR1")
save_plot("04_volcano.pdf", p_volc, 8, 7)

# ---- 11. Heatmap of top DEGs ------------------------------------------------
top_a <- head(up_siADAR1$symbol[order(up_siADAR1$padj)], 20)
top_b <- head(up_Scramble$symbol[order(up_Scramble$padj)], 20)
top   <- intersect(c(top_a, top_b), rownames(vsd))

pheatmap(assay(vsd)[top, , drop = FALSE],
         scale          = "row",
         annotation_col = as.data.frame(colData(dds)["treatment"]),
         fontsize_row   = 7,
         main           = "Top DEGs (row z-score of VST counts)",
         filename       = out("05_heatmap_top40.pdf"),
         width          = 8, 
         height         = 9)

# ---- 12. GO over-representation ---------------------------------------------
universe_ids <- bitr(rownames(dds), "SYMBOL", "ENTREZID", org.Hs.eg.db)

run_go <- function(symbols, label, prefix) {
  ids <- bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
  if (nrow(ids) < 5) {
    message(sprintf("GO [%s]: fewer than 5 mapped genes, skipping.", label))
    return(invisible(NULL))
  }
  ego <- enrichGO(gene = ids$ENTREZID, universe = universe_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(sprintf("GO [%s]: no significant terms (padj < 0.05).", label))
    return(invisible(NULL))
  }
  ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  write.csv(as.data.frame(ego), out(paste0(prefix, "_results.csv")), row.names = FALSE)
  save_plot(paste0(prefix, ".pdf"),
            dotplot(ego, showCategory = 15, title = paste("GO BP:", label)))
  invisible(ego)
}

ego_siADAR1  <- run_go(up_siADAR1$symbol,  "up in siADAR1",  "06_GO_up_siADAR1")
ego_Scramble <- run_go(up_Scramble$symbol, "up in Scramble", "07_GO_up_Scramble")

# ---- 13. Ranked GSEA (GO BP) ------------------------------------------------
# Rank by Wald stat; positive = higher in siADAR1
res_rank <- res_df[!is.na(res_df$stat), ]
sym2ent  <- bitr(res_rank$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)
ranked   <- merge(sym2ent, res_rank[, c("symbol", "stat")],
                  by.x = "SYMBOL", by.y = "symbol", sort = FALSE)
ranked   <- ranked[!duplicated(ranked$ENTREZID), ]
gene_list <- sort(setNames(ranked$stat, ranked$ENTREZID), decreasing = TRUE)

gse <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", minGSSize = 10, maxGSSize = 500,
             pvalueCutoff = 0.05, eps = 0, verbose = FALSE)

if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
  gse <- simplify(gse)
  write.csv(as.data.frame(gse), out("08_GSEA_GO_results.csv"), row.names = FALSE)
  save_plot("08_GSEA_GO.pdf",
            dotplot(gse, showCategory = 10, split = ".sign", title = "GSEA GO BP") +
              facet_grid(. ~ .sign),
            10, 8)
} else {
  message("gseGO returned nothing, skipping GSEA GO plot.")
}

# ---- 14. Hallmark GSEA (interferon response etc.) ---------------------------
if (requireNamespace("msigdbr", quietly = TRUE)) {
  h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "ncbi_gene")]
  h$ncbi_gene <- as.character(h$ncbi_gene)
  gsea_h <- GSEA(gene_list, TERM2GENE = h, eps = 0, pvalueCutoff = 0.05, verbose = FALSE)
  if (nrow(as.data.frame(gsea_h)) > 0) {
    write.csv(as.data.frame(gsea_h), out("10_GSEA_Hallmark_results.csv"), row.names = FALSE)
    save_plot("10_GSEA_Hallmark.pdf",
              dotplot(gsea_h, showCategory = 15, split = ".sign", title = "GSEA Hallmark") +
                facet_grid(. ~ .sign),
              10, 8)
  } else {
    message("No Hallmark sets significant at padj < 0.05.")
  }
} else {
  message("msigdbr not installed, skipping Hallmark GSEA.")
}

# ---- 15. KEGG over-representation -------------------------------------------
sig_ids  <- bitr(sig$symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db)
kegg_res <- enrichKEGG(gene = sig_ids$ENTREZID, organism = "hsa",
                       universe = universe_ids$ENTREZID,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_res <- setReadable(kegg_res, org.Hs.eg.db, keyType = "ENTREZID")  # gene symbols
  kegg_res@result <- kegg_res@result[kegg_res@result$Count >= 3, ]
  if (nrow(kegg_res@result) > 0) {
    write.csv(as.data.frame(kegg_res), out("09_KEGG_results.csv"), row.names = FALSE)
    save_plot("09_KEGG_Pathways.pdf",
              dotplot(kegg_res, showCategory = 15, title = "KEGG (significant DEGs)"), 9, 7)
  }
} else {
  message("KEGG returned nothing, skipping.")
}

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. All outputs are in: ", normalizePath(outdir))
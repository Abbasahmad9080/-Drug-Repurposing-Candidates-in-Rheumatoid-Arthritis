# Differential Expression Analysis - GSE294225 Dataset
# Comparing Healthy Controls (HC) vs Rheumatoid Arthritis (RA) patients

# Install required packages if not already installed
packages <- c("DESeq2", "dplyr", "ggplot2", "tidyverse", "pheatmap")
new_packages <- packages[!(packages %in% rownames(installed.packages()))]
if (length(new_packages)) install.packages(new_packages)

library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(pheatmap)

# Set working directory
setwd("e:/R/abbas/GSE294225 drug repurposing")

# Load count data
counts_data <- read.csv("GSE294225 data.csv", row.names = 1)
head(counts_data)
cat("\nCount data dimensions:", dim(counts_data), "\n")
cat("Sample IDs:", colnames(counts_data), "\n\n")

# Load metadata
metadata <- read.csv("GSE294225 metadata.csv", row.names = 1)
head(metadata)
cat("\nMetadata dimensions:", dim(metadata), "\n")

# Ensure sample names match between data and metadata
cat("\nChecking sample alignment...\n")
print(all(colnames(counts_data) == rownames(metadata)))

# Convert count data to matrix and round to integers
# (data appears to be normalized, so we round to nearest integer)
counts_matrix <- as.matrix(round(counts_data))

# Convert metadata condition to factor
metadata$condition <- factor(metadata$condition, levels = c("HC", "RA"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)

# Pre-filtering: remove genes with very low counts
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("\nGenes retained after filtering:", nrow(dds), "\n")

# Run DESeq2 analysis
cat("\nRunning DESeq2 pipeline...\n")
dds <- DESeq(dds)

# Get results (RA vs HC comparison)
res <- results(dds, contrast = c("condition", "RA", "HC"))
res_sorted <- res[order(res$padj), ]




# Summary of results
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# Extract significant genes (padj < 0.05)
sig_genes <- res_sorted[!is.na(res_sorted$padj) & res_sorted$padj < 0.05, ]
cat("\nNumber of significant genes (padj < 0.05):", nrow(sig_genes), "\n")

# Further filtering: log2FoldChange > 1 or < -1
sig_genes_lfc <- sig_genes[abs(sig_genes$log2FoldChange) > 1, ]
cat("Number of significant genes with |log2FC| > 1:", nrow(sig_genes_lfc), "\n\n")

# Show top 20 upregulated genes in RA
cat("=== Top 20 Upregulated in RA ===\n")
top_up <- sig_genes_lfc[sig_genes_lfc$log2FoldChange > 0, ]
top_up <- top_up[order(top_up$log2FoldChange, decreasing = TRUE), ]
print(head(top_up, 20))

# Show top 20 downregulated genes in RA
cat("\n=== Top 20 Downregulated in RA ===\n")
top_down <- sig_genes_lfc[sig_genes_lfc$log2FoldChange < 0, ]
top_down <- top_down[order(top_down$log2FoldChange), ]
print(head(top_down, 20))

# Create results dataframe
results_df <- as.data.frame(res_sorted)
results_df$Gene <- rownames(results_df)
results_df <- results_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Save all results
write.csv(results_df, "results/DEA_all_genes.csv", row.names = FALSE)
cat("\n✓ All DEA results saved to: results/DEA_all_genes.csv\n")

# Save significant genes
sig_results_df <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
write.csv(sig_results_df, "results/DEA_significant_genes_padj0.05.csv", row.names = FALSE)
cat("✓ Significant genes (padj < 0.05) saved to: results/DEA_significant_genes_padj0.05.csv\n")

# Save significant with fold change > 1
sig_lfc_results_df <- results_df[
  !is.na(results_df$padj) & 
  results_df$padj < 0.05 & 
  abs(results_df$log2FoldChange) > 1, 
]
write.csv(sig_lfc_results_df, "results/DEA_significant_genes_padj0.05_lfc1.csv", row.names = FALSE)
cat("✓ Significant genes (padj < 0.05, |logFC| > 1) saved to: results/DEA_significant_genes_padj0.05_lfc1.csv\n\n")

# Generate visualization: MA plot
cat("Generating MA plot...\n")
png("plots/DEG/MA_plot.png", width = 800, height = 600)
plotMA(res, ylim = c(-5, 5))
dev.off()
cat("✓ MA plot saved to: plots/DEG/MA_plot.png\n")

# Generate visualization: Volcano plot
cat("Generating Volcano plot...\n")
png("plots/DEG/Volcano_plot.png", width = 800, height = 600)
plot(res$log2FoldChange, -log10(res$padj),
     main = "Volcano Plot: RA vs HC",
     xlab = "log2(Fold Change)",
     ylab = "-log10(adjusted p-value)",
     pch = 20, col = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "red", "gray"))
abline(v = c(-1, 1), col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)
dev.off()
cat("✓ Volcano plot saved to: plots/DEG/Volcano_plot.png\n")


EnhancedVolcano::EnhancedVolcano(toptable = res, lab = res$Gene, x = 'log2FoldChange', y = 'padj',pCutoffCol = 0.05, FCcutoff = 1)

# Generate heatmap of top differentially expressed genes
cat("Generating heatmap of top DEGs...\n")
top_degs <- rownames(sig_genes_lfc)[1:min(20, nrow(sig_genes_lfc))]
heatmap_data <- counts(dds, normalized = TRUE)[top_degs, ]

png("plots/DEG/Top_DEGs_heatmap.png", width = 800, height = 600)
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = metadata,
         main = "Top 20 Differentially Expressed Genes",
         show_rownames = TRUE)
dev.off()
cat("✓ Heatmap saved to: plots/DEG/Top_DEGs_heatmap.png\n")

cat("\n=== DEA ANALYSIS COMPLETE ===\n")
cat("All results and plots have been saved to the results/ and plots/DEG/ directories.\n")



## ============================================================================
## HUB GENES FUNCTIONAL ENRICHMENT ANALYSIS
## GO (BP, MF, CC) + KEGG | Alzheimer's Disease Study
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggupset)
})

## ---- Directories ----
wdir <- getwd()
dir.create(file.path(wdir, "plots/hub_enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(wdir, "tables"),               recursive = TRUE, showWarnings = FALSE)

## ---- Thresholds ----
p_cutoff <- 0.05
q_cutoff <- 0.10
top_n    <- 15       # categories to show in plots

## ---- Tracking ----
enrich_log <- list()

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  HUB GENES ENRICHMENT ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

## ============================================================================
## 1. LOAD HUB GENES
## ============================================================================


DEGS <- read.csv('results/DEA_significant_genes_padj0.05.csv')

## Use top 20 hub gene symbols
hub_symbols <- unique(na.omit(DEGS$Gene))
hub_symbols <- hub_symbols[hub_symbols != "" & !is.na(hub_symbols)]


cat("  ✓ Top 20 hub genes (symbols) :", length(hub_symbols), "\n")
cat("  ✓ Genes:", paste(hub_symbols, collapse = ", "), "\n\n")

## ============================================================================
## 2. ID CONVERSION: SYMBOL → ENTREZ
## ============================================================================

cat("Converting gene symbols → ENTREZ IDs...\n")

id_map <- bitr(
  hub_symbols,
  fromType = "SYMBOL",
  toType   = c("ENTREZID", "ENSEMBL"),
  OrgDb    = org.Hs.eg.db
)

id_map <- id_map[!duplicated(id_map$SYMBOL), ]
entrez_ids <- unique(id_map$ENTREZID)

cat("  ✓ Successfully mapped:", length(entrez_ids), "genes\n")
cat("  ! Unmapped:", length(hub_symbols) - nrow(id_map), "genes\n\n")

enrich_log$input_genes  <- length(hub_symbols)
enrich_log$mapped_genes <- length(entrez_ids)


entrez_ids
## ============================================================================
## 3. GO ENRICHMENT — BIOLOGICAL PROCESS (BP)
## ============================================================================

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  GO Enrichment: Biological Process (BP)\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

go_bp <- tryCatch({
  clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
}, error = function(e) { cat("  ✗ GO BP failed:", e$message, "\n"); NULL })

go_bp_n <- 0
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  go_bp_n <- nrow(go_bp@result)
  cat("  ✓ Significant GO BP terms:", go_bp_n, "\n")
  
  ## Save table
  write.csv(as.data.frame(go_bp@result),
            file.path(wdir, "tables", "HubGenes_GO_BP_Enrichment.csv"),
            row.names = FALSE)
  cat("  ✓ Saved: tables/HubGenes_GO_BP_Enrichment.csv\n")
  
  show_n <- min(top_n, go_bp_n)
  
  ## Dotplot
  p1 <- enrichplot::dotplot(go_bp, showCategory = show_n,
                            font.size = 10, color = "p.adjust") +
    labs(title = "Significant DEGS: GO Biological Process",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "01_GO_BP_Dotplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p1); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "01_GO_BP_Dotplot.pdf"), width = 12, height = 8)
  print(p1); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/01_GO_BP_Dotplot.png/.pdf\n")
  
  ## Barplot
  p2 <- barplot(go_bp, showCategory = show_n, font.size = 10) +
    labs(title = "Significant DEGS: GO Biological Process",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "02_GO_BP_Barplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p2); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "02_GO_BP_Barplot.pdf"), width = 12, height = 8)
  print(p2); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/02_GO_BP_Barplot.png/.pdf\n")
  
  ## Network (cnetplot)
  p3 <- enrichplot::cnetplot(go_bp, showCategory = min(5, go_bp_n),
                             foldChange = NULL) +
    labs(title = "Significant DEGS: GO BP Gene-Concept Network") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "03_GO_BP_Cnetplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p3); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "03_GO_BP_Cnetplot.pdf"), width = 14, height = 10)
  print(p3); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/03_GO_BP_Cnetplot.png/.pdf\n")
  
  ## Emap (enrichment map)
  go_bp_sim <- enrichplot::pairwise_termsim(go_bp)
  p4 <- enrichplot::emapplot(go_bp_sim, showCategory = min(20, go_bp_n)) +
    labs(title = "Significant DEGS: GO BP Enrichment Map") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "04_GO_BP_Emapplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p4); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "04_GO_BP_Emapplot.pdf"), width = 14, height = 10)
  print(p4); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/04_GO_BP_Emapplot.png/.pdf\n\n")
  
} else {
  cat("  ! No significant GO BP terms found\n\n")
}
enrich_log$go_bp <- go_bp_n

## ============================================================================
## 4. GO ENRICHMENT — MOLECULAR FUNCTION (MF)
## ============================================================================

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  GO Enrichment: Molecular Function (MF)\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

go_mf <- tryCatch({
  clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
}, error = function(e) { cat("  ✗ GO MF failed:", e$message, "\n"); NULL })

go_mf_n <- 0
if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  go_mf_n <- nrow(go_mf@result)
  cat("  ✓ Significant GO MF terms:", go_mf_n, "\n")
  
  write.csv(as.data.frame(go_mf@result),
            file.path(wdir, "tables", "HubGenes_GO_MF_Enrichment.csv"),
            row.names = FALSE)
  cat("  ✓ Saved: tables/HubGenes_GO_MF_Enrichment.csv\n")
  
  show_n <- min(top_n, go_mf_n)
  
  ## Dotplot
  p5 <- enrichplot::dotplot(go_mf, showCategory = show_n,
                            font.size = 10, color = "p.adjust") +
    labs(title = "Significant DEGS: GO Molecular Function",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "05_GO_MF_Dotplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p5); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "05_GO_MF_Dotplot.pdf"), width = 12, height = 8)
  print(p5); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/05_GO_MF_Dotplot.png/.pdf\n")
  
  ## Barplot
  p6 <- barplot(go_mf, showCategory = show_n, font.size = 10) +
    labs(title = "Significant DEGS: GO Molecular Function",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "06_GO_MF_Barplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p6); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "06_GO_MF_Barplot.pdf"), width = 12, height = 8)
  print(p6); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/06_GO_MF_Barplot.png/.pdf\n")
  
  ## Cnetplot
  p7 <- enrichplot::cnetplot(go_mf, showCategory = min(5, go_mf_n),
  ) +
    labs(title = "Significant DEGS: GO MF Gene-Concept Network") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "07_GO_MF_Cnetplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p7); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "07_GO_MF_Cnetplot.pdf"), width = 14, height = 10)
  print(p7); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/07_GO_MF_Cnetplot.png/.pdf\n\n")
  
} else {
  cat("  ! No significant GO MF terms found\n\n")
}
enrich_log$go_mf <- go_mf_n

## ============================================================================
## 5. GO ENRICHMENT — CELLULAR COMPONENT (CC)
## ============================================================================

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  GO Enrichment: Cellular Component (CC)\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

go_cc <- tryCatch({
  clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
}, error = function(e) { cat("  ✗ GO CC failed:", e$message, "\n"); NULL })

go_cc_n <- 0
if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  go_cc_n <- nrow(go_cc@result)
  cat("  ✓ Significant GO CC terms:", go_cc_n, "\n")
  
  write.csv(as.data.frame(go_cc@result),
            file.path(wdir, "tables", "HubGenes_GO_CC_Enrichment.csv"),
            row.names = FALSE)
  cat("  ✓ Saved: tables/HubGenes_GO_CC_Enrichment.csv\n")
  
  show_n <- min(top_n, go_cc_n)
  
  ## Dotplot
  p8 <- enrichplot::dotplot(go_cc, showCategory = show_n,
                            font.size = 10, color = "p.adjust") +
    labs(title = "Significant DEGS: GO Cellular Component",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "08_GO_CC_Dotplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p8); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "08_GO_CC_Dotplot.pdf"), width = 12, height = 8)
  print(p8); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/08_GO_CC_Dotplot.png/.pdf\n")
  
  ## Barplot
  p9 <- barplot(go_cc, showCategory = show_n, font.size = 10) +
    labs(title = "Significant DEGS: GO Cellular Component",
         subtitle = paste("Top", show_n, "enriched terms")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "09_GO_CC_Barplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p9); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "09_GO_CC_Barplot.pdf"), width = 12, height = 8)
  print(p9); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/09_GO_CC_Barplot.png/.pdf\n")
  
  ## Cnetplot
  p10 <- enrichplot::cnetplot(go_cc, showCategory = min(5, go_cc_n),
  ) +
    labs(title = "Significant DEGS: GO CC Gene-Concept Network") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "10_GO_CC_Cnetplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p10); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "10_GO_CC_Cnetplot.pdf"), width = 14, height = 10)
  print(p10); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/10_GO_CC_Cnetplot.png/.pdf\n\n")
  
} else {
  cat("  ! No significant GO CC terms found\n\n")
}
enrich_log$go_cc <- go_cc_n

## ============================================================================
## 6. COMBINED GO COMPARISON PLOT (BP + MF + CC)
## ============================================================================

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  Combined GO Comparison Plot\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

tryCatch({
  go_all <- clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
  
  if (!is.null(go_all) && nrow(go_all@result) > 0) {
    write.csv(as.data.frame(go_all@result),
              file.path(wdir, "tables", "HubGenes_GO_ALL_Enrichment.csv"),
              row.names = FALSE)
    
    p_all <- enrichplot::dotplot(go_all,
                                 showCategory = 10,
                                 split        = "ONTOLOGY",
                                 font.size    = 9) +
      facet_grid(ONTOLOGY ~ ., scale = "free") +
      labs(title    = "Significant DEGS: GO Enrichment (BP + MF + CC)",
           subtitle = "Top 10 terms per ontology") +
      theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            strip.text    = element_text(face = "bold", size = 11))
    
    png(file.path(wdir, "plots/hub_enrichment", "11_GO_ALL_Combined_Dotplot.png"),
        width = 14, height = 14, units = "in", res = 300)
    print(p_all); dev.off()
    pdf(file.path(wdir, "plots/hub_enrichment", "11_GO_ALL_Combined_Dotplot.pdf"),
        width = 14, height = 14)
    print(p_all); dev.off()
    cat("  ✓ Saved: plots/hub_enrichment/11_GO_ALL_Combined_Dotplot.png/.pdf\n\n")
  }
}, error = function(e) cat("  ! Combined GO plot failed:", e$message, "\n\n"))

## ============================================================================
## 7. KEGG PATHWAY ENRICHMENT
## ============================================================================

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("  KEGG Pathway Enrichment\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

options(timeout = 600)
kegg_n <- 0

kegg_enrich <- tryCatch({
  clusterProfiler::enrichKEGG(
    gene          = entrez_ids,
    organism      = "hsa",
    pvalueCutoff  = p_cutoff,
    pAdjustMethod = "BH"
  )
}, error = function(e) { cat("  ✗ KEGG failed:", e$message, "\n"); NULL })

if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  kegg_n <- nrow(kegg_enrich@result)
  cat("  ✓ Significant KEGG pathways:", kegg_n, "\n")
  
  ## Make gene names readable
  kegg_readable <- clusterProfiler::setReadable(kegg_enrich,
                                                OrgDb    = org.Hs.eg.db,
                                                keyType  = "ENTREZID")
  
  write.csv(as.data.frame(kegg_readable@result),
            file.path(wdir, "tables", "HubGenes_KEGG_Enrichment.csv"),
            row.names = FALSE)
  cat("  ✓ Saved: tables/HubGenes_KEGG_Enrichment.csv\n")
  
  show_n <- min(top_n, kegg_n)
  
  ## Dotplot
  p11 <- enrichplot::dotplot(kegg_readable, showCategory = show_n,
                             font.size = 10, color = "p.adjust") +
    labs(title    = "Significant DEGS: KEGG Pathway Enrichment",
         subtitle = paste("Top", show_n, "pathways")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "12_KEGG_Dotplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p11); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "12_KEGG_Dotplot.pdf"), width = 12, height = 8)
  print(p11); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/12_KEGG_Dotplot.png/.pdf\n")
  
  ## Barplot
  p12 <- barplot(kegg_readable, showCategory = show_n, font.size = 10) +
    labs(title    = "Significant DEGS: KEGG Pathway Enrichment",
         subtitle = paste("Top", show_n, "pathways")) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  png(file.path(wdir, "plots/hub_enrichment", "13_KEGG_Barplot.png"),
      width = 12, height = 8, units = "in", res = 300)
  print(p12); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "13_KEGG_Barplot.pdf"), width = 12, height = 8)
  print(p12); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/13_KEGG_Barplot.png/.pdf\n")
  
  ## Cnetplot
  p13 <- enrichplot::cnetplot(kegg_readable, showCategory = min(5, kegg_n)
  ) +
    labs(title = "Significant DEGS: KEGG Gene-Concept Network") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "14_KEGG_Cnetplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p13); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "14_KEGG_Cnetplot.pdf"), width = 14, height = 10)
  print(p13); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/14_KEGG_Cnetplot.png/.pdf\n")
  
  ## Enrichment map
  kegg_sim <- enrichplot::pairwise_termsim(kegg_readable)
  p14 <- enrichplot::emapplot(kegg_sim, showCategory = min(20, kegg_n)) +
    labs(title = "Significant DEGS: KEGG Enrichment Map") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  
  png(file.path(wdir, "plots/hub_enrichment", "15_KEGG_Emapplot.png"),
      width = 14, height = 10, units = "in", res = 300)
  print(p14); dev.off()
  pdf(file.path(wdir, "plots/hub_enrichment", "15_KEGG_Emapplot.pdf"), width = 14, height = 10)
  print(p14); dev.off()
  cat("  ✓ Saved: plots/hub_enrichment/15_KEGG_Emapplot.png/.pdf\n\n")
  
} else {
  cat("  ! No significant KEGG pathways found\n\n")
}
enrich_log$kegg <- kegg_n

## ============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  STEP 7: PPI NETWORK ANALYSIS (CLEAN + ROBUST)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(STRINGdb)
library(igraph)
library(ggraph)
library(dplyr)

options(timeout = 16000)

## ------------------------------------------
## 1. DEG genes (SYMBOLS assumed)
## ------------------------------------------

deg_genes <- rownames(sig_lfc_results_df)

## ------------------------------------------
## 2. STRING mapping (FIXED)
## ------------------------------------------

gene_df <- data.frame(gene = deg_genes)

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 400
)

mapped <- string_db$map(
  gene_df,
  "gene",
  removeUnmappedRows = TRUE
)

cat("✓ Mapped genes:", nrow(mapped), "\n\n")

## ------------------------------------------
## 3. PPI network
## ------------------------------------------

ppi <- string_db$get_interactions(mapped$STRING_id)

## remove NA + self loops
ppi <- ppi[complete.cases(ppi[, c("from","to")]), ]
ppi <- ppi[ppi$from != ppi$to, ]

g <- graph_from_data_frame(ppi[, c("from","to")], directed = FALSE)

## ------------------------------------------
## 4. Hub metrics
## ------------------------------------------

deg <- degree(g)
bet <- betweenness(g)

hub_table <- data.frame(
  STRING_ID   = names(deg),
  Degree      = deg,
  Betweenness = bet
)

hub_table <- hub_table %>%
  arrange(desc(Degree))

## attach gene names
id_map <- mapped[, c("STRING_id", "gene")]

hub_table <- merge(
  hub_table,
  id_map,
  by.x = "STRING_ID",
  by.y = "STRING_id",
  all.x = TRUE
)

## ------------------------------------------
## 5. TOP 15% HUBS (Degree + Betweenness overlap)
## ------------------------------------------

deg_cut <- quantile(hub_table$Degree, 0.85, na.rm = TRUE)
bet_cut <- quantile(hub_table$Betweenness, 0.85, na.rm = TRUE)

top15_hubs <- hub_table %>%
  filter(Degree >= deg_cut & Betweenness >= bet_cut)

cat("✓ Top 15% hubs:", nrow(top15_hubs), "\n")

top15_ids <- top15_hubs$STRING_ID

g_top15 <- induced_subgraph(g, vids = V(g)[name %in% top15_ids])

V(g_top15)$gene <- id_map$gene[match(V(g_top15)$name, id_map$STRING_id)]
V(g_top15)$gene[is.na(V(g_top15)$gene)] <- V(g_top15)$name[is.na(V(g_top15)$gene)]

## ------------------------------------------
## 6. CORE HUBS (Top 2% strict)
## ------------------------------------------

deg_cut2 <- quantile(hub_table$Degree, 0.98, na.rm = TRUE)
bet_cut2 <- quantile(hub_table$Betweenness, 0.98, na.rm = TRUE)

core_hubs <- hub_table %>%
  filter(Degree >= deg_cut2 & Betweenness >= bet_cut2)

cat("✓ Core hubs (2%):", nrow(core_hubs), "\n")

core_ids <- core_hubs$STRING_ID

g_core <- induced_subgraph(g, vids = V(g)[name %in% core_ids])

V(g_core)$gene <- id_map$gene[match(V(g_core)$name, id_map$STRING_id)]
V(g_core)$gene[is.na(V(g_core)$gene)] <- V(g_core)$name[is.na(V(g_core)$gene)]

## ------------------------------------------
## 7. SAVE TABLES
## ------------------------------------------

write.csv(hub_table,
          file.path(wdir, "csvs", "all_hubs.csv"),
          row.names = FALSE)

write.csv(top15_hubs,
          file.path(wdir, "csvs", "top15_percent_hubs.csv"),
          row.names = FALSE)

write.csv(core_hubs,
          file.path(wdir, "csvs", "core_2_percent_hubs.csv"),
          row.names = FALSE)

cat("✓ CSVs saved\n\n")

## ------------------------------------------
## 8. PLOTS
## ------------------------------------------

## TOP 15% NETWORK
png(file.path(wdir, "plots/network", "Top15_Hub_Network.png"),
    width = 14, height = 10, res = 300, units = "in")

ggraph(g_top15, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(size = 4, color = "steelblue") +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Top 15% Hub Genes Network")

dev.off()

pdf(file.path(wdir, "plots/network", "Top15_Hub_Network.pdf"),
    width = 14, height = 10)

ggraph(g_top15, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(size = 4, color = "steelblue") +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Top 15% Hub Genes Network")

dev.off()

## ------------------------------------------
## CORE 2% NETWORK
## ------------------------------------------

png(file.path(wdir, "plots/network", "Core_2percent_Hub_Network.png"),
    width = 14, height = 10, res = 300, units = "in")

ggraph(g_core, layout = "stress") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(size = 12, color = "skyblue") +
  geom_node_text(aes(label = gene), repel = TRUE, size = 4) +
  theme_void() +
  ggtitle("Core Hub Genes Network (Top 2%)")

dev.off()

pdf(file.path(wdir, "plots/network", "Core_2percent_Hub_Network.pdf"),
    width = 10, height = 10)

ggraph(g_core, layout = "kk") +
  geom_edge_link(alpha = 0.4, linewidth = 0.5, color = "black") +
  geom_node_point(size = 12, color = "skyblue") +
  geom_node_text(aes(label = gene), repel = TRUE, size = 4) +
  theme_void() +
  ggtitle("Core Hub Genes Network (Top 2%)")

dev.off()

## ------------------------------------------
## LOG
## ------------------------------------------

pipeline_log$step7 <- list(
  total_hubs = nrow(hub_table),
  top15_hubs = nrow(top15_hubs),
  core_hubs  = nrow(core_hubs)
)

cat("════════ STEP 7 COMPLETE ════════\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  STEP 8: DRUG REPURPOSING\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

hub_data   <- read.csv('tables/ML_SVM_RFE_Selected_Genes.csv')
core_genes <- unique(na.omit(hub_data$gene))[1:20]

write.table(core_genes,
            file.path(wdir, "csvs", "DGIdb_input_genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("  ✓ Saved: csvs/DGIdb_input_genes.txt (", length(core_genes), "genes)\n")

## Load DGIdb results
drug_data <- read.csv("gene_interaction_results-23_05_2026.tsv")

drug_clean <- drug_data %>%
  tidyr::separate(col  = 1,
                  into = c("gene", "drug", "status", "approval", "score"),
                  sep  = "\t",
                  fill = "right")
drug_clean$score <- as.numeric(drug_clean$score)

write.csv(drug_clean,
          file.path(wdir, "csvs", "drug_interactions_clean.csv"),
          row.names = FALSE)
cat("  ✓ Saved: csvs/drug_interactions_clean.csv (", nrow(drug_clean), "rows)\n")

## Rank drugs
drug_rank <- drug_clean %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(
    targets   = n(),
    avg_score = mean(score, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(desc(targets))

write.csv(drug_rank,
          file.path(wdir, "csvs", "drug_rank_raw.csv"),
          row.names = FALSE)

## Clean drug names
drug_rank1 <- drug_rank %>%
  dplyr::mutate(
    drug = as.character(drug),
    drug = na_if(drug, "NA"),
    drug = na_if(drug, ""),
    drug = gsub("[0-9]",      "",  drug),
    drug = gsub("[^A-Za-z ]", " ", drug),
    drug = trimws(drug),
    drug = gsub("\\s+", " ",  drug)
  ) %>%
  dplyr::filter(!is.na(drug) & drug != "")



write.csv(drug_rank1,
          file.path(wdir, "csvs", "All_final_drugs_list.csv"),
          row.names = FALSE)
cat("  ✓ Saved: csvs/All_final_drugs_list.csv (", nrow(drug_rank1), "drugs after cleaning)\n")

drug_rank2 <- read.csv(file.path(wdir, "csvs", "All_final_drugs_list.csv"))


drug_rank2 <- drug_rank2 %>%
  filter(!is.na(avg_score) & !is.nan(avg_score))

drug_rank2 <- drug_rank2 %>%
  filter(!is.na(avg_score) & !is.nan(avg_score))


drug_rank2 <- drug_rank2 %>%
  arrange(desc(targets), desc(avg_score))

summary(drug_rank2$avg_score)
head(drug_rank2)

## Drug bar plot
p_drug <- head(drug_rank2, 12) %>%
  ggplot(aes(reorder(drug, targets), targets)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title    = "Top 12 Drug Candidates by Hub Gene Targets",
    subtitle = "RA – Drug Repurposing via DGIdb",
    x        = "Drug",
    y        = "Number of Hub Gene Targets"
  )

pdf(file.path(wdir, "plots/drug_repurposing", "01_Top12_Drug_Candidates.pdf"),
    width = 7, height = 5)
print(p_drug); dev.off()
png(file.path(wdir, "plots/drug_repurposing", "01_Top12_Drug_Candidates.png"),
    width = 10, height = 6, units = "in", res = 300)
print(p_drug); dev.off()
cat("  ✓ Saved: plots/drug_repurposing/01_Top12_Drug_Candidates.png/.pdf\n\n")

final12topdrugs <- head(drug_rank2, 12)


write.csv(final12topdrugs, 'csvs/final Drug targets.csv')


## ------------------------------------------
## Save Top 12 Drug Candidates
## ------------------------------------------
library(igraph)
library(ggraph)
library(dplyr)

## ------------------------------------------
## 1. Top drug-target interactions
## ------------------------------------------

top12_drugs <- head(drug_rank2$drug, 12)

drug_target_net <- drug_clean %>%
  filter(drug %in% top12_drugs) %>%
  filter(!is.na(gene),
         !is.na(drug)) %>%
  distinct(gene, drug)

## Save interaction table
write.csv(
  drug_target_net,
  file.path(
    wdir,
    "csvs/drug_repurposing",
    "Top12_Drug_Target_Interactions.csv"
  ),
  row.names = FALSE
)

cat("✓ Drug-target interactions saved\n")

## ------------------------------------------
## 2. Build network
## ------------------------------------------

g_drug <- graph_from_data_frame(
  drug_target_net,
  directed = FALSE
)

## Node type
V(g_drug)$type <- ifelse(
  V(g_drug)$name %in% drug_target_net$drug,
  "Drug",
  "Gene"
)

## ------------------------------------------
## 3. Publication-ready network
## ------------------------------------------

p_drug_net <- ggraph(g_drug, layout = "stress") +
  
  geom_edge_link(
    color = "grey75",
    alpha = 0.5,
    linewidth = 0.6
  ) +
  
  geom_node_point(
    aes(color = type,
        size = ifelse(type == "Drug", 10, 7))
  ) +
  
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3.5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c(
      Drug = "#E63946",
      Gene = "#457B9D"
    )
  ) +
  
  theme_void(base_size = 12) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.position = "right"
  ) +
  
  labs(
    title = "Drug–Target Interaction Network",
    subtitle = "Top Drug Candidates and Hub Gene Targets",
    color = "Node Type"
  )

## ------------------------------------------
## 4. Save plot
## ------------------------------------------

pdf(
  file.path(
    wdir,
    "plots/drug_repurposing",
    "Drug_Target_Network.pdf"
  ),
  width = 10,
  height = 8
)
print(p_drug_net)
dev.off()

png(
  file.path(
    wdir,
    "plots/drug_repurposing",
    "Drug_Target_Network.png"
  ),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)
print(p_drug_net)
dev.off()

cat("✓ Saved Drug_Target_Network.png/.pdf\n")







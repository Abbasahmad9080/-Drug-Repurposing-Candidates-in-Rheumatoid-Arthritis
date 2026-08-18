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

cat("Loading hub genes...\n")
hub_data   <- read.csv(file.path(wdir, "csvs", "all hub genes1.csv"))
core_hubs  <- read.csv(file.path(wdir, "csvs", "core_hub_genes.csv"))
top20_hubs <- read.csv(file.path(wdir, "csvs", "top_20_hub_genes.csv"))

## Use top 20 hub gene symbols
hub_symbols <- unique(na.omit(top20_hubs$gene))
hub_symbols <- hub_symbols[hub_symbols != "" & !is.na(hub_symbols)]

cat("  ✓ Total hub genes loaded     :", nrow(hub_data), "\n")
cat("  ✓ Core hub genes             :", nrow(core_hubs), "\n")
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
    labs(title = "Hub Genes: GO Biological Process",
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
    labs(title = "Hub Genes: GO Biological Process",
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
                              foldChange = NULL, colorEdge = TRUE) +
    labs(title = "Hub Genes: GO BP Gene-Concept Network") +
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
    labs(title = "Hub Genes: GO BP Enrichment Map") +
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
    labs(title = "Hub Genes: GO Molecular Function",
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
    labs(title = "Hub Genes: GO Molecular Function",
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
                              colorEdge = TRUE) +
    labs(title = "Hub Genes: GO MF Gene-Concept Network") +
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
    labs(title = "Hub Genes: GO Cellular Component",
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
    labs(title = "Hub Genes: GO Cellular Component",
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
                               colorEdge = TRUE) +
    labs(title = "Hub Genes: GO CC Gene-Concept Network") +
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
      labs(title    = "Hub Genes: GO Enrichment (BP + MF + CC)",
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
    labs(title    = "Hub Genes: KEGG Pathway Enrichment",
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
    labs(title    = "Hub Genes: KEGG Pathway Enrichment",
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
  p13 <- enrichplot::cnetplot(kegg_readable, showCategory = min(5, kegg_n),
                               colorEdge = TRUE) +
    labs(title = "Hub Genes: KEGG Gene-Concept Network") +
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
    labs(title = "Hub Genes: KEGG Enrichment Map") +
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
## FINAL SUMMARY
## ============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║         HUB GENES ENRICHMENT ANALYSIS — SUMMARY                      ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Input hub genes              : %-6d                              ║\n", enrich_log$input_genes))
cat(sprintf("║  Successfully mapped (ENTREZ) : %-6d                              ║\n", enrich_log$mapped_genes))
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  GO Biological Process (BP)   : %-6d significant terms             ║\n", enrich_log$go_bp))
cat(sprintf("║  GO Molecular Function (MF)   : %-6d significant terms             ║\n", enrich_log$go_mf))
cat(sprintf("║  GO Cellular Component (CC)   : %-6d significant terms             ║\n", enrich_log$go_cc))
cat(sprintf("║  KEGG Pathways                : %-6d significant pathways          ║\n", enrich_log$kegg))
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  TABLES SAVED                                                         ║\n")
cat("║    tables/HubGenes_GO_BP_Enrichment.csv                               ║\n")
cat("║    tables/HubGenes_GO_MF_Enrichment.csv                               ║\n")
cat("║    tables/HubGenes_GO_CC_Enrichment.csv                               ║\n")
cat("║    tables/HubGenes_GO_ALL_Enrichment.csv                              ║\n")
cat("║    tables/HubGenes_KEGG_Enrichment.csv                                ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  PLOTS SAVED  (plots/hub_enrichment/)                                 ║\n")
cat("║    01_GO_BP_Dotplot       02_GO_BP_Barplot    03_GO_BP_Cnetplot       ║\n")
cat("║    04_GO_BP_Emapplot      05_GO_MF_Dotplot    06_GO_MF_Barplot        ║\n")
cat("║    07_GO_MF_Cnetplot      08_GO_CC_Dotplot    09_GO_CC_Barplot        ║\n")
cat("║    10_GO_CC_Cnetplot      11_GO_ALL_Combined  12_KEGG_Dotplot         ║\n")
cat("║    13_KEGG_Barplot        14_KEGG_Cnetplot    15_KEGG_Emapplot        ║\n")
cat("║    (each saved as .png + .pdf)                                        ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("  Completed: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

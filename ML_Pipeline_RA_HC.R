## ============================================================================
## MACHINE LEARNING PIPELINE — BIOMARKER DISCOVERY
## Classes: RA vs HC | LOOCV | LASSO + RF + SVM-RFE + XGBoost
## ============================================================================

## ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(caret)
  library(randomForest)
  library(e1071)
  library(pROC)
  library(pheatmap)
  library(reshape2)
  library(xgboost)
  library(MLmetrics)
  library(gridExtra)
  library(ggrepel)
  library(RColorBrewer)
  library(viridis)
})

## ---- Directories ----
for (d in c("plots/machine_learning", "tables", "csvs")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(123)
RANDOM_STATE <- 123

## ---- Colors ----
col_RA  <- "#E63946"
col_HC  <- "#2E86AB"
col_acc <- "#A23B72"

## ---- Tracking ----
pipeline_log <- list()

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  MACHINE LEARNING PIPELINE — RA vs HC BIOMARKER DISCOVERY\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

## ============================================================================
## STEP 1: LOAD & PREPARE DATA
## ============================================================================

cat("── STEP 1: DATA LOADING ──────────────────────────────────────────────────\n\n")

input_ML <- read.csv("csvs/final input for ML.csv")

cat("  ✓ ML input shape  :", nrow(input_ML), "samples ×", ncol(input_ML), "columns\n")
cat("  ✓ Class distribution:\n")
print(table(input_ML$class))
cat("\n")

## X and Y
X     <- input_ML[, colnames(input_ML) != "class"]
Y_raw <- as.factor(input_ML$class)

## Scale X
X_scaled        <- as.data.frame(scale(X))
x_matrix        <- as.matrix(X_scaled)

## Recode: 1 → "RA", 0 → "HC"
Y <- ifelse(Y_raw == 1, "RA", "HC")
Y <- factor(Y, levels = c("HC", "RA"))

cat("  ✓ Features (genes) :", ncol(X_scaled), "\n")
cat("  ✓ RA samples       :", sum(Y == "RA"), "\n")
cat("  ✓ HC samples       :", sum(Y == "HC"), "\n\n")

pipeline_log$n_samples  <- nrow(input_ML)
pipeline_log$n_features <- ncol(X_scaled)
pipeline_log$n_RA       <- sum(Y == "RA")
pipeline_log$n_HC       <- sum(Y == "HC")

## ============================================================================
## STEP 2: LASSO FEATURE SELECTION
## ============================================================================

cat("── STEP 2: LASSO FEATURE SELECTION ──────────────────────────────────────\n\n")

set.seed(RANDOM_STATE)
lasso_model <- cv.glmnet(
  x_matrix, Y,
  family       = "binomial",
  alpha        = 1,
  type.measure = "class"
)

best_lambda <- lasso_model$lambda.min
lasso_coef  <- coef(lasso_model, s = best_lambda)
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]

cat("  ✓ Best lambda      :", round(best_lambda, 6), "\n")
cat("  ✓ LASSO selected genes (", length(lasso_genes), "):", paste(lasso_genes, collapse = ", "), "\n\n")

## Save
lasso_coef_df <- data.frame(
  Gene        = rownames(lasso_coef),
  Coefficient = as.numeric(lasso_coef)
) %>%
  dplyr::filter(Gene != "(Intercept)", Coefficient != 0) %>%
  dplyr::arrange(desc(abs(Coefficient)))

write.csv(lasso_coef_df, "tables/ML_LASSO_Selected_Genes.csv", row.names = FALSE)
cat("  ✓ Saved: tables/ML_LASSO_Selected_Genes.csv\n\n")

pipeline_log$lasso_genes <- lasso_genes

## ---- Plot: LASSO cross-validation curve
pdf("plots/machine_learning/01_LASSO_CV_Curve.pdf", width = 8, height = 5)
plot(lasso_model,
     main = "LASSO Cross-Validation (Binomial)",
     xlab = "Log Lambda",
     ylab = "Misclassification Error")
abline(v = log(lasso_model$lambda.min), col = col_RA, lty = 2, lwd = 2)
legend("topright", legend = paste("lambda.min =", round(best_lambda, 5)),
       col = col_RA, lty = 2, bty = "n")
dev.off()

png("plots/machine_learning/01_LASSO_CV_Curve.png", width = 8, height = 5,
    units = "in", res = 300)
plot(lasso_model,
     main = "LASSO Cross-Validation (Binomial)",
     xlab = "Log Lambda",
     ylab = "Misclassification Error")
abline(v = log(lasso_model$lambda.min), col = col_RA, lty = 2, lwd = 2)
legend("topright", legend = paste("lambda.min =", round(best_lambda, 5)),
       col = col_RA, lty = 2, bty = "n")
dev.off()
cat("  ✓ Saved: plots/machine_learning/01_LASSO_CV_Curve.png/.pdf\n")

## ---- Plot: LASSO coefficients
if (length(lasso_genes) > 0) {
  p <- ggplot(lasso_coef_df,
              aes(x = reorder(Gene, abs(Coefficient)),
                  y = Coefficient,
                  fill = Coefficient > 0)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = col_RA, "FALSE" = col_HC)) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          axis.title    = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(title    = "LASSO Selected Gene Coefficients",
         subtitle = paste("lambda.min =", round(best_lambda, 5)),
         x        = "Gene",
         y        = "Coefficient")

  pdf("plots/machine_learning/02_LASSO_Coefficients.pdf", width = 8, height = 5)
  print(p); dev.off()
  png("plots/machine_learning/02_LASSO_Coefficients.png",
      width = 8, height = 5, units = "in", res = 300)
  print(p); dev.off()
  cat("  ✓ Saved: plots/machine_learning/02_LASSO_Coefficients.png/.pdf\n\n")
}

## ============================================================================
## STEP 3: RANDOM FOREST FEATURE IMPORTANCE
## ============================================================================

cat("── STEP 3: RANDOM FOREST FEATURE IMPORTANCE ──────────────────────────────\n\n")

set.seed(RANDOM_STATE)
rf_model <- randomForest(
  x          = X_scaled,
  y          = Y,
  ntree      = 500,
  importance = TRUE
)

rf_importance <- importance(rf_model)
rf_genes_df   <- data.frame(
  Gene            = rownames(rf_importance),
  MeanDecreaseAcc = rf_importance[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini= rf_importance[, "MeanDecreaseGini"]
) %>%
  dplyr::arrange(desc(MeanDecreaseGini))

top_rf_genes <- rf_genes_df$Gene[1:8]

cat("  ✓ RF OOB error rate:", round(rf_model$err.rate[500, "OOB"], 3), "\n")
cat("  ✓ Top 8 RF genes   :", paste(top_rf_genes, collapse = ", "), "\n\n")

write.csv(rf_genes_df, "tables/ML_RF_Feature_Importance.csv", row.names = FALSE)
cat("  ✓ Saved: tables/ML_RF_Feature_Importance.csv\n\n")

pipeline_log$top_rf_genes <- top_rf_genes

## ---- Plot: RF importance (MeanDecreaseGini)
top_n   <- min(15, nrow(rf_genes_df))
rf_plot <- rf_genes_df[1:top_n, ]
rf_plot$Highlight <- rf_plot$Gene %in% top_rf_genes

p <- ggplot(rf_plot,
            aes(x = reorder(Gene, MeanDecreaseGini),
                y = MeanDecreaseGini,
                fill = Highlight)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = col_RA, "FALSE" = col_HC)) +
  coord_flip() +
  geom_hline(yintercept = mean(rf_plot$MeanDecreaseGini),
             linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        axis.title    = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(title    = "Random Forest — Feature Importance",
       subtitle = "Red = top 8 selected genes | Dashed = mean",
       x        = "Gene",
       y        = "Mean Decrease Gini")

pdf("plots/machine_learning/03_RF_Feature_Importance.pdf", width = 9, height = 6)
print(p); dev.off()
png("plots/machine_learning/03_RF_Feature_Importance.png",
    width = 9, height = 6, units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/03_RF_Feature_Importance.png/.pdf\n")

## ---- Plot: RF error rate convergence
rf_err_df <- data.frame(
  Trees = 1:500,
  OOB   = rf_model$err.rate[, "OOB"],
  RA    = rf_model$err.rate[, "RA"],
  HC    = rf_model$err.rate[, "HC"]
) %>%
  reshape2::melt(id.vars = "Trees", variable.name = "Class", value.name = "Error")

p2 <- ggplot(rf_err_df, aes(x = Trees, y = Error, color = Class)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("OOB" = "black", "RA" = col_RA, "HC" = col_HC)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")) +
  labs(title = "Random Forest — OOB Error Rate",
       x     = "Number of Trees",
       y     = "Error Rate",
       color = "Class")

pdf("plots/machine_learning/04_RF_OOB_Error.pdf", width = 8, height = 5)
print(p2); dev.off()
png("plots/machine_learning/04_RF_OOB_Error.png",
    width = 8, height = 5, units = "in", res = 300)
print(p2); dev.off()
cat("  ✓ Saved: plots/machine_learning/04_RF_OOB_Error.png/.pdf\n\n")

## ============================================================================
## STEP 4: SVM-RFE FEATURE SELECTION
## ============================================================================

cat("── STEP 4: SVM-RFE FEATURE SELECTION ─────────────────────────────────────\n\n")

set.seed(RANDOM_STATE)
svm_profile <- rfe(
  x          = X_scaled,
  y          = Y,
  sizes      = c(1:ncol(X_scaled)),
  rfeControl = rfeControl(functions = caretFuncs, method = "LOOCV"),
  method     = "svmLinear"
)

svm_genes   <- predictors(svm_profile)
svm_optimal <- svm_profile$optsize

cat("  ✓ SVM-RFE optimal features :", svm_optimal, "\n")
cat("  ✓ SVM-RFE selected genes   :", paste(svm_genes, collapse = ", "), "\n\n")

write.csv(data.frame(Gene = svm_genes),
          "tables/ML_SVM_RFE_Selected_Genes.csv",
          row.names = FALSE)
cat("  ✓ Saved: tables/ML_SVM_RFE_Selected_Genes.csv\n\n")

pipeline_log$svm_genes <- svm_genes

## ---- Plot: SVM-RFE accuracy vs number of variables
svm_res_df <- svm_profile$results

p <- ggplot(svm_res_df, aes(x = Variables, y = Accuracy)) +
  geom_line(color = col_HC, linewidth = 1) +
  geom_point(size = 3,
             color = ifelse(svm_res_df$Variables == svm_optimal, col_RA, col_HC)) +
  geom_vline(xintercept = svm_optimal, linetype = "dashed",
             color = col_RA, linewidth = 0.8) +
  annotate("text", x = svm_optimal + 0.3,
           y = min(svm_res_df$Accuracy) + 0.02,
           label = paste("Optimal =", svm_optimal),
           color = col_RA, hjust = 0, size = 4) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")) +
  labs(title = "SVM-RFE — LOOCV Accuracy vs Number of Features",
       x     = "Number of Variables",
       y     = "LOOCV Accuracy")

pdf("plots/machine_learning/05_SVM_RFE_Accuracy.pdf", width = 8, height = 5)
print(p); dev.off()
png("plots/machine_learning/05_SVM_RFE_Accuracy.png",
    width = 8, height = 5, units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/05_SVM_RFE_Accuracy.png/.pdf\n\n")

## ============================================================================
## STEP 5: FINAL BIOMARKERS — INTERSECTION RF ∩ SVM
## ============================================================================

cat("── STEP 5: FINAL BIOMARKERS ──────────────────────────────────────────────\n\n")

final_biomarkers <- Reduce(intersect, list(top_rf_genes, svm_genes))

if (length(final_biomarkers) == 0) {
  cat("  ! No intersection — using union of RF + SVM\n")
  final_biomarkers <- union(top_rf_genes, svm_genes)
}

cat("  ✓ Final biomarkers (", length(final_biomarkers), "):",
    paste(final_biomarkers, collapse = ", "), "\n\n")

write.csv(data.frame(Gene = final_biomarkers),
          "csvs/final_biomarkers_RA.csv",
          row.names = FALSE)
cat("  ✓ Saved: csvs/final_biomarkers_RA.csv\n\n")

pipeline_log$final_biomarkers <- final_biomarkers

## ---- Plot: Venn-style feature overlap summary
selection_summary <- data.frame(
  Gene      = unique(c(lasso_genes, top_rf_genes, svm_genes)),
  LASSO     = unique(c(lasso_genes, top_rf_genes, svm_genes)) %in% lasso_genes,
  RF        = unique(c(lasso_genes, top_rf_genes, svm_genes)) %in% top_rf_genes,
  SVM_RFE   = unique(c(lasso_genes, top_rf_genes, svm_genes)) %in% svm_genes
)
selection_summary$Methods <- rowSums(selection_summary[, c("LASSO", "RF", "SVM_RFE")])
selection_summary$Final   <- selection_summary$Gene %in% final_biomarkers

sel_long <- selection_summary %>%
  dplyr::select(Gene, LASSO, RF, SVM_RFE) %>%
  reshape2::melt(id.vars = "Gene", variable.name = "Method", value.name = "Selected") %>%
  dplyr::mutate(Selected = ifelse(Selected, "Yes", "No"))

p <- ggplot(sel_long, aes(x = Method, y = Gene, fill = Selected)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("Yes" = col_RA, "No" = "gray90")) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        axis.title    = element_text(face = "bold"),
        axis.text.x   = element_text(face = "bold"),
        legend.position = "right") +
  labs(title = "Feature Selection — Method Overlap",
       x     = "Selection Method",
       y     = "Gene",
       fill  = "Selected")

pdf("plots/machine_learning/06_Feature_Selection_Overlap.pdf", width = 7, height = max(5, nrow(selection_summary) * 0.4))
print(p); dev.off()
png("plots/machine_learning/06_Feature_Selection_Overlap.png",
    width = 7, height = max(5, nrow(selection_summary) * 0.4),
    units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/06_Feature_Selection_Overlap.png/.pdf\n\n")

## ============================================================================
## STEP 6: INDIVIDUAL GENE ROC ANALYSIS
## ============================================================================

cat("── STEP 6: INDIVIDUAL GENE ROC ANALYSIS ──────────────────────────────────\n\n")

roc_results <- list()
auc_table   <- data.frame()

for (gene in final_biomarkers) {
  roc_obj <- pROC::roc(Y, X_scaled[, gene],
                        quiet = TRUE, direction = "auto")
  auc_val          <- as.numeric(pROC::auc(roc_obj))
  roc_results[[gene]] <- list(roc = roc_obj, auc = auc_val)
  auc_table <- rbind(auc_table,
                     data.frame(Gene = gene, AUC = round(auc_val, 3)))
}

auc_table <- auc_table %>% dplyr::arrange(desc(AUC))
best_gene  <- auc_table$Gene[1]
best_auc   <- auc_table$AUC[1]

cat("  Gene AUC scores:\n")
print(auc_table, row.names = FALSE)

cat("\n  ✓ Best single gene:", best_gene, "(AUC =", best_auc, ")\n\n")

write.csv(auc_table, "tables/ML_Individual_Gene_AUC.csv", row.names = FALSE)
cat("  ✓ Saved: tables/ML_Individual_Gene_AUC.csv\n\n")

pipeline_log$best_gene     <- best_gene
pipeline_log$best_gene_auc <- best_auc

## ---- Plot: Individual gene AUC bar
p <- ggplot(auc_table,
            aes(x = reorder(Gene, AUC), y = AUC,
                fill = AUC)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  scale_fill_gradientn(colors = c(col_HC, "white", col_RA),
                       values = c(0, 0.5, 1)) +
  geom_hline(yintercept = 0.8, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "orange", linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "gray50", linewidth = 0.7) +
  geom_text(aes(label = round(AUC, 3)), hjust = -0.15, size = 4) +
  coord_flip(ylim = c(0, 1.08)) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        axis.title    = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(title    = "Individual Gene Diagnostic Performance",
       subtitle = "Dashed lines: AUC = 0.8 (red), 0.7 (orange), 0.5 (gray)",
       x        = "Gene",
       y        = "AUC")

pdf("plots/machine_learning/07_Individual_Gene_AUC.pdf", width = 8, height = 5)
print(p); dev.off()
png("plots/machine_learning/07_Individual_Gene_AUC.png",
    width = 8, height = 5, units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/07_Individual_Gene_AUC.png/.pdf\n")

## ---- Plot: Best gene ROC curve (matches R: plot(roc_curve, col="blue"))
best_roc <- roc_results[[best_gene]]$roc

pdf("plots/machine_learning/08_Best_Gene_ROC.pdf", width = 7, height = 6)
plot(best_roc,
     col      = "blue",
     lwd      = 2.5,
     main     = paste("ROC:", best_gene),
     xlab     = "False Positive Rate (1 - Specificity)",
     ylab     = "True Positive Rate (Sensitivity)",
     cex.main = 1.3,
     cex.lab  = 1.1)
polygon(c(rev(best_roc$specificities / 100),
          best_roc$specificities / 100),
        c(rev(best_roc$sensitivities / 100),
          rep(0, length(best_roc$sensitivities))),
        col = adjustcolor("blue", alpha.f = 0.10), border = NA)
legend("bottomright",
       legend = paste("AUC =", round(pROC::auc(best_roc), 3)),
       col = "blue", lwd = 2.5, bty = "n", cex = 1.1)
dev.off()

png("plots/machine_learning/08_Best_Gene_ROC.png",
    width = 7, height = 6, units = "in", res = 300)
plot(best_roc,
     col      = "blue",
     lwd      = 2.5,
     main     = paste("ROC:", best_gene),
     xlab     = "False Positive Rate (1 - Specificity)",
     ylab     = "True Positive Rate (Sensitivity)",
     cex.main = 1.3,
     cex.lab  = 1.1)
polygon(c(rev(best_roc$specificities / 100),
          best_roc$specificities / 100),
        c(rev(best_roc$sensitivities / 100),
          rep(0, length(best_roc$sensitivities))),
        col = adjustcolor("blue", alpha.f = 0.10), border = NA)
legend("bottomright",
       legend = paste("AUC =", round(pROC::auc(best_roc), 3)),
       col = "blue", lwd = 2.5, bty = "n", cex = 1.1)
dev.off()
cat("  ✓ Saved: plots/machine_learning/08_Best_Gene_ROC.png/.pdf\n")

## ---- Plot: All biomarker ROC curves overlaid
roc_colors <- colorRampPalette(c(col_RA, col_HC, col_acc, "#2A9D8F", "#E9C46A"))(
  length(final_biomarkers))

pdf("plots/machine_learning/09_All_Biomarkers_ROC.pdf", width = 8, height = 7)
plot(roc_results[[final_biomarkers[1]]]$roc,
     col      = roc_colors[1],
     lwd      = 2,
     main     = "ROC Curves — Final Biomarkers (RA vs HC)",
     cex.main = 1.2,
     cex.lab  = 1.1)
if (length(final_biomarkers) > 1) {
  for (i in 2:length(final_biomarkers)) {
    plot(roc_results[[final_biomarkers[i]]]$roc,
         col = roc_colors[i], lwd = 2, add = TRUE)
  }
}
abline(a = 0, b = 1, lty = 2, col = "gray60")
legend("bottomright",
       legend = paste0(auc_table$Gene, "  (AUC=", auc_table$AUC, ")"),
       col    = roc_colors[match(auc_table$Gene, final_biomarkers)],
       lwd    = 2, bty = "n", cex = 0.85)
dev.off()

png("plots/machine_learning/09_All_Biomarkers_ROC.png",
    width = 8, height = 7, units = "in", res = 300)
plot(roc_results[[final_biomarkers[1]]]$roc,
     col  = roc_colors[1], lwd = 2,
     main = "ROC Curves — Final Biomarkers (RA vs HC)",
     cex.main = 1.2, cex.lab = 1.1)
if (length(final_biomarkers) > 1) {
  for (i in 2:length(final_biomarkers)) {
    plot(roc_results[[final_biomarkers[i]]]$roc,
         col = roc_colors[i], lwd = 2, add = TRUE)
  }
}
abline(a = 0, b = 1, lty = 2, col = "gray60")
legend("bottomright",
       legend = paste0(auc_table$Gene, "  (AUC=", auc_table$AUC, ")"),
       col    = roc_colors[match(auc_table$Gene, final_biomarkers)],
       lwd    = 2, bty = "n", cex = 0.85)
dev.off()
cat("  ✓ Saved: plots/machine_learning/09_All_Biomarkers_ROC.png/.pdf\n\n")

## ============================================================================
## STEP 7: MULTI-MODEL LOOCV EVALUATION
## ============================================================================

cat("── STEP 7: MULTI-MODEL LOOCV EVALUATION ──────────────────────────────────\n\n")

ctrl_loo <- trainControl(
  method          = "LOOCV",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

feat_use <- if (length(final_biomarkers) > 0) final_biomarkers else colnames(X_scaled)
X_final  <- X_scaled[, feat_use, drop = FALSE]

model_results <- list()

## ---- Logistic Regression
cat("  Training: Logistic Regression ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  lr_model <- train(x = X_final, y = Y,
                    method    = "glm",
                    family    = "binomial",
                    trControl = ctrl_loo,
                    metric    = "ROC")
  model_results[["Logistic Regression"]] <- lr_model
  cat("    ✓ AUC:", round(max(lr_model$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

## ---- Random Forest
cat("  Training: Random Forest ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  rf_cv <- train(x = X_final, y = Y,
                 method    = "rf",
                 ntree     = 500,
                 trControl = ctrl_loo,
                 metric    = "ROC")
  model_results[["Random Forest"]] <- rf_cv
  cat("    ✓ AUC:", round(max(rf_cv$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

## ---- SVM Linear
cat("  Training: SVM (Linear) ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  svm_lin <- train(x = X_final, y = Y,
                   method    = "svmLinear",
                   trControl = ctrl_loo,
                   metric    = "ROC")
  model_results[["SVM (Linear)"]] <- svm_lin
  cat("    ✓ AUC:", round(max(svm_lin$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

## ---- SVM Radial
cat("  Training: SVM (Radial) ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  svm_rad <- train(x = X_final, y = Y,
                   method    = "svmRadial",
                   trControl = ctrl_loo,
                   metric    = "ROC")
  model_results[["SVM (Radial)"]] <- svm_rad
  cat("    ✓ AUC:", round(max(svm_rad$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

## ---- Naive Bayes
cat("  Training: Naive Bayes ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  nb_model <- train(x = X_final, y = Y,
                    method    = "nb",
                    trControl = ctrl_loo,
                    metric    = "ROC")
  model_results[["Naive Bayes"]] <- nb_model
  cat("    ✓ AUC:", round(max(nb_model$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

## ---- Gradient Boosting
cat("  Training: Gradient Boosting ...\n")
set.seed(RANDOM_STATE)
tryCatch({
  gbm_model <- train(x = X_final, y = Y,
                     method    = "gbm",
                     trControl = ctrl_loo,
                     metric    = "ROC",
                     verbose   = FALSE)
  model_results[["Gradient Boosting"]] <- gbm_model
  cat("    ✓ AUC:", round(max(gbm_model$results$ROC), 3), "\n")
}, error = function(e) cat("    ✗ Failed:", e$message, "\n"))

cat("\n")

## ---- Compile performance table
==========================================
  ## Compile model performance safely
  ## ==========================================

library(pROC)
library(caret)
library(dplyr)

perf_rows <- lapply(names(model_results), function(nm){
  
  m <- model_results[[nm]]
  
  if(is.null(m$pred)){
    cat("No predictions for:", nm, "\n")
    return(NULL)
  }
  
  preds <- m$pred
  
  ## For RF only
  if("mtry" %in% colnames(preds)){
    preds <- preds[preds$mtry == unique(preds$mtry)[1], ]
  }
  
  ## Check classes
  cat("\nModel:", nm, "\n")
  print(table(preds$obs))
  
  ## Skip if only one class
  if(length(unique(preds$obs)) < 2){
    cat("Skipping", nm, "- only one class present\n")
    return(NULL)
  }
  
  ## Determine positive class probability
  positive_class <- "RA"
  
  if(!(positive_class %in% colnames(preds))){
    cat("Skipping", nm, "- RA probability column missing\n")
    return(NULL)
  }
  
  ## Accuracy
  acc <- mean(preds$pred == preds$obs)
  
  ## ROC
  roc_obj <- pROC::roc(
    response = preds$obs,
    predictor = preds[[positive_class]],
    levels = c("HC","RA"),
    direction = "<",
    quiet = TRUE
  )
  
  auc_v <- as.numeric(pROC::auc(roc_obj))
  
  ## Confusion matrix
  cm <- confusionMatrix(
    preds$pred,
    preds$obs,
    positive = "RA"
  )
  
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  data.frame(
    Model       = nm,
    Accuracy    = round(acc,3),
    AUC         = round(auc_v,3),
    Sensitivity = round(sens,3),
    Specificity = round(spec,3)
  )
})

## Combine
perf_df <- do.call(
  rbind,
  Filter(Negate(is.null), perf_rows)
)

## Rank by AUC
perf_df <- perf_df %>%
  arrange(desc(AUC))

print(perf_df)

## Save
write.csv(
  perf_df,
  file.path(wdir, "csvs", "ML_Model_Performance.csv"),
  row.names = FALSE
)

## Store in pipeline log
pipeline_log$model_results <- perf_df
pipeline_log$best_model <- perf_df$Model[1]



## ============================================================================
## STEP 8: MODEL EVALUATION PLOTS
## ============================================================================

cat("\n── STEP 8: MODEL EVALUATION PLOTS ───────────────────────────────────────\n\n")

## ---- Plot: Multi-model ROC curves
roc_colors_m <- colorRampPalette(
  c(col_RA, col_HC, col_acc, "#2A9D8F", "#E9C46A", "#264653")
)(length(model_results))

cat("\n── FIXED ROC PLOT GENERATION ───────────────────────────────\n")

library(pROC)

roc_colors_m <- colorRampPalette(
  c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")
)(length(model_results))

pdf("plots/machine_learning/10_All_Models_ROC_FIXED.pdf", 9, 7)

plot(1, type = "n",
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC Curves — All Models (RA vs HC)\nLOOCV")

abline(0,1,lty=2,col="gray70")

legend_labels <- c()
valid_colors <- c()
i_color <- 1

for(nm in names(model_results)) {
  
  preds <- model_results[[nm]]$pred
  
  if(is.null(preds)) next
  
  ## Check class distribution
  if(length(unique(preds$obs)) < 2) {
    cat("Skipping", nm, "- only one class\n")
    next
  }
  
  ## Find probability column safely
  prob_col <- intersect(c("RA","HC"), colnames(preds))
  
  if(length(prob_col) == 0) {
    cat("Skipping", nm, "- no RA probability column\n")
    next
  }
  
  roc_obj <- tryCatch({
    pROC::roc(
      response = preds$obs,
      predictor = preds[[prob_col[1]]],
      levels = c("HC","RA"),
      direction = "<"
    )
  }, error = function(e) NULL)
  
  if(is.null(roc_obj)) next
  
  lines(1 - roc_obj$specificities,
        roc_obj$sensitivities,
        col = roc_colors_m[i_color],
        lwd = 2)
  
  legend_labels <- c(legend_labels,
                     paste0(nm, " (AUC=", round(auc(roc_obj),3), ")"))
  
  valid_colors <- c(valid_colors, roc_colors_m[i_color])
  i_color <- i_color + 1
}

legend("bottomright",
       legend = legend_labels,
       col = valid_colors,
       lwd = 2,
       bty = "n",
       cex = 0.9)

dev.off()

cat("✓ ROC plot saved successfully\n")

## ---- Plot: Model comparison (Accuracy + AUC + Sensitivity + Specificity)
perf_long <- perf_df %>%
  reshape2::melt(id.vars     = "Model",
                 variable.name = "Metric",
                 value.name    = "Score")

p <- ggplot(perf_long, aes(x = reorder(Model, -Score), y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.7, alpha = 0.9) +
  geom_text(aes(label = round(Score, 2)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("Accuracy"    = col_HC,
                                "AUC"         = col_RA,
                                "Sensitivity" = col_acc,
                                "Specificity" = "#2A9D8F")) +
  geom_hline(yintercept = 0.8, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        axis.title    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        legend.position = "top") +
  ylim(0, 1.15) +
  labs(title = "Model Comparison — LOOCV Performance (RA vs HC)",
       x     = "Model",
       y     = "Score",
       fill  = "Metric")

pdf("plots/machine_learning/11_Model_Comparison.pdf", width = 12, height = 6)
print(p); dev.off()
png("plots/machine_learning/11_Model_Comparison.png",
    width = 12, height = 6, units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/11_Model_Comparison.png/.pdf\n")

## ---- Plot: Confusion matrices (best model)
best_model_obj <- model_results[[pipeline_log$best_model]]
if (!is.null(best_model_obj$pred)) {
  preds_best <- best_model_obj$pred
  cm_best    <- confusionMatrix(preds_best$pred, preds_best$obs, positive = "RA")
  cm_mat     <- as.data.frame(cm_best$table)

  p <- ggplot(cm_mat, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = Freq), size = 10, fontface = "bold") +
    scale_fill_gradient(low = "white", high = col_RA) +
    theme_bw(base_size = 13) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title      = element_text(face = "bold"),
          legend.position = "none") +
    labs(title    = paste("Confusion Matrix —", pipeline_log$best_model),
         subtitle = paste("LOOCV | AUC =", 1),
         x        = "True Class",
         y        = "Predicted Class")

  pdf("plots/machine_learning/12_Confusion_Matrix_Best_Model.pdf",
      width = 5, height = 4.5)
  print(p); dev.off()
  png("plots/machine_learning/12_Confusion_Matrix_Best_Model.png",
      width = 5, height = 4.5, units = "in", res = 300)
  print(p); dev.off()
  cat("  ✓ Saved: plots/machine_learning/12_Confusion_Matrix_Best_Model.png/.pdf\n")
}

## ---- Plot: Expression heatmap of final biomarkers
expr_bio <- t(X_scaled[, feat_use, drop = FALSE])
ann_col  <- data.frame(Class = Y, row.names = colnames(expr_bio))

colnames(ann)

pdf("plots/machine_learning/13_Biomarker_Expression_Heatmap.pdf",
    width = 10, height = max(4, length(feat_use) * 0.6))
pheatmap::pheatmap(
  expr_bio,
  annotation_col  = ann_col,
  annotation_colors = list(Class = c("RA" = col_RA, "HC" = col_HC)),
  color           = colorRampPalette(c(col_HC, "white", col_RA))(50),
  breaks          = seq(-3, 3, length.out = 51),
  scale           = "none",
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_colnames   = FALSE,
  fontsize        = 11,
  cellheight      = 22,
  main            = "Final Biomarker Expression (Z-scored)\nRA vs HC"
)
dev.off()

png("plots/machine_learning/13_Biomarker_Expression_Heatmap.png",
    width = 10, height = max(4, length(feat_use) * 0.6),
    units = "in", res = 300)
pheatmap::pheatmap(
  expr_bio,
  annotation_col  = ann_col,
  annotation_colors = list(Class = c("RA" = col_RA, "HC" = col_HC)),
  color           = colorRampPalette(c(col_HC, "white", col_RA))(50),
  breaks          = seq(-3, 3, length.out = 51),
  scale           = "none",
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_colnames   = FALSE,
  fontsize        = 11,
  cellheight      = 22,
  main            = "Final Biomarker Expression (Z-scored)\nRA vs HC"
)
dev.off()
cat("  ✓ Saved: plots/machine_learning/13_Biomarker_Expression_Heatmap.png/.pdf\n")

## ---- Plot: Boxplots of each biomarker by class
box_df <- X_scaled[, feat_use, drop = FALSE]
box_df$Class <- Y
box_long <- reshape2::melt(box_df, id.vars = "Class",
                            variable.name = "Gene",
                            value.name    = "Expression")

p <- ggplot(box_long, aes(x = Class, y = Expression, fill = Class)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2, width = 0.6) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7,
              aes(color = Class)) +
  scale_fill_manual(values  = c("RA" = col_RA, "HC" = col_HC)) +
  scale_color_manual(values = c("RA" = col_RA, "HC" = col_HC)) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5),
        axis.title      = element_text(face = "bold"),
        strip.text      = element_text(face = "bold", size = 10),
        legend.position = "none") +
  labs(title = "Final Biomarkers — Expression by Class (RA vs HC)",
       x     = "Class",
       y     = "Scaled Expression (Z-score)")

plot_h <- max(5, ceiling(length(feat_use) / 3) * 3)
pdf("plots/machine_learning/14_Biomarker_Boxplots.pdf",
    width = 12, height = plot_h)
print(p); dev.off()
png("plots/machine_learning/14_Biomarker_Boxplots.png",
    width = 12, height = plot_h, units = "in", res = 300)
print(p); dev.off()
cat("  ✓ Saved: plots/machine_learning/14_Biomarker_Boxplots.png/.pdf\n\n")

## ============================================================================
## FINAL SUMMARY
## ============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║       ML PIPELINE SUMMARY — RA vs HC BIOMARKER DISCOVERY             ║\n")
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Input samples              : %-6d                              ║\n", pipeline_log$n_samples))
cat(sprintf("║  Input features (hub genes) : %-6d                              ║\n", pipeline_log$n_features))
cat(sprintf("║  RA samples                 : %-6d                              ║\n", pipeline_log$n_RA))
cat(sprintf("║  HC samples                 : %-6d                              ║\n", pipeline_log$n_HC))
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  FEATURE SELECTION                                                    ║\n")
cat(sprintf("║    LASSO genes              : %d → %s\n",
    length(lasso_genes), paste(lasso_genes, collapse=", ")) %>%
    sprintf("║  %-71s║\n", sub("║  ", "", .)) %>% cat())
cat(sprintf("║    RF top 8 genes           : %s\n",
    paste(top_rf_genes, collapse=", ")) %>%
    sprintf("║  %-71s║\n", sub("║  ", "", .)) %>% cat())
cat(sprintf("║    SVM-RFE genes (%d)        : %s\n",
    length(svm_genes), paste(svm_genes, collapse=", ")) %>%
    sprintf("║  %-71s║\n", sub("║  ", "", .)) %>% cat())
cat(sprintf("║    Final biomarkers (%d)     : %s\n",
    length(final_biomarkers), paste(final_biomarkers, collapse=", ")) %>%
    sprintf("║  %-71s║\n", sub("║  ", "", .)) %>% cat())
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  INDIVIDUAL GENE AUC                                                  ║\n")
for (i in seq_len(nrow(auc_table))) {
  cat(sprintf("║    %-25s AUC = %-6.3f                          ║\n",
              auc_table$Gene[i], auc_table$AUC[i]))
}
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  MODEL PERFORMANCE (LOOCV)                                            ║\n")
for (i in seq_len(nrow(perf_df))) {
  cat(sprintf("║    %-22s Acc:%.3f AUC:%.3f Sens:%.3f Spec:%.3f  ║\n",
              perf_df$Model[i], perf_df$Accuracy[i],
              perf_df$AUC[i], perf_df$Sensitivity[i], perf_df$Specificity[i]))
}
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Best single gene : %-20s (AUC = %.3f)               ║\n",
    pipeline_log$best_gene, pipeline_log$best_gene_auc))
cat(sprintf("║  Best model       : %-20s (AUC = %.3f)               ║\n",
    pipeline_log$best_model, pipeline_log$best_model_auc))
cat("╠═══════════════════════════════════════════════════════════════════════╣\n")
cat("║  FILES SAVED                                                          ║\n")
cat("║  csvs/  final_biomarkers_RA.csv                                       ║\n")
cat("║  tables/ ML_LASSO_Selected_Genes.csv                                  ║\n")
cat("║           ML_RF_Feature_Importance.csv                                ║\n")
cat("║           ML_SVM_RFE_Selected_Genes.csv                               ║\n")
cat("║           ML_Individual_Gene_AUC.csv                                  ║\n")
cat("║           ML_Model_Performance_LOOCV.csv                              ║\n")
cat("║  plots/machine_learning/                                              ║\n")
cat("║    01_LASSO_CV_Curve            02_LASSO_Coefficients                 ║\n")
cat("║    03_RF_Feature_Importance     04_RF_OOB_Error                       ║\n")
cat("║    05_SVM_RFE_Accuracy          06_Feature_Selection_Overlap          ║\n")
cat("║    07_Individual_Gene_AUC       08_Best_Gene_ROC                      ║\n")
cat("║    09_All_Biomarkers_ROC        10_All_Models_ROC_LOOCV               ║\n")
cat("║    11_Model_Comparison          12_Confusion_Matrix_Best_Model        ║\n")
cat("║    13_Biomarker_Expression_Heatmap  14_Biomarker_Boxplots             ║\n")
cat("║    (all saved as .png + .pdf)                                         ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("  Completed: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# --------------------------------------------------
# Script: Differential Expression Analysis (limma-voom)
# Description: Performs differential gene expression analysis across bacterial infection datasets using limma and edgeR
# ----------------------------------------------------
# Load required libraries
library(limma)
library(edgeR)
# Define condition vector
condition <- factor(c("control", "control", "control", "infected", "infected", "infected"))
#Design matrix
design <- model.matrix(~ condition)
#Function for DE analysis
run_limma <- function (count_matrix, output_name) {
  #Create DGEList
  dge <- DGEList(counts)
  # Normalize counts
  dge <- calcNormFactors(dge)
  # Voom transformation
  v <- voom(dge, design, plot=TRUE)
  #Linear modeling
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  # Extract results
  results <- topTable(fit, coef=2, number=Inf)
  # Filter significant DEGs
  deg <- results[results$adj.P.Val < 0.05 &abs(results$logFC) > 1,]
  #Save results
  write.csv(deg, paste0("DEG_results_", output_name,"_limma.csv"))
}


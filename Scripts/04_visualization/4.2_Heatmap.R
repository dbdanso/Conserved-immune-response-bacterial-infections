# -----------------------------------------------
# Script: Heatmap of Conserved Genes
#Description: Visualizes expression patterns (logFC) of conserved genes across bacterial infections.
# -----------------------------------------------
# Load required library
library(pheatmap)
# Rename columns of the merged data obtained from the conserved_gene Script with infections
colnames(merged)[2:5] <- c(
  "A. veronii AS1",
  "E. coli K12",
  "L. pneumophila 130b",
  "M. tuberculosis H37Rv"
)
# Select top variable genes based on mean absolute logFC across infections
top30 <- merged[
  order(rowMeans(abs(merged[, -1])), decreasing = TRUE),
]
[1:30, ]
# Prepare matrix for heatmap
mat_top <- as.matrix(top30[, -1]) # remove gene_symbol column
rownames(mat_top) <- top30$gene_symbol
# Scale data (Z-score per gene)
mat_scaled <- t(scale(t(mat_top)))
# Plot heatmap
pheatmap(
  mat_scaled,show_rownames = TRUE,fontsize_row = 7,
  fontsize_col = 10,main = "Top 30 Conserved Genes"
)

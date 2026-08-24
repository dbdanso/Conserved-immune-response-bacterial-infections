# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)
# Run Helper function 3 (KEGG enrichment) and extract results
res_kegg <- Kegg_up_enrichment(merged)
kegg_up <- res_kegg$kegg_up
kegg_up_read <- res_kegg$kegg_up_read
# Compute average fold change across all infections
gene_fc_values <- rowMeans(merged[, 2:5], na.rm = TRUE)
# Assign gene symbols as names to the fold change vector
names(gene_fc_values) <- merged$gene_symbol
# Convert gene symbols to ENTREZ IDs for compatibility with KEGG analysis
gene_df <- bitr(
  merged$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Align fold change values with the converted ENTREZ IDs
gene_fc <- gene_fc_values[gene_df$SYMBOL]
# Replace gene symbol names with ENTRZ IDs
names(gene_fc) <- gene_df$ENTREZID
# Extract all genes involved in enriched KEGG pathways
kegg_genes <- unique(unlist(strsplit(kegg_up@result$geneID, "/")))
# Filter fold change vector to include only genes present in KEGG pathways
gene_fc <- gene_fc[names(gene_fc) %in% kegg_genes]
# Generate cnetplot
p <- enrichplot::cnetplot(
  x = kegg_up_read, showCategory = 5,foldChange = gene_fc
)
p + scale_color_gradient2(low = "blue",high = "red")




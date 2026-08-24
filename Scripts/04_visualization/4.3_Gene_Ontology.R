# -----------------------------------------
# Script: Gene Ontology (GO) Enrichment Analysis
# Description: Identifies enriched biological processes among conserved genes across infections.
# ----------------------------------------------
# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)
# Helper function 2: Preparation of upregulated genes from merged conserved genes data
prepare_upregulated_genes <- function(merged_df) {
  # Rename columns
  colnames(merged_df)[2:5] <- c(
  "A. veronii AS1",
  "E. coli K12",
  "L. pneumophila 130b",
  "M. tuberculosis H37Rv"
)
  # Extract conserved upregulated genes
  up_all4 <- merged_df[rowSums(merged_df[, 2:5] > 0) == 4,]
  # Extract gene symbols
  gene_symbols <- up_all4$gene_symbol
  # Convert to ENTREZ IDs
  gene_df <- bitr(
    gene_symbols, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene_list <- gene_df$ENTREZID
  return(
    list(upregulated_table = up_all4, gene_symbols = gene_symbols, entrez_id = gene_list)
)
}
gene_list <- prepare_upregulated_genes(merged)

go_up <- enrichGO(gene= gene_list,OrgDb = org.Hs.eg.db,
  ont = "BP",         
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
go_up <- simplify(
  go_up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min)
go_up_df <- as.data.frame (go_up)
write.csv(go_up_df, "Upregulated_GeneOntology")
df_top <- as.data.frame(go_up)
df_top <- df_top %>%slice_min(order_by = p.adjust, n = 11, with_ties = FALSE)
df_top <- df_top%>%arrange(p.adjust) %>%slice(1:11)
df_top$GeneRatio <- sapply(df_top$GeneRatio, function(x) {
  eval(parse(text = x))
})

df_top$Description <- factor(
  df_top$Description,
  levels = rev(df_top$Description))
  ggplot(df_top,
       aes(x = GeneRatio,
           y = Description,
           size = Count,
           color = p.adjust)) +
  
  geom_point() +
  
  scale_x_continuous(limits = c(0, max(df_top$GeneRatio))) +
  
  scale_size(range = c(2, 5)) +
  
  scale_color_gradientn(colors = c(high = "blue",low = "red"),
    trans = "reverse") +
  
  labs(
    x = "Gene Ratio",
    y = "",
    title = "Top 11 Activated Biological Processes") +
  
  theme_minimal() +
  
  theme(axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5, face = "bold"))

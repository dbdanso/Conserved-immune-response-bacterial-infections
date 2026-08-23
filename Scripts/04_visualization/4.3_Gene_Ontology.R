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

colnames(merged)[2:5] <- c(
  "A. veronii AS1",
  "E. coli K12",
  "L. pneumophila 130b",
  "M. tuberculosis H37Rv")
up_all4 <- merged[rowSums(merged[, c("A. veronii AS1","E. coli K12","L. pneumophila 130b","M. tuberculosis H37Rv")] > 0) == 4,]
down_all4 <- merged[rowSums(merged[, c("A. veronii AS1","E. coli K12","L. pneumophila 130b","M. tuberculosis H37Rv")] < 0) == 4,]
write.csv(up_all4[,1], "conserved_upregulated_genes.csv")
write.csv(down_all4[,1], "conserved_downregulated_genes.csv")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)   
library(enrichplot)
library(ggplot2)
library(dplyr)
conserved_df_up <- up_all4$gene_symbol
gene_df_up <- bitr(conserved_df_up, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene_list_up <- gene_df_up$ENTREZID

go_up <- enrichGO(gene= gene_list_up,OrgDb = org.Hs.eg.db,
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

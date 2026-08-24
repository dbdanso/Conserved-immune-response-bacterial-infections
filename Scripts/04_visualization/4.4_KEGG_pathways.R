# ---------------------------------------------
# Script: KEGG Pathway Enrichment Analysis
# Description: Identifies enriched KEGG pathways among conserved genes
# ---------------------------------------------
# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.bd)
library(enrichplot)
library(ggplot2)
library(dplyr)
# Helper function 3: Run KEGG enrichment (upregulated genes)
Kegg_up_enrichment <- function(merged_df) {
  # Prepare gene list using Helper function 2
  gene_list <- prepare_upregulated_genes(merged)
  # KEGG enrichment
  kegg_up <- enrichKEGG(
  gene = gene_list, organism = "hsa",
  pvalueCutoff = 0.05
)
  # Convert ENTREZ IDs to gene symbols
  kegg_up <- setReadable(kegg_up,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
  # Convert to dataframe
  df_kegg <- as.data.frame(kegg_up)
  return(df_kegg)
}
df_kegg <- Kegg_up_enrichment(merged)
# Save full results
write.csv(df_kegg, "KEGG_upregulated.csv", row.names = FALSE)
# Select top pathways
df_kegg_top <- df_kegg %>%arrange(p.adjust) %>%slice(1:11)
# Convert GeneRatio to numeric
df_kegg_top$GeneRatio <- sapply(
  df_kegg_top$GeneRatio, function(x) 
{
  parts <-strsplit(x, "/")[[1]]
  as.numeric(parts[1])/ as.numeric(parts[2])
}
)
# Order for plotting
df_kegg_top$Description <- factor(
  df_kegg_top$Description,
  levels = rev(df_kegg_top$Description)
)
# Plot KEGG dotplot
ggplot(df_kegg_top,
       aes(x = GeneRatio, y = Description,
           size = Count,color = p.adjust)
) +
geom_point(
) +
scale_x_continuous(limits = c(0, max(df_kegg_top$GeneRatio))
) +
scale_size(range = c(2, 5)
) +
scale_color_gradientn(
    colors = c(high = "blue",low = "red"),
    trans = "reverse"
) +
labs(
    x = "Gene Ratio", y = "",
    title = "Top 11 Activated KEGG Pathways"
) +
  
theme_minimal(
) +
theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
)




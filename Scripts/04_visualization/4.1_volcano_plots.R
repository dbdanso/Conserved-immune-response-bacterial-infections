# ----------------------------------------------
# Script: Volcano Plot Visualization
#Description: Generates volcano plots for gene expression across bacterial infections
# ----------------------------------------------
# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
# Volcano plot function
plot_volcano <- function(results_df, title_name) {
  # Classify genes
  results_df$regulation <- "Not Significant"
  results_df$regulation[
    results_df$adj.P.Val < 0.05 & results_df$logFC > 1
] <- "Upregulated"
  results_df$regulation[
    results_df$adj.P.Val < 0.05 & results_df$logFC < -1
  ] <- "Downregulated"
  # Count categories
  counts <- table(results_df$regulation)
  up_count   <- counts["Upregulated"]
  down_count <- counts["Downregulated"]
  ns_count   <- counts["Not Significant"]
  # Replace NA counts with 0
  up_count <- ifelse(is.na(up_count), 0, up_count)
  down_count <- ifelse(is.na(down_count), 0, down_count)
  ns_count <- ifelse(is.na(ns_count), 0, ns_count)
  # create a column named X in results_df
  results_df$X <- rownames(results_df)
  # Convert gene symbols (Use helper function in conserved_gene script)
  results_df <- convert_symbols(results_df)
  # Select top genes for labeling
  top_up <- results_df %>% filter(logFC > 1 & adj.P.Val < 0.05)
  %>%arrange(desc(logFC)) %>%head(4)
  top_down <- results_df %>% filter(logFC < -1 & adj.P.Val < 0.05) %>%
  arrange(logFC) %>%head(4)
  top_ns <- results_df %>% filter(abs(logFC) < 1 & adj.P.Val > 0.05)
  %>%arrange(desc(logFC)) %>%head(4)
  top_genes <- rbind(top_up, top_down, top_ns)
  #Generate plot
  p <- ggplot(
    results_df, aes(x = logFC, y = -log10(adj.P.Val))
    )+
  geom_point(
    aes(color = regulation), alpha = 0.7
      )+
  geom_vline(
    xintercept = c(-1, 1), linetype = "dashed"
        )+
  geom_hline(
    yintercept = -log10(0.05), linetype = "dashed"
          )+ 
  theme_minimal(
            )+
  scale_color_manual(
    values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "green")
          )+ 
  geom_label_repel(
    data = top_genes,aes(label = gene_symbol),size = 1.5,box.padding = 0.5, point.padding = 0.1, fill = "white", segment.color = "black", segment.size = 0.5, max.overlaps = 15
        )+
  labs(
    title = title_name, x = "Log 2 Fold Change",y = "-log10(adjusted p
    value)",color = "Gene Regulation"
      )+ 
  theme_minimal()
  # Add annotation text
  p <- p + annotate(
    "text", x = max(results_df$logFC), y = 8,
           label= paste("Upregulated\n(", up_count,")"),
           hjust = 1, color = "red", size = 2.5
   )+
  annotate(
    "text", x = min(results_df$logFC), y = 8, 
           label = paste("Downregulated\n(", down_count,")"),
           hjust = 0, color = "blue", size = 2.5
)+
  annotate(
    "text", x = 0, y = 8,
           label = paste("NS\n(", ns_count,")"),
           hjust = 0.5, color = "green", size = 2.5
)
  return(p)
}


  
                        

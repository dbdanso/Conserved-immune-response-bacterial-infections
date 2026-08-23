# ----------------------------------------------
# Script: Conserved Gene Identification
# Description: Identifies differentially expressed genes common across four bacterial infections and merges their logFC values
# -----------------------------------------------
# Load required libraries
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
#Load DEG results
deg_AV <- read.csv("DEG_results_AV_limma.csv")
deg_EC <- read.csv("DEG_results_EC_limma.csv")
deg_LP <- read.csv("DEG_results_LP_limma.csv")
deg_MT <- read.csv("DEG_results_MT_limma.csv")
# Helper function 1: Convert Ensembl IDs to gene symbols
convert_symbols <- function(df) {
  #Ensure gene IDs are characters
  df$X <- as.character(df$X)
  # Remove version numbers (e.g., ENSG000001.1 -> ENSG000001)
  df$X <- sub("\\..*", "", df$X)
  # Map to gene symbols
  symbols <- mapIds(
    org.Hs.eg.db,keys = df$X, column = "SYMBOL",keytype =
    "ENSEMBL",multiVals = "first")
  df$gene_symbol <- symbols
  # Remove unmapped genes
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", ]
  return(df)
}
# Apply conversion
deg_AV <- convert_symbols(deg_AV)
deg_EC <- convert_symbols(deg_EC)
deg_LP <- convert_symbols(deg_LP)
# MT already contains gene symbols in the first column
colnames(deg_MT)[1] <- "gene_symbol"
# Extract unique gene lists
genes_AV <- unique(deg_AV$gene_symbol)
genes_EC <- unique(deg_EC$gene_symbol)
genes_LG <- unique(deg_LP$gene_symbol)
genes_MT <- unique(deg_MT$gene_symbol)
# Identify conserved genes
common_genes <- Reduce(
  intersect, list(genes_AV, genes_EC, genes_LP, genes_MT)
)
# Number of conserved genes
length(common_genes)
#Save list
write.csv(common_genes, "conserved_genes.csv", row.names = FALSE)
# Subset DEG tables to conserved genes
common_AV <- deg_AV[deg_AV$gene_symbol %in% common_genes, ]
common_EC <- deg_EC[deg_EC$gene_symbol %in% common_genes, ]
common_LG <- deg_LP[deg_LP$gene_symbol %in% common_genes, ]
common_MT <- deg_MT[deg_MT$gene_symbol %in% common_genes, ]
# Rename logFC columns
colnames(common_AV)[colnames(common_AV) == "logFC"] <- "logFC_AV"
colnames(common_EC)[colnames(common_EC) == "logFC"] <- "logFC_EC"
colnames(common_LP)[colnames(common_LP) == "logFC"] <- "logFC_LP"
colnames(common_MT)[colnames(common_MT) == "logFC"] <- "logFC_MT"
# Merge logFC values across infections
common1 <- common_AV[, c("gene_symbol", "logFC_AV")]
common2 <- common_EC[, c("gene_symbol", "logFC_EC")]
common3 <- common_LG[, c("gene_symbol", "logFC_LP")]
common4 <- common_MT[, c("gene_symbol", "logFC_MT")]
merged <- Reduce(
  function(x, y) merge(x, y, by="gene_symbol"),list(common1,common2,common3,   common4)
)
# Check for duplicates
any(duplicated(merged$gene_symbol))
                 

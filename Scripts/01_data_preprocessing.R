# -----------------------------------------------
# Script: Data Preprocessing (E. coli)
# Description: Loads raw count data and extracts THP-1 macrophage samples for control and E. coli-infected conditions
# -----------------------------------------------
# Load required library
library(readr)
# Load raw count dataset
raw_data <- read_tsv (data/"GSE273835_Raw_counts_for_THP.tsv")
# Convert to data frame (optional, for compatibility)
raw_data <- as.data.frame (raw_data)
# ------------------------------------------------
# Extract relevant samples
# THP_M = Control macrophages
# THP_E = E. coli-infected macrophages
counts_EC <- raw_data[, c("THP_M_rpt1","THP_M_rpt2","THP_M_rpt3", "THP_E_rpt1","THP_E_rpt2","THP_E_rpt3")]
# -----------------------------------------------
# Assign gene IDs as rownames
rownames(counts_EC) <- raw_data$gene_id

# ---------------------------------------------------
# Script: Data Preprocessing (M. tuberculosis)
# Description: Loads raw count data and extracts THP-1 macrophage samples for control and M. tuberculosis-infected conditions
# -------------------------------------------------
#Load required library
library(readxl)
#Load raw count dataset

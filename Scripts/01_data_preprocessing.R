# -----------------------------------------------
# Script: Data Preprocessing (E. coli)
# Description: Loads raw count data and extracts THP-1 macrophage samples for control and E. coli-infected conditions
# -----------------------------------------------
# Load required library
library(readr)
# Load raw count dataset
raw_data <- read_tsv ("GSE273835_Raw_counts_for_THP.tsv")
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
# Script: Data Preprocessing (A. veronii)
# Description: Loads raw count data and extracts THP-1 macrophage samples for control and A. veronii-infected conditions
# -----------------------------------------------------
library(readr)
raw_data <- read_tsv("GSE273835_Raw_counts_for_THP.tsv")
raw_data <- as.data.frame(raw_data)
counts_AV <- raw_data[, c("THP_M_rpt1","THP_M_rpt2","THP_M_rpt3", "THP_AV_rpt1","THP_AV_rpt2","THP_AV_rpt3")]
rownames(counts_AV) <- raw_data$gene_id

# ---------------------------------------------------
# Script: Data Preprocessing (M. tuberculosis)
# Description: Loads raw count data and extracts THP-1 macrophage samples for control and M. tuberculosis-infected conditions
# -------------------------------------------------
library(readxl)
raw_data <- read_excel("GSE275580_Complete_Read_Count_Data.xlsx")
raw_data <- as.data.frame(raw_data)
counts_MT <- raw_data[, c("Non_Treat_1","Non_Treat_3","Non_Treat_4", "Infected_2","Infected_3","Infected_4")]
rownames(counts_MT) <- make.unique(as.character(raw_data$gene_id))

# ------------------------------------------------------------------
# Script: Data Preprocessing (L. pneumophila)
# Description: Loads count data from separate control and L. pneumophila-infected samples, and merges into single count matrix
# --------------------------------------------------------------------
# Load required library
library(readr)
# Load datasets
ctrl1 <- read_tsv("GSM8750874_UI2.tabular")
ctrl2 <- read_tsv("GSM8750875_UI3.tabular")
ctrl3 <- read_tsv("GSM8750876_UI4.tabular")
inf1  <- read_tsv("GSM8750865_WT2.tabular")
inf2  <- read_tsv("GSM8750866_WT3.tabular")
inf3  <- read_tsv("GSM8750867_WT4.tabular")
# Rename count columns
# Second column contains read counts
colnames(ctrl1)[2] <- "ctrl1"
colnames(ctrl2)[2] <- "ctrl2"
colnames(ctrl3)[2] <- "ctrl3"
colnames(inf1)[2]  <- "inf1"
colnames(inf2)[2]  <- "inf2"
colnames(inf3)[2]  <- "inf3"
# Merge datasets by Gene ID
merged <- ctrl1
merged <- merge(merged, ctrl2, by = "Geneid")
merged <- merge(merged, ctrl3, by = "Geneid")
merged <- merge(merged, inf1, by = "Geneid")
merged <- merge(merged, inf2, by = "Geneid")
merged <- merge(merged, inf3, by = "Geneid")
# Create count matrix
rownames(merged) <- merged$Geneid
counts_LP <- merged[, -1] # remove Geneid column
counts_LP <- as.matrix(counts_LP)
# Preview
head(counts_LP)

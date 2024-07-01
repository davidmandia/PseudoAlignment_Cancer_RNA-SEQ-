# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load necessary libraries
library(biomaRt)

# Read in differential expression results
results <- read.csv("differential_expression_results.csv")

# Print summary of the data frame
summary(as.data.frame(results))

# Filter for significantly differentially expressed genes (optional)
sig_genes <- subset(results, padj < 0.05)

# Print significant genes
print(sig_genes)

# Extract gene names (assuming 'X' contains gene names)
gene_list <- sig_genes$X

# Initialize an empty data frame to store gene IDs
gene_conversion <- data.frame()

# Convert Ensembl Transcript IDs to Entrez Gene IDs using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(attributes = c('ensembl_transcript_id_version', 'entrezgene_id'),
                         filters = 'ensembl_transcript_id_version',
                         values = gene_list,
                         mart = ensembl)

# Merge gene conversion data with results based on ensembl_transcript_id_version
results_with_ids <- merge(results, gene_conversion, by.x = "X", by.y = "ensembl_transcript_id_version", all.x = TRUE)

# Save the updated data frame with Entrez Gene IDs as CSV
write.csv(results_with_ids, "differential_expression_results_with_ids.csv", row.names = FALSE)

# Print the first few rows of the updated data frame
head(results_with_ids)

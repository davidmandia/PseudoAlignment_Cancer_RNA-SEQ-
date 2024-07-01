# Set up a personal library path
personal_lib <- "~/R/libs"
if (!dir.exists(personal_lib)) dir.create(personal_lib, recursive = TRUE)
.libPaths(c(personal_lib, .libPaths()))

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = personal_lib)

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", lib = personal_lib)

if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr", lib = personal_lib)

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")



# Load necessary libraries

library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(readr)


# Define sample names and output directory
sample_names <- c("SRR5458657", "SRR5458658", "SRR5458661", "SRR5458662")
output_dir <- "kallisto_output_combined"

# Create the metadata dataframe
metadata <- data.frame(
  Sample = sample_names,
  Condition = c('Control', 'Control', 'Treatment', 'Treatment')
)
rownames(metadata) <- metadata$Sample

# Initialize an empty list to store dataframes
dfs <- list()

# Load data from each sample's abundance.tsv file
for (sample in sample_names) {
  file_path <- file.path(output_dir, sample, 'abundance.tsv')
  if (file.exists(file_path)) {
    df <- read_tsv(file_path, col_types = cols())
    df <- df[, c("target_id", "est_counts")]
    colnames(df) <- c("target_id", sample)
    dfs[[sample]] <- df
  } else {
    warning(paste(file_path, "does not exist"))
  }
}

# Combine all dataframes into a single dataframe
combined_df <- Reduce(function(x, y) merge(x, y, by = "target_id"), dfs)
rownames(combined_df) <- combined_df$target_id
combined_df <- combined_df[, -1]


# Round the estimated counts to the nearest integers
combined_df <- round(combined_df)


# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = combined_df, colData = metadata, design = ~ Condition)

# Perform the differential expression analysis
dds <- DESeq(dds)

# Get the results
results <- results(dds)

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

results_with_ids <- na.omit(results_with_ids)

# Save the updated data frame with Entrez Gene IDs as CSV
write.csv(results_with_ids, "differential_expression_results_with_ids.csv", row.names = FALSE)



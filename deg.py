import os
import pandas as pd
import scanpy as sc

# Define the sample names and the output directory
sample_names = ["SRR5458657", "SRR5458658", "SRR5458661", "SRR5458662"]
output_dir = "kallisto_output_combined"

# Initialize an empty list to store dataframes
dfs = []

# Load data from each sample's abundance.tsv file
for sample in sample_names:
    file_path = os.path.join(output_dir, sample, 'abundance.tsv')
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        df = df[['est_counts']]  # Use estimated counts
        df.columns = [sample]  # Rename column to the sample name
        dfs.append(df)
    else:
        print(f"Warning: {file_path} does not exist")

# Combine all dataframes into a single dataframe
combined_df = pd.concat(dfs, axis=1)

# Create the metadata
metadata = pd.DataFrame({
    'Sample': sample_names,
    'Condition': ['Control', 'Control', 'Treatment', 'Treatment']
})

# Set Sample column as index
metadata.set_index('Sample', inplace=True)

# Create an AnnData object
adata = sc.AnnData(X=combined_df.T, obs=metadata)

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform the data
sc.pp.log1p(adata)

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata, 'Condition', method='wilcoxon')

# Retrieve and format the results
de_results = adata.uns['rank_genes_groups']

# Extract results for the first group
group_name = de_results['names'].dtype.names[0]
result_df = pd.DataFrame({
    'names': de_results['names'][group_name],
    'logfoldchanges': de_results['logfoldchanges'][group_name],
    'pvals': de_results['pvals'][group_name],
})

# Show the top differentially expressed genes
print(result_df.head())

# Save the results
result_df.to_csv('differential_expression_results.csv', index=False)

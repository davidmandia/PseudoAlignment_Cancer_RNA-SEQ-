import pandas as pd
import numpy as np

# Load the CSV file
df = pd.read_csv('differential_expression_results_with_ids.csv')  # Adjust the separator if necessary
# Calculate -log10(padj)
df['-log10(padj)'] = -np.log10(df['padj'])

# Calculate the combined score: log2FoldChange * -log10(padj)
df['combined_score'] = df['log2FoldChange'] * df['-log10(padj)']

# Sort by the combined score in descending order
df_sorted = df.sort_values(by='combined_score', ascending=False)

# Save as .rnk file
df_sorted[['entrezgene_id', 'combined_score']].to_csv('output.rnk', sep='\t', header=False, index=False)
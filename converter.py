import pandas as pd
import numpy as np

# Load the CSV file
df = pd.read_csv('data.csv', sep='\t')  # Adjust the separator if necessary

# Calculate -log10(padj)
df['-log10(padj)'] = -np.log10(df['padj'])

# Sort by -log10(padj) descending
df_sorted = df.sort_values(by='-log10(padj)', ascending=False)

# Save as .rnk file
df_sorted[['entrezgene_id', '-log10(padj)']].to_csv('output.rnk', sep='\t', header=False, index=False)

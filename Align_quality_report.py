
import os
import pandas as pd
import matplotlib.pyplot as plt


OUTPUT_DIR = "kallisto_output_combined"


# Function to parse Kallisto output and gather metrics
def parse_kallisto_output(output_dir):
    abundance_file = os.path.join(output_dir, "abundance.tsv")
    df = pd.read_csv(abundance_file, sep='\t')
    mapping_rate = (df['est_counts'] > 0).sum() / len(df) * 100
    return {
        'Mapping Rate (%)': mapping_rate,
        'Mean Fragment Length': 200,  # Replace with actual values if available
        'Standard Deviation': 20     # Replace with actual values if available
    }

# Define sample names
sample_names = ["SRR5458657", "SRR5458658", "SRR5458661", "SRR5458662"]

# Create a DataFrame to store metrics
metrics = pd.DataFrame(columns=['Sample', 'Mapping Rate (%)', 'Mean Fragment Length', 'Standard Deviation'])


# Populate DataFrame with metrics
for sample_name in sample_names:
    output_dir = os.path.join(OUTPUT_DIR, sample_name)
    metrics.loc[len(metrics)] = [sample_name] + list(parse_kallisto_output(output_dir).values())

# Print metrics
print("\nAlignment Quality Metrics:")
print(metrics)

# Plotting example (adjust as needed)
plt.figure(figsize=(10, 6))
plt.bar(metrics['Sample'], metrics['Mapping Rate (%)'], color='blue')
plt.xlabel('Sample')
plt.ylabel('Mapping Rate (%)')
plt.title('Mapping Rate for Samples')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('mapping_rate.png')
plt.show()

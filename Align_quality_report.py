import os
import pandas as pd
import matplotlib.pyplot as plt
import json

OUTPUT_DIR = "kallisto_output_combined"

# Function to parse Kallisto output and gather metrics
def parse_kallisto_output(output_dir):
    run_info_file = os.path.join(output_dir, "run_info.json")
    with open(run_info_file, 'r') as f:
        run_info = json.load(f)
    mapping_rate = (run_info['n_pseudoaligned'] / run_info['n_processed']) * 100
    return {
        'Mapping Rate (%)': mapping_rate,
        'Mean Fragment Length': run_info.get('mean_fragment_length', 200),
        'Standard Deviation': run_info.get('sd_fragment_length', 20)
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

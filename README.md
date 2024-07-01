# RNA-Seq Analysis Pipeline

This repository contains scripts for performing RNA-Seq data analysis using Kallisto for pseudo-alignment, quality assessment of alignment, generating an RNK file, and DESeq2 for differential expression analysis.

## Contents

- `pseudo_align.sh`: A Bash script to download and install Kallisto, download the transcriptome reference, create a Kallisto index, download RNA-Seq data, and perform pseudo-alignment.
- `Download_input.sh`: A Bash script to download RNA-Seq FASTQ files from specified URLs.
- `Align_quality_report.py`: A Python script to assess the quality of the alignment performed by Kallisto and generate a mapping rate plot.
- `rnk_file.py`: A Python script to generate an RNK file from differential expression results.
- `deg.r`: An R script to perform differential expression analysis using DESeq2 and annotate results with Entrez Gene IDs.

## Prerequisites

- Unix-like operating system (Linux or macOS)
- `wget` for downloading files
- `tar` for extracting compressed files
- R with the following packages: `BiocManager`, `DESeq2`, `readr`, `biomaRt`, `clusterProfiler`, `org.Hs.eg.db`
- Python with the following packages: `pandas`, `matplotlib`

## Installation

### Kallisto

The `pseudo_align.sh` script downloads and installs Kallisto. No separate installation steps are needed.

### R Packages

The `deg.r` script installs the required R packages if they are not already installed. Ensure you have a working installation of R.

### Python Packages

Install the required Python packages using `pip`:

```bash
pip install pandas matplotlib
```

## Usage

#### Step 1: Pseudo-Alignment with Kallisto

This script performs the following steps:

- Downloads and installs Kallisto
- Downloads the Ensembl transcriptome reference and creates a Kallisto index
- Downloads RNA-Seq FASTQ files using Download_input.sh
- Runs Kallisto to perform pseudo-alignment on the downloaded FASTQ files
#### Step 2: Alignment Quality Assessment
Run the Align_quality_report.py script to assess the quality of the Kallisto alignment:

```bash
python Align_quality_report.py
```
This script performs the following steps:

- Parses the Kallisto output to gather alignment metrics
- Prints the alignment quality metrics
- Generates a plot of the mapping rate for each sample and saves it as mapping_rate.png



### Step 3: Differential Expression Analysis with DESeq2
Run the deg.r script in R:

** Run in Rstudio deg.r **

This script performs the following steps:
1. Loads necessary R packages
2. Reads abundance estimates from Kallisto output
3. Combines data into a single matrix and performs differential expression analysis using DESeq2
4. Annotates results with Entrez Gene IDs and saves the results to differential_expression_results_with_ids.csv

Alternatively, I have created a deg.py file that uses scanpy.
However, the script and results obtained with R fit better with the rest of the anlysis


### Step 4: Gene Set enrichment Analysis

###### Generating RNK File for GSEA
Run the rnk_file.py script to generate an RNK file from differential expression results:

```bash
python rnk_file.py
```
This script performs the following steps:

1. Loads differential expression results from differential_expression_results_with_ids.csv
2. Calculates -log10(padj) and combined_score
3. Sorts genes based on combined_score in descending order
4. Saves the gene Entrez IDs and combined scores as output.rnk

##### Perform GSEA on Webgestalt 2019 using output.rnk as input file

1. Obtain results for up-regulated and downregulated pathways
2. Combine results 

### Step 5: Report

Please see my report with the visualization of the results involving signalling 
File called "GSEA report.docx"

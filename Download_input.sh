#!/bin/bash

# Create a directory to store the input FASTQ files if it doesn't exist
mkdir -p input_fastq

# Define URLs for the FASTQ files
URLS=(
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR545/007/SRR5458657/SRR5458657.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR545/008/SRR5458658/SRR5458658.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR545/001/SRR5458661/SRR5458661.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR545/002/SRR5458662/SRR5458662.fastq.gz"
)

# Loop through each URL and download the corresponding data using wget if not already present
for URL in "${URLS[@]}"
do
    # Extract the filename from the URL
    FILENAME=$(basename "$URL")
    
    # Check if the file already exists in the input_fastq directory
    if [ ! -f "input_fastq/$FILENAME" ]; then
        echo "Downloading data from $URL..."
        wget -nc "$URL" -P input_fastq/
    else
        echo "File $FILENAME already exists in input_fastq directory. Skipping download."
    fi
done

echo "Download complete."




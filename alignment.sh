#!/bin/bash
sudo apt-get update
sudo apt-get install bowtie

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

# Download genome reference FASTA file
echo "Downloading genome reference FASTA file..."
wget -O genome_annotation.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

# Unzip genome annotation archive
echo "Extracting genome annotation archive..."
unzip genome_annotation.zip -d genome_annotation

# Locate the reference genome FASTA file
REFERENCE_GENOME="genome_annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

# Perform alignment using Bowtie
# Align each downloaded FASTQ file
for FILE in input_fastq/*.fastq.gz
do
    # Extract the filename without extension
    BASENAME=$(basename "$FILE" .fastq.gz)
    
    # Perform alignment
    echo "Aligning $FILE to $REFERENCE_GENOME..."
    bowtie -q "$REFERENCE_GENOME" "$FILE" > "${BASENAME}_alignments.sam"
done

echo "Alignment complete."


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

# Perform alignment using Bowtie
# Assuming your reference genome is named reference_genome.fasta and located in the current directory
REFERENCE_GENOME="reference_genome.fasta"

# Align each downloaded FASTQ file
for FILE in input_fastq/*.fastq.gz
do
    # Extract the filename without extension
    BASENAME=$(basename "$FILE" .fastq.gz)
    
    # Perform alignment
    echo "Aligning $FILE to $REFERENCE_GENOME..."
    bowtie -q "$REFERENCE_GENOME" "$FILE" > "${BASENAME}_alignments.sam"
done

echo "Alignment complete."




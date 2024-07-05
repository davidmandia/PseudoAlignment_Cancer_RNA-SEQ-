
# Download transciptome reference and GFF
wget -O genome_annotation.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
unzip genome_annotation.zip
##Remove file due to size 
rm genome_annotation.zip

# Check if STAR is already installed
if command -v STAR &> /dev/null; then
    echo "STAR is already installed. Skipping installation."
else
    # Install STAR
    echo "Downloading STAR..."
    wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
    tar -xzf 2.7.11b.tar.gz
    cd STAR-2.7.11b/source
    
    # Build STAR
    echo "Building STAR..."
    make STAR
    
    # Move STAR executable to /usr/local/bin/ (requires sudo)
    echo "Moving STAR executable to /usr/local/bin/..."
    sudo mv STAR /usr/local/bin/
fi

# Verify STAR installation
echo "Checking STAR version..."
STAR --version



# Prepare Genome Index with STAR
# Create directory for genome index if it doesn't exist
if [ ! -d "ncbi_dataset/data/human_genome_index" ]; then
    mkdir -p ncbi_dataset/data/human_genome_index
fi



STAR --runMode genomeGenerate \
     --genomeDir ncbi_dataset/data/human_genome_index/ \
     --genomeFastaFiles ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
     --sjdbGTFfile ncbi_dataset/data/GCF_000001405.40/genomic.gff \
     --runThreadN 8



# Align 

STAR --runThreadN 8 \
     --genomeDir /path/to/human_genome_index \
     --readFilesIn /path/to/trimmed_output/SRR5458657_trimmed.fastq.gz /path/to/trimmed_output/SRR5458658_trimmed.fastq.gz /path/to/trimmed_output/SRR5458661_trimmed.fastq.gz /path/to/trimmed_output/SRR5458662_trimmed.fastq.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix /path/to/alignment_output/prefix_



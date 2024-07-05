
#!/bin/bash
# Obtain input file by running download script
bash Download_input 

## Trim reads

### install package 
sudo apt-get install sra-toolkit

pip install cutadapt





# Set the version of Kallisto 
KALLISTO_VERSION="0.46.2"

# Set the download URL
KALLISTO_URL="https://github.com/pachterlab/kallisto/releases/download/v$KALLISTO_VERSION/kallisto_linux-v$KALLISTO_VERSION.tar.gz"

# Download the tar.gz archive
wget $KALLISTO_URL -O kallisto.tar.gz

# Extract the archive
tar -xzf kallisto.tar.gz

# Move the Kallisto binary to a PATH (so can be called using kallisto on the command line)
# Note: This might require sudo permissions
sudo mv kallisto/kallisto /usr/local/bin/

# Clean up the downloaded and extracted files
rm -rf kallisto.tar.gz kallisto

# Verify the installation
kallisto version

# Download ensembl Transcriptome reference

wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# Create Kallisto index of homo sapiens

kallisto index -i human_index.idx Homo_sapiens.GRCh38.cdna.all.fa


# Create output directory 
mkdir kallisto_output_combined



# Run alignment 
kallisto quant -i human_index.idx -o SRR5458657 --single -l 200 -s 20 SRR5458657.fastq.gz 
kallisto quant -i human_index.idx -o SRR5458658 --single -l 200 -s 20 SRR5458658.fastq.gz 
kallisto quant -i human_index.idx -o SRR5458661 --single -l 200 -s 20 SRR5458661.fastq.gz
kallisto quant -i human_index.idx -o SRR5458662 --single -l 200 -s 20 SRR5458662.fastq.gz 




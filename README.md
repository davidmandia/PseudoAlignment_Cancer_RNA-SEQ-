# Learning
#sudo apt update



## Trying different index 
pip install kb-python
conda create -n kallisto_env kallisto=0.46.1
conda init
conda activate kallisto_env

# sudo apt install kallisto
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz

kallisto index -i transcriptome.idx --gtf gene_annotation.gtf transcripts.fasta

kallisto index -i human_index.idx Homo_sapiens.GRCh38.cdna.all.fa

mkdir kallisto_output_combined
kallisto quant -i human_index.idx -o SRR5458661 --single -l 200 -s 20 SRR5458661.fastq.gz  # Per each file 

pip install pandas

pip install matplotlib

pip install scanpy pandas


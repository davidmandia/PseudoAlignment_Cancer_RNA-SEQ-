# Learning
#sudo apt update
# sudo apt install kallisto
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# kallisto index -i transcriptome_index.idx Homo_sapiens.GRCh38.cdna.all.fa
wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz
tar -xf human_index_nac.tar.xz

## Trying different index 
pip install kb-python
conda create -n kallisto_env kallisto=0.46.1
conda init
conda activate kallisto_env



wget https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
tar -xzf v0.46.1.tar.gz
cd kallisto-0.46.1
mkdir build
cd build
cmake ..
make
sudo make install
kallisto quant -i index.idx -o kallisto_output/SRR5458657 --single -l 200 -s 20 SRR5458657.fastq.gz




kallisto quant -i index.idx -o kallisto_output/SRR5458657 -b 100 -t 4 SRR5458657.fastq.gz
kallisto quant -i index.idx -o kallisto_output/SRR5458658 -b 100 -t 4 SRR5458658.fastq.gz
kallisto quant -i index.idx -o kallisto_output/SRR5458661 -b 100 -t 4 SRR5458661.fastq.gz
kallisto quant -i index.idx -o kallisto_output/SRR5458662 -b 100 -t 4 SRR5458662.fastq.gz
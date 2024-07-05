mkdir input_fastq_trimmed -p

## AGATCGGAAGAGC  is the Illumina adapter sequence

cutadapt \
  -a AGATCGGAAGAGC \
  -o input_fastq_trimmed/SRR5458657_trimmed.fastq.gz \
  input_fastq/SRR5458657.fastq.gz > input_fastq_trimmed/trimming_report_SRR5458657.txt

cutadapt \
  -a AGATCGGAAGAGC \
  -o input_fastq_trimmed/SRR5458658_trimmed.fastq.gz \
  input_fastq/SRR5458658.fastq.gz > input_fastq_trimmed/trimming_report_SRR5458658.txt

cutadapt \
  -a AGATCGGAAGAGC \
  -o input_fastq_trimmed/SRR5458661_trimmed.fastq.gz \
  input_fastq/SRR5458661.fastq.gz > input_fastq_trimmed/trimming_report_SRR5458661.txt

cutadapt \
  -a AGATCGGAAGAGC \
  -o input_fastq_trimmed/SRR5458662_trimmed.fastq.gz \
  input_fastq/SRR5458662.fastq.gz > input_fastq_trimmed/trimming_report_SRR5458662.txt


echo "If successful, please delete the input_fastq directory given its large size
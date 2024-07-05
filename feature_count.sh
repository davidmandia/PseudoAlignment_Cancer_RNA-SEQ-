# sudo apt update
# sudo apt install samtools
# conda install -c bioconda subread


# Assuming you have aligned and generated SAM files already
# Convert SAM to BAM and sort by coordinate
#!/bin/bash

# Convert SAM to BAM, Sort by Coordinate, and Index
#!/bin/bash

# Convert SAM to BAM, Sort by Coordinate, and Index
for SAM_FILE in *_alignments.sam
do
    BASENAME=$(basename "$SAM_FILE" .sam)
    if [ -f "$SAM_FILE" ]; then
        echo "Converting $SAM_FILE to BAM..."
        samtools view -bS -o "${BASENAME}.bam" "$SAM_FILE"
        
        echo "Sorting BAM file ${BASENAME}.bam..."
        samtools sort -o "${BASENAME}_sorted.bam" "${BASENAME}.bam"
        
        echo "Indexing sorted BAM file ${BASENAME}_sorted.bam..."
        samtools index "${BASENAME}_sorted.bam"
    else
        echo "ERROR: SAM file $SAM_FILE not found."
    fi
done

# Perform feature counting
GTF_FILE="genome_annotation/ncbi_dataset/data/GCF_000001405.40/genomic.gff"  
for SORTED_BAM_FILE in *_sorted.bam
do
    BASENAME=$(basename "$SORTED_BAM_FILE" _sorted.bam)
    if [ -f "$SORTED_BAM_FILE" ]; then
        echo "Counting features in $SORTED_BAM_FILE..."
        featureCounts -T 8 -a "$GTF_FILE" -o "${BASENAME}_counts.txt" "$SORTED_BAM_FILE"
    else
        echo "ERROR: Sorted BAM file $SORTED_BAM_FILE not found."
    fi
done

echo "Feature counting complete."

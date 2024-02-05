#!/bin/sh

## Download reads
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/050/SRR26705750/SRR26705750.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/046/SRR26705746/SRR26705746.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/044/SRR26705744/SRR26705744.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/040/SRR26705740/SRR26705740.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/048/SRR26705748/SRR26705748.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/045/SRR26705745/SRR26705745.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/049/SRR26705749/SRR26705749.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/041/SRR26705741/SRR26705741.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/042/SRR26705742/SRR26705742.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/039/SRR26705739/SRR26705739.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/043/SRR26705743/SRR26705743.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR267/047/SRR26705747/SRR26705747.fastq.gz

## Trimming
module load cutadapt

DIRECTORY="/ibex/scratch/tawfiqre/workflow_class"
QUALITY_CUTOFF=30

# Loop through all .fastq.gz files in the directory
for FILE in "$DIRECTORY"/*.fastq.gz; do
    OUTPUT_FILE="${FILE%.fastq.gz}_trimmed.fastq.gz"
    cutadapt -j 4 -q "$QUALITY_CUTOFF" -o "$OUTPUT_FILE" "$FILE"
done

## STAR alignment to CHM13
module load star
# Download reference genome (CHM13) and GTF annotation file from NCBI
wget https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF

# Index the genome
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir /ibex/scratch/tawfiqre/workflow_class/ncbi-genomes-2024-01-29 \
     --genomeFastaFiles /ibex/scratch/tawfiqre/workflow_class/ncbi-genomes-2024-01-29/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
     --sjdbGTFfile /ibex/scratch/tawfiqre/workflow_class/ncbi-genomes-2024-01-29/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf


# Run STAR on all reads
STAR_INDEX="/ibex/scratch/tawfiqre/workflow_class/ncbi-genomes-2024-01-29"
READS_DIR="/ibex/scratch/tawfiqre/workflow_class"
OUTPUT_DIR="/ibex/scratch/tawfiqre/workflow_class"
THREADS=16

# Loop through all trimmed fastq files in the directory
for FILE in "$READS_DIR"/*_trimmed.fastq.gz; do
    OUTPUT_PREFIX="$OUTPUT_DIR"/$(basename "$FILE" _trimmed.fastq.gz)
    STAR --genomeDir "$STAR_INDEX" \
         --readFilesIn "$FILE" \
         --runThreadN $THREADS \
         --outFileNamePrefix "$OUTPUT_PREFIX" \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand zcat
done

## Extract features
module load samtools
module load subread

BAM_DIR="/ibex/scratch/tawfiqre/workflow_class"
ANNOTATION_GTF="/ibex/scratch/tawfiqre/workflow_class/ncbi-genomes-2024-01-29/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
OUTPUT_DIR="/ibex/scratch/tawfiqre/workflow_class"
ALL_COUNTS_FILE="$OUTPUT_DIR/all_counts.tsv"

# Indexing all BAM files
for BAM_FILE in "$BAM_DIR"/*Aligned.sortedByCoord.out.bam; do
    samtools index "$BAM_FILE"
done

# Running featureCounts for all BAM files and combining the results
featureCounts -a "$ANNOTATION_GTF" -o "$ALL_COUNTS_FILE" "$BAM_DIR"/*Aligned.sortedByCoord.out.bam

# Post-processing to format the file for loading into R
# Skipping the header and merging all count files into a single file
tail -n +2 "$ALL_COUNTS_FILE" | cut -f 1,7- > "$OUTPUT_DIR/combined_counts.tsv"

echo "Count table is ready at $OUTPUT_DIR/combined_counts.tsv"

#!/bin/bash

# Define input files and reference genome
REF="/path/to/your/reference_genome.fa"
FASTQ1="/path/to/your/sample_name_1.fastq.gz"
FASTQ2="/path/to/your/sample_name_2.fastq.gz"
KNOWN_SITES="/path/to/your/known_variants.vcf.gz"
BED="/path/to/your/sample_name.bed"

# Convert FASTQ to BAM
bwa mem -t 3 -R "@RG\tID:1\tSM:sample_name\tPL:ILLUMINA\tLB:lib1\tPU:unit1" $REF $FASTQ1 $FASTQ2 > sample_name.sam
samtools view -bS sample_name.sam | samtools sort -o sample_name.bam
rm sample_name.sam

# Mark duplicates
gatk MarkDuplicates -I sample_name.bam -O sample_name.markdup.bam -M sample_name.metrics.txt

# Index BAM file
samtools index sample_name.markdup.bam

# Base Quality Score Recalibration
gatk BaseRecalibrator -R $REF -I sample_name.markdup.bam --known-sites $KNOWN_SITES -O sample_name.recal_data.table
gatk ApplyBQSR -R $REF -I sample_name.markdup.bam --bqsr-recal-file sample_name.recal_data.table -O sample_name_final.bam

# Validate BAM file
gatk ValidateSamFile -I sample_name_final.bam -MODE SUMMARY

# Convert BAM to VCF using Mutect2
gatk Mutect2 -R $REF -I sample_name_final.bam -O sample_name_unfiltered_variants.vcf -L $BED 
gatk FilterMutectCalls -R $REF -V sample_name_unfiltered_variants.vcf -O sample_name_filtered_variants.vcf -L $BED

# Completion message
echo "Pipeline completed successfully. Output VCF: sample_name_filtered_variants.vcf"

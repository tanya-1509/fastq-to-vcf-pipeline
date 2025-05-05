#!/bin/bash

# Define input files and reference genome
REF="/path/to/reference_genome.fa/"
FASTQ1="/home/tanya/PROJECTS/GATK_best_practices/TEST/SRR7890852_1.fastq.gz"
FASTQ2="/home/tanya/PROJECTS/GATK_best_practices/TEST/SRR7890852_2.fastq.gz"
KNOWN_SITES="/home/tanya/PROJECTS/GATK_best_practices/parabricks_sample/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz"
BED="/home/tanya/PROJECTS/GATK_best_practices/TEST/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38_Mod.bed"

# Convert FASTQ to BAM
bwa mem -t 3 -R "@RG\tID:1\tSM:SRR7890852\tPL:ILLUMINA\tLB:lib1\tPU:unit1" $REF $FASTQ1 $FASTQ2 > SRR7890852.sam
samtools view -bS SRR7890852.sam | samtools sort -o SRR7890852.bam
rm SRR7890852.sam

# Mark duplicates
gatk MarkDuplicates -I SRR7890852.bam -O SRR7890852.markdup.bam -M SRR7890852.metrics.txt

# Index BAM file
samtools index SRR7890852.markdup.bam

# Base Quality Score Recalibration
gatk BaseRecalibrator -R $REF -I SRR7890852.markdup.bam --known-sites $KNOWN_SITES -O SRR7890852.recal_data.table
gatk ApplyBQSR -R $REF -I SRR7890852.markdup.bam --bqsr-recal-file SRR7890852.recal_data.table -O SRR7890852_final.bam

# Validate BAM file
gatk ValidateSamFile -I SRR7890852_final.bam -MODE SUMMARY

# Convert BAM to VCF using Mutect2
gatk Mutect2 -R $REF -I SRR7890852_final.bam -O SRR7890852_unfiltered_variants.vcf -L $BED 
gatk FilterMutectCalls -R $REF -V SRR7890852_unfiltered_variants.vcf -O SRR7890852_filtered_variants.vcf -L $BED

# Completion message
echo "Pipeline completed successfully. Output VCF: SRR7890852_filtered_variants.vcf"

#!/usr/bin/env bash

# Exit on errors
set -e

# Reference genome
REF="./ref_seq/hg38.fa"

# Input files 
NORMAL_R1="PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz"
NORMAL_R2="PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz"
TUMOR_R1="PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz"
TUMOR_R2="PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz"

# Read Group Information
RG_ID_NORM="Norm"
RG_SM_NORM="PA221MH-Norm"
RG_PL="ILLUMINA"
RG_LB="lib09-Norm"
RG_PU="unit1"

RG_ID_TUMOR="Tumor"
RG_SM_TUMOR="PA220KH-Tumor"
RG_LB="lib09-Tumor"

# Align reads using BWA MEM with Read Group info included
echo "Aligning normal tissue reads..."
bwa mem -M -Y -t 2 -R "@RG\tID:$RG_ID_NORM\tSM:$RG_SM_NORM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU" "$REF" "$NORMAL_R1" "$NORMAL_R2" | \
    samtools view -b -S - > normal_aligned.bam

echo "Aligning tumor tissue reads..."
bwa mem -M -Y -t 2 -R "@RG\tID:$RG_ID_TUMOR\tSM:$RG_SM_TUMOR\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU" "$REF" "$TUMOR_R1" "$TUMOR_R2" | \
    samtools view -b -S - > tumor_aligned.bam

# Sort and index bam files
samtools sort normal_aligned.bam -o normal_sorted.bam
samtools sort tumor_aligned.bam -o tumor_sorted.bam

samtools index normal_sorted.bam
samtools index tumor_sorted.bam


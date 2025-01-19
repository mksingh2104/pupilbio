#!/usr/bin/env bash

# Exit on errors
set -e

# For normal sample
bcftools mpileup -f ../ref_seq/hg38.fa -Q 20 -B -Ou ../normal_sorted.bam | bcftools call -c -v -Oz -o normal.vcf.gz
bcftools index normal.vcf.gz

# For tumor sample 
bcftools mpileup -f ../ref_seq/hg38.fa -Q 20 -B -Ou ../tumor_sorted.bam | bcftools call -c -v -Oz -o tumor.vcf.gz
bcftools index tumor.vcf.gz


#!/usr/bin/env bash

# Exit on errors
set -e

# Reference genome
REF="ref_seq/hg38.fa"

# Sorted BAM files from alignment
NORMAL_BAM="normal_sorted.bam"
TUMOR_BAM="tumor_sorted.bam"

# Output VCFs
RAW_SOMATIC_VCF="somatic.vcf.gz"
FILTERED_SOMATIC_VCF="somatic_filtered.vcf.gz"

# Germline resource. Couldn't find for hg19 and therefore using hg38.
GERMLINE_RESOURCE="af-only-gnomad.hg38.vcf.gz"

echo "Calling somatic mutations with Mutect2"
./gatk-4.2.5.0/gatk Mutect2 \
    -R "$REF" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal "PA221MH-Norm" \
    --germline-resource "$GERMLINE_RESOURCE" \
    -O "$RAW_SOMATIC_VCF"

echo "Filtering Mutect2 calls"
./gatk-4.2.5.0/gatk FilterMutectCalls \
    -V "$RAW_SOMATIC_VCF" \
    -R "$REF" \
    -O "$FILTERED_SOMATIC_VCF"

echo "Somatic mutation calling complete. Filtered VCF: $FILTERED_SOMATIC_VCF"

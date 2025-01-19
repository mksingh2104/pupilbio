#!/usr/bin/env bash

# Exit on errors
set -e

bcftools filter -i 'QUAL>20 && DP>10' -Oz -o normal.filtered.vcf.gz normal.vcf.gz
bcftools filter -i 'QUAL>20 && DP>10' -Oz -o tumor.filtered.vcf.gz tumor.vcf.gz

bcftools index normal.filtered.vcf.gz
bcftools index tumor.filtered.vcf.gz

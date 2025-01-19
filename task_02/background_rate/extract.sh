#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' ../mpileup/normal.filtered.vcf.gz > normal_variants.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' ../mpileup/normal.vcf.gz > normal_variants.txt

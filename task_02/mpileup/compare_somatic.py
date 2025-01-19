#!/usr/bin/env python3

import sys
import cyvcf2

# Open the VCFs
normal_vcf = cyvcf2.VCF("normal.filtered.vcf.gz")
tumor_vcf = cyvcf2.VCF("tumor.filtered.vcf.gz")

# Storing normal variants in a dictionary for comparison with tumor: (chrom, pos, ref, alt)
normal_dict = {}

def get_alt_fraction(variant):
    #print( variant.format('AD') )
    ad = variant.format('AD')[0]  # can hold for multiple record. we want first
    if ad is None:
        return None
    ref_depth, alt_depths = ad[0], ad[1:]  # splits AD into ref and alt
    total_depth = sum(ad)
    if total_depth == 0:
        return None
    alt_depth = sum(alt_depths)  # sums up all alternate allele depths if multiple alts
    return alt_depth / total_depth if total_depth > 0 else None

# 1) Read normal VCF, store variant info
for variant in normal_vcf:
    #print(variant)
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alts = variant.ALT  # list of ALT alleles
    for alt in alts:
        key = (chrom, pos, ref, alt)
        alt_fraction = get_alt_fraction(variant)
        if alt_fraction is not None:
            normal_dict[key] = alt_fraction

# 2) Compare with tumor VCF
print("CHROM\tPOS\tREF\tALT\tNormal_AF\tTumor_AF\tSomatic?")
for variant in tumor_vcf:
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alts = variant.ALT  # list of ALT alleles
    for alt in alts:
        key = (chrom, pos, ref, alt)
        tumor_af = get_alt_fraction(variant)
        normal_af = normal_dict.get(key, 0.0)  # default 0 if not found
        if tumor_af is None or tumor_af == 0:
            continue
        is_somatic = (tumor_af >= 0.1 and normal_af <= 0.01)
        if is_somatic:
            print(f"{chrom}\t{pos}\t{ref}\t{alt}\t{normal_af:.3f}\t{tumor_af:.3f}\t{is_somatic}")


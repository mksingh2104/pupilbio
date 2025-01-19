#!/usr/bin/env python3

import statistics

input_file = "normal_variants.txt"
ALT_FRACTION_THRESHOLD = 0.30  # excluding anything >=30% alt fraction. arbitrary. can be changed

alt_fractions = []

with open(input_file, "r") as f:
    for line in f:
        # CHROM POS REF ALT DP AD
        cols = line.strip().split()
        #if len(cols) < 6:
        #    continue
        chrom, pos, ref, alt = cols[0], cols[1], cols[2], cols[3]
        dp = cols[4]
        ad = cols[5]

        # Convert dp to integer
        try:
            dp_val = int(dp)
        except ValueError:
            continue

        # AD could be comma separated: e.g. "30,2"
        ad_values = ad.split(',')
        if len(ad_values) < 2:
            continue

        try:
            ref_depth = int(ad_values[0])
            alt_depth = int(ad_values[1])
        except ValueError:
            continue

        # Calculate alt fraction
        if dp_val > 0:
            alt_fraction = alt_depth / dp_val
        else:
            continue

        # filter out high-frequency alleles (likely real variants)
        if alt_fraction < ALT_FRACTION_THRESHOLD:
            alt_fractions.append(alt_fraction)

print( alt_fractions )
background_rate = statistics.median(alt_fractions)
print(f"Median background mutation rate: {background_rate:.6e}")

# Simple approach
# RPM_needed = (1 / background_rate) * 1e6
rpm_needed = (1 / background_rate) * 1e6
print(f"Reads per million required for confident call: {rpm_needed:.2f}")

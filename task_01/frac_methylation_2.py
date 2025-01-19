import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

###############################
# 1) compute_pattern_enrichment_fisher
###############################
def compute_pattern_enrichment_fisher(
    input_csv,
    tissue1="cfDNA",
    tissue2="Islet",
    specificity_threshold=0.01, # minimizing false positive in tissue #1.
    alpha=0.05,
    output_csv=None
):

    # input data
    df = pd.read_csv(input_csv)

    # formatting
    old_patterns = ["`000","`001","`010","`011","`100","`101","`110","`111"]
    new_patterns = [p.replace("`", "") for p in old_patterns]
    rename_map    = dict(zip(old_patterns, new_patterns))
    df = df.rename(columns=rename_map)

    # Patterns
    pattern_cols = ["000","001","010","011","100","101","110","111"]

    # 1) Group by (CpG_Coordinates, Tissue) and sum pattern counts
    df_grouped = (
        df.groupby(["CpG_Coordinates", "Tissue"])[pattern_cols]
          .sum()
          .reset_index()
    )

    print( df_grouped.head(5) )

    # total reads for that locus + tissue
    df_grouped["total_reads"] = df_grouped[pattern_cols].sum(axis=1)

    # 2) Melt => each row has (CpG_Coordinates, Tissue, pattern, read_count)
    df_melt = df_grouped.melt(
        id_vars=["CpG_Coordinates","Tissue","total_reads"],
        value_vars=pattern_cols,
        var_name="pattern",
        value_name="read_count"
    )

    print( df_melt.head(5) )

    # 3) Pivot to get Tissue1 vs Tissue2 side-by-side
    df_pivot = df_melt.pivot_table(
        index=["CpG_Coordinates","pattern"],
        columns="Tissue",
        values=["read_count","total_reads"],
        fill_value=0
    ).reset_index()

    print( df_pivot.head(5) )

    # Flatten multi-level columns
    df_pivot.columns = [
        "_".join([str(c) for c in col if c]) for col in df_pivot.columns.values
    ]
    df_pivot.rename(columns={
        "CpG_Coordinates_":"CpG_Coordinates",
        "pattern_":"pattern"
    }, inplace=True)

    print( df_pivot.head(5) )

    # 4) Compute fraction for Tissue1, Tissue2
    df_pivot["fraction_Tissue1"] = (
        df_pivot[f"read_count_{tissue1}"] / df_pivot[f"total_reads_{tissue1}"].replace(0, np.nan)
    ).fillna(0)
    df_pivot["fraction_Tissue2"] = (
        df_pivot[f"read_count_{tissue2}"] / df_pivot[f"total_reads_{tissue2}"].replace(0, np.nan)
    ).fillna(0)

    print( df_pivot.head(5) )

    # 5) Fisher's exact test => Tissue2 fraction > Tissue1 fraction
    pvals = []
    for idx, row in df_pivot.iterrows():
        c2 = row.get(f"read_count_{tissue2}", 0)
        t2 = row.get(f"total_reads_{tissue2}", 0)
        c1 = row.get(f"read_count_{tissue1}", 0)
        t1 = row.get(f"total_reads_{tissue1}", 0)

        if (t2 == 0 and t1 == 0):
            pvals.append(1.0)
            continue

        table = [[c2, t2 - c2],
                 [c1, t1 - c1]]
        _, pval = fisher_exact(table, alternative="greater")
        pvals.append(pval)

    df_pivot["p_value"] = pvals

    # 6) Filter by Tissue1 fraction <= specificity_threshold
    #    => ensures minimal false positives in Tissue1
    cond_specific = (df_pivot["fraction_Tissue1"] <= specificity_threshold)
    df_spec = df_pivot[cond_specific].copy()

    # 7) Filter p-value < alpha
    df_signif = df_spec[df_spec["p_value"] < alpha].copy()

    # 8) Save & return
    if output_csv:
        df_signif.to_csv(output_csv, index=False)

    return df_signif


##############################
# 4) Putting It All Together
##############################
if __name__ == "__main__":

    final_patterns = compute_pattern_enrichment_fisher(
        input_csv="PupilBioTest_PMP_revA.csv",
        tissue1="cfDNA",
        tissue2="Islet",
        specificity_threshold=0.01,
        alpha=0.05,
        output_csv="pmp_enrichment.csv"
    )

    print("\nSignificant PMPs with high specificity:\n")
    print(final_patterns.head(20))
    
    # we pick "top 10" by smallest p-value:
    final_patterns.sort_values("p_value", inplace=True)
    top10_pmps = final_patterns.head(10).copy()
    print("\nTop 10 PMPs:\n", top10_pmps[["CpG_Coordinates","pattern","fraction_Tissue1","fraction_Tissue2","p_value"]])


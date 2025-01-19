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
    """
    """

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
# 2) Part (b) - Coverage threshold using Fisher's test
##############################
from scipy.stats import fisher_exact

def estimate_coverage_fisher(
    f1, f2,
    n1=1000000,           # coverage for Tissue #1
    alpha=0.05,
    power=0.8,
    max_coverage=1000000,
    step=100
):
    """
    """
    np.random.seed(42)
    num_sims = 1000

    # If there's no difference or fraction is reversed, we won't find coverage
    if f2 <= f1:
        return None

    for n2 in range(10, max_coverage+1, step):
        count_significant = 0
        for _ in range(num_sims):
            k1 = np.random.binomial(n1, f1)
            k2 = np.random.binomial(n2, f2)

            table = [[k2, n2 - k2],
                     [k1, n1 - k1]]
            _, pval = fisher_exact(table, alternative="greater")
            if pval < alpha:
                count_significant += 1

        current_power = count_significant / num_sims
        if current_power >= power:
            return n2

    return None

def estimate_thresholds_for_top10_pmps(
    top10_df,  # must contain columns fraction_Tissue1, fraction_Tissue2
    coverage_tissue1=1000000, # or any known coverage for Tissue #1
    alpha=0.05,
    power=0.8
):
    """
    """
    results = []
    for idx, row in top10_df.iterrows():
        cpg_coord = row["CpG_Coordinates"]
        patt = row["pattern"]
        f1   = row["fraction_Tissue1"]
        f2   = row["fraction_Tissue2"]

        n_required = estimate_coverage_fisher(
            f1=f1, f2=f2,
            n1=coverage_tissue1,
            alpha=alpha, power=power,
            max_coverage=100000,
            step=100
        )
        results.append({
            "CpG_Coordinates": cpg_coord,
            "pattern": patt,
            "fraction_T1": f1,
            "fraction_T2": f2,
            "coverage_Tissue1": coverage_tissue1,
            "coverage_required_T2": n_required
        })
    return pd.DataFrame(results)


##############################
# 4) Putting It All Together
##############################
if __name__ == "__main__":


    # 1) Identify significant PMPs using the Fisher-based enrichment analysis
    final_patterns = compute_pattern_enrichment_fisher(
        input_csv="PupilBioTest_PMP_revA.csv",
        tissue1="cfDNA",
        tissue2="Islet",
        specificity_threshold=0.01,
        alpha=0.05,
        output_csv=None  # or specify a filename to save
    )

    print("\nSignificant PMPs:\n", final_patterns.head())

    # 2) Sort the PMPs by smallest p-value and select top 10
    final_patterns.sort_values("p_value", inplace=True)
    top10_pmps = final_patterns.head(10).copy()

    ###############
    # (b) Coverage Threshold Estimation (Fisher-based)
    ###############
    print("\nEstimating coverage required for Tissue #2 (Fisher-based approach)...")
    coverage_df = estimate_thresholds_for_top10_pmps(
        top10_pmps,
        coverage_tissue1=1000000,  # 1M
        alpha=0.05,
        power=0.8
    )
    print("\nCoverage needed for top10 PMPs:\n")
    print(coverage_df)



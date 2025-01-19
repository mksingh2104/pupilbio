import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# -----------------------------
# PART A: Per-row coverage
# -----------------------------
pattern_map = {
    '`000': [1, 1, 1],
    '`001': [0, 0, 1],
    '`010': [0, 1, 0],
    '`011': [0, 1, 1],
    '`100': [1, 0, 0],
    '`101': [1, 0, 1],
    '`110': [1, 1, 0],
    '`111': [1, 1, 1]
}

df = pd.read_csv('PupilBioTest_PMP_revA.csv')

# Split the "CpG_Coordinates" column into three separate CpG positions
df[['CpG1', 'CpG2', 'CpG3']] = df['CpG_Coordinates'].str.split(':', expand=True)

# Initialize coverage columns
df['CpG1_Coverage'] = 0
df['CpG2_Coverage'] = 0
df['CpG3_Coverage'] = 0

# Sum coverage based on methylation patterns
for pat, mask in pattern_map.items():
    if pat not in df.columns:
        continue
    # mask = [x1, x2, x3], where x1..x3 is 0 or 1
    df['CpG1_Coverage'] += df[pat] * mask[0]
    df['CpG2_Coverage'] += df[pat] * mask[1]
    df['CpG3_Coverage'] += df[pat] * mask[2]

print(df.head(5)) #[['CpG_Coordinates','CpG1','CpG2','CpG3','CpG1_Coverage','CpG2_Coverage','CpG3_Coverage']])

df_list = []
# Group by Tissue, Sample_ID, and Replicate
for tissue, dft in df.groupby(['Tissue']):

    print(tissue)
    # 1) Melt the DataFrame so that CpG1, CpG2, CpG3 coverage appear in separate rows
    df_cp1 = dft[['strand']].copy()
    df_cp1['CpG'] = dft['CpG1']
    df_cp1['Coverage'] = dft['CpG1_Coverage']
    df_cp1['Sample_ID'] = dft['Sample_ID']
    df_cp1['Replicate'] = dft['Replicate']
    #df_cp1['Tissue'] = dft['Tissue']
    #Sample_ID Replicate
    print( df_cp1.head(5) )

    df_cp2 = dft[['strand']].copy()
    df_cp2['CpG'] = dft['CpG2']
    df_cp2['Coverage'] = dft['CpG2_Coverage']
    df_cp2['Sample_ID'] = dft['Sample_ID']
    df_cp2['Replicate'] = dft['Replicate']
    #df_cp2['Tissue'] = dft['Tissue']
    print( df_cp2.head(5) )

    df_cp3 = dft[['strand']].copy()
    df_cp3['CpG'] = dft['CpG3']
    df_cp3['Coverage'] = dft['CpG3_Coverage']
    df_cp3['Sample_ID'] = dft['Sample_ID']
    df_cp3['Replicate'] = dft['Replicate']
    #df_cp3['Tissue'] = dft['Tissue']
    print( df_cp3.head(5) )

    df_melted = pd.concat([df_cp1, df_cp2, df_cp3], ignore_index=True)
    print( df_melted.head(5) )

    # df_temp = (
    #    df_melted
    #    .groupby(['CpG','strand', 'Sample_ID', 'Replicate'], as_index=False)['Coverage']
    #    .sum()
    #)

    df_max = ( df_melted.groupby(['CpG', 'strand', 'Sample_ID', 'Replicate'], as_index=False)['Coverage'].max() )
    print( df_max.head(5) )

    df_sum = ( df_max.groupby(['CpG'], as_index=False)['Coverage'].sum())
    print( df_sum.head(5) )

    df_sum['Tissue'] = tissue[0]
    df_list.append(df_sum)
#
df_final = pd.concat(df_list, ignore_index=True)
print(df_final.head(15))

# -----------------------------
# NEW: (1a) Median and CV Calculation
# -----------------------------
# For each tissue, calculate:
#   - The median single-CpG coverage
#   - Coefficient of Variation (CV) = std / mean

stats_df = (
    df_final
    .groupby("Tissue")["Coverage"]
    .agg(["median", "mean", "std"])
    .reset_index()
)
stats_df["cv"] = stats_df["std"] / stats_df["mean"]

print("\n--- Coverage Stats by Tissue ---")
print(stats_df)


# A) BOX PLOT - Summed coverage distribution by tissue
plt.figure(figsize=(8, 6))
sns.boxplot(data=df_final, x='Tissue', y='Coverage')
plt.title("Box Plot of Summed Coverage by Tissue")
plt.xlabel("Tissue")
plt.ylabel("Summed Coverage (f + r)")
plt.tight_layout()
plt.show()

# Define total number of CpG sites for normalization
max_coverage = int(df_final['Coverage'].max())

# B) Percentage of CpGs with coverage >= X, by Tissue
plt.figure(figsize=(12, 8))
plt.title('Percentage of CpG Sites with Coverage ≥ X by Tissue')
plt.xlabel('Coverage (X)')
plt.ylabel('Percentage of CpG Sites (%)')
plt.grid(True)

tissues = df_final['Tissue'].unique()

for tissue in tissues:
    tissue_data = df_final[df_final['Tissue'] == tissue]
    max_coverage_tissue = tissue_data['Coverage'].max()
    total_cpg_sites = tissue_data['CpG'].nunique()

    coverage_range = np.arange(0, max_coverage_tissue + 1, 100)

    percentages = [
        100 * (tissue_data['Coverage'] >= cov).sum() / total_cpg_sites
        for cov in coverage_range
    ]

    plt.plot(
        coverage_range, percentages,
        label=f'Tissue {tissue}',
        marker='o', linestyle='-', markersize=4
    )

plt.legend(title='Tissue Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

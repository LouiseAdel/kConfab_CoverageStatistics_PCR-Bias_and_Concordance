# script2_run_jaccard.py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Input: converted long-format TSV files
input_folder = "Jaccardinput"
output_folder = "jaccard_heatmaps_triangle"
binary_output_folder = "binary_matrices"

os.makedirs(output_folder, exist_ok=True)
os.makedirs(binary_output_folder, exist_ok=True)

tsv_files = glob.glob(os.path.join(input_folder, "*_long.tsv"))
print(f"🔍 Found {len(tsv_files)} files in {input_folder}")

required_cols = {'Sample', 'Isoform', 'Gene', 'protocol', 'support_percent', 'sample_type'}

def compute_jaccard(df, file_name):
    df["support_percent"] = pd.to_numeric(df["support_percent"], errors="coerce")

    # Log excluded rows
    bad_rows = df[df["support_percent"].isna()]
    log_path = os.path.join(binary_output_folder, os.path.basename(file_name).replace(".tsv", "_excluded_rows.txt"))
    if not bad_rows.empty:
        bad_rows.to_csv(log_path, sep="\t", index=False)
        print(f"🧾 Logged excluded rows: {log_path}")

    df = df.dropna(subset=["support_percent"])

    df['isoform_id'] = df['Gene'].astype(str) + "|" + df['Isoform'].astype(str) + "|" + df['sample_type'].astype(str)
    df['presence'] = df['support_percent'].apply(lambda x: 1 if x > 0 else 0)

    binary_matrix = df.pivot_table(index='isoform_id',
                                   columns='protocol',
                                   values='presence',
                                   aggfunc='max',
                                   fill_value=0)

    binary_file_path = os.path.join(binary_output_folder, os.path.basename(file_name).replace(".tsv", "_binary_matrix.txt"))
    binary_matrix.to_csv(binary_file_path, sep="\t")
    print(f"📝 Saved binary matrix: {binary_file_path}")

    def jaccard_similarity(col1, col2):
        intersection = (col1 & col2).sum()
        union = (col1 | col2).sum()
        return intersection / union if union != 0 else 0

    protocols = binary_matrix.columns
    jaccard_df = pd.DataFrame(index=protocols, columns=protocols)

    for p1 in protocols:
        for p2 in protocols:
            jaccard_df.loc[p1, p2] = jaccard_similarity(binary_matrix[p1], binary_matrix[p2])

    jaccard_df = jaccard_df.astype(float)
    mask = np.triu(np.ones_like(jaccard_df, dtype=bool), k=1)

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        jaccard_df,
        annot=True,
        mask=mask,
        cmap="RdBu_r",
        vmin=0,
        vmax=1,
        fmt=".2f",
        square=True,
        linewidths=0.5,
        linecolor='white',
        cbar_kws={"label": "Jaccard Similarity"}
    )
    plt.title(f"Jaccard Similarity: {os.path.basename(file_name)}")
    plt.tight_layout()
    output_path = os.path.join(output_folder, os.path.basename(file_name).replace(".tsv", "_jaccard.png"))
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ Saved Jaccard heatmap: {output_path}")

# Process all files
for file in tsv_files:
    try:
        df = pd.read_csv(file, sep="\t")
        if required_cols.issubset(df.columns):
            compute_jaccard(df, file)
        else:
            missing = required_cols - set(df.columns)
            print(f"⚠️ Skipped {file} – missing required columns: {missing}")
    except Exception as e:
        print(f"❌ Error processing {file}: {e}")

print("🎯 Jaccard analysis complete.")

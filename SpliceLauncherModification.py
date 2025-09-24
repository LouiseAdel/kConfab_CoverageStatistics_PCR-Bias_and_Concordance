#python SpliceLauncherModification.py
import pandas as pd

# === 1. Load SpliceLauncher file ===
# Replace with your actual filename
input_file = "splicelauncher_output.xlsx"
df = pd.read_excel(input_file)


# Identify count columns (exclude P_* probability columns)
count_cols = [c for c in df.columns if c.endswith(".count") and not c.startswith("P_")]

# Keep only useful columns
df = df[["Gene", "AnnotJuncs"] + count_cols]

# === 2. Define event map (gene-aware) ===
# Each event: which gene, which inclusion junctions, which skip/alt junctions
event_map = {
    # --- BRCA2 ---
    "BRCA2_∆19":        { 'gene': "BRCA2", 'incl': ["18_19","19_20"], 'skip': ["∆19"] },
    "BRCA2_∆20":        { 'gene': "BRCA2", 'incl': ["19_20","20_21"], 'skip': ["∆20"] },
    "BRCA2_∆4_7":       { 'gene': "BRCA2", 'incl': ["3_4","7_8"], 'skip': ["∆4_7"] },
    "BRCA2_∆5":         { 'gene': "BRCA2", 'incl': ["4_5","5_6"], 'skip': ["∆5"] },
    "BRCA2_∆5_7":       { 'gene': "BRCA2", 'incl': ["4_5","7_8"], 'skip': ["∆5_7"] },
    "BRCA2_∆25":        { 'gene': "BRCA2", 'incl': ["24_25","25_26"], 'skip': ["∆25"] },
    "BRCA2_∆17_18":     { 'gene': "BRCA2", 'incl': ["16_17","18_19"], 'skip': ["∆17_18"] },
    "BRCA2_∆18":        { 'gene': "BRCA2", 'incl': ["17_18","18_19"], 'skip': ["∆18"] },

    # --- BRCA1 ---
    "BRCA1_∆4q(22)":    { 'gene': "BRCA1", 'incl': ["3_4","4_5"], 'skip': ["∆4q(22)"] },
    "BRCA1_∆4":         { 'gene': "BRCA1", 'incl': ["3_4","4_5"], 'skip': ["∆4"] },
    "BRCA1_∆3":         { 'gene': "BRCA1", 'incl': ["2_3","3_4"], 'skip': ["∆3"] },

    "BRCA1_∆8_11":      { 'gene': "BRCA1", 'incl': ["7_8","11_12"], 'skip': ["∆8_11"] },
    "BRCA1_∆10q(3309)": { 'gene': "BRCA1", 'incl': ["9_10","10_11"], 'skip': ["∆10q(3309)"] },
    "BRCA1_∆10":        { 'gene': "BRCA1", 'incl': ["9_10","10_11"], 'skip': ["∆10"] },

    "BRCA1_∆9_10":      { 'gene': "BRCA1", 'incl': ["8_9","10_11"], 'skip': ["∆9_10"] },
    "BRCA1_∆8_10":      { 'gene': "BRCA1", 'incl': ["7_8","10_11"], 'skip': ["∆8_10"] },
    "BRCA1_∆9":         { 'gene': "BRCA1", 'incl': ["8_9","9_10"], 'skip': ["∆9"] },
    "BRCA1_∆8_9":       { 'gene': "BRCA1", 'incl': ["7_8","9_10"], 'skip': ["∆8_9"] },
    "BRCA1_∆8":         { 'gene': "BRCA1", 'incl': ["7_8","8_9"], 'skip': ["∆8"] },

    "BRCA1_▼8p(21)":    { 'gene': "BRCA1", 'incl': ["7_8","8_9"], 'skip': ["▼8p(21)"] },

    "BRCA1_∆22":        { 'gene': "BRCA1", 'incl': ["21_22","22_23"], 'skip': ["∆22"] },
    "BRCA1_∆21_22":     { 'gene': "BRCA1", 'incl': ["20_21","22_23"], 'skip': ["∆21_22"] },
    "BRCA1_∆21":        { 'gene': "BRCA1", 'incl': ["20_21","21_22"], 'skip': ["∆21"] },
    "BRCA1_∆20":        { 'gene': "BRCA1", 'incl': ["19_20","20_21"], 'skip': ["∆20"] },
}

# === 3. Helper to get inclusion/skip counts ===
def get_incl_skip_counts(df, gene, incl_labels, skip_labels, count_cols):
    # Inclusion: if two junctions, take their average; otherwise sum
    incl_mask = (df["Gene"] == gene) & (df["AnnotJuncs"].isin(incl_labels))
    incl_df = df.loc[incl_mask, count_cols]
    if len(incl_labels) == 2:
        incl_counts = incl_df.sum() / 2
    else:
        incl_counts = incl_df.sum()

    # Skipping: always sum
    skip_mask = (df["Gene"] == gene) & (df["AnnotJuncs"].isin(skip_labels))
    skip_counts = df.loc[skip_mask, count_cols].sum()

    return incl_counts, skip_counts

# === 4. Calculate PSI and Percent Skipping ===
results = []
for event, meta in event_map.items():
    incl_counts, skip_counts = get_incl_skip_counts(df, meta['gene'], meta['incl'], meta['skip'], count_cols)
    denom = incl_counts + skip_counts

    # PSI and Percent Skipping as percentages
    psi = (incl_counts / denom.replace(0, pd.NA)) * 100
    psi = psi.round(2)
    percent_skipping = (skip_counts / denom.replace(0, pd.NA)) * 100
    percent_skipping = percent_skipping.round(2)

    tmp = pd.DataFrame({
        "Event": event,
        "Gene": meta['gene'],
        "Sample": psi.index,
        "Inclusion_reads": incl_counts.values,
        "Skipping_reads": skip_counts.values,
        "Denominator": denom.values,
        "PSI": psi.values,
        "Percent_Skipping": percent_skipping.values
    })

    tmp["Denominator_Junctions"] = ", ".join(meta['incl'] + meta['skip'])
    results.append(tmp)



psi_df = pd.concat(results, ignore_index=True)

# === 5. Save to file ===
#psi_df.to_csv("BRCA1_BRCA2_events_PSI.csv", index=False)
#print("✅ PSI values written to BRCA1_BRCA2_events_PSI.csv")
psi_df = pd.concat(results, ignore_index=True)
psi_df.to_excel("BRCA1_BRCA2_events_PSI.xlsx", index=False)
print("✅ PSI values written to BRCA1_BRCA2_events_PSI.xlsx")

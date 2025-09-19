# SpliceLauncher Modification
Python script for recalculating alternative splicing event percentages from SpliceLauncher junction count tables. The workflow imports junction counts, defines inclusion as the average of the two canonical flanking junctions, and defines event evidence as the direct alternative junction(s). Event percentages are then computed as:
Event %=Event reads / (Average inclusion reads+Event reads) * 100% 
This ensures values are bounded between 0–100% and directly comparable across single‑ and multi‑exon events.

Filename: SpliceLauncherModification.py

# PCR-Bias
Script for conducting log2 Fold-change and Fishers Exact on levels of the isoform Δ10q from long-read RNA-sequencing protocols compared to direct RNA long-read sequensing for comparisons on kConfab Cell lines

Within the script, printed results are available

Filename: PCRBias_Statistics.R

# Concordance
Script for conducting Jaccard Similarity calculations and plotting for long-read RNA-sequencing protocol comparisons on kConfab Cell lines.

Contains:

Python script - Run_jaccard.py

Example input data - ExampleInput_Long.tsv

Outdated files: OutdatedExampleInput_long.tsv and outdated_Jaccard.py

# Coverage Statistics
Script for conducting Kruskal-Wallis and Dunn's Test on Coverage Data from long-read RNA-sequencing protocol comparisons on kConfab Cell lines

Results are also presented within the file

File name: CoverageStatistics.R


# Overall Notes
Script development was aided by Copilot and ChatGPT. 

#Commands for Calculation of Fisher's Exact Test and Log2 Fold-Change - Isoform Rediscovery, PCR Bias

#Comparison of Δ10q levels in controls (Sample 15 and 19) 
#from Protocols (P) to levels detected by Direct RNA-seq (ONT-Direct)
#For Short-read analysis, the Full-length reads are the average of junctions supporting the exon 9-10 and 10-11 connection without Δ10q.

#Libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

#PCR Bias Statistics ----
#Individual Protocols ----
#Read data file
data_controls <- read_excel("PCRBias_delta10q_Levels_controlsonly.xlsx")
print(data_controls)

#SampleID    Protocol Delta10q_Percentage Delta10q_Reads FullLength_Reads
#Control_S15        1                4.43           2598          58646  
#Control_S19        1               22.9            3328          14558  
#Control_S15        2               17.0              10             59  
#Control_S19        2                8.51              4             47  
#Control_S15        3               21.3              36            169  
#Control_S19        3               16.7              34            204  
#Control_S15        4               44.5           17850          40139  
#Control_S19        4               43.8           11671          26670  
#Control_S15        5               44.0            1464           3331  
#Control_S19        5               46.7            2259           4834  
#Control_S15        6               13.5             577           4274  
#Control_S19        6               12.2            2344          19213  
#Control_S15        7                5                 8            160  
#Control_S19        7               12.2              14            115  
#Control_S15        8                8.22              3             36.5
#Control_S19        8                7.41              5             67.5

#Protocol 1 is PacBio-WT, 2 ONT-WT-A, 3 ONT-WT-B, 4 ONT-Amplicon-1, 5 ONT-Amplicon-2, 6 ONT-Capture, 7 ONT-Direct, and 8 Ill-SR


#Combine control counts
combined_counts <- data_controls %>%
  group_by(Protocol) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    FullLength_Reads = sum(FullLength_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + FullLength_Reads)
  )

# Separate out Protocol 7 (Direct RNA-seq)
direct_rna <- combined_counts %>%
  filter(Protocol == 7) %>%
  select(Delta10q_Reads, FullLength_Reads, Delta10q_Percentage)

#Fisher's Exact Test

fisher_results <- combined_counts %>%
  filter(Protocol != 7) %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, FullLength_Reads,
        direct_rna$Delta10q_Reads, direct_rna$FullLength_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()


#Log2 Fold-Change
logfc_data <- combined_counts %>%
  filter(Protocol != 7) %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (direct_rna$Delta10q_Percentage + 1e-6))
  )


#Comnbined Log2FC and Fisher
results_individual <- logfc_data %>%
  left_join(fisher_results %>% select(Protocol, Fisher_p), by = "Protocol")

# Optional: Print results
print(results_individual)

#Result - Individual Protocols:
#Protocol Delta10q_Reads FullLength_Reads Delta10q_Percentage  Log2FC Fisher_p
#Ill-SR                      8              104                7.14 -0.0525 1e+0          NS
#ONT-Amplicon-1          29521            66809               30.6   2.05   1.89e-22      S
#ONT-Amplicon-2           3723             8165               31.3   2.08   6.61e-23      S
#ONT-Capture              2921            23487               11.1   0.578  4.97e-2       S
#ONT-WT-A                   14              106               11.7   0.655  1.79e-1       NS
#ONT-WT-B                   70              373               15.8   1.09   6.13e-4       S
#PacBio-WT                5926            73204               7.49   0.0158 1e+0          NS

#Grouped by amount of PCR:  ----

data_controls_group <- data_controls %>%
  mutate(
    Protocol_Group = case_when(
      Protocol %in% c(1, 2, 3, 6) ~ "PCR_cDNA_only",
      Protocol %in% c(4, 5) ~ "PCR_cDNA_plus_enrichment",
      Protocol == 7 ~ "Direct_RNA",
      TRUE ~ "Other"
    )
  )

# Combine control counts per protocol group
combined_groups <- data_controls_group %>%
  filter(Protocol_Group != "Other") %>%
  group_by(Protocol_Group) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    FullLength_Reads = sum(FullLength_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + FullLength_Reads)
  )

# Extract Direct RNA group counts
direct_rna_group <- combined_groups %>%
  filter(Protocol_Group == "Direct_RNA") %>%
  select(Delta10q_Reads, FullLength_Reads, Delta10q_Percentage)

# Log2 Fold-Change for grouped protocols vs Direct RNA
logfc_data_group <- combined_groups %>%
  filter(Protocol_Group != "Direct_RNA") %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (direct_rna_group$Delta10q_Percentage + 1e-6))
  )

# Fisher's Exact Test for grouped protocols vs Direct RNA
fisher_results_group <- combined_groups %>%
  filter(Protocol_Group != "Direct_RNA") %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, FullLength_Reads,
        direct_rna_group$Delta10q_Reads, direct_rna_group$FullLength_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()

# Combine Log2FC and Fisher p-values for grouped data
results_grouped <- logfc_data_group %>%
  left_join(fisher_results_group %>% select(Protocol_Group, Fisher_p), by = "Protocol_Group")

# Print grouped protocol results
print(results_grouped)

#Result - Grouped by amount of PCR:
#Protocol_Group           Delta10q_Reads FullLength_Reads Delta10q_Percentage Log2FC Fisher_p
#PCR_cDNA_only                      8931            97170                8.42  0.184 6.01e- 1
#PCR_cDNA_plus_enrichment          33244            74974               30.7   2.05  1.79e-22






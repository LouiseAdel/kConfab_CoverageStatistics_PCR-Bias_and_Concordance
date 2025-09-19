#Commands for Calculation of Fisher's Exact Test and Log2 Fold-Change - Isoform Rediscovery, PCR Bias

#Comparison of Δ10q levels in controls (Sample 15 and 19) from Protocols (P) to levels detected by Direct RNA-seq (ONT-Direct)
#For Short-read analysis, the Full-length reads are the average of junctions supporting the exon 9-10 and 10-11 connection without Δ10q.
#Additional comparison of Δ10q levels in all controls (S10, S12, S13, S14, S15, S16, S17, S18, S19, and S20) against short-read levels (ONT-direct excluded due to lack of sequencing of 8 control lines)

#Libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

#PCR Bias Statistics Comparison to Direct RNA-seq----
#Individual Protocols ----
#Read data file
data_controls <- read_excel("PCRBias_delta10q_Levels_controlsonly.xlsx")
print(data_controls)

#SampleID    Protocol Delta10q_Percentage Delta10q_Reads NonDelta10q_Reads
#Control_S15        1                4.43           2598           56048  
#Control_S19        1               22.9            3328           11230  
#Control_S15        2               17.0              10              49  
#Control_S19        2                8.51              4             112  
#Control_S15        3               21.3              36             133  
#Control_S19        3               16.7              34             170  
#Control_S15        4               44.5           17850           22289  
#Control_S19        4               43.8           11671           14999  
#Control_S15        5               44.0            1464            1867  
#Control_S19        5               46.7            2259            2575  
#Control_S15        6               13.5             577            3697  
#Control_S19        6               12.2            2344           16869  
#Control_S15        7                5                 8             152  
#Control_S19        7               12.2              14             101  
#Control_S15        8                8.22              3              33.5
#Control_S19        8                7.41              5              62.5

#Protocol 1 is PacBio-WT, 2 ONT-WT-A, 3 ONT-WT-B, 4 ONT-Amplicon-1, 5 ONT-Amplicon-2, 6 ONT-Capture, 7 ONT-Direct, and 8 Ill-SR


#Combine control counts
combined_counts <- data_controls %>%
  group_by(Protocol) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    NonDelta10q_Reads = sum(NonDelta10q_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + NonDelta10q_Reads)
  )

# Separate out Protocol 7 (Direct RNA-seq)
direct_rna <- combined_counts %>%
  filter(Protocol == 7) %>%
  select(Delta10q_Reads, NonDelta10q_Reads, Delta10q_Percentage)

#Fisher's Exact Test

fisher_results <- combined_counts %>%
  filter(Protocol != 7) %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, NonDelta10q_Reads,
        direct_rna$Delta10q_Reads, direct_rna$NonDelta10q_Reads),
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
#Protocol Delta10q_Reads NonDelta10q_Reads Delta10q_Percentage  Log2FC Fisher_p
#Ill-SR                      8              96                7.69 -0.0566 1e+0
#ONT-Amplicon-1          29521            37288               44.2   2.47   5.44e-40
#ONT-Amplicon-2           3723             4442               45.6   2.51   1.20e-41
#ONT-Capture              2921            20566               12.4   0.637  2.66e-2
#ONT-WT-A                   14              161                8     0      1e+0
#ONT-WT-B                   70              303               18.8   1.23   9.41e-5
#PacBio-WT                5926             67278              8.10  0.0171  1e+0

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
    NonDelta10q_Reads = sum(NonDelta10q_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + NonDelta10q_Reads)
  )

# Extract Direct RNA group counts
direct_rna_group <- combined_groups %>%
  filter(Protocol_Group == "Direct_RNA") %>%
  select(Delta10q_Reads, NonDelta10q_Reads, Delta10q_Percentage)

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
      c(Delta10q_Reads, NonDelta10q_Reads,
        direct_rna_group$Delta10q_Reads, direct_rna_group$NonDelta10q_Reads),
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
#Protocol_Group           Delta10q_Reads NonDelta10q_Reads Delta10q_Percentage Log2FC Fisher_p
#PCR_cDNA_only                      8931            88308               9.18   0.199 6.00e-1
#PCR_cDNA_plus_enrichment          33244            41730               44.3   2.47  2.5e-40


#PCR Bias Statistics Comparison to Short-read RNA-seq----
#Individual Protocols ----
#Read data file
data_controls_short <- read_excel("PCRBias_delta10q_Levels_controlsonly_All.xlsx")

#Combine control counts
combined_counts_short <- data_controls_short %>%
  group_by(Protocol) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    NonDelta10q_Reads = sum(NonDelta10q_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + NonDelta10q_Reads)
  )

# Separate out Protocol 7 (Direct RNA-seq)
short_rna <- combined_counts_short %>%
  filter(Protocol == 8) %>%
  select(Delta10q_Reads, NonDelta10q_Reads, Delta10q_Percentage)

#Fisher's Exact Test

fisher_results_short <- combined_counts_short %>%
  filter(Protocol != 8) %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, NonDelta10q_Reads,
        short_rna$Delta10q_Reads, short_rna$NonDelta10q_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()


#Log2 Fold-Change
logfc_data_short <- combined_counts_short %>%
  filter(Protocol != 8) %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (short_rna$Delta10q_Percentage + 1e-6))
  )


#Comnbined Log2FC and Fisher
results_individual_short <- logfc_data_short %>%
  left_join(fisher_results_short %>% select(Protocol, Fisher_p), by = "Protocol")

# Optional: Print results
print(results_individual_short)

#Result - Individual Protocols compared to short-read:
#Protocol Delta10q_Reads NonDelta10q_Reads Delta10q_Percentage  Log2FC Fisher_p
#PacBio-WT         60722            378723            13.8     0.568   2.29e-4
#ONT-WT-A            49               391             11.1     0.257   3.17e-1
#ONT-WT-B           173              1052             14.1     0.600   1.78e-3
#ONT-Amp-1         75160            89392             45.7     2.29    2.41e-103
#ONT-Amp-2         15011            18381             45.0     2.27    1.61e-98
#ONT-Capture       10149            56652             15.2     0.705   3.20e-6

#Grouped by amount of PCR against Short-read:  ----

data_controls_group_short <- data_controls_short %>%
  mutate(
    Protocol_Group = case_when(
      Protocol %in% c(1, 2, 3, 6) ~ "PCR_cDNA_only",
      Protocol %in% c(4, 5) ~ "PCR_cDNA_plus_enrichment",
      Protocol == 8 ~ "Short_RNA",
      TRUE ~ "Other"
    )
  )

# Combine control counts per protocol group
combined_groups_short <- data_controls_group_short %>%
  filter(Protocol_Group != "Other") %>%
  group_by(Protocol_Group) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    NonDelta10q_Reads = sum(NonDelta10q_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + NonDelta10q_Reads)
  )

# Extract Direct RNA group counts
Short_rna_group <- combined_groups_short %>%
  filter(Protocol_Group == "Short_RNA") %>%
  select(Delta10q_Reads, NonDelta10q_Reads, Delta10q_Percentage)

# Log2 Fold-Change for grouped protocols vs Direct RNA
logfc_data_group_short <- combined_groups_short %>%
  filter(Protocol_Group != "Short_RNA") %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (direct_rna_group$Delta10q_Percentage + 1e-6))
  )

# Fisher's Exact Test for grouped protocols vs Direct RNA
fisher_results_group_short <- combined_groups_short %>%
  filter(Protocol_Group != "Short_RNA") %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, NonDelta10q_Reads,
        Short_rna_group$Delta10q_Reads, Short_rna_group$NonDelta10q_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()

# Combine Log2FC and Fisher p-values for grouped data
results_grouped_short <- logfc_data_group_short %>%
  left_join(fisher_results_group_short %>% select(Protocol_Group, Fisher_p), by = "Protocol_Group")

# Print grouped protocol results
print(results_grouped_short)
#Protocol_Group           Delta10q_Reads NonDelta10q_Reads Delta10q_Percentage Log2FC  Fisher_p
#PCR_cDNA_only                     71093            436818                14.0  0.807 1.32e-  4
#PCR_cDNA_plus_enrichment          90171            107773                45.6  2.51  7.27e-103





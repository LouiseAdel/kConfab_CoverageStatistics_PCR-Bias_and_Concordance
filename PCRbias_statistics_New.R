
setwd("C:/Users/louis/Documents/Phd/ENIGMA/Udkast Artikel")

# Libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)

# -----------------------------
# File1: replicate-level data (Ill-SR comparisons)
# -----------------------------
df <- read_excel("delta10q_levels_barplot.xlsx")

# Drop duplicate summary column
df <- df %>% select(-Delta10q_Percentage)

# Reshape controls to long format
controls_long <- df %>%
  filter(SampleID == "Controls") %>%
  select(Protocol, starts_with("Control")) %>%
  pivot_longer(cols = starts_with("Control"),
               names_to = "Control", values_to = "Value") %>%
  filter(!is.na(Value))

# Assign groups
controls_long <- controls_long %>%
  mutate(Group = case_when(
    Protocol %in% c("ONT-Amp.-1", "ONT-Amp.-2") ~ "Amplicon",
    Protocol %in% c("PacBio-WT", "ONT-WT-A", "ONT-WT-B", "ONT-Capture") ~ "Capture+WT",
    Protocol == "ONT-Direct" ~ "Direct",
    Protocol == "Ill-SR" ~ "Short-read"
  ))

# Extract Ill-SR replicate values
ill_sr_values <- controls_long %>%
  filter(Group == "Short-read") %>%
  pull(Value)
ill_sr_median <- median(ill_sr_values, na.rm = TRUE)

# ---- Wilcoxon vs Ill-SR (protocol-level)
compare_protocol_to_sr <- function(df, protocol_name) {
  g_values <- df %>% filter(Protocol == protocol_name) %>% pull(Value)
  test <- wilcox.test(g_values, ill_sr_values, exact = TRUE, correct = FALSE)
  tibble(
    Level = "Protocol",
    Name = protocol_name,
    N_test = length(g_values),
    N_ref = length(ill_sr_values),
    Median_test = median(g_values),
    Median_ref = ill_sr_median,
    log2FC_vs_SR = log2(median(g_values) / ill_sr_median),
    Test = "Wilcoxon_vs_SR",
    Statistic = test$statistic,
    p_value = test$p.value
  )
}

wilcox_protocol_vs_sr <- controls_long %>%
  distinct(Protocol) %>%
  filter(Protocol != "Ill-SR") %>%
  pull(Protocol) %>%
  map_dfr(~compare_protocol_to_sr(controls_long, .x))

# ---- Wilcoxon vs Ill-SR (group-level)
compare_group_to_sr <- function(df, group_name) {
  g_values <- df %>% filter(Group == group_name) %>% pull(Value)
  test <- wilcox.test(g_values, ill_sr_values, exact = TRUE, correct = FALSE)
  tibble(
    Level = "Group",
    Name = group_name,
    N_test = length(g_values),
    N_ref = length(ill_sr_values),
    Median_test = median(g_values),
    Median_ref = ill_sr_median,
    log2FC_vs_SR = log2(median(g_values) / ill_sr_median),
    Test = "Wilcoxon_vs_SR",
    Statistic = test$statistic,
    p_value = test$p.value
  )
}

wilcox_group_vs_sr <- bind_rows(
  compare_group_to_sr(controls_long, "Amplicon"),
  compare_group_to_sr(controls_long, "Capture+WT"),
  compare_group_to_sr(controls_long, "Direct")
)

# -----------------------------
# File2: matched controls with read counts (Direct comparisons)
# -----------------------------
file2 <- read_excel("PCRBias_delta10q_Levels_Direct.xlsx")

# Assign groups
file2 <- file2 %>%
  mutate(Group = case_when(
    Protocol %in% c("ONT-Amp.-1", "ONT-Amp.-2") ~ "Amplicon",
    Protocol %in% c("PacBio-WT", "ONT-WT-A", "ONT-WT-B", "ONT-Capture") ~ "Capture+WT",
    Protocol == "ONT-Direct" ~ "Direct",
    Protocol == "Ill-SR" ~ "Short-read"
  ))

# Extract Direct replicate percentages
direct_values <- file2 %>%
  filter(Group == "Direct") %>%
  pull(Delta10q_Percentage)
direct_median <- median(direct_values, na.rm = TRUE)

# ---- Fisher vs Direct (protocol-level, log2FC from medians)
compare_protocol_to_direct <- function(df, protocol_name) {
  g_values <- df %>% filter(Protocol == protocol_name) %>% pull(Delta10q_Percentage)
  sub <- df %>% filter(Protocol == protocol_name)
  mat <- matrix(c(sum(sub$Delta10q_Reads),
                  sum(sub$NonDelta10q_Reads),
                  sum(file2$Delta10q_Reads[file2$Group == "Direct"]),
                  sum(file2$NonDelta10q_Reads[file2$Group == "Direct"])),
                nrow = 2, byrow = TRUE)
  test <- fisher.test(mat)
  tibble(
    Level = "Protocol",
    Name = protocol_name,
    N_test = length(g_values),
    N_ref = length(direct_values),
    Median_test = median(g_values),
    Median_ref = direct_median,
    log2FC_vs_Direct = log2(median(g_values) / direct_median),
    Test = "Fisher_vs_Direct",
    Statistic = NA,
    p_value = test$p.value
  )
}

fisher_protocol_vs_direct <- file2 %>%
  distinct(Protocol) %>%
  filter(Protocol != "ONT-Direct") %>%
  pull(Protocol) %>%
  map_dfr(~compare_protocol_to_direct(file2, .x))

# ---- Fisher vs Direct (group-level, log2FC from medians)
compare_group_to_direct <- function(df, group_name) {
  g_values <- df %>% filter(Group == group_name) %>% pull(Delta10q_Percentage)
  sub <- df %>% filter(Group == group_name)
  mat <- matrix(c(sum(sub$Delta10q_Reads),
                  sum(sub$NonDelta10q_Reads),
                  sum(file2$Delta10q_Reads[file2$Group == "Direct"]),
                  sum(file2$NonDelta10q_Reads[file2$Group == "Direct"])),
                nrow = 2, byrow = TRUE)
  test <- fisher.test(mat)
  tibble(
    Level = "Group",
    Name = group_name,
    N_test = length(g_values),
    N_ref = length(direct_values),
    Median_test = median(g_values),
    Median_ref = direct_median,
    log2FC_vs_Direct = log2(median(g_values) / direct_median),
    Test = "Fisher_vs_Direct",
    Statistic = NA,
    p_value = test$p.value
  )
}

fisher_group_vs_direct <- bind_rows(
  compare_group_to_direct(file2, "Amplicon"),
  compare_group_to_direct(file2, "Capture+WT"),
  compare_group_to_direct(file2, "Short-read")
)

# -----------------------------
# Master sheet: combine all results
# -----------------------------
master_results <- bind_rows(
  wilcox_protocol_vs_sr,
  wilcox_group_vs_sr,
  fisher_protocol_vs_direct,
  fisher_group_vs_direct
)

# -----------------------------
# Export results
# -----------------------------
write_xlsx(
  list(
    "Wilcoxon_vs_ShortRead_PROTOCOL" = wilcox_protocol_vs_sr,
    "Wilcoxon_vs_ShortRead_GROUP"    = wilcox_group_vs_sr,
    "Fisher_vs_Direct_PROTOCOL"      = fisher_protocol_vs_direct,
    "Fisher_vs_Direct_GROUP"         = fisher_group_vs_direct,
    "Master_All"                     = master_results
  ),
  "Delta10q_Comparisons.xlsx"
)


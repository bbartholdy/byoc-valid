# Compare species abundance for biofilm samples and modern plaque/calculus
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(phyloseq)
library(ANCOMBC)
taxa_table <- readr::read_tsv("04-analysis/decontam/post-decontam_taxatable.tsv")
analysis_metadata <- readr::read_tsv("01-documentation/analysis-metadata.tsv")
experiment_metadata <- readr::read_tsv("01-documentation/experiment-metadata.tsv")
source("02-scripts/functions.R")

# Data prep ---------------------------------------------------------------

taxa_table_long <- taxa_table %>%
  pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "count")

# convert otu tables to phyloseq object
byoc_matrix <- taxa_table_long %>%
  filter(str_detect(sample, "SYN")) %>%
  pivot_wider(names_from = "sample", values_from = "count") %>%
  column_to_rownames(var = "#OTU ID") 

byoc_otu <- byoc_matrix %>%
  otu_table(taxa_are_rows = T)

byoc_meta <- experiment_metadata %>%
  filter(`#SampleID` %in% analysis_metadata$`#SampleID`) %>%
  column_to_rownames(var = "#SampleID") %>% 
  sample_data()

byoc_phyloseq <- phyloseq(byoc_otu, byoc_meta)

plaque_matrix <- taxa_table_long %>%
  filter(
    sample %in% filter(
      analysis_metadata, str_detect(Env, "calculus|plaque") # comparative samples and model calculus
      )$`#SampleID`
    ) %>%
  pivot_wider(names_from = "sample", values_from = "count") %>%
  column_to_rownames(var = "#OTU ID") %>%
  otu_table(taxa_are_rows = T)

plaque_meta <- analysis_metadata %>%
  filter(str_detect(Env, "calculus|plaque")) %>%
  column_to_rownames(var = "#SampleID") %>%
  sample_data()

plaque_phyloseq <- phyloseq(plaque_matrix, plaque_meta)


# Differential abundance --------------------------------------------------

# within experiment
byoc_da <- ANCOMBC::ancombc(
  byoc_phyloseq, 
  formula = "Env",
  group = "Env",
  global = F,
  p_adj_method = "fdr"
  )

# plaque and calculus
plaque_da <- ANCOMBC::ancombc(
  plaque_phyloseq, 
  formula = "Env",
  group = "Env",
  global = F,
  p_adj_method = "fdr"
  )

byoc_log_abund <- bias_correct(byoc_da, byoc_otu)

byoc_logf_full <- as_tibble(byoc_da$res$lfc, rownames = "species") %>%
  pivot_longer(-species, values_to = "lfc") %>%
  full_join(
    as_tibble(byoc_da$res$se, rownames = "species") %>%
      pivot_longer(-species, values_to = "se"),
    by = c("species", "name")
  ) %>%
  full_join(
    as_tibble(byoc_da$res$q, rownames = "species") %>%
      pivot_longer(-species, values_to = "q_value"),
    by = c("species", "name")
  ) %>%
  mutate(upper = lfc + se,
         lower = lfc - se,
         name = str_remove(name, "^Env"),
         abn = case_when(sign(lfc) == -1 ~ "model_calculus",
                         TRUE ~ name)) %>%
  rename(env = name)


# table with log-fold change statistics between artificial calculus and others
# negative lfc means higher abundance in artificial calculus
# positive lfc means lower abundance in art calculus
plaque_logf_change <- as_tibble(plaque_da$res$lfc, rownames = "species") #%>%
plaque_logf_se <- as_tibble(plaque_da$res$se, rownames = "species") #%>%
plaque_logf_q <- as_tibble(plaque_da$res$q_val, rownames = "species")

plaque_logf_full <- plaque_logf_change %>%
  pivot_longer(-species, values_to = "lfc") %>%
  full_join(
    plaque_logf_se %>%
      pivot_longer(-species, values_to = "se"),
    by = c("species", "name")
  ) %>%
  full_join(
    plaque_logf_q %>%
      pivot_longer(-species, values_to = "q_value"),
    by = c("species", "name")
  ) %>%
  mutate(lower = lfc - se,
         upper = lfc + se,
         name = str_remove(name, "^Env"),
         abn = case_when(sign(lfc) == -1 ~ "model_calculus",
                         TRUE ~ name)) %>%
  rename(env = name)

write_tsv(byoc_logf_full, "04-analysis/diff-abund/byoc_logf-full.tsv")
write_tsv(byoc_logf_full, "05-results/byoc_logf-full.tsv")
write_tsv(plaque_logf_full, "04-analysis/diff-abund/plaque_logf-full.tsv")
write_tsv(plaque_logf_full, "05-results/plaque_logf-full.tsv")




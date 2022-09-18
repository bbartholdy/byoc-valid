# Compare species abundance for biofilm samples and modern plaque/calculus
library(tidyverse)
library(phyloseq)
library(ANCOMBC)
taxa_table <- readr::read_tsv("05-results/post-decontam_taxatable.tsv")
analysis_metadata <- readr::read_tsv("01-documentation/analysis-metadata.tsv")
experiment_metadata <- readr::read_tsv("01-documentation/experiment-metadata.tsv")
#bac_properties <- read_tsv("01-documentation/species-properties.tsv")
# function to calculate bias-corrected log-observed abundances
# from vignette("ANCOMBC")

#' function to calculate bias-corrected log-observed abundances
#' from `vignette("ANCOMBC")`
#' 
#' @param x `ancombc` object.
bias_correct <- function(x, otu_table) {
  require(microbiome)
  samp_frac <- x$samp_frac
  samp_frac[is.na(samp_frac)] <- 0 # replace NAs with 0
  log_otu <- log(microbiome::abundances(otu_table) + 1)
  log_otu_adj <- exp(t(t(log_otu) - samp_frac)) %>% # bias corrected abund
    as_tibble(rownames = "species")
  return(log_otu_adj)
}


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
         abn = case_when(sign(lfc) == -1 ~ "byoc_calculus",
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
         abn = case_when(sign(lfc) == -1 ~ "byoc_calculus",
                         TRUE ~ name)) %>%
  rename(env = name)

write_tsv(byoc_logf_full, "05-results/byoc_logf-full.tsv")
write_tsv(plaque_logf_full, "05-results/plaque_logf-full.tsv")




























# function to get pallette for top species/genus

function(x, pallette, n) {
  col_pallette <- palette(n = n)
  names(col_pallette) <- levels(x)
  col_scale <- scale_fill_manual(name = "", values = col_pallette)
}
#my_pallette <- c(viridisLite::viridis(n = length(top_species_names)), "grey") 
#names(my_pallette) <- levels(top_abundance$top_species)
#my_scale <- scale_fill_manual(name = "top_species", values = my_pallette)




# Diff abundance ----------------------------------------------------------

# exp_day <- ANCOMBC::ancombc(
#   experiment_phyloseq, 
#   formula = "Env",
#   group = "Env",
#   global = T,
#   p_adj_method = "fdr")
# 
# ANCOMBC::ancom(experiment_phyloseq, main_var = "Env")
# 
# # from vignette("ANCOMBC")
# samp_frac <- exp_day$samp_frac
# # Replace NA with 0
# samp_frac[is.na(samp_frac)] = 0 
# # add 1 to counts for log-transformation
# log_exp_otu <- log(microbiome::abundances(experiment_otu) + 1)
# 
# # Adjust the log observed abundances
# log_exp_otu_adj = t(t(log_exp_otu))
# # Bias-corrected log observed abundances
#   # Show the first 6 samples
# round(log_exp_otu_adj[, 1:6], 2) %>% 
#   as_tibble(rownames = "species") %>%
#   .$SYN001.A0101
# 
# class(otu_table(taxatable_matrix, taxa_are_rows = T))



#da_env_taxa <- ANCOMBC::ancombc(otu_phyloseq, formula = "Env")

# taxatable_abundance_long <- taxatable %>%
#   pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "count") %>% 
#   rename(species = `#OTU ID`) %>% 
#   group_by(sample) %>%
#   mutate(rel_abund = count / sum(count))
# 
# species_abundance_long <- taxatable_abundance_long %>%
#   #bind_rows(taxatable_abundance) %>%
#   mutate( # fix species names
#     species = str_remove(species, "\\["), # remove square brackets from names
#     species = str_remove(species, "\\]") # probably not the most elegant solution
#     ) 
#   
# 
# top_species_names <- species_abundance_long %>% 
#   group_by(sample) %>% 
#   arrange(desc(rel_abund)) %>% 
#   slice_head(n = 4) %>% # take the top 4 species from each sample
#   .$species %>% 
#   unique()
# 
# genus_abundance_long <- species_abundance_long %>%
#   mutate(genus = str_extract(species, "^([\\w\\-]+)")) %>% 
#   dplyr::select(!c(species,rel_abund)) %>% 
#   group_by(sample, genus) %>%
#   summarise(count = sum(count, na.rm = T)) %>% 
#   mutate(rel_abund = count / sum(count))
# 
# top_genus_names <- genus_abundance_long %>% 
#   group_by(sample) %>% 
#   arrange(desc(rel_abund)) %>% 
#   slice_head(n = 4) %>% 
#   .$genus %>% 
#   unique()

# Export data -------------------------------------------------------------

#write_tsv(species_abundance_long, "05-results/species-abundance_long.tsv")
#write_tsv(genus_abundance_long, "05-results/genus-abundance_long.tsv")


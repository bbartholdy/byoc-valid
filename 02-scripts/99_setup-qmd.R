## Script to load all required objects into qmd documents

# Metagenomic analysis ----------------------------------------------------


# Helper objects
env_controls <- c("indoor_air", "sediment", "stool", "skin") # vector to remove environmental controls

day_order <- dna_experiment_metadata %>%
  filter(`#SampleID` %in% dna_analysis_metadata$`#SampleID`) %>% 
  mutate(Env = factor(
    Env, 
    levels = c("saliva", "medium", "byoc_calculus") # force level order so it doesn't order alphabetically
    ) 
  ) %>% 
  group_by(Env) %>%
  arrange(day, .by_group = T) %>%
  filter(`#SampleID` %in% colnames(sourcetracker2)) %>%
  mutate(
    rm = if_else(`#SampleID` %in% dna_analysis_metadata$`#SampleID`, F, T),
    col = case_when(rm == T ~ "red", # excluded samples coloured red
                    rm == F ~ "black"))



# convert species table to long format for easier wrangling
species_counts_long <- otu_decontam %>%
  pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "count") %>%
  rename(species = `#OTU ID`) #%>%
  #filter(sample %in% dna_analysis_metadata$`#SampleID`) # does not change data frame size

# collapse species counts into genus counts
genus_counts_long <- species_counts_long %>%
  mutate(genus = stri_extract(species, regex = "\\w+")) %>%
  group_by(sample, genus) %>%
  summarise(count = sum(count)) %>% 
  group_by(sample) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup()


# convert alpha diversity data to long format

alpha_div_long <- alpha_div %>%
  pivot_longer(cols = where(is.numeric), names_to = "index") %>%
  left_join(dna_metadata, by = c("sample" = "#SampleID")) %>%
  left_join(dna_experiment_metadata, by = c("sample" = "#SampleID", "Env")) %>%
  mutate(day_grouped = case_when( # group days to increase sample size
    day < 6 ~ "inoc",
    day > 6 & day < 24 ~ "treatm",
    day == 24 ~ "model"),
    day_grouped = factor(day_grouped, levels = c("inoc", "treatm", "model"))
  )

# summary stats for alpha diversity in experiment
alpha_summ_byoc <- alpha_div_long %>%
  filter(
    stri_detect(sample, fixed = "SYN") # isolate samples from this study
  ) %>%
  group_by(index, day_grouped) %>% # group by sample type
  summarise(mean = mean(value),
    sd = sd(value))


# summary stats for grouped samples by Env
alpha_summ_all <- alpha_div_long %>%
  filter(
    Env != "skin",
    Env != "sediment",
    Env != "stool",
    Env != "indoor_air") %>%
  group_by(index, Env) %>%
  summarise(mean = mean(value),
    sd = sd(value))


# table containing the oxygen tolerance of species
species_properties <- bac_properties %>%
  right_join(species_counts_long, by = "species") %>%
  mutate(genus = stri_extract(species, regex = "\\w+")) %>%
  left_join(
    rename(genus_oxytol, genus_oxytol = `Oxygen tolerance`), by = "genus"
  ) %>%
  mutate(`Oxygen tolerance` = case_when(
    is.na(`Oxygen tolerance`) ~ genus_oxytol, # if no species-level info available
    TRUE ~ `Oxygen tolerance`)
  ) %>%
  group_by(sample) %>%
  mutate(
    abs = if_else(is.na(abs), FALSE, abs), # is ABS species? TRUE/FALSE
    rel_abund = count / sum(count)
  ) %>%
  ungroup()

# sPCA analyses

byoc_explain_var <- spca_byoc$prop_expl_var$X # experiment samples
# projection of samples onto PCs
byoc_princomp <- spca_byoc$x %>% # experiment samples
  as_tibble(rownames = "sample") %>%
  left_join(dna_experiment_metadata, by = c("sample" = "#SampleID"))
# loadings of species
byoc_pca_loadings <- spca_byoc$loadings$X %>% # experiment samples
  as_tibble(rownames = "species") %>%
  #left_join(bac_properties, by = "species") %>% 
  left_join(species_properties, by = "species") %>% 
  dplyr::select(species,PC1,PC2, `Oxygen tolerance`) %>%
  arrange(desc(abs(PC1))) %>%
  distinct(species, .keep_all = T)

compar_explain_var <- spca_species$prop_expl_var$X # comparative samples
# projection of samples onto PCs
compar_princomp <- spca_species$x %>% # comparative samples
  as_tibble(rownames = "sample") %>%
  left_join(dna_metadata, by = c("sample" = "#SampleID"))
compar_pca_loadings <- spca_species$rotation %>% # comparative samples
  as_tibble(rownames = "species") %>%
  left_join(species_properties, by = "species") %>% 
  dplyr::select(species,PC1,PC2,PC3,`Oxygen tolerance`) %>%
  arrange(desc(abs(PC1))) %>% 
  distinct(species, .keep_all = T)


# FTIR analysis -----------------------------------------------------------



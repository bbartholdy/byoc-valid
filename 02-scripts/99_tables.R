# Metagenomics ------------------------------------------------------------

dna_samples <- dna_experiment_metadata %>%
  #filter(Study == "this_study") %>%
  dplyr::select(`#SampleID`, Env, day) #%>%
# dna_samples %>%
#   knitr::kable(col.names = c("Sample ID", "Sample type", "Sampling day"))

dna_samples_table <- dna_experiment_metadata %>%
  mutate(Env = factor(Env, levels = c("saliva", "medium", "model_calculus"))) %>% 
  group_by(Env, day) %>%
  count() %>%
  arrange(day, Env)


# FTIR --------------------------------------------------------------------

ftir_samples_table <- ftir_metadata %>%
  filter(source == "Artificial") %>%
  mutate(sample_type = as.factor(sample_type)) %>% 
  group_by(day) %>% 
  summarise(
    sample_type = sample_type,
    n = n(),
    weight = mean(weight_mg, na.rm = T)
  ) %>%
  distinct(day, .keep_all = T) %>% 
  relocate(sample_type, .before = day)
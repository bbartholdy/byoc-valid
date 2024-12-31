# Metagenomics ------------------------------------------------------------

dna_samples <- dna_experiment_metadata %>%
  dplyr::select(`#SampleID`, Env, day)

dna_samples_table <- dna_experiment_metadata %>%
  mutate(Env = factor(Env, levels = c("saliva", "medium", "model_calculus"))) %>% 
  group_by(Env, day) %>%
  count() %>%
  arrange(day, Env)


# FTIR --------------------------------------------------------------------

ftir_samples_table <- ftir_metadata %>%
  filter(
    !is.na(day) &
      !is.na(analysis_id)
  ) %>%
  group_by(day) %>% 
  summarise(
    n = n(),
    weight = mean(weight_mg, na.rm = T)
  ) %>%
  distinct(day, .keep_all = T)

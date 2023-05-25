# Figures generated from the metagenomic analysis

# clean labels

oral_source_order <- c(
  "saliva",
  "buccal_mucosa",
  "medium",
  "plaque",
  "vitro_biofilm",
  "model_calculus",
  "modern_calculus"
)

oral_source_labels <- c(
  "buccal_mucosa" = "buccal mucosa",
  "model_calculus" = "model calculus",
  "modern_calculus" = "modern calculus",
  "saliva" = "saliva",
  "subgingival_plaque" = "sub plaque",
  "supragingival_plaque" = "supra plaque",
  "vitro_biofilm" = "in vitro biofilm"
)

byoc_source_labels <- c(
  "saliva" = "saliva",
  "medium" = "medium",
  "model_calculus" = "model calculus"
)

# helper objects to set x-axis limits
min_byoc_PC1 <- min(byoc_pca_loadings$PC1)
min_byoc_PC2 <- min(byoc_pca_loadings$PC2)
max_byoc_PC1 <- max(byoc_pca_loadings$PC1)
max_byoc_PC2 <- max(byoc_pca_loadings$PC2)
min_comp_PC1 <- min(compar_pca_loadings$PC1)
min_comp_PC2 <- min(compar_pca_loadings$PC2)
max_comp_PC1 <- max(compar_pca_loadings$PC1)
max_comp_PC2 <- max(compar_pca_loadings$PC2)

## Alpha-diversity

# facet labels

indices <- c(
  "shannon" = "Shannon Index",
  "pilou_even" = "Pielou Evenness",
  "richness" = "Number of species"
)

div_byoc_fig <- alpha_div_long %>%
  filter(
    str_detect(sample, "SYN"),
    index == "shannon" | index == "pilou_even" | index == "richness"
  ) %>%
  mutate(index = factor(index, levels = c("shannon", "pilou_even", "richness"))) %>%
  ggplot(aes(x = day_grouped, y = value)) +
  geom_violin(aes(col = day_grouped, fill = day_grouped), alpha = 0.5) +
  geom_boxplot(width = 0.12) +
  facet_wrap(~ index, scales = "free_y", labeller = as_labeller(indices))

div_compar_fig <- alpha_div_long %>%
  filter(!Env %in% env_controls,
    index %in% c("shannon", "richness", "pilou_even")
  ) %>%
  mutate(index = factor(index, levels = c("shannon", "pilou_even", "richness"))) %>%
  mutate(Env = case_when(str_detect(Env, "plaque") ~ "plaque",
    TRUE ~ Env),
    Study = case_when(
      is.na(Study) ~ "other",
      TRUE ~ Study
    ),
    Env = factor(Env, levels = oral_source_order)) %>%
  ggplot(aes(x = Env, y = value)) +
    geom_violin(aes(col = Env, fill = Env), alpha = 0.5) +
    #geom_boxplot(width = 0.1, fill = "transparent") +
    geom_jitter(aes(shape = Study), width = 0.1, alpha = 0.6) +
    facet_wrap(~ index, scales = "free_y", ncol = 1, labeller = as_labeller(indices))

## Beta-diversity

byoc_spca_base <- byoc_princomp %>%
  ggplot(aes(x = PC1, y = PC2, col = as_factor(day), shape = Env)) +
  geom_point(size = 4, stroke = 1) +
  geom_vline(xintercept = 0, size = 0.2) +
  geom_hline(yintercept = 0, size = 0.2)

byoc_comp_1 <- byoc_pca_loadings %>%
  arrange(desc(PC1)) %>% 
  mutate(species = fct_reorder(species, desc(PC1))) %>%
  slice(c(
    1:20,
    seq(from = nrow(byoc_pca_loadings)-19, to = nrow(byoc_pca_loadings))
  )) %>%
  ggplot(aes(x = PC1, y = species, fill = `Oxygen tolerance`)) +
  geom_bar(stat = "identity", width = 0.8) +
  #scale_x_continuous(limits = c(min_byoc_PC1,max_byoc_PC1)) +
  #scale_y_discrete(position = "right") +
  scale_fill_manual(
    values = viridisLite::turbo(
      n = length(levels(as.factor(bac_properties$`Oxygen tolerance`)))
    )
  )

byoc_comp_2 <- byoc_pca_loadings %>%
  arrange(desc(PC2)) %>%
  mutate(species = fct_reorder(species, desc(PC2))) %>%
  slice(c(
    1:20,
    seq(from = nrow(byoc_pca_loadings)-19, to = nrow(byoc_pca_loadings))
  )) %>%
  ggplot(aes(y = PC2, x = species, fill = `Oxygen tolerance`)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = viridisLite::turbo(
      n = length(levels(as.factor(bac_properties$`Oxygen tolerance`)))
    )
  )

compar_spca_base <- compar_princomp %>%
  mutate(Study = case_when(
    Study != "this_study" | is.na(Study) ~ "other",
    TRUE ~ Study
  )) %>% 
  ggplot(aes(x = PC1, y = PC2, col = Study, shape = Env)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, size = 0.2) +
  geom_hline(yintercept = 0, size = 0.2)

compar_comp_1 <- compar_pca_loadings %>%
  arrange(desc(PC1)) %>% 
  mutate(species = fct_reorder(species, desc(PC1))) %>%
  slice(c(
    1:20,
    seq(from = nrow(compar_pca_loadings)-19, to = nrow(compar_pca_loadings))
  )) %>%
  ggplot(aes(x = PC1, y = species, fill = `Oxygen tolerance`)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_x_continuous(position = "bottom", limits = c(min_comp_PC1,max_comp_PC1)) +
  scale_fill_manual(
    values = viridisLite::turbo(
      n = length(levels(as.factor(bac_properties$`Oxygen tolerance`)))
    )
  )

compar_comp_2 <- compar_pca_loadings %>%
  arrange(desc(PC2)) %>%
  mutate(species = fct_reorder(species, desc(PC2))) %>%
  slice(c(
    1:20,
    seq(from = nrow(compar_pca_loadings)-19, to = nrow(compar_pca_loadings))
  )) %>%
  ggplot(aes(y = PC2, x = species, fill = `Oxygen tolerance`)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(limits = c(min_comp_PC2,max_comp_PC2), position = "right") +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values = viridisLite::turbo(
      n = length(levels(as.factor(bac_properties$`Oxygen tolerance`)))
    )
  )

## Core genera

genus_counts_short <- genus_counts_long %>%
  left_join(dna_metadata, by = c("sample" = "#SampleID")) %>%
  filter(str_detect(Env, "plaque|calculus")) %>%
  group_by(Env, genus) %>%
  summarise(rel_abund = mean(rel_abund)) %>% # mean relative abundance in each sample type
  mutate(
    genus = case_when(rel_abund <= 0.05 ~ NA_character_,
      TRUE ~ genus),
    label = str_extract(genus, "\\w{2}") # create 2-letter abbreviation for labels
  ) #%>% 
#mutate(genus = as.factor(genus))
genus_names <- sort(unique(genus_counts_short$genus)[-1])
#genus_cols <- viridisLite::viridis(length(genus_names))
#names(genus_cols) <- genus_names
#filter(rel_abund > 0.05) %>% 
core_genera_fig <- genus_counts_short %>% 
  ggplot(aes(x = "", y = rel_abund, fill = genus, label = label)) +
  geom_col(col = "grey50") +
  coord_polar(theta = "y", start = 0) +
  geom_label(position = position_stack(vjust = 0.5), col = "white", show.legend = F) +
  facet_wrap(~ Env, labeller = labeller(Env = c(
    "model_calculus" = "model calculus",
    "modern_calculus" = "modern calculus",
    "subgingival_plaque" = "subgingival plaque",
    "supragingival_plaque" = "supragingival plaque"
  )))

## Log fold change

byoc_logf_base <- byoc_logf_full %>% 
  group_by(env) %>% 
  arrange(desc(abs(lfc))) %>%
  slice_head(n = 30) %>% # causes loss of some env values
  ungroup() %>%
  arrange(desc(abs(lfc))) %>%
  slice_head(n = 30) %>% # extract top 20 lfc from modern_calculus
  bind_rows(filter( # recombine with other env values so all are included in plot
    byoc_logf_full, 
    env != "saliva",
    species %in% .$species
  )
  ) %>%
  left_join(bac_properties, by = "species") %>% 
  ggplot(aes(x = lfc, y = reorder(species, lfc), col = `Oxygen tolerance`, shape = abn)) +
  geom_point(size = 1.5) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ env)

compar_da_pca <- plaque_logf_full %>%
  inner_join(compar_pca_loadings, by = "species")

compar_logf_base <- compar_da_pca %>%
  arrange(desc(abs(lfc))) %>%
  slice_head(n = 90) %>%
  ggplot(aes(x = lfc, y = reorder(species, lfc), col = `Oxygen tolerance`, shape = abn)) +
  geom_point() +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ env)

compar_logf_pc1 <- compar_da_pca %>%
  arrange(desc(abs(PC1))) %>% # order by absolute value of loadings
  slice_head(n = 90) %>% # show top 30 species (3 groups * 30 = 90)
  ggplot(aes(x = lfc, y = reorder(species, lfc), col = `Oxygen tolerance`, shape = abn)) +
  geom_point(size = 1.5) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ env, labeller = labeller(
    env = oral_source_labels
  ))

compar_logf_pc2<- compar_da_pca %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n = 90) %>% 
  mutate(`Oxygen tolerance` = factor(
    `Oxygen tolerance`, levels = c("aerobe", "anaerobe", "facultative anaerobe"))
  ) %>%
  ggplot(aes(x = lfc, y = reorder(species, lfc), col = `Oxygen tolerance`, shape = abn)) +
  geom_point(size = 1.5) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~ env, labeller = labeller(
    env = oral_source_labels
  ))


# Supplementary figures ---------------------------------------------------

# species composition

  # Top positive loadings on PC1
species_pos_pc1 <- clr_compar_long %>%
  left_join(pca_loadings, by = "species") %>%
  left_join(as_tibble(spca_species$x, rownames = "sample"), by = "sample") %>%
  left_join(select(dna_metadata, `#SampleID`, Env), by = c("sample" = "#SampleID")) %>% 
  #left_join(analysis_colours, by = c("sample" = "#SampleID")) %>%
  filter(species %in% arrange(pca_loadings, desc(PC1))$species[1:200]) %>%
  #arrange(sample, PC1.y) %>%
  #ggplot(aes(x = reorder(sample, PC1.y), y = reorder(species, PC1.x), fill = clr_count)) +
  ggplot(aes(x = fct_reorder(sample, Env), y = fct_reorder(species, PC1.x), fill = clr_count)) +
  geom_tile()

# Top negative loadings on PC1
species_neg_pc1 <- clr_compar_long %>%
  left_join(pca_loadings, by = "species") %>%
  left_join(as_tibble(spca_species$x, rownames = "sample"), by = "sample") %>%
  left_join(select(dna_metadata, `#SampleID`, Env), by = c("sample" = "#SampleID")) %>% 
  #left_join(analysis_colours, by = c("sample" = "#SampleID")) %>%
  #arrange(species, PC1.y) %>% # arrange by species loading (PC1.x)
  filter(species %in% arrange(pca_loadings, PC1)$species[1:200]) %>%
  arrange(sample, PC1.x) %>% # arrange by sample loading (PC1.y)
  ggplot(aes(x = fct_reorder(sample, Env), y = fct_reorder(species, desc(PC1.x)), fill = clr_count)) +
  geom_tile()

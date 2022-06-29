library(mixOmics)
library(vegan)
library(tidyverse)
library(here)

otu_table <- read_tsv(here("05-results/post-decontam_taxatable.tsv"))
metadata <- read_tsv(here("01-documentation/metadata.tsv"))
analysis_metadata <- read_tsv(here("01-documentation/analysis-metadata.tsv"))
experiment_metadata <- read_tsv(here("01-documentation/experiment-metadata.tsv")) %>%
  filter(`#SampleID` %in% analysis_metadata$`#SampleID`)


# Data prep ---------------------------------------------------------------

species_otu_table <- otu_table %>%
  column_to_rownames(var = "#OTU ID") %>%
  t() %>%
  as_tibble(rownames = "sample")

comp_metadata <- analysis_metadata %>%
  filter(
    Env != "medium",
    Env != "sediment", 
    Env != "skin",  
    Env != "stool", 
    Env != "indoor_air",
    #Env != "modern_calculus"
  )

sample_metadata <- experiment_metadata %>%
  filter(`#SampleID` %in% analysis_metadata$`#SampleID`)

species_otu_matrix <- species_otu_table %>%
  #dplyr::select(!contains("Enterococcus")) %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

byoc_otu_matrix <- species_otu_matrix %>%
  subset(rownames(species_otu_matrix) %in% experiment_metadata$`#SampleID`)

comp_otu_matrix <- species_otu_matrix %>%
  subset(rownames(species_otu_matrix) %in% comp_metadata$`#SampleID`)

# Alpha diversity ---------------------------------------------------------


alpha_div_shan <- vegan::diversity(species_otu_matrix)
alpha_div_inv <- vegan::diversity(species_otu_matrix, index = "invsimpson")
alpha_div_unb <- simpson.unb(species_otu_matrix,inverse = T)

alpha_div <- 
  as_tibble(
    alpha_div_shan, rownames = "sample"
    ) %>%
  full_join(
    as_tibble(alpha_div_inv, rownames = "sample"),
    by = "sample"
    ) %>%
  full_join(
    as_tibble(alpha_div_unb, rownames = "sample"),
    by = "sample"
  ) %>%
  rename(shannon = value.x,
         simp_inv = value.y,
         simp_unb = value)
  

write_tsv(alpha_div, here("05-results/alpha-diversity.tsv"))


# Beta diversity ----------------------------------------------------------

# Experiment samples
clr_byoc <- logratio.transfo(byoc_otu_matrix, "CLR", 1)

spca_byoc <- spca(clr_byoc, ncomp = 10, scale = F)
byoc_explain_var <- spca_byoc$prop_expl_var$X

exp_pca_loadings <- spca_byoc$loadings$X %>%
  as_tibble(rownames = "species") %>%
  dplyr::select(species,PC1,PC2) %>%
  arrange(desc(abs(PC1)))

save(spca_byoc, file = here("05-results/spca_byoc.rda"))
write_tsv(exp_pca_loadings, here("05-results/experiment-pca-loadings.tsv"))

# Comparative samples

clr_species <- logratio.transfo(comp_otu_matrix, "CLR", 1)

spca_species <- spca(clr_species, ncomp = 10, scale = F)

spca_compar <- spca_species$x %>%
  as_tibble(rownames = "sample")
compar_explain_var <- spca_species$prop_expl_var$X

pca_loadings <- spca_species$rotation %>%
  as_tibble(rownames = "species") %>% 
  dplyr::select(species, PC1, PC2, PC3) %>%
  arrange(desc(abs(PC1)))

pca_loadings

pca_loadings %>%
  arrange(desc(abs(PC2)))

save(spca_species, file = here("05-results/spca_species.rda"))
write_tsv(pca_loadings, here("05-results/all-pca-loadings.tsv"))

Y <- tibble("#SampleID" = rownames(clr_species)) %>%
  left_join(
    comp_metadata,
    by = "#SampleID"
    ) %>%
  .$Env
species_splsda <- splsda(clr_species, Y)

plotIndiv(species_splsda, group = Y, ind.names = F, legend = T)


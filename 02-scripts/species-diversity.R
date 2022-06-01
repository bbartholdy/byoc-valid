library(vegan)
library(mixOmics)
library(tidyverse)

taxatable <- readr::read_tsv("05-results/species-abundance_long.tsv")
metadata <- readr::read_tsv("01-documentation/metadata.tsv")


# Data prep ---------------------------------------------------------------

otu_table <- taxatable %>% 
  column_to_rownames(var = "species") %>%
  t() #%>%
  as.matrix()

  # OTU table with one row per sample and one column per taxon
ref_taxatable <- ref_taxatable_long %>% 
  dplyr::select(sample, species, count) %>%
  pivot_wider(names_from = species, values_from = count) 

ref_otu_table <- ref_taxatable %>%
  column_to_rownames(var = "sample") %>%
  mutate(across(everything(), replace_na, 0)) #%>%
  as.matrix()

sample_metadata <- metadata %>%
  filter(SourceSink == "sink",
         `#SampleID` != "LIB030.A0117")

# need to combine OTU tables with missing taxa given NA, and then converted to 0
otu_table %>%
  bind_rows(ref_otu_table)

# PCA ---------------------------------------------------------------------

logratio.transfo(logratio = "CLR", offset = 1)

spca

# mixOmics

init_pca <- pca(sample_counts_offset, ncomp = 10, logratio = "CLR")

plot(init_pca, type = "l")

# PCA of samples grouped by day

plotIndiv(init_pca, group = group_env)

sample_pcs <- init_pca$variates$X

# combine metadata with PC coordinate values
sample_pcs_long <- sample_pcs %>%
  as_tibble(rownames = "#SampleID") %>% 
  inner_join(sample_metadata, by = "#SampleID") %>% 
  pivot_longer(cols = PC1:PC10, names_to = "PC", values_to = "coordinate")


# Alpha-diversity ---------------------------------------------------------

# Shannon index for samples

diversity(otu_table)
diversity(ref_otu_table)

alpha_diversity <- diversity(sample_counts) #%>% # calculate Shannon index
  na.omit


# Beta-diversity ----------------------------------------------------------

sample_counts_offset <- sample_counts %>%
    relocate(sample_metadata$`#SampleID`) %>% # relocate columns to match metadata
    as.matrix() %>%
    t() %>%
    + 1
group_day <- as.factor(sample_metadata$day)
group_treatment <- sample_metadata$treatment
group_treatment[is.na(group_treatment)] <- 0
group_plate <- as.factor(sample_metadata$plate)
group_env <- as.factor(sample_metadata$Env)



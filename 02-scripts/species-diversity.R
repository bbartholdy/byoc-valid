library(vegan)
library(mixOmics)
library(tidyverse)

taxatable <- readr::read_csv("03-data/taxatable_decontam.csv")

sample_counts <- taxatable %>% 
  select(!LIB030.A0117) %>%
  column_to_rownames("species")


# Alpha-diversity ---------------------------------------------------------

alpha_diversity <- diversity(sample_counts) %>% # calculate Shannon index
  na.omit

# combine Shannon index values with species counts (for which index could be calculated)



# Beta-diversity ----------------------------------------------------------



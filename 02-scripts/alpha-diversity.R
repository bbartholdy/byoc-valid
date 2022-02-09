library(vegan)
library(tidyverse)

species_counts <- readr::read_csv("03-data/taxatable_decontam.csv")

species_counts %>% 
  select(!LIB030.A0117) %>%
  column_to_rownames("species") %>%
  #diversity() %>% 
  na.omit() %>%
  view



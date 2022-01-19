library(vegan)

species_counts <- readr::read_csv("03-data/taxatable_decontam.csv")

species_counts %>% 
  column_to_rownames("species") %>%
  diversity() %>% 
  view



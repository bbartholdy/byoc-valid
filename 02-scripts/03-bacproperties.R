# Bacterial properties

library(BacDive)
library(tidyverse)


# Upload data -------------------------------------------------------------

metadata <- readr::read_tsv("01-documentation/metadata.tsv")
species_abundance_long <- readr::read_tsv("05-results/species-abundance_long.tsv")
bacdive_oxytol <- readr::read_csv(
  "03-data/2022-05-27_bacdive_oxytol-search.csv",
  skip = 2
  )

# Functions ---------------------------------------------------------------






# access_credentials <- Sys.getenv(c("BACDIVE_USER", "BACDIVE_PW"))
# bacdive_token <- open_bacdive(access_credentials[[1]], access_credentials[[2]])

# get list of bacterial species from all samples

all_species_names <- unique(species_abundance_long$species) %>%
  as_tibble()
write_tsv(all_species_names, "species-list-for-bacdive.txt", col_names = F)

# filter species that are included in this study and the comparative samples
sample_oxytol <- bacdive_oxytol %>%
  mutate("Oxygen tolerance" = case_when(
    `Oxygen tolerance` == "aerobe" | 
      `Oxygen tolerance` == "obligate aerobe" ~ "aerobe",
    `Oxygen tolerance` == "aerotolerant" | 
      `Oxygen tolerance` == "facultative anaerobe" |
      `Oxygen tolerance` == "microaerophile" |
      `Oxygen tolerance` == "facultative aerobe" ~ "facultative anaerobe",
    `Oxygen tolerance` == "anaerobe" |
      `Oxygen tolerance` == "obligate anaerobe" ~ "anaerobe")) %>% 
  filter(
    #str_match(species, all_species_names$species),
    species %in% all_species_names$value, # create less strict matching criterion to incorporate strains
    is_type_strain_header == 1) %>%
  select(species, `Oxygen tolerance`) #%>%
#mutate(species = str_remove(species, "\\[")) %>%
#mutate(species = str_remove(species, "\\]"))

# add abs species column

# list of ABS species from Nikitkova et al. 2013
abs_vector <- c("Streptococcus mitis", "Streptococcus oralis",
                "Streptococcus gordonii", "Streptococcus sanguinis",
                "Streptococcus cristatus", "Streptococcus anginosus",
                "Streptococcus salivarius", "Streptococcus vestibularis",
                "Streptococcus mutans")

all_species_names %>%
  filter(str_detect(value, regex("streptococcus", ignore_case = T))) %>% view()

bac_properties <- sample_oxytol %>%
  mutate(abs = case_when(species %in% abs_vector ~ TRUE,
                         TRUE ~ FALSE))

write_tsv(bac_properties, "01-documentation/species-properties.tsv")


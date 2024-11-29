# Bacterial properties

library(dplyr)
library(tibble)
library(stringr)
library(readr)

# Upload data -------------------------------------------------------------

metadata <- readr::read_tsv("01-documentation/metadata.tsv")
taxatable <- readr::read_tsv("04-analysis/decontam/post-decontam_taxatable.tsv")
bacdive_oxytol <- readr::read_csv(
  "03-data/2022-12-12_bacdive-oxytol-search.csv",
  skip = 2
)


# Bacterial properties ----------------------------------------------------

# get list of bacterial species from all samples
all_species_names <- taxatable %>% 
  mutate(`#OTU ID` = str_remove(`#OTU ID`, "\\]"),
         `#OTU ID` = str_remove(`#OTU ID`, "\\[")) %>%
    .$`#OTU ID`
all_species_names[which(str_detect(all_species_names, coll("]")))]
length(all_species_names)

# combine similar oxygen tolerances
bacdive_oxytol_comb <- bacdive_oxytol %>%
  filter(
    #str_match(species, all_species_names$species),
    species %in% all_species_names,
  ) %>%
  mutate("Oxygen tolerance" = case_when(
    `Oxygen tolerance` == "aerobe" | 
      `Oxygen tolerance` == "obligate aerobe" ~ "aerobe",
    `Oxygen tolerance` == "aerotolerant" | 
      `Oxygen tolerance` == "facultative anaerobe" |
      `Oxygen tolerance` == "microaerophile" |
      `Oxygen tolerance` == "facultative aerobe" ~ "facultative anaerobe",
    `Oxygen tolerance` == "anaerobe" |
      `Oxygen tolerance` == "obligate anaerobe" ~ "anaerobe")) #%>% 

# distribution of oxygen tolerance within genera

genus_oxytol <- bacdive_oxytol_comb %>%
  mutate(genus = str_extract(species, "\\w+")) %>%
  group_by(genus) %>% 
  count(`Oxygen tolerance`, sort = T) %>% # arrange by highest counts
  distinct(genus,.keep_all = T) %>% # keep only the oxytol with highest count for each genus
  dplyr::select(!n)

# filter species that are included in this study and the comparative samples
sample_oxytol <- bacdive_oxytol_comb %>%
  group_by(species) %>%
  arrange(desc(is_type_strain_header), .by_group = T) %>%
  dplyr::select(species, `Oxygen tolerance`) %>%
  distinct(species, .keep_all = T) # retains only first row of duplicates; preference to type strains (not ideal)

  
# add abs species column

# list of ABS species from Nikitkova et al. 2013
abs_vector <- c("Streptococcus mitis", "Streptococcus oralis",
                "Streptococcus gordonii", "Streptococcus sanguinis",
                "Streptococcus cristatus", "Streptococcus anginosus",
                "Streptococcus salivarius", "Streptococcus vestibularis",
                "Streptococcus mutans")

bac_properties <- sample_oxytol %>%
  mutate(abs = case_when(
    species %in% abs_vector ~ TRUE,
    TRUE ~ FALSE
  ))

# produce a species list to use with the bacdive API (future analysis)
write_tsv(as_tibble(all_species_names), "04-analysis/bacdive/species-list-for-bacdive.txt", col_names = F)

write_tsv(bac_properties, "04-analysis/bacdive/species-properties.tsv")
write_tsv(genus_oxytol, "04-analysis/bacdive/genus-O2tolerance.tsv")

# what genera are missing in the bacdive search?
genus_in_bacdive <- bacdive_oxytol %>%
  mutate(genus = str_extract(species, "\\w+")) %>%
  .$genus %>%
  unique()
genus_in_samples <- taxatable %>% 
  mutate(genus = str_extract(`#OTU ID`, "\\w+")) %>%
  .$genus %>%
  unique()
genus_in_samples[!genus_in_samples %in% genus_in_bacdive]

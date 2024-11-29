# Bacterial properties

library(dplyr)
library(tibble)
library(stringi)
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
  mutate(`#OTU ID` = stri_replace(`#OTU ID`, "", regex =  "\\]"),
         `#OTU ID` = stri_replace(`#OTU ID`, "", regex = "\\[")) %>%
    .$`#OTU ID`

# see if previous operation removed all names with the ']' symbol
  # expect: character(0)
all_species_names[which(stri_detect(all_species_names, fixed = "]"))]
# make sure no species names were removed
stopifnot(length(all_species_names) == nrow(taxatable))

# combine similar oxygen tolerances
bacdive_oxytol_comb <- bacdive_oxytol %>%
  filter(
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
  mutate(genus = stri_extract(species, regex = "\\w+")) %>%
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
  mutate(genus = stri_extract(species, regex = "\\w+")) %>%
  .$genus %>%
  unique()
genus_in_samples <- taxatable %>% 
  mutate(genus = stri_extract(`#OTU ID`, regex = "\\w+")) %>%
  .$genus %>%
  unique()
genus_in_samples[!genus_in_samples %in% genus_in_bacdive]

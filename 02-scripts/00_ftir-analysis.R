library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringi)
library(here)
library(patchwork)

# set a consistent theme across plots
options(ggplot2.discrete.colour = function() scale_colour_viridis_d(),
        ggplot2.continuous.colour = function() scale_colour_viridis_c())


# Functions ---------------------------------------------------------------

clean_ftir <- function(x, normalise = FALSE) {
  x_clean <- x %>%
    filter(abs != 0) %>%
    arrange(desc(wavenumber)) # unnecessary - can be reversed in ggplot
  
  if(normalise == TRUE){
    x_norm <- x_clean %>%
      mutate(abs = abs / max(abs) * 100)
    return(x_norm)
  } else {
    return(x_clean) 
  }
}


# Upload data -------------------------------------------------------------

byoc_ftir <- readr::read_csv(here("03-data/FTIR/byoc-ftir.csv")) %>%
  mutate(sample_id = tube)
grind_data_raw <- readr::read_csv(here("03-data/FTIR/grind-curve_archDC.csv"))
ftir_files <- list.files(here("04-analysis/FTIR"), "(?i).CSV", full.names = T)
ftir_data_list <- lapply(ftir_files, read_csv, col_names = c("wavenumber", "abs"))
sample_names <- list.files(here("04-analysis/FTIR/")) %>%
  stri_replace("", regex = "(?i).CSV")
names(ftir_data_list) <- sample_names

# generate metadata for samples

ftir_metadata <- sample_names %>%
  as_tibble() %>%
  mutate(sample = value,
         sample_id = sample) %>% # keep original sample ID
  separate_rows(sample, sep = "\\+") %>% # separate combined samples 
  mutate(
    sample_id = stri_extract(sample, regex = "^[A-Z0-9]+.[A-Z0-9]+|[A-Za-z0-9\\-]+"),
    day = stri_extract(sample_id, regex = "(?<=F)[0-9]+"),
    well = stri_extract(sample_id, regex = "^(?<=A-Z0-9).[A-Z0-9]+"),
    source = case_when(
      stri_detect(sample_id, fixed = "F") ~ "Artificial", 
      stri_detect(sample_id, fixed = "Arch") ~ "Archaeological",
      stri_detect(sample_id, fixed = "modern") ~ "Modern"),
    sample_id = case_when(source == "Modern" ~ sample_id,
                        TRUE ~ sample_id),
    comb = case_when(stri_detect(value, regex = "\\+") ~ value,
                   TRUE ~ NA_character_),
    grind = case_when(stri_detect(value, regex = "(?<=_grind_)[a-f]$") ~ TRUE, # collapse grind samples
                    TRUE ~ FALSE)) %>%
  select(!c(value, sample)) %>%
  mutate(day = as.numeric(day),
         comb = stri_remove(comb, regex = "_grind_[a-z]")) %>% 
  full_join(byoc_ftir, by = c("day", "sample_id")) %>%
  distinct(sample_id, .keep_all = T) %>%
  select(!c(tube, sample, experiment))

# Prepare data ------------------------------------------------------------

grind_sample_order <- c(
  "Archaeological calculus",
  "Artificial calculus day 16",
  "Artificial calculus day 20",
  "Artificial calculus day 24",
  #"Archaeological bone",
  #"Bone-Dentine",
  #"Bone-Dentine_2",
  "Enamel",
  "Enamel_2",
  "Enamel_3")

grind_data <- grind_data_raw %>%
  filter(Sample != "Synthetic") %>% # not sure what 'Synthetic' refers to...
  mutate(day = stri_extract(Sample, regex = "(?<=F)[0-9]+"),
         grind = stri_extract(Sample, regex = "[a-f]$"),
         Sample = case_when(
           #Sample == "Bone +Dentine" ~ "Bone-Dentine",
           #Sample == "Bone +Dentine_2" ~ "Bone-Dentine_2",
           stri_detect(Sample, fixed = "F") ~ "Artificial calculus", 
           stri_detect(Sample, fixed = "MB11") ~ "Archaeological calculus",
           TRUE ~ Sample),
         Sample_day = if_else(!is.na(day), paste(Sample, "day", day), Sample),
         Sample_day = factor(Sample_day, levels = grind_sample_order)
  )

for(i in 1:length(ftir_data_list)){
  sample_name <- names(ftir_data_list)[i]
  ftir_data_list[[i]] <- ftir_data_list[[i]] %>%
    mutate(sample = sample_name)
}

# remove 0.0 absorbance and reverse wavenumber
  # normalise?
ftir_data_list_cleaned <- lapply(ftir_data_list, clean_ftir)

# combine samples with their metadata

ftir_data <- do.call(bind_rows, ftir_data_list_cleaned)

ftir_metadata <- ftir_data %>%
  distinct(sample) %>%
  mutate(
    #analysis_id = sample,
    analysis_id = sample
    ) %>% # keep original sample ID
  separate_rows(sample, sep = "\\+") %>% # separate combined samples
  mutate(
    sample_id = stri_extract(sample, regex = "^[A-Z0-9]+.[A-Z0-9]+|[A-Za-z0-9\\-]+"),
    day = stri_extract(sample_id, regex = "(?<=F)[0-9]+"),
    well = stri_extract(sample_id, regex = "^(?<=A-Z0-9).[A-Z0-9]+"),
    source = case_when(
      stri_detect(sample_id, fixed = "F") ~ "Artificial",
      stri_detect(sample_id, fixed = "Arch") ~ "Archaeological",
      stri_detect(sample_id, fixed = "modern") ~ "Modern"),
    grind = case_when(stri_detect(sample, regex = "(?<=_grind_)[a-f]$") ~ TRUE, # collapse grind samples
                      TRUE ~ FALSE)) %>%
  mutate(
    sample_id = case_when(
      stri_detect(analysis_id, fixed = "modern-ref") ~ analysis_id, # re-separate modern
      TRUE ~ sample_id
      )
  ) %>%
  mutate(
    day = as.numeric(day),
    #comb = str_remove(comb, "_grind_[a-z]")
  ) %>%
  select(!c(sample)) %>%
  inner_join(byoc_ftir, by = c("day", "sample_id")) %>% 
  group_by(sample_id) %>%
  mutate(analysis_id = paste0(analysis_id, collapse = ";")) %>%
  distinct(sample_id, .keep_all = T) %>%
  select(!c(tube, sample, experiment))

# need a better way to combine FTIR data with metadata
  # currently doesn't work because grind extensions (_a:f) are removed from metadata

write_csv(grind_data, "05-results/grind-data_cleaned.csv")
write_csv(ftir_data, "05-results/ftir-data.csv")
#write_csv(ftir_data_long, "05-results/ftir-data_long.csv")
write_tsv(ftir_metadata, "01-documentation/ftir-metadata.tsv")

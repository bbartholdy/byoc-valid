library(readr)
library(dplyr)
library(stringr)
library(here)
source("02-scripts/functions.R")

ftir_sample_metadata <- read_csv(here("01-documentation/ftir_model-metadata.csv"))
grind_data_raw <- read_csv(here("03-data/FTIR/ftir_calculus-grind-curve.csv"))
ftir_files <- list.files(here("03-data/FTIR"), "(?i).CSV", full.names = T)[-87] # remove grind curve calculations
ftir_data_list <- lapply(ftir_files, read_csv, col_names = c("wavenumber", "abs"))
analysis_id <- list.files(here("03-data/FTIR/"), pattern = "*.CSV") %>%
  stri_replace("", regex = "(?i).CSV", opts_regex = list(case_insensitive = T))
analysis_id <- c(analysis_id, "modern-ref_1", "modern-ref_2")
names(ftir_data_list) <- analysis_id

# fix typos from file names
analysis_id[analysis_id == "F12.1A5+F12.B1_B"] <- "F12.1A5+F12.1B1_B"
analysis_id[analysis_id == "F12.1A5+F12.B1"] <- "F12.1A5+F12.1B1"

analysis_id[analysis_id == "F7"] <- NA # unknown context
analysis_id <- analysis_id[!is.na(analysis_id)]

# generate additional metadata for samples --------------------------------

ftir_analysis_metadata <- analysis_id %>%
  as_tibble() %>%
  mutate(analysis = value) %>% # keep original sample ID
  separate_rows(value, sep = "\\+") %>% # separate combined samples 
  
  mutate(sample_id = stri_extract(value, regex = "^[A-Z0-9]+.[A-Z0-9]+|[A-Za-z0-9\\-]+")) |>
  group_by(sample_id) |>
  mutate(
    analysis_id = paste0(analysis, collapse = ";"),
  ) |>
  ungroup() |>
  distinct(sample_id, .keep_all = T) |>
  select(!c(value, analysis))

ftir_metadata <- ftir_analysis_metadata |>
  full_join(ftir_sample_metadata, by = c("sample_id" = "tube")) |>
  mutate(source = case_when(
      str_detect(sample_id, "F") ~ "Artificial", 
      str_detect(sample_id, "Arch") ~ "Archaeological",
      str_detect(sample_id, "modern") ~ "Modern"),
      sample = "dental_calculus"
  )

write_tsv(ftir_metadata, "01-documentation/ftir-metadata.tsv")


# Combine spectra into single file ---------------------------------------

for(i in 1:length(ftir_data_list)){
  sample_name <- names(ftir_data_list)[i]
  ftir_data_list[[i]] <- ftir_data_list[[i]] %>%
    mutate(analysis_id = sample_name)
}

# remove 0.0 absorbance and reverse wavenumber
  # normalise?
ftir_data_list_cleaned <- lapply(ftir_data_list, clean_ftir)

# combine samples with their metadata

ftir_data <- do.call(bind_rows, ftir_data_list_cleaned)

write_csv(ftir_data, "04-analysis/FTIR/ftir_full-data.csv")
write_csv(ftir_data, "05-results/ftir_full-data.csv")


# Prepare grind data -----------------------------------------------------

grind_sample_order <- c(
  "Archaeological calculus",
  "Artificial calculus day 16",
  "Artificial calculus day 20",
  "Artificial calculus day 24",
  "Enamel",
  "Enamel_2",
  "Enamel_3")

grind_data <- grind_data_raw %>%
  filter(Sample != "Synthetic") %>% # not sure what 'Synthetic' refers to...
  mutate(day = stri_extract(Sample, regex = "(?<=F)[0-9]+"),
         grind = stri_extract(Sample, regex = "[a-f]$"),
         Sample = case_when(
           stri_detect(Sample, fixed = "F") ~ "Artificial calculus", 
           stri_detect(Sample, fixed = "MB11") ~ "Archaeological calculus",
           TRUE ~ Sample),
         Sample_day = if_else(!is.na(day), paste(Sample, "day", day), Sample),
         Sample_day = factor(Sample_day, levels = grind_sample_order)
  )

write_csv(grind_data, "04-analysis/FTIR/ftir_grind-data.csv")
write_csv(grind_data, "05-results/ftir_grind-data.csv")

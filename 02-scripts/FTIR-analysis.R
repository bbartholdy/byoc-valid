library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
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

# function to isolate locate peaks (maxima) and label them on plots

test_spectrum <- readr::read_csv("04-analysis/FTIR/ArchDC_MB11_grind_c.CSV", col_names = c("wavenumber", "abs"))
x <- clean_ftir(test_spectrum)
x <- filter(x, wavenumber < 700 & wavenumber >500)
# need to ignore region 2000-1700

#' @param n integer. If given, only returns the highest n peaks.
#' @param abs_cutoff numeric. peak height cutoff. A value is given as the minimum
#' height from a peak to the adjacent valley, for the peak to be considered as a true peak.
#' @param range vector or list of vectors containing two values to indicate a range to identify peaks. 
peak_finder <- function(x, n, abs_cutoff = 0.1, range){
  # could pmax function be useful?
  # useful ranges:
  range1 <- c(3800,2800)
  range2 <- c(1700,1300)
  range3 <- c(1300,1100)
  range4 <- c(1100, 900)
  range5 <- c(900, 700)
  range6 <- c(700,500)
  ignore_range <- c(2000,1700)
  
  abs <- x$abs
  # abs
  # diff(sign(diff(round(abs,3))), differences = 2)
  # sign(diff(abs, differences = 1))
  # peaks <- which(diff(sign(diff(abs))) == -2)
  # valleys <- which(diff(sign(diff(abs))) == 2)
  # which(diff(sign(diff(abs, differences = 1)), differences = 4) == -6)
  
  abs_signs <- c(0, sign(diff(abs))) # vector of signs associated with each abs value
  peaks <- which(diff(sign(diff(abs))) == -2)
  valleys <- which(diff(sign(diff(abs))) == 2)
    
  x$signs <- abs_signs
  x$peak <- x$abs
  x$peak[peaks] <- TRUE
  x$peak <- ifelse(x$peak == TRUE, TRUE, FALSE)
  x$valley <- x$abs
  x$valley[valleys] <- TRUE
  x$valley <- ifelse(x$valley == TRUE, TRUE, FALSE)
  x$wavenumber[peaks]
  peak_abs <- x$abs[peaks]
  
  peak_define <- x %>%
    filter(peaks == T | valleys == T) %>%
    mutate(
      diff = c(
        0, sapply(2:nrow(.), 
                  function(x) .$abs[x] + .$abs[x-1]
               )
        )
      ) %>%
    filter(diff > abs_cutoff)
  
  if (!is.null(n)){
    n_peaks <- peak_define %>%
      arrange(desc(abs)) %>% 
      slice_head(n)
  }
  
  return(peak_define)
}

# Upload data -------------------------------------------------------------

byoc_ftir <- readr::read_csv(here("03-data/FTIR/byoc-ftir.csv")) %>%
  mutate(sample_id = tube)
grind_data_raw <- readr::read_csv(here("03-data/FTIR/grind-curve_archDC.csv"))
ftir_files <- list.files(here("04-analysis/FTIR"), "(?i).CSV", full.names = T)
ftir_data_list <- lapply(ftir_files, read_csv, col_names = c("wavenumber", "abs"))
sample_names <- list.files(here("04-analysis/FTIR/")) %>%
  str_remove("(?i).CSV")
names(ftir_data_list) <- sample_names

# generate metadata for samples

ftir_metadata <- sample_names %>%
  as_tibble() %>%
  mutate(sample = value,
         analysis_id = sample) %>% # keep original sample ID
  separate_rows(sample, sep = "\\+") %>% 
  mutate(
    sample_id = str_extract(sample, "^[A-Z0-9]+.[A-Z0-9]+|[A-Za-z0-9\\-]+"),
    day = str_extract(sample_id, "(?<=F)[0-9]+"),
    well = str_extract(sample_id, "^(?<=A-Z0-9).[A-Z0-9]+"),
    source = case_when(
      str_detect(sample_id, "F") ~ "Artificial", 
      str_detect(sample_id, "Arch") ~ "Archaeological",
      str_detect(sample_id, "modern") ~ "Modern"),
    sample_id = case_when(source == "Modern" ~ sample,
                        TRUE ~ sample_id),
    comb = case_when(str_detect(value, "\\+") ~ value,
                   TRUE ~ NA_character_),
    grind = case_when(str_detect(value, "(?<=_grind_)[a-f]$") ~ TRUE, # collapse grind samples
                    TRUE ~ FALSE)) %>%
  select(!c(value, sample)) %>%
  mutate(day = as.numeric(day),
         comb = str_remove(comb, "_grind_[a-z]")) %>% 
  full_join(byoc_ftir, by = c("day", "sample_id")) %>%
  distinct(sample_id, .keep_all = T) %>%
  select(!c(tube, sample, experiment))

# Prepare data ------------------------------------------------------------

grind_sample_order <- c(
  "Archaeological calculus",
  "Artificial calculus day 16",
  "Artificial calculus day 20",
  "Artificial calculus day 24",
  "Archaeological bone",
  "Bone-Dentine",
  "Bone-Dentine_2",
  "Enamel",
  "Enamel_2",
  "Enamel_3")

grind_data <- grind_data_raw %>%
  filter(Sample != "Synthetic") %>%
  mutate(day = stringr::str_extract(Sample, "(?<=F)[0-9]+"),
         grind = stringr::str_extract(Sample, "[a-f]$"),
         Sample = case_when(
           Sample == "Bone +Dentine" ~ "Bone-Dentine",
           Sample == "Bone +Dentine_2" ~ "Bone-Dentine_2",
           str_detect(Sample, "F") ~ "Artificial calculus", 
           str_detect(Sample, "MB11") ~ "Archaeological calculus",
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

ftir_data_long <- inner_join(ftir_metadata, ftir_data, by = c("analysis_id" = "sample"))

write_csv(grind_data, "05-results/grind-data_cleaned.csv")
write_csv(ftir_data_long, "05-results/ftir-data_long.csv")
write_tsv(ftir_metadata, "01-documentation/ftir-metadata.tsv")

library(tidyverse)
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

grind_data_raw <- readr::read_csv(here("03-data/FTIR/grind-curve_archDC.csv"))
ftir_files <- list.files(here("04-analysis/FTIR"), "(?i).CSV", full.names = T)
ftir_data_list <- lapply(ftir_files, read_csv, col_names = c("wavenumber", "abs"))
sample_names <- list.files(here("04-analysis/FTIR/")) %>%
  str_remove("(?i).CSV")
names(ftir_data_list) <- sample_names

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

# generate metadata for samples

ftir_metadata <- sample_names %>%
  as_tibble() %>%
  rename(sample = value) %>%
  mutate(
    sample_name = str_extract(sample, "^[A-Z0-9]+.[A-Z0-9]+|[A-Za-z0-9]+"),
    day = str_extract(sample_name, "(?<=F)[0-9]+"),
    grind = str_extract(sample, "(?<=_grind_)[a-f]$"),
    well = str_extract(sample_name, "^(?<=A-Z0-9).[A-Z0-9]+"),
    source = case_when(
      str_detect(sample_name, "F") ~ "Artificial", 
      str_detect(sample_name, "Arch") ~ "Archaeological",
      str_detect(sample_name, "modern") ~ "Modern"))

# combine samples with their metadata

ftir_data <- do.call(bind_rows, ftir_data_list_cleaned)

ftir_data_long <- inner_join(ftir_metadata, ftir_data, by = "sample")

write_csv(grind_data, "05-results/grind-data_cleaned.csv")
write_csv(ftir_data_long, "05-results/ftir-data_long.csv")
write_tsv(ftir_metadata, "01-documentation/ftir-metadata.tsv")

# Visualise spectra -------------------------------------------------------

# byoc calculus spectra

ftir_data_long %>% 
  filter(day == 24, # only byoc calculus
         is.na(grind)) %>% # remove grind samples
  ggplot(aes(x = wavenumber, y = abs, col = sample)) +
    geom_line()

# change over course of experiment

ftir_data_long %>%
  filter(source == "Artificial",
         is.na(grind)) %>% # remove grind samples
  ggplot(aes(x = wavenumber, y = abs))

# compare byoc to arch calculus

ftir_data_long %>%
  filter(str_detect(sample, "(?<=MB11)_grind_[dc]") |
           sample == "F24.2A4" | sample == "F24.1D3") %>%
  ggplot(aes(x = wavenumber, y = abs)) +
    #geom_line(aes(col = source), alpha = 0.4) +
    geom_line(aes(col = sample)) +
    theme_classic() +
    scale_x_reverse()

# range of spectra in arch calculus grind curve

ftir_data_long %>%
  filter(sample_name == "ArchDC") %>% 
  ggplot(aes(x = wavenumber, y = abs)) +
    geom_line(alpha = 0.4, col = "grey") +
    geom_line(aes(col = sample)) +
    theme_classic()

# Plots of grind curves (MOVED TO FTIR-analysis.Rmd)

grind_all_plot <- grind_data %>%
  #group_by(day, Sample) %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = F) +
  labs(col = "Sample", shape = "Sample",
       x = "FWHM of the 1035 peak",
       y = "Splitting factor")

# isolate calculus samples to see diffs between days and the 'real deal'
grind_calc_plot <- grind_data %>%
  filter(Sample == "Artificial calculus" | Sample == "Archaeological calculus") %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = F) +
  labs(col = "Sample", shape = "Sample",
       x = "FWHM of the 1035 peak",
       y = "Splitting factor") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_colour_viridis_d(end = 0.4)
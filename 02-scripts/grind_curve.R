library(tidyverse)
library(here)

options(ggplot2.discrete.colour = function() scale_colour_viridis_d(),
        ggplot2.continuous.colour = function() scale_colour_viridis_c())

# set a consistent theme across plots

grind_data_raw <- readr::read_csv(here("analyses/FTIR/grind-curve_archDC.csv"))

grind_data <- grind_data_raw %>%
  filter(Sample != "Synthetic") %>%
  mutate(day = stringr::str_extract(Sample, "(?<=F)[0-9]+"),
         grind = stringr::str_extract(Sample, "[a-f]$"),
         Sample = case_when(
           str_detect(Sample, "F") ~ "Artificial calculus", 
           str_detect(Sample, "MB11") ~ "Archaeological calculus",
           TRUE ~ Sample),
         Sample_day = if_else(!is.na(day), paste(Sample, "day", day), Sample
                              )
         )

# Plot of grind curves
grind_data %>%
  group_by(day, Sample) %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  labs(col = "Sample", shape = "Sample",
       x = "FWHM of the 1035 peak",
       y = "Splitting factor")

# isolate calculus samples to see diffs between days and the 'real deal'
grind_data %>%
  filter(Sample == "Artificial calculus" | Sample == "Archaeological calculus") %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  labs(col = "Sample", shape = "Sample",
       x = "FWHM of the 1035 peak",
       y = "Splitting factor")

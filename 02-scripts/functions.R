# Easy access to re-used objects
  # especially long-format data frames that are to big to export

library(tidyr)
library(dplyr)

# converts SourceTracker2 output to long format
sourcetracker2_longer <- function(){
  sourcetracker2 %>%
  pivot_longer(cols = where(is.numeric), 
    values_to = "proportion", 
    names_to = "SampleID") %>%
  rename(source = ...1)
}

ftir_spect_plot <- function(.data, sample_name){
  .data %>%
    dplyr::filter(
      analysis_id == {{ sample_name }}
    ) %>%
    ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
      geom_line() +
      theme_classic() +
      scale_x_reverse()
}

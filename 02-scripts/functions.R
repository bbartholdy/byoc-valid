# Easy access to re-used objects
  # especially long-format data frames that are to big to export

library(tidyr)
library(dplyr)

# converts SourceTracker2 output to long format
sourcetracker2_longer <- function(data){
  data %>%
    pivot_longer(
      cols = where(is.numeric), 
      values_to = "proportion", 
      names_to = "SampleID"
    ) %>%
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

#' function to calculate bias-corrected log-observed abundances
#' from `vignette("ANCOMBC")`
#' 
#' @param x `ancombc` object.
bias_correct <- function(x, otu_table) {
  require(microbiome)
  samp_frac <- x$samp_frac
  samp_frac[is.na(samp_frac)] <- 0 # replace NAs with 0
  log_otu <- log(microbiome::abundances(otu_table) + 1)
  log_otu_adj <- exp(t(t(log_otu) - samp_frac)) %>% # bias corrected abund
    as_tibble(rownames = "species")
  return(log_otu_adj)
}

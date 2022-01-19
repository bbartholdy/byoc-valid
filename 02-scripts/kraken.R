# kraken
library(tidyverse)
# Upload all kraken reports

file_names <- list.files("03-data/kraken/", "_report")
file_paths <- paste0("03-data/kraken/", file_names)
sample_names <- gsub(".unmapped.*", "", file_names)

kraken_reports <- vector(mode = "list", length = length(file_paths))
kraken_reports <- lapply(file_paths, readr::read_tsv, col_names = F)
names(kraken_reports) <- file_names

# extract columns 2 and 8 for all cases where col 6 has "S"
kraken_out <- lapply(kraken_reports, function(x) subset(x, x$X6 == "S")[,c(2,8)])

# merge all *_report files
  # add Sample name column to each data frame, then merge all data frames

for(i in 1:length(kraken_out)){
  kraken_out[[i]]$sample <- sample_names[i]
}

kraken_otu <- kraken_out %>%
  reduce(full_join, by = c("X2", "X8", "sample")) %>%
  rename(count = X2,
         species = X8)

# Species in the rows and samples in the column

kraken_seqtab <- kraken_otu %>%
  pivot_wider(names_from = species, values_from = count)

kraken_taxatable <- kraken_otu %>%
  pivot_wider(names_from = sample, values_from = count)

write_csv(kraken_otu, "03-data/species_counts.csv")
write_csv(kraken_seqtab, "03-data/sequence_table.csv")
write_csv(kraken_taxatable, "03-data/taxa_table.csv")

# quick search for ABS
kraken_otu %>%
  filter(sample != "LIB030.A0117") %>%
  filter(grepl("Rothia", species))

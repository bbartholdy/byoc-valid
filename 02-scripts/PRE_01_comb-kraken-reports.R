# Script to create species tables from Kraken reports
library(dplyr)
library(purrr)
library(here)

# Upload data -------------------------------------------------------------

metadata <- readr::read_tsv("01-documentation/dna-metadata.tsv")

# Upload all kraken reports

file_names <- list.files("04-analysis/kraken/", "_report")
file_paths <- paste0("04-analysis/kraken/", file_names)
sample_names <- gsub(".unmapped.*", "", file_names)

kraken_reports <- vector(mode = "list", length = length(file_paths))
kraken_reports <- lapply(file_paths, readr::read_tsv, col_names = F)
names(kraken_reports) <- file_names

# extract columns 2 and 8 for all cases where col 6 has "S"
kraken_data <- lapply(kraken_reports, function(x) subset(x, x$X6 == "S")[,c(2,8)])

# merge all *_report files
# add Sample name column to each data frame, then merge all data frames

for(i in 1:length(kraken_data)){
  kraken_data[[i]]$sample <- sample_names[i]
}

kraken_otu_long <- kraken_data %>%
  reduce(full_join, by = c("X2", "X8", "sample")) %>%
  rename(count = X2,
         species = X8)

readr::write_csv(kraken_otu_long, here("03-data/kraken-OTU_long.csv"))
# Script to create species tables from Kraken reports
library(tidyverse)

sample_metadata <- readr::read_csv(here("03-data/sample_metadata.csv"))
source_samples <- readr::read_csv(here("03-data/source_samples_key.csv"))

# Upload all kraken reports

file_names <- list.files("04-analysis/kraken/", "_report")
file_paths <- paste0("04-analysis/kraken/", file_names)
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

kraken_otu_filtered <- kraken_otu #%>%
  #mutate(count = replace_na(count, 0)) %>%
  #filter(count > 1000)

# Species in the rows and samples in the column

kraken_seqtab <- kraken_otu %>%
  pivot_wider(names_from = species, values_from = count)

kraken_otufilter_table <- kraken_otu_filtered %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  mutate(across(where(is.numeric), replace_na, 0)) #%>%
  #rename("OTU ID" = species)

# quick search for ABS
kraken_otu %>%
  filter(sample != "LIB030.A0117") %>%
  filter(grepl("Rothia", species))

# Prepare metadata --------------------------------------------------------

sample_metadata$sample[-1] <- paste0(sample_metadata$sample[-1], "0101") # match metadata names to sequence names

colnames(source_samples) <- c("sample", "source") # match column names to sample_metadata

source_samples <- source_samples %>%
  mutate(SourceSink = "source")

metadata <- sample_metadata %>% # combine source metadata with samples
  mutate(SourceSink = "sink") %>%
  bind_rows(source_samples) %>%
  rename(Env = source,
         "#SampleID" = sample) %>%
  filter(`#SampleID` %in% sample_names) # subset successful sequences

#metadata <- subset(metadata, metadata$sample %in% sample_names) # subset successful sequences (?)

# mapping for SourceTracker

st_map <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>% 
  write_tsv(here("04-analysis/ST-map.txt"))

write_csv(kraken_otu, "04-analysis/species_counts.csv")
write_csv(kraken_seqtab, "04-analysis/sequence_table.csv")
write_tsv(kraken_otufilter_table, "04-analysis/OTUfilter_table.tsv")

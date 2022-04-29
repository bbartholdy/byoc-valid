library(tidyverse)
library(here)

# Upload data -------------------------------------------------------------

#sample_metadata <- readr::read_csv(here("03-data/sample_metadata.csv"))
kraken_otu <- readr::read_csv(here("03-data/kraken-OTU_long.csv"))
metadata <- readr::read_tsv("03-data/metadata.tsv")
#source_samples <- readr::read_csv(here("03-data/source_samples_key.csv"))

file_names <- list.files("04-analysis/kraken/", "_report")
#file_paths <- paste0("04-analysis/kraken/", file_names)
sample_names <- gsub(".unmapped.*", "", file_names)

# Filter OTU table --------------------------------------------------------

# filter species by relative abundance
  # species with lower than 0.001% relative abundance removed
kraken_otu_filtered <- kraken_otu %>%
  mutate(count = replace_na(count, 0)) %>%
  group_by(sample) %>%
  mutate(rel_abund = count / sum(count),
         sum_abund = sum(rel_abund)) %>%  # check that each sample == 1
  filter(rel_abund >= 0.00001) 

# Species in the rows and samples in the column

kraken_otufilter_table <- kraken_otu_filtered %>%
  select(!c(rel_abund, sum_abund)) %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  mutate(across(where(is.numeric), replace_na, 0)) #%>%
#rename("OTU ID" = species)

#write_csv(kraken_otu, here("04-analysis/species_counts.csv"))
#write_csv(kraken_seqtab, here("04-analysis/sequence_table.csv"))
write_tsv(kraken_otufilter_table, here("04-analysis/OTUfilter_table.tsv"))


# Prepare metadata --------------------------------------------------------

# sample_metadata$sample[-1] <- paste0(sample_metadata$sample[-1], "0101") # match metadata names to sequence names
# 
# colnames(source_samples) <- c("sample", "source") # match column names to sample_metadata
# 
# source_metadata <- source_samples %>%
#   mutate(SourceSink = "source")
# 
# metadata <- sample_metadata %>% # combine source metadata with samples
#   mutate(SourceSink = "sink") %>%
#   bind_rows(source_metadata) %>%
#   rename(Env = source,
#          "#SampleID" = sample) %>%
#   filter(`#SampleID` %in% sample_names) # subset successful sequences

#metadata <- subset(metadata, metadata$sample %in% sample_names) # subset successful sequences (?)

# mapping for SourceTracker

st_map <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  write_tsv(here("04-analysis/sourcetracker/ST-map.txt"))

st_map_plaque_comb <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  mutate(Env = if_else(Env == "supragingival_plaque" | Env == "subgingival_plaque",
                        "plaque", Env)) %>%
  write_tsv(here("04-analysis/sourcetracker/ST_comb-plaque-map.txt"))

# st_noplaque_map <- st_map %>%
#   filter(Env != "supragingival_plaque",
#          Env != "subgingival_plaque") %>%
#   write_tsv(here("04-analysis/ST-noplaque-map.txt"))
# 
# kraken_otu_noplaque <- kraken_otu %>%
#   filter(sample %in% st_noplaque_map$`#SampleID`,
#          count > 1000) %>%
#   pivot_wider(names_from = sample, values_from = count) %>%
#   mutate(across(where(is.numeric), replace_na, 0))

#write_tsv(kraken_otu_noplaque, here("04-analysis/OTUnoplaque_table.tsv"))
#write_tsv(metadata, here("03-data/metadata.tsv"))

# Script to create species tables from Kraken reports
library(tidyverse)
library(here)

# Upload data -------------------------------------------------------------

metadata <- readr::read_tsv("01-documentation/metadata.tsv")

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

#write_csv(kraken_otu_long, here("03-data/kraken-OTU_long.csv"))


file_names <- list.files("04-analysis/kraken/", "_report")
#file_paths <- paste0("04-analysis/kraken/", file_names)
sample_names <- gsub(".unmapped.*", "", file_names)

# Filter OTU table --------------------------------------------------------

# filter species by relative abundance
  # species with lower than 0.001% relative abundance removed
kraken_otu_filtered <- kraken_otu_long %>%
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

# names of BYOC samples
which_samples <- metadata %>%
  filter(SourceSink == "sink") %>%
  .$`#SampleID`

# isolate the BYOC samples
sample_taxatable <- kraken_otu_filtered %>% 
  select(!c(rel_abund, sum_abund)) %>% 
  filter(sample %in% which_samples) %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  mutate(across(where(is.numeric), replace_na, 0))

write_tsv(sample_taxatable, here("04-analysis/pre-decontam_sample_taxatable.tsv"))
write_tsv(kraken_otufilter_table, here("04-analysis/OTUfilter_table.tsv"))

# mapping for SourceTracker

st_map <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  write_tsv(here("04-analysis/sourcetracker/ST-map.txt"))

st_map_plaque_comb <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  mutate(Env = if_else(Env == "supragingival_plaque" | Env == "subgingival_plaque",
                        "plaque", Env)) %>%
  write_tsv(here("04-analysis/sourcetracker/ST_comb-plaque-map.txt"))

# map with combined plaque sources and no sediments
st_map_nosedi <- st_map_plaque_comb %>% 
  filter(Env != "sediment") %>% 
  write_tsv(here("04-analysis/sourcetracker/ST_no-sedi.txt"))


library(here)
library(dplyr)
library(tidyr)

# Data upload -------------------------------------------------------------

kraken_otu_long <- readr::read_csv("04-analysis/kraken/kraken-OTU_long.csv")
metadata <- readr::read_tsv("01-documentation/dna-metadata.tsv")

# isolate library

lib_sample <- kraken_otu_long %>%
  filter(sample == "LIB030.A0117")

# Filter OTU table --------------------------------------------------------

# filter species by relative abundance
  # species with lower than 0.001% relative abundance removed
kraken_otu_filtered <- kraken_otu_long %>%
  filter(sample != "LIB030.A0117") %>% 
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

# mapping for SourceTracker

st_map <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  write_tsv(here("04-analysis/sourcetracker/ST-map.txt"))

st_map_plaque_comb <- metadata %>%
  select(`#SampleID`, Env, SourceSink) %>%
  mutate(Env = if_else(Env == "supragingival_plaque" | Env == "subgingival_plaque",
                        "plaque", Env)) %>%
  filter(Env != "buccal_mucosa",
         Env != "vitro_biofilm") %>% 
  write_tsv(here("04-analysis/sourcetracker/ST_comb-plaque-map.txt"))

# map with combined plaque sources and no sediments
st_map_nosedi <- st_map_plaque_comb %>% 
  filter(Env != "sediment",
         Env != "buccal_mucosa") %>% 
  write_tsv(here("04-analysis/sourcetracker/ST_no-sedi.txt"))

write_tsv(lib_sample, here("04-analysis/lib_sample.tsv"))
write_tsv(kraken_otufilter_table, here("04-analysis/OTUfilter_table.tsv"))

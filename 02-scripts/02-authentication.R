# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)
library(tidyverse)
library(here)

# upload data
lib_conc <- readr::read_tsv(here("01-documentation/SYN_library_quant.tsv")) # library concentrations
otu_filtered_table <- readr::read_tsv(here("04-analysis/pre-decontam_sample_taxatable.tsv"))
sourcetracker2 <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions.txt"))
sourcetracker2_stdevs <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions_stds.txt"))
sourcetracker2.2 <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker2_output2/mixing_proportions.txt"))
sourcetracker2.2_stdevs <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker2_output2/mixing_proportions_stds.txt"))
sample_taxatable <- readr::read_tsv(here("04-analysis/pre-decontam_sample_taxatable.tsv"))

metadata <- readr::read_tsv(here("01-documentation/metadata.tsv"))
file_names <- list.files(here("04-analysis/kraken/"), "_report")
sample_names <- gsub(".unmapped.*", "", file_names)


# SourceTracker -----------------------------------------------------------

# sourcetracker with plaque source samples included

sourcetracker2_long <- sourcetracker2 %>%
  pivot_longer(cols = where(is.numeric), 
               values_to = "proportion", 
               names_to = "SampleID") %>%
  rename(source = ...1)

sourcetracker2.2_long <- sourcetracker2.2 %>%
  pivot_longer(cols = where(is.numeric), 
               values_to = "proportion", 
               names_to = "SampleID") %>%
  rename(source = ...1)


# Arrange samples by age
day_order <- metadata %>%
  filter(SourceSink == "sink") %>%
  arrange(day) %>%
  filter(`#SampleID` %in% colnames(sourcetracker2.2)) %>%
  .$`#SampleID`

# day_order <- metadata %>%
#   filter(SourceSink == "sink") %>%
#   arrange(day) %>%
#   .$`#SampleID`

# day_order <- day_order[day_order %in% sourcetracker2$SampleID]

sourcetracker2.2_long %>%
  #mutate(SampleID = fct_reorder(SampleID, day_order)) %>%
  ggplot(aes(x = SampleID, proportion, fill = source)) +
    geom_col() +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(limits = day_order)

#ggsave(here("sourcetracker_plot3.png"), width = 10, height = 7, units = "in")

# sourcetracker with plaque source samples removed


# contributions in problematic samples

problem_samples <- sourcetracker2_long %>%
  filter(source == "Unknown",
         proportion > 0.5) %>%
  .$SampleID

#file_names_problems <- paste0("04-analysis/sourcetracker/sourcetracker2_output/", problem_samples, ".feature_table.txt")
file_names_problems <- paste0("04-analysis/sourcetracker/sourcetracker2_output/", problem_samples, ".feature_table.txt")
problem_data <- lapply(as.list(file_names_problems), read_tsv)

file_names_all <- paste0(
  "04-analysis/sourcetracker/sourcetracker2_output/", 
  unique(sourcetracker2_long$SampleID), 
  ".feature_table.txt")
all_data <- lapply(as.list(file_names_all), read_tsv)

names(problem_data) <- problem_samples
names(all_data) <- unique(sourcetracker2_long$SampleID)

problem_list_long <- lapply(problem_data, function(x) x %>% 
                              pivot_longer(cols = where(is.numeric),
                                           values_to = "count",
                                           names_to = "taxon") %>%
                              rename(source = ...1) %>%
                              filter(source == "Unknown")
       )

problem_data_long <- problem_list_long %>% 
  map_df(~ as_tibble(.x), .id = "SampleID")

all_list_long <- lapply(all_data, function(x) x %>% 
                              pivot_longer(cols = where(is.numeric),
                                           values_to = "count",
                                           names_to = "taxon") %>%
                              rename(source = ...1) %>%
                              filter(source == "Unknown")
)

all_data_long <- all_list_long %>% 
  map_df(~ as_tibble(.x), .id = "SampleID")


# How many unknowns are actually oral taxa

iso_database <- load_database(cuperdec_database_ex, target = "oral")
oral_taxa <- iso_database %>%
  filter(Isolation_Source == TRUE)

#ggsave(here("unknowns_plot1.png"), width = 10, height = 7, units = "in")

all_data_long %>% 
  filter(count > 0) %>% 
  mutate(oral_source = if_else(taxon %in% oral_taxa$Taxon, "oral", "other")) %>%
  #group_by(SampleID) %>%
  #mutate(proportion = ) # proportion of taxon count in sample
  ggplot(aes(x = SampleID, y = count, fill = oral_source)) +
  geom_col(position = "fill") +
  theme(axis.text.x = element_text(angle = 90))

# cuperdec ----------------------------------------------------------------

taxa_table <- load_taxa_table(sample_taxatable)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(metadata,
                           sample_col = "#SampleID",
                           source_col = "Env"
                           )

curves <- calculate_curve(taxa_table, iso_database)
filter_result <- simple_filter(curves, 60)
plot_cuperdec(curves, metadata_table, filter_result)


# decontam ----------------------------------------------------------------

lib_conc$`Full Library Id` %in% names(otu_filtered_table)

#kraken_seqtab[is.na(kraken_seqtab)] <- 0 # convert NAs to 0

filtered_seqtab <- sample_taxatable %>%
  column_to_rownames(var = "species") %>%
  as.matrix() %>%
  t()

# test for contaminants
which_contaminants <- isContaminant(filtered_seqtab, 
                              conc = lib_conc$`Quantification post-Indexing total`)

contaminants <- filter(which_contaminants, contaminant == TRUE)

# filter out contaminant species
  # make sure no oral taxa are on the list
contaminant_species <- row.names(contaminants)

true_contaminants <- contaminants %>%
  filter(!contaminant_species %in% oral_taxa$Taxon)

taxatable_decontam <- sample_taxatable %>%
  filter(!species %in% true_contaminants$contaminant)

# how many contaminants per sample?
sample_taxatable %>%
  select(!c(LIB030.A0117, SYN001.A0101, SYN002.A0101, SYN003.A0101)) %>%
  pivot_longer(cols = where(is.numeric), names_to = "SampleID", values_to = "count") %>%
  filter(species %in% contaminant_species) %>%
  #mutate(SampleID = fct_reorder(SampleID, day_order)) %>%
  group_by(SampleID) %>%
  summarise(contam = sum(count)) %>%
  ggplot(aes(x = SampleID, y = contam)) +
    geom_col()

write_tsv(taxatable_decontam, here("05-results/taxatable_decontam.tsv"))
write_tsv(contaminants, here("04-analysis/list-of-contaminants.tsv"))

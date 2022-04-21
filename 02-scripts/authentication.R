# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)
library(tidyverse)
library(here)

# upload data
lib_conc <- readr::read_tsv(here("03-data/SYN_library_quant.tsv")) # library concentrations
#kraken_seqtab <- readr::read_csv(here("04-analysis/sequence_table.csv"))
kraken_taxatab <- readr::read_csv(here("04-analysis/OTUfilter_table.tsv"))
sourcetracker <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker_output/mixing_proportions.txt"))
sourcetracker_stdevs <- readr::read_tsv(here("04-analysis/sourcetracker/sourcetracker_output/mixing_proportions_stds.txt"))

metadata <- readr::read_csv(here("03-data/sample_metadata.csv"))
file_names <- list.files(here("04-analysis/kraken/"), "_report")
sample_names <- gsub(".unmapped.*", "", file_names)


# SourceTracker -----------------------------------------------------------

# sourcetracker with plaque source samples included
sourcetracker_long <- sourcetracker %>%
  pivot_longer(cols = where(is.numeric), 
               values_to = "proportion", 
               names_to = "source")

# Arrange samples by age
day_order <- metadata %>%
  filter(SourceSink == "sink") %>%
  arrange(day) %>%
  .$`#SampleID`

day_order <- day_order[day_order %in% sourcetracker$SampleID]
#day_order <- rep(day_order, each = 6)

sourcetracker_long %>% 
  #mutate(SampleID = fct_reorder(SampleID, day_order)) %>%
  ggplot(aes(x = SampleID, proportion, fill = source)) +
    geom_col() +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(limits = day_order)

ggsave(here("sourcetracker_plot3.png"), width = 10, height = 7, units = "in")

# sourcetracker with plaque source samples removed
sourcetracker2_long <- sourcetracker2 %>%
  pivot_longer(cols = where(is.numeric), 
               values_to = "proportion", 
               names_to = "source")

sourcetracker2_long %>% 
  ggplot(aes(SampleID, proportion, fill = source)) +
    geom_col() +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(limits = day_order)

# cuperdec ----------------------------------------------------------------

taxa_table <- load_taxa_table(kraken_taxatab)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(metadata,
                           sample_col = "#SampleID",
                           source_col = "Env"
                           )

curves <- calculate_curve(taxa_table, iso_database)
filter_result <- simple_filter(curves, 60)
plot_cuperdec(curves, metadata_table, filter_result)


# decontam ----------------------------------------------------------------

lib_conc$`Full Library Id` %in% names(kraken_seqtab)

kraken_seqtab[is.na(kraken_seqtab)] <- 0 # convert NAs to 0

kraken_seqtab <- kraken_seqtab %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

# test for contaminants
contaminants <- isContaminant(kraken_seqtab, 
                              conc = lib_conc$`Quantification post-Indexing total`)

# filter out contaminant species

taxatable_decontam <- kraken_taxatab %>%
  filter(!contaminants$contaminant)

write_csv(taxatable_decontam, here("04-analysis/taxatable_decontam.csv"))

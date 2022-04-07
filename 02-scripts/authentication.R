# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)
library(tidyverse)
library(here)

# upload data
lib_conc <- readr::read_tsv(here("03-data/SYN_library_quant.tsv")) # library concentrations
kraken_seqtab <- readr::read_csv(here("03-data/sequence_table.csv"))
kraken_taxatab <- readr::read_csv(here("03-data/taxa_table.csv"))
file_names <- list.files(here("03-data/kraken/"), "_report")
sample_names <- gsub(".unmapped.*", "", file_names)

# cuperdec ----------------------------------------------------------------

taxa_table <- load_taxa_table(kraken_taxatab)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(metadata,
                           sample_col = "sample",
                           source_col = "source"
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

write_csv(taxatable_decontam, here("03-data/taxatable_decontam.csv"))

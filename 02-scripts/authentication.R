# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)

# upload library concentrations and data
lib_conc <- readr::read_tsv("03-data/SYN_library_quant.tsv")
kraken_seqtab <- readr::read_csv("03-data/sequence_table.csv")

# cuperdec

cuperdec::cuperdec_taxatable_ex
cuperdec::cuperdec_database_ex

taxa_table <- load_taxa_table(kraken_taxatable)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(cuperdec_metadata_ex,
                           sample_col = "#SampleID",
                           source_col = "Env"
)
curves <- calculate_curve(taxa_table, iso_database)
plot_cuperdec(curves, metadata_table)




data(cuperdec_taxatable_ex)
data(cuperdec_database_ex)
data(cuperdec_metadata_ex)

taxa_table <- load_taxa_table(cuperdec_taxatable_ex)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(cuperdec_metadata_ex,
                           sample_col = "#SampleID",
                           source_col = "Env"
)

curves <- calculate_curve(taxa_table, iso_database)
plot_cuperdec(curves, metadata_table)



# decontam

lib_conc$`Full Library Id` %in% names(kraken_seqtab)

kraken_seqtab[is.na(kraken_seqtab)] <- 0 # convert NAs to 0

kraken_seqtab <- kraken_seqtab %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

# test for contaminants
isContaminant(kraken_seqtab, conc = lib_conc$`Quantification post-Indexing total`)

# filter out contaminant species




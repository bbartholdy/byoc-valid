# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)

# upload library concentrations and data
lib_conc <- readr::read_tsv("03-data/SYN_library_quant.tsv")
kraken_seqtab <- readr::read_csv("03-data/sequence_table.csv")
kraken_taxatab <- readr::read_csv("03-data/taxa_table.csv")
metadata <- readr::read_csv("03-data/sample_metadata.csv")

metadata$sample[-1] <- paste0(metadata$sample[-1], "0101") # match metadata names to sequence names
metadata <- subset(metadata, metadata$sample %in% sample_names) # subset successful sequences (?)


# cuperdec ----------------------------------------------------------------

cuperdec::cuperdec_metadata_ex$Env

taxa_table <- load_taxa_table(kraken_taxatab)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
metadata_table <- load_map(metadata,
                           sample_col = "sample",
                           source_col = "source"
                           )

curves <- calculate_curve(taxa_table, iso_database)
filter_result <- simple_filter(curves, 50)
mean(filter_result$Passed) # all passed
plot_cuperdec(curves, metadata_table, filter_result)


# decontam ----------------------------------------------------------------

lib_conc$`Full Library Id` %in% names(kraken_seqtab)

kraken_seqtab[is.na(kraken_seqtab)] <- 0 # convert NAs to 0

kraken_seqtab <- kraken_seqtab %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

# test for contaminants
isContaminant(kraken_seqtab, conc = lib_conc$`Quantification post-Indexing total`)

# filter out contaminant species




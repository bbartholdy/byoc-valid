# authentication of samples

# see https://github.com/ZandraFagernas/unified_protocol/blob/master/UP_supp_figures_cleaned.Rmd
# https://benjjneb.github.io/DecontamManuscript/Analyses/oral_contamination.html

library(decontam)
library(cuperdec)
library(tidyverse)
library(here)

# upload data
metadata <- readr::read_tsv("01-documentation/metadata.tsv")
lib_conc <- readr::read_tsv("01-documentation/SYN_library_quant.tsv") # library concentrations
otu_filtered_table <- readr::read_tsv("04-analysis/decontam/pre-decontam_OTUfiltered-table_from-biom.tsv", skip = 1)
#comp_taxatable <- readr::read_tsv("04-analysis/pre-decontam_comparative_taxatable.tsv")
sourcetracker2 <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions.txt")
sourcetracker2_stdevs <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions_stds.txt")
sourcetracker2.2 <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output2/mixing_proportions.txt")
sourcetracker2.2_stdevs <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output2/mixing_proportions_stds.txt")
env_controls <- readr::read_csv("03-data/environmental-controls.csv")

file_names <- list.files(here("04-analysis/kraken/"), "_report")
sample_names <- gsub(".unmapped.*", "", file_names)


# SourceTracker -----------------------------------------------------------

sourcetracker2_long <- sourcetracker2 %>%
  pivot_longer(cols = where(is.numeric), 
               values_to = "proportion", 
               names_to = "SampleID") %>%
  rename(source = ...1)

# filter out samples where the proportion of indoor_air exceeds oral source
remove_samples <- sourcetracker2_long %>% 
  pivot_wider(names_from = "source", values_from = "proportion") %>% 
  filter(indoor_air / (`modern calculus`+ plaque + saliva) > 1) %>%
  .$SampleID

analysis_metadata <- metadata %>%
  filter(!`#SampleID` %in% remove_samples)

write_tsv(analysis_metadata, "01-documentation/analysis-metadata.tsv")

# cuperdec ----------------------------------------------------------------

# taxa_table <- load_taxa_table(sample_taxatable)
iso_database <- load_database(cuperdec_database_ex, target = "oral")
# metadata_table <- load_map(metadata,
#                            sample_col = "#SampleID",
#                            source_col = "Env"
#                            )
# 
# curves <- calculate_curve(taxa_table, iso_database)
# filter_result <- simple_filter(curves, 60)
# plot_cuperdec(curves, metadata_table, filter_result)


# decontam ----------------------------------------------------------------

lib_conc_select <- filter(lib_conc, `Full Library Id` %in% rownames(otu_filtered_matrix))

#kraken_seqtab[is.na(kraken_seqtab)] <- 0 # convert NAs to 0
otu_filtered_matrix <- otu_filtered_table %>%
  column_to_rownames(var = "#OTU ID") %>%
  select(which(colnames(.) %in% analysis_metadata$`#SampleID`)) %>%
  t()

# negative controls
neg_controls <- analysis_metadata %>%
  filter(`#SampleID` %in% rownames(otu_filtered_matrix)) %>% # need to filter this out in dataprep or kraken script
  mutate(neg = case_when(
    Env %in% c("sediment", "skin", "indoor_air") ~ TRUE,
    TRUE ~ FALSE)) %>%
  rename(sample = `#SampleID`) %>%
  dplyr::select(sample, neg)

# Need to order neg_controls according to otu_filtered_matrix

analysis_metadata %>%
  filter(!`#SampleID` %in% rownames(otu_filtered_matrix)) # these samples are missing in otu_filtered_matrix

decontam_neg <- neg_controls %>%
  arrange(rownames(otu_filtered_matrix)) %>%
  .$neg
# comp_matrix <- comp_taxatable %>%
#   column_to_rownames(var = "sample") #%>%
#   as.matrix()

# combine with samples in same order as env_controls$neg

# filtered_seqtab <- sample_taxatable %>%
#   column_to_rownames(var = "species") %>%
#   as.matrix() %>%
#   t()
# 
# select_seqtab <- otu_select_samples %>%
#   column_to_rownames(var = "species") %>%
#   as.matrix() %>%
#   t()

# test for contaminants using prevalence method with environmental samples

# test for contaminants
which_contaminants <- isContaminant(
  otu_filtered_matrix,
  #conc = lib_conc$`Quantification post-Indexing total`,
  neg = decontam_neg,
  #threshold = 0.1
  )

# which_contaminants <- isContaminant(
#   otu_filtered_matrix, 
#   #conc = lib_conc_select$`Quantification post-Indexing total`,
#   neg = decontam_neg,
#   method = "prevalence"
#   )

contaminants <- filter(which_contaminants, contaminant == TRUE)

# number of contaminants identified
sum(contaminants$contaminant)

# filter out contaminant species
  # make sure no oral taxa are on the list
contaminant_species <- row.names(contaminants)

oral_taxa <- iso_database
true_contaminants <- contaminants %>%
  filter(!contaminant_species %in% oral_taxa$Taxon)
list_contaminants <- filter(true_contaminants, contaminant == TRUE)

# filter out putative contaminants from OTU table
otu_decontam <- otu_filtered_table %>%
  filter(!`#OTU ID` %in% rownames(list_contaminants))

# how many contaminants per sample?
otu_decontam %>%
  dplyr::select(!c(SYN001.A0101, SYN002.A0101, SYN003.A0101)) %>%
  pivot_longer(
    cols = where(is.numeric), 
    names_to = "SampleID", 
    values_to = "count") %>%
  filter(species %in% contaminant_species) %>%
  #mutate(SampleID = fct_reorder(SampleID, day_order)) %>%
  group_by(SampleID) %>%
  summarise(contam = sum(count)) %>%
  ggplot(aes(x = SampleID, y = contam)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90))

write_tsv(otu_decontam, here("05-results/post-decontam_taxatable.tsv"))
write_tsv(list_contaminants, here("04-analysis/decontam/list-of-contaminants.tsv"))

# authentication of samples

library(decontam)
library(cuperdec)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(tibble)
source("02-scripts/functions.R")

# upload data
metadata <- readr::read_tsv("01-documentation/dna-metadata.tsv")
lib_sample <- readr::read_tsv("05-results/lib_sample.tsv") %>%
  rename("#OTU ID" = species)
lib_conc <- readr::read_tsv("05-results/SYN_DNA_concentrations.tsv") # library concentrations
#otu_table <- readr::read_tsv("04-analysis/pre-decontam_sample_taxatable.tsv")
otu_filtered_table <- readr::read_tsv("04-analysis/qiime/pre-decontam_OTUfiltered-table_from-biom.tsv", skip = 1)
#comp_taxatable <- readr::read_tsv("04-analysis/pre-decontam_comparative_taxatable.tsv")
sourcetracker2 <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions.txt")
sourcetracker2_stdevs <- readr::read_tsv("04-analysis/sourcetracker/sourcetracker2_output/mixing_proportions_stds.txt")

file_names <- list.files(here("04-analysis/kraken/"), "_report")
sample_names <- gsub(".unmapped.*", "", file_names)


# Data prep ---------------------------------------------------------------

# combine library sample with the filtered table
otu_comb_long <- otu_filtered_table %>% 
  pivot_longer(
    cols = where(is.numeric), 
    names_to = "sample", 
    values_to = "count"
    ) %>%
  bind_rows(lib_sample)

# SourceTracker -----------------------------------------------------------

sourcetracker2_long <- sourcetracker2_longer(sourcetracker2)

# contributions in problematic samples
problem_samples <- sourcetracker2_long %>%
  filter(source == "indoor_air",
         proportion > 0.2) %>%
  .$SampleID

# pivot_wider(names_from = "source", values_from = "proportion")

file_names_problems <- paste0("04-analysis/sourcetracker/sourcetracker2_output/", problem_samples, ".feature_table.txt")

file_names_all <- paste0(
  "04-analysis/sourcetracker/sourcetracker2_output/", 
  unique(sourcetracker2_long$SampleID), 
  ".feature_table.txt")
all_data <- lapply(as.list(file_names_all), read_tsv)

names(all_data) <- unique(sourcetracker2_long$SampleID)

all_list_long <- lapply(all_data, function(x) x %>% 
                          pivot_longer(cols = where(is.numeric),
                                       values_to = "count",
                                       names_to = "taxon") %>%
                          rename(source = ...1)
)

all_data_long <- all_list_long %>% 
  map_df(~ as_tibble(.x), .id = "SampleID") %>%
  filter(count > 0)

# How many unknowns are actually oral taxa

iso_database <- load_database(cuperdec_database_ex, target = "oral")
oral_taxa <- iso_database %>%
  filter(Isolation_Source == TRUE)


# filter out samples where the proportion of indoor_air exceeds oral source
remove_samples1 <- sourcetracker2_long %>% 
  pivot_wider(names_from = "source", values_from = "proportion") %>% 
  filter(indoor_air / (`modern calculus`+ plaque + saliva) > 1) %>%
  .$SampleID

remove_samples2 <- all_data_long %>%
  filter(count > 0) %>% 
  mutate(oral_source = if_else(taxon %in% oral_taxa$Taxon, "oral", "other")) %>%
  group_by(SampleID, oral_source) %>%
  summarise(count = sum(count)) %>% 
  mutate(prop = count / sum(count)) %>% 
  filter(oral_source == "oral" & prop < 0.7) %>% 
  .$SampleID

# remove samples present in both remove_samples1 and remove_samples2
remove_samples <- remove_samples1[remove_samples1 %in% remove_samples2]

analysis_metadata <- metadata %>%
  filter(!`#SampleID` %in% remove_samples) %>%
  dplyr::select(!c(Project, download_link))

otu_removed_table <- otu_filtered_table %>%
  dplyr::select(which(!colnames(.) %in% remove_samples))
 
otu_removed_long <- otu_comb_long %>%
  filter(!sample %in% remove_samples)

# decontam ----------------------------------------------------------------

otu_prev_matrix <- otu_removed_table %>%
  column_to_rownames(var = "#OTU ID") %>%
  select(which(colnames(.) %in% analysis_metadata$`#SampleID`)) %>%
  t()

otu_freq_matrix <- otu_removed_long %>% 
  filter(str_detect(sample, "SYN")) %>%
  bind_rows(lib_sample) %>%
  pivot_wider(names_from = "#OTU ID", values_from = "count") %>%
  column_to_rownames(var = "sample") %>% 
  mutate(across(everything(), replace_na, 0)) %>% 
  as.matrix()

lib_conc_select <- filter(lib_conc, Library_Id %in% rownames(otu_freq_matrix))

# negative controls
neg_controls <- analysis_metadata %>%
  filter(`#SampleID` %in% rownames(otu_prev_matrix)) %>% # need to filter this out in dataprep or kraken script
  mutate(neg = case_when(
    Env %in% c("sediment", "skin", "indoor_air") ~ TRUE,
    TRUE ~ FALSE)) %>%
  rename(sample = `#SampleID`) %>%
  dplyr::select(sample, neg)

# Need to order neg_controls according to otu_filtered_matrix

analysis_metadata %>%
  filter(!`#SampleID` %in% colnames(otu_filtered_table)) # these samples are missing in otu_filtered_matrix

decontam_neg <- neg_controls %>%
  arrange(rownames(otu_prev_matrix)) %>%
  .$neg

  # frequency method with library concentrations
which_contaminants_freq <- isContaminant(
  otu_freq_matrix,
  conc = lib_conc_select$`Quantification_post-Indexing_total`,
  #neg = decontam_neg,
  threshold = 0.99
  )
  # prevalence method with environmental samples
which_notcontaminants_prev  <- isNotContaminant(
  otu_prev_matrix,
  #conc = lib_conc_select$`Quantification_post-Indexing_total`,
  neg = decontam_neg,
  detailed=TRUE,
  threshold = 0.01
)

which_contaminants_prev  <- isContaminant(
  otu_prev_matrix,
  #conc = lib_conc_select$`Quantification_post-Indexing_total`,
  neg = decontam_neg,
  threshold = 0.99
)

freq_contaminants <- filter(which_contaminants_freq, contaminant == TRUE)
prev_contaminants <- filter(which_contaminants_prev, contaminant == TRUE)
prev_notcontaminants <- filter(which_notcontaminants_prev, not.contaminant == FALSE)

contaminants <- prev_notcontaminants %>%
  mutate(contaminant = if_else(not.contaminant == TRUE, FALSE, TRUE)) %>%
  select(!not.contaminant) %>%
  bind_rows(
    prev_contaminants,
    freq_contaminants
  )

# number of contaminants identified by each method
sum(which_contaminants_freq$contaminant)
sum(which_contaminants_prev$contaminant)
sum(which_notcontaminants_prev$not.contaminant) # number of non-contaminants

# filter out contaminant species
contaminant_species <- unique(
  c(
    row.names(prev_notcontaminants),
    row.names(prev_contaminants), 
    row.names(freq_contaminants)
    )
  )
  # make sure no oral taxa are on the list
    # import list of oral taxa from cuperdec
oral_taxa <- filter(cuperdec_database_ex, isolation_source == "oral")
true_contaminants <- contaminant_species[!contaminant_species %in% oral_taxa$species]

# filter out putative contaminants from OTU table
otu_decontam <- otu_removed_table %>%
  filter(!(`#OTU ID` %in% true_contaminants))

# export sourcetracker and metadata
write_tsv(all_data_long, "04-analysis/sourcetracker/source-comb_long.tsv")
write_tsv(analysis_metadata, "01-documentation/dna-analysis-metadata.tsv")

# export decontam results
write_tsv(neg_controls, here("04-analysis/decontam/negative-controls.tsv"))
write_tsv(otu_decontam, here("04-analysis/decontam/post-decontam_taxatable.tsv"))
file.symlink("04-analysis/decontam/post-decontam_taxatable.tsv", "05-results/metagenomics/post-decontam_taxatable.tsv")
#write_tsv(otu_decontam, here("05-results/post-decontam_taxatable.tsv")) # copy for results dir
write_tsv(as_tibble(true_contaminants), here("04-analysis/decontam/list-of-contaminants.txt"), col_names = F)
file.symlink("04-analysis/decontam/list-of-contaminants.txt", "05-results/metagenomics/list-of-contaminants.txt")

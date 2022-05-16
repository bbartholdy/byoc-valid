library(tidyverse)
library(here)

# Upload data -------------------------------------------------------------

hmasm_data <- readr::read_csv(here("01-documentation/HMASM.csv"))
hospital_surfaces <- readr::read_tsv(
  here("01-documentation/hospital_surfaces_url.txt"), 
  col_names = F)

# HMP samples -------------------------------------------------------------

names(hmasm_data)
# select relevant samples
source_names <- c("supragingival_plaque", "subgingival_plaque", "saliva")

set.seed(42)
supra <- hmasm_data %>%
  filter(`Body Site` == "supragingival_plaque") %>%
  filter(`Reads File Size` < 2000000000) %>%
  slice_sample(n = 10, replace = F)

# select buccal_mucosa sites
mucosa <- hmasm_data %>%
  filter(`Body Site` == "buccal_mucosa") %>%
  filter(`Reads File Size` < 2000000000) %>%
  slice_sample(n = 10, replace = F)

select_data <- hmasm_data %>%
  filter(`Body Site` == "subgingival_plaque" | `Body Site` == "saliva") %>%
  bind_rows(supra, mucosa)

readr::write_tsv(select_data, here("01-documentation/HMP_source_samples.tsv"))

urls <- paste0("http://downloads.hmpdacc.org", select_data$`Reads file location`) %>%
  as_tibble()

mucosa_urls <- paste0("http://downloads.hmpdacc.org", mucosa$`Reads file location`) %>%
  as_tibble()

saliva_source <- select_data %>%
  filter(`Body Site` == "saliva")
urls <- paste0("http://downloads.hmpdacc.org", select_data$`Reads file location`) %>%
  as_tibble()

readr::write_delim(urls, here("01-documentation/saliva_sample_urls.txt"), col_names = F)
readr::write_delim(mucosa_urls, here("01-documentation/buccal-mucosa_sample_urls.txt"), col_names = F)


# Other samples -----------------------------------------------------------

hospital_meta <- hospital_surfaces %>%
  rename(url = X1) %>%
  mutate("#SampleID" = str_extract(url, "[SRR0-9_1-2.fastq.gz]+$")) %>%
  mutate(`#SampleID` = str_remove(`#SampleID`, "(?<=_)[1-2.fastq.gz]+$")) %>%
  mutate(`#SampleID` = str_remove(`#SampleID`, "[[:punct:]]"),
         Env = "indoor",
         Project = "PRJNA376580",
         SourceSink = "source",
         source = "hospital") %>%
  select(!url)
         #`#SampleID` = str_remove(`#SampleID`, "[2.fastq.gz]")) %>% # need better solution
  view()
  
mucosa_meta <- mucosa %>%
  rename(Env = `Body Site`,
         "#SampleID" = `SRS ID`) %>%
  mutate(SourceSink = "source",
         Project = "HMP") %>%
  select(Env, `#SampleID`, SourceSink)

metadata <- metadata %>%
  bind_rows(hospital_meta, mucosa_meta)


# Input TSV ---------------------------------------------------------------

metadata <- readr::read_tsv("01-documentation/metadata.tsv")
indoor_samples <- readr::read_delim("indoor-samples.txt", delim = "\t", col_names = F)
tsv_template <- readr::read_tsv("01-documentation/TSV_template.tsv")

# Metadata

sample_names <- unique(str_extract(indoor_samples$X1, "^[SRR0-9]+"))
file_names <- indoor_samples

indoor_meta <- tibble("#SampleID" = sample_names, "Env" = "indoor_air", "SourceSink" = "source",
                      "Project" = "PRJNA422794")

metadata <- indoor_meta %>%
  bind_rows(metadata)

# Input TSV

r1 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", file_names$X1[str_detect(file_names$X1, "_1")])
r2 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", file_names$X1[str_detect(file_names$X1, "_2")])

indoor_input <- tibble(
  "Sample_Name" = sample_names,
  "Library_ID" = sample_names,
  "Lane" = 1,
  "Colour_Chemistry" = 4,
  "SeqType" = "PE",
  "Organism" = NA,
  "Strandedness" = "double",
  "UDG_Treatment" = "none",
  "R1" = r1,
  "R2" = r2,
  "BAM" = NA
)

write_tsv(indoor_input, "01-documentation/indoor_input.tsv")

metadata <- metadata %>%
  filter(Env != "indoor_air")


mucosa_names <- mucosa_meta$`#SampleID`
indoor_names <- unique(hospital_meta$`#SampleID`)

mucosa_r1 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", mucosa_names, ".denovo_duplicates_marked.trimmed.1.fastq")
mucosa_r2 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", mucosa_names, ".denovo_duplicates_marked.trimmed.2.fastq")

indoor_r1 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", indoor_names, "_1.fastq.gz")
indoor_r2 <- paste0("/home/bartholdybp/data1/databases/refseq/source/", indoor_names, "_2.fastq.gz")


mucosa_input <- tibble(
  "Sample_Name" = mucosa_names,
  "Library_ID" = mucosa_names,
  "Lane" = 1,
  "Colour_Chemistry" = 4,
  "SeqType" = "PE",
  "Organism" = "Human",
  "Strandedness" = "double",
  "UDG_Treatment" = "none",
  "R1" = mucosa_r1,
  "R2" = mucosa_r2,
  "BAM" = NA
)

hospital_input <- tibble(
  "Sample_Name" = indoor_names,
  "Library_ID" = indoor_names,
  "Lane" = 1,
  "Colour_Chemistry" = 4,
  "SeqType" = "PE",
  "Organism" = NA,
  "Strandedness" = "double",
  "UDG_Treatment" = "none",
  "R1" = indoor_r1,
  "R2" = indoor_r2,
  "BAM" = NA
)

mucosa_indoor_input <- mucosa_input %>%
  bind_rows(hospital_input)

write_tsv(mucosa_indoor_input, here("04-analysis/eager/mucosa-indoor_input.tsv"))

calculus_input <- readr::read_tsv(here("04-analysis/eager/calculus_input.tsv"))

source_input <- read_tsv(here("04-analysis/eager/source_input.tsv"))

source_input <- bind_rows(mucosa_indoor_input, calculus_input, source_input)
view(source_input)

source_input <- source_input %>%
  mutate(Organism = if_else(str_detect(Sample_Name, "ERR"), NA_character_, Organism))

source_input %>%
  mutate(test = str_detect(Sample_Name, "[ERR]+")) %>% view

write_tsv(source_input, here("04-analysis/eager/source_input.tsv"))

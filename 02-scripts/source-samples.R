library(dplyr)
library(magrittr)
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

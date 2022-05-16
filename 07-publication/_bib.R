library(here)
library(rbbt)

file_names <- list.files(here("07-publication"), pattern = ".qmd")

if(file.exists(here("07-publication/references.bib"))) {
  bbt_update_bib(here("07-publication", file_names), 
                 here("07-publication/references.bib"),
                 overwrite = T)
} else {
  refs <- rbbt::bbt_detect_citations(file_names)
  bbt_write_bib(here("07-publication/references.bib"), 
                keys = refs, 
                ignore = bbtignore, 
                overwrite = T, 
                translator = "bibtex")
}


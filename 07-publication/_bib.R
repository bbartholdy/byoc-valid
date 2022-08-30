library(here)
library(rbbt)

file_names <- list.files(here("07-publication"), pattern = ".qmd")


try(
if(file.exists(here("07-publication/references.bib"))) {
  bbt_update_bib(here("07-publication", file_names), 
                 here("07-publication/references.bib"),
                 overwrite = T)
} else {
  bbt_write_bib(here("07-publication/references.bib"), 
                keys = rbbt::bbt_detect_citations(file_names), 
                ignore = bbtignore, 
                overwrite = T, 
                translator = "bibtex")
}
)

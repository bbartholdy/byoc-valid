library(here)
library(rbbt)
library(stringr)

file_names <- list.files(here("07-publication"), pattern = ".qmd", full.names = T)
file_names <- file_names[str_detect(file_names, "_", negate = T)]
bbtkeys <- bbt_detect_citations(file_names)
bbtignore <- bbtkeys[str_detect(bbt_detect_citations(file_names), "^fig-")]

try(
if(file.exists(here("07-publication/references.bib"))) {
  bbt_update_bib(
    file_names, 
    here("07-publication/references.bib"),
    overwrite = T,
    ignore = bbtignore
    )
} else {
  bbt_write_bib(
    here("07-publication/references.bib"), 
    keys = bbtkeys, 
    ignore = bbtignore, 
    overwrite = T, 
    translator = "bibtex"
    )
}
)

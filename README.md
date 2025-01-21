# Assessing the validity of a calcifying oral biofilm model as a suitable proxy for dental calculus

<!-- using the Opinionated Bioinformatics Project Directory Structure
https://github.com/paleobiotechnology/analysis-project-structure -->

![GitHub Release](https://img.shields.io/github/v/release/bbartholdy/byoc-valid)

## Table of Contents

<!-- TOC depthfrom:2 depthto:2 -->

- [Table of Contents](#table-of-contents)
- [Preamble](#preamble)
- [General Organisation](#general-organisation)
- [Structure](#structure)
- [Directory Descriptions](#directory-descriptions)

<!-- /TOC -->

## Structure

### Overall

The overall structure is as in the tree structure below, with the main top-level folders listed in numerical order, alongside a main repository README and a couple of additional useful files.

Brief summary descriptions of the main folders are as follows ([more details here](#directory-descriptions)):

- `01-documentation`: Contains initial metadata about samples and data files that are used for downstream analysis. This includes locations of comparative data.
- `02-scripts`: Contains all scripts and code notebooks used in the day-to-day analysis during the project. Can optionally include sub-directories for each language (e.g., R, python).
- `03-data`: Contains all the large raw, or common-starting point files for all downstream analyses. For (meta)genomics these are normally files such as BAM, SAM, FASTQ, FASTA etc.
- `04-analysis`: Contains all the output from software, tools, and notebooks of all analyses. This is the main 'working' directory of the project.
- `05-results`: Contains copies of all final output from all `04-analysis` (i.e., without intermediate files). These will be used for the bare-minimal reproducible results for reports and publication.
- `06-reports`: Contains presentations, summary notebooks of particular stages or packages of the project. Used for informing the final publication.
- `07-publication`: Contains main text, figures, supplementary files and data. Optionally can formatted with bookdown for pretty online rendering with direct links to intermediate files in `04-analysis`.


### `renv.lock`

This is a file that captures the R environment using the [**renv**](https://rstudio.github.io/renv/articles/renv.html), including packages and R version.

When reproducing the R code, you can use the **renv** package and the function `renv::restore()` to restore the R packages that are needed to
run the R code.


### `.conda_environment.yml`

This is a file that utilises the [`conda`](https://docs.conda.io/en/latest/) (mostly) portable packaging system. You use this, or multiple files, to define all software, and specific versions of said software, used in the project or analysis packages.

After installation of conda, the environment can be created and activated as follows

```bash
conda env create -f conda_environment.yml
conda activate <NAME_OF_PROJECT>
```


## Directory Descriptions

### `01-documentation/`

This directory contains primarily files that you need before you can do any analyses. This often includes files that contain information of wet lab processing of new samples, or metadata files of publicly available datasets used as comparative datasets.

In addition, we recommend – upon final publication of the project – placing in this directory a final table that contains the location(s) of all novel data generated in this project on public archives. For sequencing data this will normally be the table that can be exported from the ENA or SRA that has the FTP or Aspera URLs of the uploaded FASTQ files.

Most of these files will be simple text files in tabular format, such as CSV, TSV, XLSX, or TXT files.

### `02-scripts/`

This contains all scripts and notebooks used in the 'day-to-day' analysis of the project. All of these scripts/notebooks produce both intermediate and final files used in the analysis of the project.

These can be things such as simple shell/bash scripts (`.sh`) or R (`.R`) scripts, or notebooks such as RMarkdown/Notebooks (`.Rmd`) or Jupyter notebooks (`.ipynb`). It is optional how this is structured internally, whether by analysis component prefixes, directories per analysis component, or by programming language.

These scripts and notebooks should only use relative links to refer to input and output files present in the directory, and not to point to anything else present on your given machine or infrastructure.

File names should be descriptive so readers can find relevant files. Abbreviations or acronyms are not recommended, as these are often difficult to understand for people not intimately involved in the project.


### `04-analysis/`

This directory contains the bulk of the analysis carried out in the project.

Most scripts in `02-scripts/` will refer to this directory, in terms of input and output files, and working directories. This directory will NOT be uploaded to the Git repository as the output of many analyses will be very large (see `03-data/`). However, it should be made sure that the internal directory structure of this folder can be reconstructed based on the notebooks and scripts that list the actual analyses performed.

Due to this directory being the main working directory for most of the project, the structure of this directory can be more flexible, and structured in a way that fits the project best. However, we generally recommend structuring this directory on a tool-by-tool basis, i.e., each tool will have one directory. Within this, there will be different analysis runs as required.

Either symbolic links or relative paths should be used to link the output of one tool as input to another tool (i.e., do NOT copy these files). _However_, once an analysis is finalised and will not be run again - the final and relevant output files should be _copied_ into `05-results/`, as these _will_ be uploaded to the online Git repository.

### `05-results/`

This directory contains all the final files generated in `04-analysis` (copied over), and will be the main files that can be referred to in the final publication supplementary information. These should consist of tables, text files, and **small** binary files that are important for understanding the interpretations made in the final publication.

This should NOT include intermediate, temporary or working/scratch files.

We recommend structuring this section either by sub-directories, or file prefixes, that correspond to each 'logical' analysis section of the publication (this can be thought as each section of the main text of the publication itself). This allows someone entering the repository to find a specific file required to answer their specific question, and other possible files.

For example, there could be a quality control set of files indicated with the key `QUAL` which consists of all relevant log files and metric tables that describe the outcome of preprocessing of sequencing data. A second analysis section could be for phylogenetic analyses, which would be indicated with `PHYLO`.

### `06-reports/`

This is a recommended directory that contains documents or files that can be useful for summarising main results, including interpretation, and ultimately can be used to help inform the writing of the final publication. It will be uploaded to the online Git repository.

This directory can include things such as presentation slides, or finalised notebooks that summarise the outcomes of each different analysis section. For the latter, these files should gather and aggregate various output files from `05-results`, perform different summary analyses, and generate summary plots and figures in a reproducible manner. In this case we again recommend to use software environment managers or containers to allow for such notebooks to be executed in away that allows for such reproducibility.

It is also recommended to organise these in analysis-batch specific directories or prefixes that correspond to the input files in `05-results`.

These reports or notebooks can go into more technical detail that any file for the final publication, including describing explorations and/or failures for future prosperity. Not all the contents of these reports will necessarily go into the `07-publication` directory.

### `07-publication/`

This directory contains all the final files for publication. Manuscript files
are written in Quarto (*.qmd*). The *_extension* folder contains the Quarto [arXiv
extension](https://github.com/mikemahoney218/quarto-arxiv) for formatting.

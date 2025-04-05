# Assessing the validity of a calcifying oral biofilm model as a suitable proxy for dental calculus

<!-- using the Opinionated Bioinformatics Project Directory Structure
https://github.com/paleobiotechnology/analysis-project-structure -->

[![GitHub Release](https://img.shields.io/github/v/release/bbartholdy/byoc-valid)](https://github.com/bbartholdy/byoc-valid/releases/tag/v23.05.0) [![DOI](https://img.shields.io/badge/DOI-10.12688/openreseurope.19129.1-blue)](https://doi.org/10.12688/openreseurope.19129.1)


## Table of Contents

<!-- TOC depthfrom:2 depthto:2 -->

- [Table of Contents](#table-of-contents)
- [Structure](#structure)
- [Computational reproducibility](#computational-reproducibility)
- [Directory Descriptions](#directory-descriptions)

<!-- /TOC -->

## Structure

### Overall

The overall structure is as in the tree structure below, with the main top-level folders listed in numerical order, alongside a main repository README and a couple of additional useful files.

Brief summary descriptions of the main folders are as follows ([more details here](#directory-descriptions)):

- `01-documentation`: Contains initial metadata about samples and data files that are used for downstream analysis. This includes locations of comparative data.
- `02-scripts`: Contains all scripts used in the analysis during the project.
- `03-data`: Contains all the large raw, or common-starting point files for all downstream analyses. **These are too big to include in the GitHub repository**.
- `04-analysis`: Contains all the output from software, tools, and notebooks of all analyses. This is the main 'working' directory of the project.
- `05-results`: Contains copies of all final output from all `04-analysis` (i.e., without intermediate files). These will be used for the bare-minimal reproducible results for reports and publication.
- `06-reports`: Contains presentations, summary notebooks of particular stages or packages of the project. Used for informing the final publication.
- `07-publication`: Contains main text, figures, supplementary files and data. Optionally can formatted with bookdown for pretty online rendering with direct links to intermediate files in `04-analysis`.


## Computational reproducibility

To the best of my ability, I have provided all the code needed to reproduce the analysis. I welcome
all attempts to reproduce/replicate this study, and I am happy to hear about any issues that were
encountered by the reproducer and any ways that I can improve this work.

### `renv.lock`

This is a file that captures the R environment using the [**renv**](https://rstudio.github.io/renv/articles/renv.html), including packages and R version.

When reproducing the R code, you can use the **renv** package and the function `renv::restore()` to restore
the R packages that are needed to run the R code.

### `01-documentation/software_versions.csv`

This file contains the software and versions that were used to run the DNA preprocessing in [EAGER](https://nf-co.re/eager/) (using
[Kraken2](https://github.com/DerrickWood/kraken2) for metagenomic classification).

### `01-documentation/conda_versions.tsv`

This file contains the dependencies for QIIME2 and SourceTracker2.


### Analysis

Pre-processing of the DNA was done using EAGER, with the SLURM scripts
in `02-scripts/` prefixed with 'PRE'. The output files from EAGER and Kraken were
combined in R (`02-scripts/01-comb-kraken-reports.R`).
OTU table was filtered for relative abundance. Percent abundance of each taxon
across all samples was calculated and then taxa with lower than 0.001% abundance
were filtered out (`02-scripts/01-dataprep.R`).

Authentication was done using QIIME2 (and `filter_samples_from_otu_table.py` from QIIME) and SourceTracker2.
The steps are presented in `02-scripts/AUTH_01_ST2.md`.

Oxygen tolerance for bacterial species was retrieved from [BacDive](https://bacdive.dsmz.de) on
2022-08-26. Alpha and Beta diversity were calculated using the [**vegan**](https://vegandevs.github.io/vegan/) and
[mixOmics](https://mixomics.org/) R packages (`02-scripts/DIV_01_alpha-beta.R). More details can
be found in `06-reports/metagen-diversity.qmd`.

Differential abundance was calculated using the [ANCOMBC](https://bioconductor.org/packages/release/bioc/html/ANCOMBC.html) R package
(`02-scripts/DIFF_01_lfc.R`). More details can be found in `06-reports/metagen-diffabund.qmd`.

FTIR data were cleaned and processed in R (`02-scripts/FTIR_00_data-prep.R`). More details can be found in `06-reports/FTIR-analysis.qmd`.

Preparation of analysis, figures, and tables for the manuscript can be found in the scripts with the prefix 'OUT' (`02-scripts/OUT_*),
and the manuscript source file is `07-publication/index.qmd`.


## Directory Descriptions

### `01-documentation/`

This directory contains primarily files that you need before you can do any analyses. This often includes files that contain information of wet lab processing of new samples, or metadata files of publicly available datasets used as comparative datasets.

In addition, we recommend – upon final publication of the project – placing in this directory a final table that contains the location(s) of all novel data generated in this project on public archives. For sequencing data this will normally be the table that can be exported from the ENA or SRA that has the FTP or Aspera URLs of the uploaded FASTQ files.

Most of these files will be simple text files in tabular format, such as CSV, TSV, XLSX, or TXT files.

### `02-scripts/`

This contains all scripts and notebooks used in the 'day-to-day' analysis of the project. All of these scripts/notebooks produce both intermediate and final files used in the analysis of the project. These contain .slurm, .R, and .md files.
Naming of the scripts is according to stage of analysis. PRE =
Pre-processing; AUTH = authentication; DIV = diversity; DIFF = differential abundance; FTIR = FTIR analysis; OUT = output-related files (used in publication).
Most of the outputs from these scripts can be found in `04-analysis/`.

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

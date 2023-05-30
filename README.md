# BYOC Validation

<!-- using the Opinionated Bioinformatics Project Directory Structure
https://github.com/paleobiotechnology/analysis-project-structure -->

![Version Badge](https://img.shields.io/badge/github/v/release/bbartholdy/byoc-valid?display_name=tag)

## Table of Contents

<!-- TOC depthfrom:2 depthto:2 -->

- [Table of Contents](#table-of-contents)
- [Preamble](#preamble)
- [General Organisation](#general-organisation)
- [Structure](#structure)
- [Directory Descriptions](#directory-descriptions)

<!-- /TOC -->

## Preamble

To ensure reproducibility and readability within bioinformatics projects, we propose the following project repository structure. We hope that it will allow other researchers to easily follow how data was processed and analyses were performed.

This structure has been developed out of (meta)genomics projects, however it should be generic enough for other fields to adopt.

## General Organisation

The three main concepts for this structure is as follows:

1. Everything should be self-contained (as far as possible) within a single git repository
2. Minimise hard-disk footprint: link, do not copy
3. Every component should follow a logical flow, being easily be linked together, and followable by an outsider

For **points 1 and 2**, this means that all links should be relative to the repository itself. By doing, so, this prevents ballooning hard-disk space usage that occurs by repeatedly copying the same large files used in input for different analyses.

Scripts, notebooks, and similar should utilise relative paths where the repository folder is the top-level part of the directory tree, and refer to files or directories that traverse between different internal folders (such as from raw data to analysis, or between analyses folders) should be with shortcuts such as `../../`.

> :open_book: Example: an analysis command for analysis B, which uses output from analysis A as input, should be written as follows:
>
>```bash
>tool_b --input ../analysis_A/output_A.txt --output output_b
>```

Alternatively, symbolic links (symlinks) can also be alternatively utilised, e.g. generated with `ln -s`, to refer to files within a given directory, even if stored in a different one. Symlinks should still be made relative within the repository.

> :open_book: Example:
>
>```bash
>mkdir analysis_b/ && cd analysis_b/
>mkdir input_b/ output_b/
>cd input_b/
>ln -s ../../analysis_a/output/*.txt .
>cd ..
>tool_b --input input_b/*.txt --output output_b
>```

Likewise, the use of a single conda environment or container (docker, singularity, etc.), is recommended. This means that all software versions are also exactly the same as used by researchers aiming to reproduce the analysis.

For **point 3**, this is actualised here by following a semi-strict directory structure that follows the major steps of a typical academic project of approximately raw data → processing and cleaning → analysis → publication, and heavy use of documentation via README files and reproducible notebooks (e.g., RNotebook, jupyter notebooks etc.).

Ultimately, write and structure everything with the vision that you are another researcher (or yourself!) trying to re-run your analyses in the future.

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


```bash
tree -I renv
```

```bash
.
├── 01-documentation
│   ├── analysis-metadata.tsv
│   ├── buccal-mucosa_sample_urls.txt
│   ├── conda_versions.tsv
│   ├── experiment-metadata.tsv
│   ├── ftir-metadata.tsv
│   ├── genus-O2tolerance.tsv
│   ├── metadata.tsv
│   ├── nonHMP_source.txt
│   ├── software_versions.csv
│   ├── source_samples_key.csv
│   ├── species-properties.tsv
│   ├── SYN_DNA_concentrations.tsv
│   └── TSV_template.tsv
├── 02-scripts
│   ├── 00-comb-kraken-reports.R
│   ├── 00_ftir-analysis.R
│   ├── 01-dataprep.R
│   ├── 02-authentication.R
│   ├── 03-bacproperties.R
│   ├── 04-diversity.R
│   ├── 05-diff-abundance.R
│   ├── 99_ftir-figures.R
│   ├── 99_metagen-figures.R
│   ├── 99_setup-qmd.R
│   ├── 99_tables.R
│   ├── functions.R
│   └── _source-samples.R
├── 03-data
│   ├── 2022-05-27_bacdive_oxytol-search.csv
│   ├── 2022-12-12_bacdive-oxytol-search.csv
│   ├── 2023-03-24_bacdive_oxytol-search.csv
│   ├── FTIR
│   ├── kraken-OTU_long.csv
│   ├── raw_data
│   ├── reference_genomes
│   └── sample_metadata.csv
├── 04-analysis
│   ├── decontam
│   ├── eager
│   ├── FTIR
│   ├── kraken
│   ├── OTUfilter_table.tsv
│   └── sourcetracker
├── 05-results
│   ├── all-pca-loadings_ext.tsv
│   ├── all-pca-loadings.tsv
│   ├── alpha-diversity.tsv
│   ├── byoc_logf-full.tsv
│   ├── clr-byoc.tsv
│   ├── clr-compar-extended.tsv
│   ├── clr-compar.tsv
│   ├── experiment-pca-loadings.tsv
│   ├── ftir
│   ├── ftir-data.csv
│   ├── ftir-data_long.csv
│   ├── grind-data_cleaned.csv
│   ├── lib_sample.tsv
│   ├── plaque_logf-full.tsv
│   ├── post-decontam_taxatable.tsv
│   ├── README
│   ├── spca_byoc.rda
│   ├── spca_species_ext.rda
│   └── spca_species.rda
├── 06-reports
│   ├── diff-abund.nb.html
│   ├── FTIR-analysis.docx
│   ├── FTIR-analysis.html
│   ├── FTIR-analysis.nb.html
│   ├── FTIR-analysis.Rmd
│   ├── metagen-authentication.nb.html
│   ├── metagen-authentication.Rmd
│   ├── metagen-diffabund.nb.html
│   ├── metagen-diffabund.Rmd
│   ├── metagen-diversity.nb.html
│   ├── metagen-diversity.Rmd
│   ├── metagen-explore.nb.html
│   ├── metagen-explore.Rmd
│   ├── metagen-summary.html
│   ├── metagen-summary.nb.html
│   └── metagen-summary.Rmd
├── 07-publication
│   ├── 01-intro.qmd
│   ├── 02-matmet.qmd
│   ├── 03-results.qmd
│   ├── 04-disc.qmd
│   ├── apa.csl
│   ├── _archive
│   ├── arxiv.sty
│   ├── _bib.R
│   ├── _extensions
│   ├── figures
│   ├── _freeze
│   ├── index_cache
│   ├── index_files
│   ├── index.qmd
│   ├── index.tex
│   ├── orcidlink.sty
│   ├── _output
│   ├── _quarto.yml
│   ├── references.bib
│   ├── _supp-mat.qmd
├── byoc-valid.Rproj
├── CHANGELOG_files
│   └── mediabag
├── CHANGELOG.md
├── Dockerfile
├── README.md
├── renv
└── renv.lock
```

An important component of this structure to note is within  scripts, results, and reports, that files are prefixed with short three or four letter codes in all capitals (e.g. QUAL, ANA). This indicates batches of files that are all related to a single analysis component or question that is addressed in the project. These prefixes can also be replaced with specific directories instead. The purpose of this structure is to make it easier for readers to rapidly find scripts of interest, and possible scripts of interest. This system is preferable over, e.g., numbering as it allows more flexibility in case analysis components change in order in the final publication, or morph and merge.

### Repository main `README.md`

This markdown (`.md`) file acts as the main description of the repository as a whole.

This should be the main port-of-call for anyone who wants to explore, replicate, or try and find any analysis or results file. Generally you can describe the contents of each top-level directory, point to the original publication and other related metadata.

### Repository `.gitignore`

This top-level repository file tells Git which files or folders should not be uploaded to online Git repositories. This is particularly important not just for data-sensitive files (login tokens, etc.), but also for ensuring very large files are not uploaded. For example, GitHub has a strict limit of a single file can be no more than 100MB.

In particular, the directories that should be added to the .gitignore`:

- `03-data`
- `04-analysis`

For safety, we also recommend adding the following file formats to `.gitignore`, when working on (meta)genomics projects.

- `*.bam`
- `*.sam.gz`
- `*.sam`
- `*.fastq.gz`
- `*.fq.gz`
- `*.fastq`
- `*.fq`
- `*.fasta.gz`
- `*.fna.gz`
- `*.ffn.gz`
- `*.fa.gz`
- `*.fasta`
- `*.fna`
- `*.ffn`
- `*fa`

It is critical, however, that all files can be reconstructed. This comes through having highly descriptive, easy-to-follow, and reproducible notebooks and scripts.

Furthermore, we recommend providing sufficient information or automated scripts available at publication that allows a user to download the raw data uploaded to public repositories (such as the [`ENA`](https://www.ebi.ac.uk/ena/browser/home) or [`SRA`](https://www.ncbi.nlm.nih.gov/sra) for (meta)genomics) and replicates the repository structure of ignored directories that was on the original machine used for analysis.

### `.conda_environment.yml`

This is a file that utilises the [`conda`](https://docs.conda.io/en/latest/) (mostly) portable packaging system. You use this, or multiple files, to define all software, and specific versions of said software, used in the project or analysis packages.

> 💡 This can alternatively be replaced with a [Docker](https://docker.com) description file, for even more robust reproducibility.

This allows visitors of the repository who wish to recreate steps or all your analysis but using the same software configuration. Generally, use of the [`bioconda`](https://bioconda.github.io) software channel is recommended for defining software and versions, as this is designed specifically for bioinformatics, and will likely contain the majority of tools will need for projects particularly in (meta)genomics.

After installation of conda, the environment can be created and activated as follows

```bash
conda env create -f conda_environment.yml
conda activate <NAME_OF_PROJECT>
```

Note that you may need to update this overtime, as you add more analyses. If this happens often, it is recommended to make multiple conda environmental files for different analyses components.

> ⚠️ In large projects with many software, you may end up getting dependency conflicts, where two software require two different versions of the same software. In such a case we recommend splitting up the conda environment into multiple for each component.

## Directory Descriptions

### `01-documentation/`

This directory contains primarily files that you need before you can do any analyses. This often includes files that contain information of wet lab processing of new samples, or metadata files of publicly available datasets used as comparative datasets.

In addition, we recommend – upon final publication of the project – placing in this directory a final table that contains the location(s) of all novel data generated in this project on public archives. For sequencing data this will normally be the table that can be exported from the ENA or SRA that has the FTP or Aspera URLs of the uploaded FASTQ files.

Most of these files will be simple text files in tabular format, such as CSV, TSV, XLSX, or TXT files.

### `02-scripts/`

This contains all scripts and notebooks used in the 'day-to-day' analysis of the project. All of these scripts/notebooks produce both intermediate and final files used in the analysis of the project.

These can be things such as simple shell/bash scripts (`.sh`), python (`.py`), or R (`.R`) scripts, or notebooks such as RMarkdown/Notebooks (`.Rmd`) or Jupyter notebooks (`.ipynb`). It is optional how this is structured internally, whether by analysis component prefixes, directories per analysis component, or by programming language.

These scripts and notebooks should only use relative links to refer to input and output files present in the directory, and not to point to anything else present on your given machine or infrastructure.

File names should be descriptive so readers can find relevant files. Abbreviations or acronyms are not recommended, as these are often difficult to understand for people not intimately involved in the project.

### `03-data/`

This directory directory contains all raw data, and/or initial pre-processed data that would be repeatedly used as input into downstream analyses. This can also include databases and other very large, unmodifiable, and reusable files.

Importantly, this directory is **not** uploaded to the Git repository, due to the very large size of raw sequencing files for (meta)genomic datasets (i.e., GBs in size). Therefore, it is critical that it must be easily reconstructible, i.e., that data can be easily downloaded (e.g. from public databases from the ENA/SRA, in the context of DNA sequencing data), as well as structured in the way expected by all downstream steps.

We recommend to split such a directory such as you separate newly generated data (which would either be symlinked from a centralised institutional raw-data dump, or downloaded from your sequencing facilities' data transfer system), from already public data.

Recommended subdirectories with descriptions are described below.

#### `raw_data/`

Includes unpublished raw data generated for the project. For (meta)genomics, we recommend organising it in the form one sub-directory within `raw_data` per sequencing library (e.g. all FASTQs of all lanes of the same library from the same run, both R1 and R2, would all be in one sub-directory). Sequencing data in here will often still have adapters and other artefacts.

#### `published_data/`

Includes previously published data, such as comparative datasets that are downloaded from public repositories (such as ENA/SRA for (meta)genomics). These should be organised similarly as `raw_data`, however may contain partially processed data (e.g. already adapter clipped).

It is important that the exact version of the data, date, and command that was used to download such data is recorded in `01-documentation`, for example in a metadata sheet. Tables that can be downloaded from respective databases, such as from the ENA or SRA web browser(s), can be useful for this.

#### `reference_genomes/`

For (meta)genomics, often analysts will compare DNA sequencing reads against a reference genome. Reference genome sequences are often stored in text formats such as `fasta`, alongside additional auxiliary files such as annotation files. These should also be placed here, again with a suggestion of one per organism. Often these are downloaded from the [NCBI Genome database](https://www.ncbi.nlm.nih.gov/genome/).

It is important that the exact version of the data, date, and command that was used to download such data is recorded in `01-documentation`, for example in a metadata sheet. Particularly important is the specific assembly or accession number of a given genome is recorded, as genome assemblies can be updated over time (and sometimes removed).

#### `databases/`

In addition to large sequencing files, particularly for metagenomics, tools require large database files consisting of hundreds and thousands of reference genomes. Such database directories and/or files should be stored in this directory, with a recommendation of having one sub-directory per tool, and then within this, the different types databases that maybe used in downstream analysis.

In some cases you may have two variants of this, one with 'raw' files, and another in 'indexed'. If you index the database yourself, you should also consistently organise this across all subdirectories.

In `01-documentation/` it should be recorded where/how each database was downloaded, and if manually indexed, and the corresponding scripts should be recorded in `02-scripts`.

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
are written in Quarto (*.qmd*). The *_extension* folder contains the Quarto arXiv
extension for formatting.

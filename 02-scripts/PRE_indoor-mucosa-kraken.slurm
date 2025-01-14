#!/bin/bash
#SBATCH --job-name=mucosa-indoor-source_kraken
#SBATCH --output=%x_%j.out
#SBATCH --partition=cpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=b.p.bartholdy@arch.leidenuniv.nl
#SBATCH --mail-type="ALL"

# load modules
module load Java/11.0.2
module load Singularity/3.6.1-Go-1.14
module load Nextflow/21.03.0
module load Kraken2/2.1.1-gompi-2020a-Perl-5.30.2

echo "## Module loaded"
module list

# current directory
CWD=$(pwd)
echo "## Starting off from $CWD"

# transfer job files to a job folder
mkdir jobs/$SLURM_JOB_ID data1/analysis/$SLURM_JOB_ID
cp byoc.config source_kraken.slurm jobs/$SLURM_JOB_ID/

ANALYSIS_DIR=~/data1/analysis/$SLURM_JOB_ID/

SDIR=/scratchdata/$SLURM_JOB_USER/$SLURM_JOB_ID

# kraken database
#kraken2-build --standard --db "$SDIR"/krakendb

# move to local scratch
echo "## Move to $SDIR"
cd $SDIR
echo "## Moved to $(pwd)"

cp "$CWD"/byoc.config "$SDIR"/

#export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzOTc1fS40NmZhNzQwODIwYWQ4NmExNTA2YjQyMDE0NjgyMzZlNzljY2Y1ZjZi
#export NXF_VER=21.03.0

# run nextflow
echo "## Running nextflow"
nextflow run nf-core/eager -r 2.4.4 \
-c byoc.config \
-profile alice,byoc \
-name mucosa-indoor-source_kraken \
--outdir "$ANALYSIS_DIR"/results \
-work-dir "$ANALYSIS_DIR"/work \
--input "$CWD"/byoc_analysis/eager/mucosa-indoor_input.tsv \
--fasta "$CWD"/data/databases/refseq/human/GRCh38_latest_genomic.fna \
--metagenomic_tool 'kraken' \
--database "$CWD"/data/databases/kraken/k2_standard_20210517.tar.gz
echo "## Running nextflow ... Done"

# cleaning up intermediate files
echo "## Cleaning up intermediate files"
cp ./.nextflow.log "$ANALYSIS_DIR"
nextflow clean -f -k
echo "## Cleaning up intermediate files ... Done"

# check directory content
echo "## Check directory content"
ls ./
echo "## Check directory content ... Done"

echo "## eager_test done. Have a nice day."

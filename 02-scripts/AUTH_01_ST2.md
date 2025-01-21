OTU table was filtered for relative abundance. Percent abundance of each taxon
across all samples was calculated and then taxa with lower than 0.001% abundance
were filtered out. See */02-scripts/01-dataprep.R*.

Installation of Qiime2 and dev version of SourceTracker2 via conda:

```sh
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-linux-conda.yml
rm qiime2-2022.2-py38-linux-conda.yml # cleanup
conda activate qiime2-2022.2
pip install https://github.com/biota/sourcetracker2/archive/master.zip
```

Installation of Qiime1 to access `filter_samples_from_otu_table.py` script

```sh
conda create -n qiime1 python=2.7 qiime -c bioconda
```

Convert OTU table from *.tsv* to *.biom*.

```sh
biom convert -i 04-analysis/OTUfilter_table.tsv -o 04-analysis/sourcetracker/OTUfilter-table-from-tsv_json.biom --table-type="OTU table" --to-json
```

Filter OTU table

```sh
filter_samples_from_otu_table.py \
-i 04-analysis/sourcetracker/OTUfilter-table-from-tsv_json.biom \
-o 04-analysis/sourcetracker/OTUs1000filter_table.biom \
-n 1000
```

Table summary

```sh
biom summarize-table -i 04-analysis/sourcetracker/OTUs1000filter_table.biom > 04-analysis/sourcetracker/summary_OTUs1000filter-table.txt
```

Convert to TSV for use with the **decontam** package.

```sh
biom convert -i 04-analysis/sourcetracker/OTUs1000filter_table.biom -o 04-analysis/decontam/pre-decontam_OTUfiltered-table_from-biom.tsv --to-tsv
```

Run SourceTracker2 on the filtered OTU table with rarefaction depth of 1000
for both source and samples. Samples and sources were mapped in the
[ST_comb-plaque-map.txt](../04-analysis/sourcetracker/ST_comb-plaque-map.txt).
The plaque source is a combination of supragingival and subgingival plaque.

```sh
conda activate qiime2-2022.2
```

Run with indoor_air sources included and sediment removed

```sh
sourcetracker2 \
    -i 04-analysis/sourcetracker/OTUs1000filter_table.biom \
    -m 04-analysis/sourcetracker/ST_comb-plaque-map.txt \
    --source_sink_column SourceSink  \
    --source_column_value source \
    --source_rarefaction_depth 1000 \
    --sink_rarefaction_depth 1000 \
    --sink_column_value sink \
    --source_category_column Env \
    -o 04-analysis/sourcetracker/sourcetracker2_output \
    --jobs 2 \
    --per_sink_feature_assignments
```

The full results can be found in *04-analysis/sourcetracker/sourcetracker2_output*

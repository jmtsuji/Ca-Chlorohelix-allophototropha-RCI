# 16S ribosomal RNA gene amplicon analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

**NOTE: for each code section provided below, the code ought to be run from within this `16S_rRNA_gene_amplicon_analysis` directory.**

## Install AXIOME3
[AXIOME3](https://github.com/neufeld/axiome3) is a convenient workflow for analyzing gene amplicon sequencing data using QIIME2.  

Here, we'll install and use commit 1ec1ea6 based on QIIME version 2019.10.

Note that you will need to have pre-installed [miniconda (e.g., miniconda3)](https://docs.conda.io/en/latest/miniconda.html).

```bash
git clone https://github.com/neufeld/AXIOME3.git

cd AXIOME3

git checkout 1ec1ea6

# This step (creating the environment) could take time...
conda env create --name axiome3_1ec1ea6 --file conda_env_file/qiime2-2019.10_AXIOME3.yml
```

Done. Before running AXIOME3, make sure to activate the environment by running `conda activate axiome3_1ec1ea6`.

## Analyze the data
Analysis is done within the AXIOME3 folder cloned from Github.

For classification, you will need to download and train the Silva classifier, which is described 
[here](https://docs.qiime2.org/2019.10/tutorials/feature-classifier/) but not shown in this tutorial.

Make the config file
```bash
# Actiavte the environment by running:
# conda activate axiome3_Jun2020_1ec1ea6

# Generate the run config using the provided manifest file in this repo
python luigi_config_generator.py \
	--manifest ../qiime2_manifest.tsv \
	--sample-type SampleData[PairedEndSequencesWithQuality] \
	--input-format PairedEndFastqManifestPhred33V2
```
Note that, in reality, I ran these two samples in a larger data analysis run with other 16S rRNA gene amplicon samples. 
This makes little difference on the output ASVs, however.

Configure `configuration/luigi.cfg` by editing it in a text editor:
- Denoise: n_threads = 10
- Taxonomic Classification: n_jobs = 10
- Make sure to point AXIOME3 to the location of the classifier file you made
- Left everything else as DEFAULTS

Then run AXIOME3
```bash
python 16S_pipeline.py Core_Analysis --local-scheduler 2>&1 | tee ../axiome3.log
```
Once finished, the ASV table will be available at `AXIOME3/output/exported/ASV_table_combined.tsv`.


**In lieu of uploading all intermediate analysis files, a summary of the non-rarefied amplicon sequencing variant (ASV) table is provided at `ASV_table_non_rarefied.tsv`. **
**ASV table data was used to generate Supplementary Data 6 in the manuscript.**

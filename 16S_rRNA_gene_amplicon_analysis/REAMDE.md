# 16S ribosomal RNA gene amplicon analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020


## Data download
You can download the data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="16S_rRNA_gene_data_accessions.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
r1_urls=($(cat "${accession_data_filepath}" | cut -f 6 | tail -n +2))
r2_urls=($(cat "${accession_data_filepath}" | cut -f 7 | tail -n +2))

# Iternatively download the FastQ files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  r1_url=${r1_urls[${j}]}
  r2_url=${r2_urls[${j}]}

  # Download
  echo "[$(date -u)]: Downloading '${sample_id}'"
  wget -O "${output_dir}/${sample_id}_R1.fastq.gz" "${r1_url}"
  wget -O "${output_dir}/${sample_id}_R2.fastq.gz" "${r2_url}"
done

echo "[$(date -u)]: Finished."
```

You now have FastQ files for the two datasets.

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


**In lieu of uploading all intermediate analysis files, a summary of the non-rarefied amplicon sequencing variant (ASV) table is provided at `ASV_table_non_rarefied.tsv`. 
ASV table data was used to generate Supplementary Data 1 in the manuscript.**

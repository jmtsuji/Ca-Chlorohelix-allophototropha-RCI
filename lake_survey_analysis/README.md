# Lake survey metagenome/metatranscriptome analysis
Part of the larger '*Ca.* Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

**NOTE: for each code section provided below, the code ought to be run from within this `lake_survey` directory.**

## Data download
You can download the metagenome data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="metagenome_data_accessions_ncbi.tsv"
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

You now have FastQ files for several of the metagenome datasets.  

Unfortunately, some of the metagenomes are stored in the JGI database, and there is not a simple programmatic way to download these. 
See `metagenome_data_accessions_jgi.tsv` for the IMG Genome IDs of each of these metagenomes, if you wish to manually download them.

To download the metatranscriptome data from NCBI, repeat the above code block using the `metatranscriptome_data_accessions.tsv` file 
in place of `metagenome_data_accessions_ncbi.tsv`.

## Lake metagenome analysis

### Install ATLAS
[ATLAS](https://github.com/metagenome-atlas/atlas) is a metagenome QC, assembly, and binning workflow.  

We'll install and use version 2.1.4 for metagenome data.

Note that you will need to have pre-installed [miniconda (e.g., miniconda3)](https://docs.conda.io/en/latest/miniconda.html).

```bash
# This step (creating the environment) could take time...
conda create -n atlas_2.1.4 -c bioconda -c conda-forge metagenome-atlas=2.1.4
```

Done. Before running ATLAS, make sure to activate the environment by running `conda activate atlas_2.1.4`.

If it is your first time using ATLAS, note that a large amount of database files will be auto-downloaded during the run.

### Process the metagenome data
A config file and sampel ID file for the run have already been created at `lake_metagenomes/config.yaml` and `lake_metagenomes/samples.tsv`.   
If you want to create config files for yourself, you can use the `atlas init` command as documented in the ATLAS repo.  

Notes:
- make sure you modify the `config.yaml` file in the `lake_metagenomes` folder so that `database_dir` is a real directory on your machine.
- you might have to modify the filepaths for the samples in the `samples.tsv` file so that they are correct on your system.

Start the ATLAS run
```bash
# Activate the environment by running:
# conda activate atlas_2.1.4

cd lake_metagenomes
# Notice that samples.tsv and config.yaml files are provided here

atlas run -w . -c config.yaml -j 50 all --reason 2>&1 | tee atlas_run.log

cd ..
```
Note this might take several weeks of computational runtime and considerable (e.g., 200 GB) RAM!! 
If you just want the genome bins, you can skip running ATLAS yourself and just directly download the bins (see below).

Run the GTDB taxonomy classifier on the samples after completing the ATLAS run:
```bash
TODO
```


### Notes from my ATLAS runs
The above command (i.e., a single `atlas run` command) is a simplified version of how ATLAS was actually run when analyzing the data for this paper. 
When working on this paper, the ATLAS analysis was actually spread across several different ATLAS run attempts rather than run in a single ATLAS 
command. For transparency's sake, here are a few major notes from those ATLAS runs:

```
- In actuality, I ran the QC and assembly modules of ATLAS in several batches as different metagenomes became available. Then combined all runs to perform the subsequent steps in the pipeline such as genome binning.
- Ran JGI metagenomes in `interleaved` mode for the QC portion of the pipeline
- Had to manually generate the read counts of each metagenome file following QC. Sometimes the read counts file failed to generate because of how I specified the ATLAS rules.
- ATLAS 2.1.4 had an error where the 'genomes/genomes' output folder could not be recognized as input to subsequent post-dereplication steps in the pipeine. I changed the name to 'genome/Genomes' and updated the genome directory specifier in the config file to get around this issue.

CONFIRM:
- I wonder if I actually fixed the read counts issue??
- Were some assembled using a different assembler?
```

### Output summary
If you'd like to see the set of 756 MAGs generated by the ATLAS run, you can find them, along with their gene annotations, at the [Zenodo repo associated with this code repository](https://doi.org/10.5281/zenodo.3930110). Specifically, look for `lake_survey_MAGs.tar.gz` and `lake_survey_MAGs_eggnog_annotations.tar.gz`.

I also included a few raw output files in the `lake_metagenomes/summary` folder for reference and to show how the run stats were calculated (below):
- `read_counts.tsv`: Read counts for metagenomes after QC
- `combined_contig_stats.tsv`: Assembly statistics for each metagenome
- `completeness_checkm.tsv`: CheckM stats for dereplicated metagenome-assembled genome (MAG)
- `raw_counts_genomes.tsv`: reads mapped to each MAG
- `gtdbtk.bac120.summary.tsv` and `gtdbtk.ar122.summary.tsv`: GTDB taxonomy for the bacterial and archaeal MAGs, respectively

### Summary statistics
To generate Supplementary Data 3 with the relative abundances of the MAGs, I ran a simple script to combine various ATLAS output files and calculate relative abundances:

```bash
generate_MAG_table.py \
# TODO
```


## Lake metatranscriptome analysis

### Install ATLAS
Here, we'll use ATLAS commit 59da38f to perform sample QC. This commit was made about two weeks before the full release of ATLAS version 2.2.0

Note that you will need to have pre-installed [miniconda (e.g., miniconda3)](https://docs.conda.io/en/latest/miniconda.html).
```bash
# This step (creating the environment) could take time...
conda create -n atlas_2.2.0-dev1 -c conda-forge -c bioconda python=3.6 snakemake pandas bbmap=37.78 click=7 ruamel.yaml biopython

# Then, clone the ATLAS repo
git clone https://github.com/jmtsuji/atlas.git
cd atlas
git checkout 59da38f
# This commit is also available from the main ATLAS repo, but here I use my fork because it will come in handy in a few steps (see below)

# And install ATLAS via that repo into the conda env. Now, editing that repo will adjust ATLAS itself
conda activate atlas_2.2.0-dev1
pip install --editable .
cd ..
```

One minor typo had to be fixed in the code in order for it to run. I had to edit line 39 of `atlas/rules/cat_taxonomy.smk` from `java_mem` to `mem`. This rule is not even used in this ATLAS run, but it prevents ATLAS from starting due to the typo.

### Analyze metatranscriptome data
A config file for the run and sample ID guide file are included at `lake_metatranscriptomes/config.yaml` and `lake_metatranscriptomes/samples.tsv`

Perform QC on the metatranscriptome reads using ATLAS `59da38f`:
```bash
# Activate the environment by running:
# conda activate atlas_2.2.0-dev1

cd lake_metatranscriptomes
# Notice that samples.tsv and config.yaml files are provided here

atlas run -w . -c config.yaml -j 50 qc --reason 2>&1 | tee atlas_run_QC.log

cd ..
```

Then, switch the atlas repo over to the custom `maprna` branch, commit `96e47df`:
```bash
cd atlas
git checkout 96e47df
cd ..
```

Then, perform RNA read mapping. Note that we will use the same config.yaml and samples.tsv files as the first step. As one important point, one variable in the config file, `genome_dir`, points to the location of the MAGs generated in the lake metagenome analysis above. It has been set in the config file to `../lake_metagenomes/genomes/Genomes`, but you might need to change this if you've put the MAGs in a different directory. I am also not sure if the relative path will work here in the config file; I used an absolute path in reality when running this code, but I changed it to relative here so that filepaths are kept contained within the Github repo.
```bash
# Activate the environment by running:
# conda activate atlas_2.2.0-dev1

cd lake_metatranscriptomes

# Final pre-run prep - the gene annotation files from the metagenome run are also needed for featureCounts to work
mkdir -p genomes
ln -s ../lake_metagenomes/genomes/annotations genomes

atlas run -w . -c config.yaml -j 40 genomes --reason --until summarize_gene_counts 2>&1 | tee atlas_run_mapping.log

cd ..
```
What's going on under the hood here:
- The typical `genomes` module of ATLAS is run (without customization) up until the QC-processed metatranscriptome reads are mapped to the genome bins
- From there, the custom `quantify.smk` snakefile that I wrote is run until the last step, rule `summarize_gene_counts`. Here is specifically what is done:
  1. Rule `combine_gff_annotations`: the GFF-format annotation files from all MAGs (in `genomes/annotations/genes`) are concatenated into a single huge GFF file. This works because all contig and ORF IDs are unique for each MAG generated by ATLAS.
  2. Rule `run_feature_counts`: featureCounts is run a single time using the huge GFF file as input, along with the BAM files (for each mapping to genome bins) for all input metatranscriptomes. The output is a large, raw tab-separated table of the number of RNA read counts to each gene across all MAGs.
  3. Rule `parse_feature_counts_output`: parse the ORF IDs to get the ID of each MAG; clean up the table format.
  4. Rule `split_feature_counts_table`: split the featureCounts table to make a single table for each MAG, in case helpful.
  5. Rule `summarize_gene_counts`: touch an empty output file to signify that all steps are complete.

### Run notes
The above ATLAS code is a streamlined version of the code that was used to analyze the data. Here I note a few important points about the actual data analysis that was performed:

```
- The QC read counts to each metatranscriptome failed to generate during the ATLAS run, like was mentioned above for the metagenome ATLAS run. I counted the number of QC reads for each metatranscriptome after the ATLAS run and then created a `read_counts_QC.tsv` file with format matching the typical report file output by ATLAS
- The `raw_counts_genomes.tsv` file also did not generate because of how I specified the maprna rule. I manually generated that by summarizing the mapped reads against each contig using the BAM files output by ATLAS and then summing those counts for the MAGs those contigs belonged to.
```

### Output summary
I included a few raw output files in the `lake_metatranscriptomes/summary` folder for reference and to show how the run stats were calculated (below):
- `read_counts_QC.tsv`: Read counts for metatranscriptomes after QC
- `raw_counts_genomes.tsv`: number of mapped reads on each MAG
- **One file, `lake_survey_MAGs_featureCounts.tsv.gz`, was put in the [Zenodo repo associated with this code repository](https://doi.org/10.5281/zenodo.3930110) due to its large size**
  - This file contains the raw gene counts for each metatranscriptome against each protein-coding gene in each of the 756 MAGs.

### Summary statistics
To generate Supplementary Data 4 with the relative expressions of the MAGs, I ran a simple script to combine various ATLAS output files and calculate relative abundances:

```bash
generate_MAG_table.py \
# TODO

# CheckM and GTDB files are the same as used for metagenomes
```

I then re-normalized this table afterwards to be normalized by the total read counts to all MAGs instead of to the total number of metatranscriptome reads, as described in the methods in the paper. Then, I performed some simple mean/standard deviation stats on these output tables in Excel to make Supplementary Data 4.

## Generation of Figure 3

### Figure 3a: Bubble plot of relative expression of MAGs
The Jupyter notebook, `data_viz_Fig3a/Fig_3a_plotter.ipynb`, summarizes read mapping stats and generates Fig. 3a. 
Input files are in `data_viz_Fig3a/input_data`.


### Figure 3b (and Extended Data Fig. 7): up/down regulation of genes in the two "_Ca._ Chloroheliales" MAGs
The Jupyter notebook, `data_viz_Fig3b/Fig_3b_plotter.ipynb`, summarizes read count normalization methods and plotting of both figures.  

Input files need to be downloaded (as a `.tar.gz` file) from the Zenodo repo [using this link](https://zenodo.org/record/5131685/files/lake_survey_Ca_Chloroheliales_MAGs_info.tar.gz?download=1). The tarball should be extracted, and the contents should be saved in `data_viz_Fig3b/input_data`. You're then ready to run the Jupyter notebook.

Data from these analyses was also summarized to make Supplementary Data 5.

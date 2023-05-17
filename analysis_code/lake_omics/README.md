# Lake survey metagenome/metatranscriptome analysis
Part of the larger '*Ca.* Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

**NOTE: for each code section provided below, the code ought to be run from within this `lake_omics` directory.**

This README describes how lake metagenome/metatranscriptome survey data were analyzed, including generation of selected figures and supplementary data files.

## Data download
For info on the datasets that were analyzed, see the `source_data` folder, specifically the `README` and the following info files:
- `lake_metagenome_data_accessions_jgi.tsv`
- `lake_metagenome_data_accessions_ncbi.tsv`
- `lake_metatranscriptome_data_accessions.tsv`

## Lake metagenome analysis
In the manuscript, the main purpose of the lake metagenome analyses was to:
- Recover metagenome-assembled genomes (MAGs) from the lakes to search for RCI-encoding relatives of "_Ca_. Chx. allophototropha"
  - Those MAGs will then be used as the basis for RNA read mapping work with metatranscriptome data
- Determine the relative abundances of "_Ca_. Chx. allophototropha" relatives in Boreal Shield lakes
- Detect relative abundances of phototrophic "_Ca_. Chx. allophototropha" relatives at the unassembled read level to cross check MAG-based data

To analyze the metagenomes, I ran them through the [ATLAS](https://github.com/metagenome-atlas/atlas) pipeline and then did 
a few post-run analyses.

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

Note that I had to modify two of the ATLAS code files in order to run the samples:
- `tree.py` had to be modified to work (parses the CheckM tree - I do not actually use the results of this in the paper, I think, but this bug had to be fixed for ATLAS to run.)
  - Specifically, `tree.py` could not find `parsers_checkm.py` for import. This appears to be a relative filepath issue. I did a somewhat hacky fix:
    - Copied the contents of `parsers_checkm.py` to the top of the `tree.py` file so that `parsers_checkm.py` no longer needed to be called.
    - Removed the import statement of `parsers_checkm`
    - Removed the `parsers_checkm` module call (and just directly called the `parsers_checkm` function) in `if __name__ == "__main__"`
- CAT taxonomy is slow -- it would take weeks to taxonomically classify the genome bins from this study. I decided to use GTDB taxonomy instead.
  - Thus, I disabled CAT taxonomy classification within ATLAS
  - I made a single edit to the `Snakefile` by simply commenting out `"genomes/taxonomy/taxonomy.tsv",` under `rule genomes`.

(Both of these issues were addressed in future releases of ATLAS.)

### Process the metagenome data
A config file and sampel ID file for the run have already been created at `lake_metagenomes/config.yaml` and `lake_metagenomes/samples.tsv`.   
If you want to create config files for yourself, you can use the `atlas init` command as documented in the ATLAS repo.  

Notes:
- make sure you modify the `config.yaml` file in the `lake_metagenomes` folder so that `database_dir` is a real directory on your machine.
  - You will also want to check that the memory and thread counts match what you hope to use on your server.
- you might also have to modify the filepaths for the samples in the `samples.tsv` file so that they are correct on your system.

Start the ATLAS run
```bash
# Activate the environment by running:
# conda activate atlas_2.1.4

cd lake_metagenomes
# Notice that samples.tsv and config.yaml files are provided here

atlas run -w . -c config.yaml -j 50 all --reason --keep-going 2>&1 | tee atlas_run_part1.log

cd ..
```
This might take several weeks of computational runtime and considerable (e.g., 200 GB) RAM!! 
If you just want the genome bins, you can skip running ATLAS yourself and just directly download the bins from 
the [Zenodo repo associated with this code repository](https://doi.org/10.5281/zenodo.3930110). Specifically, look for `lake_survey_MAGs.tar.gz`.

Note that the above ATLAS run will get you close to the end of the pipeline but will run into two errors and stop. 
You then have to work around those errors and keep going. (These errors were fixed in subsequent versions of ATLAS.)

#### Fix 1: move the genomes directory
```
# Issue is documented at https://github.com/metagenome-atlas/atlas/issues/231
# ATLAS summarizes the genome bins (after dereplication) in genomes/genomes, but then it cannot recognize genomes/genomes for downstream steps
# Thus you need to rename the genomes/genomes folder, specify the new name as a config parameter, and keep going

cd lake_metagenomes
mv genomes/genomes genomes/Genomes

atlas run -w . -c config.yaml -j 50 genomes --reason --config genome_dir=genomes/Genomes 2>&1 | tee atlas_run_part2.log

# genomes module is now done. Continue to the genecatalog module
# Note that I changed to 20 threads here because there were other jobs running on our server, but you do not need to do this. (I had to temporarily change the thread requirements in the config file to 20 threads as well.)
atlas run -w . -c config.yaml -j 20 genecatalog --reason 2>&1 | tee atlas_run_part3.log
# you'll then run into the error fixed by "Fix 2" below.

cd ..
```

#### Fix 2: Work around an issue in Genecatalog's `rule cluster_genes` (and related rules).  
The rule looks for an output file that does not exist. 
Rather than fixing the snakemake code, I instead just ran the exact commands run by ATLAS manually on the command line, then kept going.
```
# Activate the MMSeqs conda env created by Snakemake. Yours might be in a different location.
# You can find it by going to the `conda_envs` folder within your ATLAS database directory and searching for which of the .yaml files mentions MMSeqs
# Then, to get the path of that conda env, just use the path to the .yaml file without the .yaml suffix.
# i.e., in my case, the correct environment was specified by /Data/reference_databases/atlas/2.1.x/conda_envs/d5bf8789.yaml, so I ran:
conda activate /Data/reference_databases/atlas/2.1.x/conda_envs/d5bf8789

cd lake_metagenomes

### Step 1: manually finish `rule cluster_genes`
mmseqs createdb Genecatalog/all_genes/predicted_genes.faa Genecatalog/all_genes/predicted_genes.db \
  > >(tee  logs/Genecatalog/clustering/cluster_proteins.log)

mkdir -p atlas_tmp/mmseqs

mmseqs linclust -c 0.9 \
  --min-seq-id 0.99 \
  --threads 20 Genecatalog/all_genes/predicted_genes.db Genecatalog/clustering/protein_clusters.db \
  atlas_tmp \
  > >(tee -a  logs/Genecatalog/clustering/cluster_proteins.log)

rm -fr  atlas_tmp \
  > >(tee -a  logs/Genecatalog/clustering/cluster_proteins.log)


### Step 2: manually finish `rule get_rep_proteins`
mmseqs createtsv Genecatalog/all_genes/predicted_genes.db \
  Genecatalog/all_genes/predicted_genes.db \
  Genecatalog/clustering/protein_clusters.db \
  Genecatalog/orf2gene_oldnames.tsv  > >(tee   logs/Genecatalog/clustering/get_rep_proteins.log)

mmseqs result2repseq Genecatalog/all_genes/predicted_genes.db Genecatalog/clustering/protein_clusters.db \
  Genecatalog/protein_catalog.db  > >(tee -a  logs/Genecatalog/clustering/get_rep_proteins.log)

mmseqs result2flat Genecatalog/all_genes/predicted_genes.db Genecatalog/all_genes/predicted_genes.db \
  Genecatalog/protein_catalog.db Genecatalog/representatives_of_clusters.fasta  \
  > >(tee -a  logs/Genecatalog/clustering/get_rep_proteins.log)


### Step 3: manually finish `rule rename_protein_catalog`
# Switch back to the main conda env for ATLAS
conda activate atlas_2.1.4

## Then start python
python

#####################
## Within python, run:
import pandas as pd

# From utils.py
def gen_names_for_range(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]

gene2proteins= pd.read_csv("Genecatalog/orf2gene_oldnames.tsv", index_col=1, header=None,sep='\t')
protein_clusters_old_names= gene2proteins[0].unique()
map_names = dict(zip(protein_clusters_old_names, gen_names_for_range(len(protein_clusters_old_names),'Gene')))
gene2proteins['Gene'] = gene2proteins[0].map(map_names)
gene2proteins.index.name='ORF'
gene2proteins['Gene'].to_csv("Genecatalog/clustering/orf2gene.tsv",sep='\t',header=True)
exit()
#####################

### Step 4: Final post-run cleanup
rm Genecatalog/orf2gene_oldnames.tsv
rm Genecatalog/all_genes/predicted_genes.db*
rm Genecatalog/clustering/protein_clusters.db.*
rm Genecatalog/protein_catalog.db.*

### Step 5: Resume ATLAS and finish the analysis
atlas run -w . -c config.yaml -j 20 genecatalog --reason 2>&1 | tee atlas_run_part4.log

# Final cleanup
rm Genecatalog/representatives_of_clusters.fasta
```
All done! (Whew.) You should have a complete ATLAS run now.

#### Taxonomy classification
Run the GTDB taxonomy classifier on the samples after completing the ATLAS run:

I used GTDBTk version `0.3.2` with database release 89, which can be installed as follows:
```bash
# Install the GTDB classifier via conda
# Using GTDBTk `v0.3.2`, database release 89
# Note that the DB download will take some time!

conda create -n gtdbtk_0.3.2 -c bioconda gtdbtk=0.3.2
conda activate gtdbtk_0.3.2

# Download the GTDBTk database, release 89
download-db.sh
```
Note that the GTDBTk version `0.3.2` and/or the DB release 89 might no longer be active due to rapid development of the GTDBTK, however. 
For future analyses, I recommend installing the latest version of the GTDBTk with the newest database.

Run the classifier:
```bash
atlas_dir="lake_metagenomes"
genome_dir="${atlas_dir}/genomes/Genomes"
out_dir="${atlas_dir}/genomes/taxonomy_gtdbtk"
genome_extension=fasta
output_prefix=gtdbtk
min_percent_aa=10 # 10 is the default
threads=40

# Activate the conda env by running:
# conda activate gtdbtk_0.3.2

mkdir -p "${out_dir}"
gtdbtk classify_wf \
   --genome_dir "${genome_dir}" \
   --out_dir "${out_dir}" \
   -x ${genome_extension} \
   --min_perc_aa ${min_percent_aa} \
   --prefix ${output_prefix} \
   --cpus ${threads}
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
To generate Supplementary Data 3 with the relative abundances of the MAGs, I ran a simple personal script to combine various ATLAS output files and calculate relative abundances:

```bash
cd lake_metagenomes/summary`
git clone https://github.com/jmtsuji/atlas2-helpers
cd atlas2-helpers
git checkout f461ba3
cd ..

conda create pandas python=3.6 pandas=0.24
conda activate pandas

# Normalized to assembled reads, using GTDB taxonomy
atlas2-helpers/scripts/generate_MAG_table.py \
  -o MAG_table_DNA_to_assembled.tsv \
  -g raw_counts_genomes.tsv \
  -R combined_contig_stats.tsv \
  -t gtdbtk.bac120.summary.tsv gtdbtk.ar122.summary.tsv \
  -c completeness_checkm.tsv \
  2>&1 | tee MAG_table_DNA_to_assembled.tsv
  
# Normalized to total metagenome reads (not used in this paper, but provided here for reference)
atlas2-helpers/scripts/generate_MAG_table.py \
  -o MAG_table_DNA_to_unassembled.tsv \
  -g raw_counts_genomes.tsv \
  -r read_counts.tsv \
  -t gtdbtk.bac120.summary.tsv gtdbtk.ar122.summary.tsv \
  -c completeness_checkm.tsv \
  2>&1 | tee MAG_table_DNA_to_unassembled.tsv

cd ../..
```
`MAG_table_DNA_to_assembled.tsv` was then opened in Excel, and the visual appearance of the table was cleaned up to make Supplementary Data 3. 
This table nicely summarizes the taxonomy of all MAGs and their percent relative abundances in all metagenome samples.

### Unassembled read search for "_Ca_. Chx. allophototropha"-like RCI genes
As a cross-comparison to the MAG-based data above, the QC'ed metagenome data were searched at the unassembled read level for the "_Ca_. Chx. allophototropha" like RCI gene.

#### Predict short amino acids
Install software - here, I create a custom conda environment for a modified version of FragGeneScanPlusPlus
```bash
# Install bbmap and build dependencies for FGSpp
conda create -n FGSpp_471fdf7 -c conda-forge -c bioconda bbmap=38.91 meson ninja
conda activate FGSpp_471fdf7
# bbmap version 38.91

# Now install FGSpp
git clone https://github.com/LeeBergstrand/FragGeneScanPlusPlus.git
cd FragGeneScanPlusPlus
git checkout 471fdf7 # after fixing the segmentation violation
meson build
ninja -C build

cp -r train build/FGS++ "${CONDA_PREFIX}/bin"

cd ..
rm -rf FragGeneScanPlusPlus # if desired

conda deactivate
```

Install personal code wrapper for FGS++/bbmap, `predict_short_orfs.sh `, v1.3.4
```bash
git clone https://github.com/jmtsuji/basic-sequence-analysis.git
cd basic-sequence-analysis
git checkout v1.3.4

conda activate FGSpp_471fdf7

cp -r scripts/* "${CONDA_PREFIX}/bin"
cd .. && rm -rf basic-sequence-analysis

conda deactivate
```

Gather the QC'ed reads
```bash
mkdir -p "unassembled_read_search" && cd "$_"
mkdir -p 01_metagenome_R1_files && cd "$_"

# TODO - get atlas dir name
find "../../lake_metagenomes" -maxdepth 3 -mindepth 3 -type f -name "*_QC_R1.fastq.gz" | \
  xargs -I {} ln -s {} .

find "../../lake_metagenomes" -maxdepth 3 -mindepth 3 -type f -name "*_QC_R1.fastq.gz" | \
  sort -h | \
  sed 's\../../lake_metagenomes/.*/sequence_quality_control/\\g' | \
  sed 's/_QC_R1.fastq.gz//g' \
  > samples.list

cd ..
```
Note that the IDs here might be a bit different than in the MAG table, e.g., underscores are used instead of dashes, and depth numbers might have leading zeros.

Predict amino acids - this could take time, so it is good to run inside a `screen`:
```bash
mkdir -p "02_orfs" && cd "$_"
mkdir -p logs

conda activate FGSpp_471fdf7

r1_filepaths=($(find -L "../01_metagenome_R1_files" -type f -name "*_QC_R1.fastq.gz" | sort -h))

for r1_filepath in ${r1_filepaths[@]}; do

  sample_basename="${r1_filepath##*/}"
  sample_basename="${sample_basename%_QC_R1.fastq.gz}"
  echo "[ $(date -u) ]: ${sample_basename}"
  logfile="logs/${sample_basename}.log"
  
  predict_short_orfs.sh -p 50 -t illumina_10 -l "${logfile}" "${r1_filepath}" > "${sample_basename}.faa"
  
  if [ $? != 0 ]; then
    echo "[ $(date -u) ]: ${sample_basename}: run FAILED. Moving on... but keeping output"
  fi

done

echo "[ $(date -u) ]: Done."

cd ..
```

#### Scan for functional genes

Download HMMs of interest
```bash
mkdir -p "03_hmms" && cd "$_"

### Download HMMs of interest
curl -LOJ http://fungene.cme.msu.edu/hmm_download.spr?hmm_id=31 # rpoB.hmm - June 2009, FunGene
# v2.0.1 of Github repo
wget https://raw.githubusercontent.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/v2.0.1/genome_bin_analysis/hidden_markov_models/fmoA_Chloroheliales.hmm
wget https://raw.githubusercontent.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/v2.0.1/genome_bin_analysis/hidden_markov_models/pscA_Chloroheliales.hmm

# Lengths (based on the LEN flag at the top of the flags)
# fmoA:  373
# pscA:  596
# rpoB: 2842

cd ..
```

Scan using hmmsearch
```bash
mkdir -p 04_hmmsearch && cd "$_"
conda activate general # install hmmsearch 3.3.2 in this env

# Variables
hmmsearch_outfile=hmmsearch.tsv
evalue="1e-10"

##### FUNCTION #####
hmmsearch_beautify() {
  faa_file=$1
  hmm_file=$2
  evalue=$3
  cpu=4
  
  sample_basename="${faa_file##*/}"
  sample_basename="${sample_basename%.faa}"

  echo "[ $(date -u) ]: '${sample_basename}.faa': starting search using HMM '${hmm_file##*/}'" >&2
  
  # hmmsearch 3.3.2
  hmmsearch -o /dev/null --tblout /dev/stdout -E "${evalue}" --cpu ${cpu} "${hmm_file}" "${faa_file}" | \
    grep -v "^#" | \
    tr -s " " "\t" | \
    cut -f 1,3,5,6,7 | \
    sed "s/^/${sample_basename}\t/g"
    
  echo "[ $(date -u) ]: '${sample_basename}.faa': finished search using HMM '${hmm_file##*/}'" >&2
}
export -f hmmsearch_beautify
#####################

find "../02_orfs" -type f -name "*.faa" | sort -h > faa_files.list
find "../03_hmms" -type f -name "*.hmm" | sort -h > hmm_files.list

# This command is super cool! It is the equivalent of a nested for loop.
printf "sample\tquery\thmm\tevalue\tscore\tbias\n" > "${hmmsearch_outfile}"
parallel -j 12 hmmsearch_beautify {1} {2} ${evalue} :::: faa_files.list hmm_files.list :::  \
  >> "${hmmsearch_outfile}" 2> hmmsearch.log

# Check totals
# cat hmmsearch.tsv | cut -f 1,3 | tail -n +2 | uniq -c | column -t | sort -h

# Made hmm lengths table
printf "hmm\tlength\n" > hmm_lengths.tsv
printf "fmoA_Chx\t373\n" >> hmm_lengths.tsv
printf "pscA_Chx\t596\n" >> hmm_lengths.tsv
printf "rpoB\t2842\n" >> hmm_lengths.tsv
```

To filter to reliable "_Ca_. Chx. allophototropha"-like RCI genes, run BLASTP. I will check fmoA at the same time:
```bash
blastp_check() {
  faa_file=$1
  hmmsearch_table=$2
  query_id=$3
  query_file=$4

  sample_id="${faa_file##*/}"
  sample_id="${sample_id%.faa}"
  
  echo "[ $(date -u) ]: '${sample_id}.faa': '${query_id}': starting search." >&2
  
  grep "${query_id}" "${hmmsearch_table}" | grep ${sample_id} | cut -f 2 > ${sample_id}_${query_id}.list

  if [[ $(cat ${sample_id}_${query_id}.list | wc -l) -ge 1 ]]; then
  
    # seqtk 1.3-r106
    seqtk subseq "${faa_file}" ${sample_id}_${query_id}.list > ${sample_id}_${query_id}.faa
    rm ${sample_id}_${query_id}.list

    # blastp 2.10.1+ (to stdout)
    blastp -query "${query_file}" -subject ${sample_id}_${query_id}.faa \
      -outfmt "6 qseqid sseqid pident evalue bitscore qcovhsp qstart qend qlen sstart send slen" | \
      sed "s/^/${query_id}\t/g" | \
      sed "s/^/${sample_id}\t/g"

    rm ${sample_id}_${query_id}.faa
    echo "[ $(date -u) ]: '${sample_id}.faa': '${query_id}': finished search." >&2
  
  else

    echo "[ $(date -u) ]: '${sample_id}.faa': '${query_id}': no hits to search; finished." >&2
    rm ${sample_id}_${query_id}.list

  fi
  
}
export -f blastp_check

# Make sure you run the following commands in a conda env where seqtk 1.3-r106 and blast 2.10.1 are installed.

# Github repo v2.0.1
wget https://raw.githubusercontent.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/v2.0.1/genome_bin_analysis/sequences_of_interest/pscA.faa
wget https://raw.githubusercontent.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/v2.0.1/genome_bin_analysis/sequences_of_interest/fmoA.faa

printf "sample\thmm\tqseqid\tsseqid\tpident\tevalue\tbitscore\tqcovhsp\tqstart\tqend\tqlen\tsstart\tsend\tslen\n" > pscA_fmoA_blastp.tsv
parallel -j 50 blastp_check {1} hmmsearch.tsv {2}_Chx {2}.faa :::: faa_files.list ::: pscA fmoA \
  >> pscA_fmoA_blastp.tsv 2> pscA_fmoA_blastp.log
```
Used this to filter based on pident / evalue in Jupyter lab, in `summarize_hmmsearch_results.ipynb`

Made figure: `pscA_abundances.pdf`. This will be modified to make Figure 4b.

Two TSV-format tables are also included here that summarize some raw data:
- `counts-per-sample-blast-filtered.tsv`: raw counts of pscA and fmoA after BLAST filtration for each metagenome
- `percent-per-lake-blast-filtered-length-normalized-max.tsv`: percent relative abundances following normalization to rpoB hits and to HMM length. Max. percent value for each lake is shown. These data should be the same as shown in the Figure 4b plot.

Done! On to metatranscriptome work.

## Lake metatranscriptome analysis
Lake metatranscriptome data was mapped to the set of MAGs generated from metagenome data (above) to determine the relative gene expression 
of "_Ca_. Chx. allophototropha" relatives in Boreal Shield lakes.

### Install ATLAS
Here, we'll use ATLAS commit 59da38f to perform sample QC. This commit was made about two weeks before the full release of ATLAS version 2.2.0.

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
A config file for the run and sample ID guide file are included at `lake_metatranscriptomes/config.yaml` and `lake_metatranscriptomes/samples.tsv`. 
Like mentioned for metagenome data, you will need to set the database directory path, raw read filepaths, and thread/memory requirements to match your server setup.

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

One more typo needs to be fixed in the code in order for it to run. 
There is a hanging bracket that needs to be replaced with a comma in line 239 of `atlas/Snakefile` 
(i.e., `"feature_counts":"genomes/expression/gene_counts.tsv"}` should be changed to `"feature_counts":"genomes/expression/gene_counts.tsv",`). 
This issue is fixed in the subsequent commit `001e417`, which you are welcome to use in place of `96e47df` if you'd like.

Then, perform RNA read mapping. Note that we will use the same config.yaml and samples.tsv files as the first step. 
As one important point, one variable in the config file, `genome_dir`, points to the location of the MAGs generated in the lake metagenome analysis 
above. It has been set in the config file to `../lake_metagenomes/genomes/Genomes`, but you might need to change this if you've put the MAGs in 
a different directory. I am also not sure if the relative path will work here in the config file; I used an absolute path in reality when running 
this code, but I changed it to relative here so that filepaths are kept contained within the Github repo.
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

### Post-run analysis
In order to calculate read mapping statistics needed for Figure 3a, I had to do a couple analyses outside of ATLAS.

Firstly, in my case, the `read_counts.tsv` file failed to generate during the ATLAS run. This was a simple mistake on my end and should not happen in your case if you run the code as presented above. 
I had to manually generate that file by counting the number of entries in the QC-processed metatranscriptome files. 
I summarized these counts in `read_counts_QC.tsv` (provided in the `summary` folder), which mimics the format of the real `read_counts.tsv` file.

Secondly, the `raw_counts_genomes.tsv` file is not auto-generated as part of the maprna rule. 
I manually generated that file by counting the mapped reads against each contig, based on the BAM files output by ATLAS. 
I then summed those counts for each MAG. The code to do this analysis is actually fairly simple:
```bash
cd lake_metatranscriptomes/genomes/alignments
mkdir -p summarize_genome_counts

printf "count\tcontig-id\trep-id\n" > summarize_genome_counts/L221_L304_RNA_read_counts_contigs.tsv

alignment_files=(Jun2018_L221_05m_A.bam Jun2018_L221_05m_B.bam Jun2018_L221_05m_C.bam 
Jun2018_L304_06m_A.bam Jun2018_L304_06m_B.bam Jun2018_L304_06m_C.bam
Sep2017_L227_04m_A.bam Sep2017_L227_04m_B.bam Sep2017_L227_04m_C.bam)

for alignment_file in ${alignment_files[@]}; do
  echo "${alignment_file}"
  
  samtools view "${alignment_file}" | \
    cut -f 3 | \
    uniq -c | \
    sed "s/^ \+//g" | \
    sed "s/$/ ${alignment_file%.bam}/g" | \
    tr -s ' ' $'\t' \
    > summarize_genome_counts/${alignment_file%.bam}_contig_counts.tsv
    
  cat summarize_genome_counts/${alignment_file%.bam}_contig_counts.tsv >> summarize_genome_counts/read_counts_contigs.tsv
    
done
```

The output file `read_counts_contigs.tsv` has a format like:
```
count	contig-id	rep-id
538	MAG001_40	Jun2018_L221_05m_A
2	MAG001_73	Jun2018_L221_05m_A
12	MAG001_121	Jun2018_L221_05m_A
```
Where `count` is the number of mapped reads to a contig. Thankfully, ATLAS labels all contigs by the MAG ID, then an underscore, then the contig ID. 
This makes it easy to convert this file into total RNA counts per genome. 
I implement this conversion in the iPython notebook `generate_raw_counts_genomes.ipynb`, which is included in the `summary` folder along with `read_counts_contigs.tsv.gz` (gzipped for space savings).

(Note: in reality the L227 samples were processed independently from the 221/304 samples here, and the resulting `raw_counts_genomes.tsv` tables 
were merged, but everything is summarized a single code block here for simplicity (should not change the result).)

With these files in place, you are now ready for all downstream statistics work.

### ATLAS output summary
I included key ATLAS output files and post-analysis files in the `lake_metatranscriptomes/summary` folder, namely:
- `read_counts_QC.tsv`: Read counts for metatranscriptomes after QC
- `raw_counts_genomes.tsv`: number of mapped reads on each MAG
- One file, `lake_survey_MAGs_featureCounts.tsv.gz`, was put in the [Zenodo repo associated with this code repository](https://doi.org/10.5281/zenodo.3930110) due to its large size
  - This file contains the raw gene counts for each metatranscriptome against each protein-coding gene in each of the 756 MAGs.

### Summary statistics
To generate Supplementary Data 4 with the relative expressions of the MAGs, I ran a simple, personal script to combine various ATLAS 
output files and calculate relative abundances of MAGs based on RNA-seq data:

```bash
cd lake_metatranscriptomes/summary`
git clone https://github.com/jmtsuji/atlas2-helpers
cd atlas2-helpers
git checkout 1e08f0a # version 0.1.0
cd ..

conda create pandas python=3.6 pandas=0.24
conda activate pandas

# CheckM and GTDB files are the same as used for metagenomes
atlas2-helpers/scripts/generate_MAG_table.py \
  -o MAG_table_RNA_to_unassembled.tsv \
  -g raw_counts_genomes.tsv \
  -r read_counts_QC.tsv \
  -t ../../lake_metagenomes/summary/gtdbtk.bac120.summary.tsv ../../lake_metagenomes/summary/gtdbtk.ar122.summary.tsv \
  -c ../../lake_metagenomes/summary/completeness_checkm.tsv \
  2>&1 | tee MAG_table_RNA_to_unassembled.tsv
  
cd ../..
```
Note: in reality the L227 samples were done independently here, and the resulting MAG tables were merged, but this is shown as a single 
command for simplicity (should not change the result).

The resulting `MAG_table_RNA_to_unassembled.tsv` became the basis for Supplementary Data 4. 
Like mentioned in the paper, I re-normalized `MAG_table_RNA_to_unassembled.tsv` in Excel, such that MAG relative abundances 
summed to 100% for each sample (i.e., were normalized based on the total read counts to all MAGs instead of being 
normalized based on the total number of metatranscriptome reads for each sample).
I also performed some simple mean/standard deviation stats on these output tables in Excel when making Supplementary Data 4.

## Generation of selected panels of Figure 4

### Figure 4b: relative abundances of "_Ca_. Chx. allophototropha"-like RCI genes based on unassembled reads
See description in the lake metagenomes section above (at the end of the section). Calculations are shown in a Jupyter notebook in the `unassembled_read_search` folder, along with intermediate data.

### Figure 4d: Bubble plot of relative expression of MAGs
The Jupyter notebook, `data_viz_Fig4d/Fig_3d_plotter.ipynb`, summarizes read mapping stats and generates Fig. 3d. 
Input files are in `data_viz_Fig4d/input_data`.

After production, the resulting raw PDF `Figure_4d_raw.pdf` was edited in Inkscape (e.g., to clean up axis labels and modify the bubble colour scheme).

### Figure 4e (and Extended Data Fig. 7): up/down regulation of genes in the two "_Ca._ Chloroheliales" MAGs
The Jupyter notebook, `data_viz_Fig4e/Fig_3e_plotter.ipynb`, summarizes read count normalization methods and plotting of both figures.  

Input files need to be downloaded (as a `.tar.gz` file) from the Zenodo repo [using this link](https://zenodo.org/record/5131685/files/lake_survey_Ca_Chloroheliales_MAGs_info.tar.gz?download=1). 
The tarball should be extracted, and the contents should be saved in `data_viz_Fig4e/input_data`. You're then ready to run the Jupyter notebook.

The resulting raw PDFs, `Figure_4e_raw.pdf` and `Figure_ED7_raw.pdf`, are in the folder. These received minor edits via Inkscape to finalize.

Data from these analyses were also summarized to make Supplementary Data 5.

Done!

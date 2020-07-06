# Metagenome analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020

**NOTE: for each code section provided below, the code ought to be run from within this `metagenome_analysis` directory.**

## Data download
You can download the data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="metagenome_data_accessions.tsv"
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

You now have FastQ files for the three datasets.

## Install ATLAS
[ATLAS](https://github.com/metagenome-atlas/atlas) is a metagenome QC, assembly, and binning workflow.  

Here, we'll install and use version 2.2.0.

Note that you will need to have pre-installed [miniconda (e.g., miniconda3)](https://docs.conda.io/en/latest/miniconda.html).

```bash
# This step (creating the environment) could take time...
conda create -n atlas_2.2.0 -c bioconda -c conda-forge metagenome-atlas=2.2.0
```

Note that one bug has to be edited before using the code, as described [here](https://github.com/metagenome-atlas/atlas/issues/296):
```bash
conda activate atlas_2.2.0

# Find binning.snakefile
find $CONDA_PREFIX -name "binning.snakefile"
# E.g., in ${CONDA_PREFIX}/lib/python3.6/site-packages/atlas/rules/binning.snakefile

# **MANUALLY ENTER THE FILE LOCATION HERE IN PLACE OF THE EXAMPLE PATH**
binning_snakefile_path="${CONDA_PREFIX}/lib/python3.6/site-packages/atlas/rules/binning.snakefile"
# *******

# Then add the header=False flag to one line (should be line 262):
cp "${binning_snakefile_path}" "${binning_snakefile_path}.backup"
sed -i '/new_d\.to_csv/c\        new_d.to_csv(output[0],sep='\t',header=False)' "${binning_snakefile_path}"
```

Done. Before running ATLAS, make sure to activate the environment by running `conda activate atlas_2.2.0`.


## Process the metagenome data
If it is your first time using ATLAS, note that a large amount of database files will be auto-downloaded during the run.


### Ca. Chloroheliales bin L227-5C
Config files for the run have already been created. 
If you want to create config files for yourself, you can use the `atlas init` command as documented in the ATLAS repo.  

Notes:
- make sure you modify the config.yaml file in the `atlas_L227_5C` folder so that `database_dir` is a real directory on your machine.
- you might have to modify the filepaths for the samples in the `samples.tsv` file so that they are correct on your system.

Run the sample
```bash
# Activate the environment by running:
# conda activate atlas_2.2.0

cd atlas_L227_5C
# Notice that samples.tsv and config.yaml files are provided here

atlas run -w . -c config.yaml -j 50 all --reason 2>&1 | tee atlas_run.log

cd ..
```

### Ca. Chx. allophototropha
Config files are also provided.

This case is slightly more complicated, though. 
We will use the Tell-Read QC processor and Tell-Link assembler to assemble the data ahead of time and then put these data into ATLAS partway through.

#### Assemble read cloud data separately
Download the raw MiSeq output (not yet demultiplexed into FastQ files) from the Zenodo repo corresponding to this publication:
```bash
zenodo_url="https://zenodo.org/record/3930111/files/Capt_S15_sequencer_data_raw.tar.gz"
wget -O - "${zenodo_url}" | tar -xvzf -
# Creates a directory called `Capt_S15_sequencer_data_raw`
```

De-multiplex and perform QC on the raw read cloud sequencer data, then assemble.  
First, you must **manually install Tell-Read and Tell-Link from the [Universal Sequencing website](https://www.universalsequencing.com/analysis-tools)**. 
Here, version 0.9.7 was used for Tell-Read and version 1.0.0 was used for Tell-Link. Then process the data:
```bash
# Demutiplex using Tell-Read
data_dir="Capt_S15_sequencer_data_raw"
# No underscores in output_dir name!!
output_dir="CaptS15TellRead"
# **MANUALLY insert the path to the Tell-Read script here:**
script_path="tellread-release/run_tellread.sh"

"${script_path}" \
  -i "${data_dir}" \
  -o "${output_dir}" \
  -s "T502" \
  -g "NONE" \
  2>&1 | tee CaptS15TellRead.log


# Assemble using Tell-Link
# This was performed by staff at Universal Sequencing Technology with code like:
input_dir="CaptS15TellRead/Full"
output_dir="CaptS15TellLink"
global_kmer=65
local_kmer=35

script_path="tellink-release/run_tellink.sh"

"${script_path}" \
  -r1 "${input_dir}/CaptS15TellRead_R1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz" \
  -r2 "${input_dir}/CaptS15TellRead_R2_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz" \
  -i1 "${input_dir}/CaptS15TellRead_I1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz" \
  -k "${global_kmer}"
  -lc "${local_kmer}"
  -o "${output_dir}" \
  -p 502 \
  2>&1 | tee CaptS15TellLink.log

# In reality, several different global/local kmer lengths were tested iteratively to find the combination yielding the longest scaffolds
# Use the output file `CaptS15TellLink/502/scaffolds.full.fasta`
```

If you want to skip all this, you can just download the Tell-Link scaffolds directly from the same Zenodo data repo:
```bash
zenodo_url="https://zenodo.org/record/3930111/files/scaffold.full.fasta.gz"
wget -O - "${zenodo_url}" | gunzip > scaffolds.full.fasta
```

#### Run ATLAS
Notes:
- make sure you modify the config.yaml file in the `atlas_L227_5C` folder so that `database_dir` is a real directory on your machine.
- you might have to modify the filepaths for the samples in the `samples.tsv` file so that they are correct on your system.

```bash
# Activate the environment by running: 
# conda activate atlas_2.2.0

cd "atlas_Chx_allophototropha"

# Finish the QC step
atlas run -w . -c config.yaml -j 50 all --reason -p \
  --until finished_QC \
  2>&1 | tee atlas_run_Qc.log

# Run the assembler
atlas run -w . -c config.yaml -j 50 all --reason -p \
  --until run_spades \
  2>&1 | tee atlas_run_assembly.log

# Manually swap out the SPAdes scaffolds with the ones from Tell-Link for the S15 sample:
mv Capt_S15/assembly/scaffolds.fasta Capt_S15/assembly/scaffolds.fasta.backup
mv Capt_S15/assembly/contigs.fasta Capt_S15/assembly/contigs.fasta.backup
ln ../CaptS15TellLink/502/scaffolds.full.fasta Capt_S15/assembly/scaffolds.fasta 
ln ../CaptS15TellLink/502/scaffolds.full.fasta Capt_S15/assembly/contigs.fasta

# Now finish the run
atlas run -w . -c config.yaml -j 50 all --reason -p \
  2>&1 | tee atlas_run_rest.log

cd ..
```

You now have a set of uncurated genome bins for both enrichments under the `genomes/genomes` folder within each ATLAS run folder.


## Genome bin relative abundances
```bash
# Installed the MAG table generator
git clone https://github.com/jmtsuji/atlas2-helpers.git
cd atlas2-helpers
git checkout 1e08f0a
PATH=${PATH}:${PWD}/scripts

cd ..

# Create and activate a simple conda env containing pandas
conda create -y -n pandas -c anaconda pandas
conda activate pandas

# Generated MAG tables based on assembled read counts
source_dir_1="atlas_Chx_allophototropha"
source_dir_2="atlas_L227_5C"
output_dir="relative_abundances"
output_filepath_1="${output_dir}/Capt_MAG_table_to_assembled.tsv"
output_filepath_2="${output_dir}/L227_5C_MAG_table_to_assembled.tsv"

mkdir -p "${output_dir}"

generate_MAG_table.py -o "${output_filepath_1}" \
  -a "${source_dir_1}" \
  -R "${source_dir_1}/stats/combined_contig_stats.tsv" \
  -t "${source_dir_1}/genomes/taxonomy/gtdb/gtdbtk.bac120.summary.tsv" \
  2>&1 | tee "${output_filepath_1%.tsv}.log"

generate_MAG_table.py -o "${output_filepath_2}" \
  -a "${source_dir_2}" \
  -R "${source_dir_2}/stats/combined_contig_stats.tsv" \
  -t "${source_dir_2}/genomes/taxonomy/gtdb/gtdbtk.bac120.summary.tsv" \
  2>&1 | tee "${output_filepath_2%.tsv}.log"
```
The output MAG abundance tables are provided in that folder for reference and are used to generate Extended Data Fig. 1 (see below).

## Curation of the genome bins
I identified the two '_Ca_. Chloroheliales' genome bins from the above ATLAS run folders by looking for the bins classified to the _Chloroflexota_ phylum 
based on the `genomes/taxonomy/gtdb/gtdbtk.bac120.summary.tsv` file in each folder.

Then, I manually removed suspect scaffolds based on the methods described in the Supplementary Materials of the paper.

Curated genomes are available from NCBI (see the `Ca_Chloroheliales_genome_analysis` folder in this repo for more details).


## Unassembled read-based analyses
Install FragGeneScanPlusPlus:
```bash
conda create -n FGSpp_Aug2019_9a203d8 -c bioconda bbmap=38.75
conda activate FGSpp_9a203d8

# Make a share folder
mkdir -p "${CONDA_PREFIX}/share"
cd "${CONDA_PREFIX}/share"

# Add the repo
git clone https://github.com/unipept/FragGeneScanPlusPlus.git
cd FragGeneScanPlusPlus
# commit 9a203d8 (Aug. 22, 2019)
git checkout 9a203d8

make

# Final binary is `FGSpp` in this folder
# Add the repo to the PATH
mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
if [[ ! -f "${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh" ]]; then
echo '#!/bin/sh' > "${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh"
fi
echo "export PATH=\${PATH}:${CONDA_PREFIX}/share/FragGeneScanPlusPlus" >> \
"${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh"
chmod 755 "${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh"
# Now restart the env

conda activate FGSpp_9a203d8
```

Predicted ORFs directly from short read data:
```bash
# Activate the conda env
# conda activate FGSpp_9a203d8

set -euo pipefail # So the code will fail if any of the FGS jobs have an error; code can be a bit sensitive.

work_dir="unassembled_read_analysis"
source_dir=downloads
output_dir="${work_dir}/faa_files"
logfile="${output_dir}/FGSpp.log"
threads=20
chunk_size=100 # number of sequences in a chunk; increase for faster run but more RAM requirement

mkdir -p "${output_dir}"
cd "${output_dir}"

echo "[ $(date -u) ]: Predicting amino acid sequences of short read FastQ files" | tee "${logfile}"

# Auto-set the location of the training files (same directory as the FGSpp binary)
if [[ ! -d train ]]; then
  echo "[ $(date -u) ]: Linking the training file directory" | tee -a "${logfile}"
  train_file_dir="$(which FGSpp)"
  train_file_dir="${train_file_dir%FGSpp}train"
  ln -s "${train_file_dir}" .
fi


fastq_filepaths=($(find -L "${source_dir}" -type f -name "*_R1.fastq.gz" | sort -h))
echo "[ $(date -u) ]: Scanning ${#fastq_filepaths[@]} FastQ files:" | tee -a "${logfile}"

for fastq_filepath in ${fastq_filepaths[@]}; do
  output_filepath="${output_dir}/${fastq_filepath##*/}"
  output_filepath="${output_filepath%.fastq.gz}.frag.faa"

  echo "[ $(date -u) ]: ${fastq_filepath##*/}" | tee -a "${logfile}"

  reformat.sh in="${fastq_filepath}" out=stdout.fa t=${threads} fastawrap=0 trimreaddescription=t 2>/dev/null | \
    FGSpp -s stdin -o stdout -w 0 -t illumina_10 -c ${chunk_size} -p ${threads} > "${output_filepath}"
done

echo "[ $(date -u) ]: Done." | tee -a "${logfile}"

cd ../..
# You might need to re-start your bash session to get rid of some of the annoying effects of `set -euo pipefail` on normal Terminal work..
```

Get HMM files
```bash
mkdir -p unassembled_read_analysis/hmm_files

# Get rpoB from FunGene
curl -LOJ http://fungene.cme.msu.edu/hmm_download.spr?hmm_id=31
# Saves as 'rpoB.hmm'; Nov. 2009 version

mv rpoB.hmm unassembled_read_analysis/hmm_files
```

Run MetAnnotate v0.9.2
```bash
work_dir="unassembled_read_analysis"
# **MANUALLY provide the MetAnnotate data directory. If you haven't used MetAnnotate before, 
# you'll need to download this as described in [the README](https://github.com/MetAnnotate/MetAnnotate/tree/develop)
refseq_dir="/Analysis/metannotate/data"
orf_dir="${work_dir}/faa_files"
hmm_dir="${work_dir}/hmm_files"
output_dir="${work_dir}/metannotate"

mkdir -p "${output_dir}"
cd "${output_dir}"

# Enter the MetAnnotate Docker container (note that the Docker container will auto-install when you run this command; you just need to have docker on your server)
wget https://github.com/MetAnnotate/MetAnnotate/releases/download/0.9.2/enter-metannotate
./enter-metannotate "${refseq_dir}" "${orf_dir}" "${hmm_dir}" "${output_dir}"

# Now, inside the docker container:
threads=10
echo ${threads} > MetAnnotate/concurrency.txt
ref_UID=$(stat -c "%u" /home/linuxbrew/output)
sudo chown -R linuxbrew output
metannotate-wrapper-docker sequence orf_files hmm_files 2>&1 | tee output/metannotate_wrapper_docker.log
sudo chown -R $ref_UID /home/linuxbrew/output

exit

cd ../..
```
Will now have a rpoB taxonomy file at `unassembled_read_analysis/hmm_files/rpoB_[number]_all_annotations_[random_code].tsv`.  
In this folder, you'll find the file `relative_abundances/rpoB_all_annotations.tsv` as a reference copy of that file.

Relative abundances of taxa from the enrichment cultures using both unassembled read-based methods and read mapping to genome bins (above) are used in Extended Data Fig. 1.

## Relative abundance plot
Extended Data Fig. 1 shows the relative abundances of taxa in the enrichment cultures based on the unassembled read-based profiles and the MAG-based 
read mapping profiles generated above. To make Extended Data Fig. 1, the R script `Figure_ED1_plotter.R` was used, which is in the 
`relative_abundances` folder.

## Scanning metagenomes containing potential '_Ca._ Chloroheliales' members for photosynthesis genes
This was only briefly mentioned in the paper, but because two genomes from the GTDB were classified to the same order as '_Ca_. Chx. allophototropha', 
I examined the corresponding raw metagenome files that those two genome bins were derived from to search for signs of Type I reaction center-associated 
genes. I used the custom HMMs developed as decribed in the `genome_bin_analysis` folder.

Downloaded the metagenomes
```bash
# NOTE: these are large in size (several GB)

guide_filepath="other_metagenomes/other_metagenomes_accessions.tsv"
output_dir="other_metagenomes/downloads"

mkdir -p "${output_dir}"
cd "${output_dir}"

run_IDs=($(cut -d $'\t' -f 5 "${guide_filepath}" | tail -n +2))
directions=($(cut -d $'\t' -f 6 "${guide_filepath}" | tail -n +2))
urls=($(cut -d $'\t' -f 7 "${guide_filepath}" | tail -n +2))

echo "[ $(date -u) ]: Downloading ${#run_IDs[@]} files"
for i in $(seq 1 ${#run_IDs[@]}); do

  # Set zero ordered
  j=$((${i}-1))

  # Get variables
  run_ID=${run_IDs[${j}]}
  direction=${directions[${j}]}
  url=${urls[${j}]}

  output_filename="${run_ID}_${direction}.fastq.gz"
  echo "[ $(date -u) ]: Downloading '${output_filename}'"

  wget -O "${output_filename}" "${url}"

done
echo "[ $(date -u) ]: Done."

cd ../..
```

Then predicted short fragments for all files using a slightly modified version of the FGSpp code above (just change the output dir to other_metagenomes/faa_predictions and don't limit to R1 files in the 'find' search). But note that I had to set `-w` to 1 instead of 0 for the files from Kantor et al., 2015, possibly because 
the program wasn't able to handle the longer 2x250 bp read lengths.

Then used hmmsearch to scan for photosynthesis-associated genes:
```bash
work_dir="other_metagenomes"
orf_dir="${work_dir}/faa_predictions"
hmm_dir="../genome_bin_analysis/hidden_markov_models"
output_dir="${work_dir}/hmmsearch"
output_filepath="${output_dir}/hmm_hits.tsv"
evalue=1e-1
threads=20

# Install hmmsearch using:
# conda create -n hmmsearch_3.1b2 -c bioconda hmmer=3.1b2
# conda activate hmmsearch_3.1b2

mkdir -p "${output_dir}/raw"
cd "${output_dir}"

# Find input files
hmm_files=($(find "${hmm_dir}" -type f -iname "*.hmm" | sort -h))
orf_files=($(find "${orf_dir}" -type f -iname "*.faa" | sort -h))

# Initialize output file
printf "target\thmm_name\thmm_accession\tevalue\tscore\tbias\n" > "${output_filepath}"

for orf_file in ${orf_files[@]}; do
  orf_name="${orf_file##*/}"
  orf_name="${orf_name%*.frag.faa}"

  for hmm_file in ${hmm_files[@]}; do
    hmm_name="${hmm_file##*/}"
    hmm_name="${hmm_name%.*}"

    echo "[ $(date -u) ]: ${hmm_name}: ${orf_name}"
    # hmmsearch v3.1b2
    hmmsearch -o /dev/null --tblout /dev/stdout -E ${evalue} --cpu ${threads} "${hmm_file}" "${orf_file}" | \
      tee "raw/${hmm_name}_${orf_name}_raw.txt" | \
      grep -v "^#" | tr -s " " "\t" | cut -f 1,3-7 >> "${output_filepath}"
  done
done

echo "[ $(date -u) ]: Finished."
```
There were not many hits at all - just ~13 PscA hits. I pulled the predicted amino acid sequences for these hits and ran them against 
RefSeq via online BLASTP (Mar. 11th, 2020), and all had close (>95%) identity to PSI from chloroplasts or to other known PSI/RCI genes.

Lastly, checked the percent relative abundances of the detected genome bins within the metagenome to make sure that they are of sufficient 
relative abundance to get reliable hits of phototrophy genes (if they indeed encode them).
```bash
# Downloaded the genome bins
mkdir -p other_metagenomes/bins
wget -O other_metagenomes/bins/Chloroflexi_bin_54-19.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1/GCA_001898225.1_ASM189822v1_genomic.fna.gz
wget -O other_metagenomes/bins/Chloroflexi_RRmetagenome_bin16.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/243/865/GCA_003243865.1_ASM324386v1/GCA_003243865.1_ASM324386v1_genomic.fna.gz

# Then mapped the reads
work_dir="other_metagenomes"
output_dir="${work_dir}/read_mapping"
source_dir="other_metagenomes/downloads"
bin_dir="other_metagenomes/bins"
threads=40

# Need to use bbmap; install via:
# conda create -n bbmap_38.75 -c bioconda bbmap=38.75
# conda activate bbmap_38.75

mkdir -p "${output_dir}"
cd "${output_dir}"

bins=($(find "${bin_dir}" -type f -iname "*.fna.gz" | sort -h))
fastq_R1s=($(find "${source_dir}" -type f -iname "*_R1.fastq.gz" | sort -h))

for bin in ${bins[@]}; do
  bin_basename="${bin##*/}"
  bin_basename="${bin_basename%.fna.gz}"

  for fastq_R1 in ${fastq_R1s[@]}; do

    fastq_basename="${fastq_R1##*/}"
    fastq_basename="${fastq_basename%_R1.fastq.gz}"
    fastq_R2="${fastq_R1%_R1.fastq.gz}_R2.fastq.gz"
    outm="${bin_basename}_${fastq_basename}.sam"
    outcov="${bin_basename}_${fastq_basename}_cov.txt"

    echo "[ $(date -u) ]: ${bin_basename}: ${fastq_basepath}"

    bbmap.sh -Xmx20g ref="${bin}" in="${fastq_R1}" in2="${fastq_R2}" outm="${outm}" covstats="${outcov}" \
      minid=0.95 maxsites=5 nodisk=t threads=${threads} 2>&1 | tee "${fastq_basename}.log"
  done
done

cd ../..
```

Results for the 54-19 bin for its corresponding metagenome (SRR3901702):
```
Percent mapped:                         0.852
Percent proper pairs:                   0.830
Average coverage:                       15.475
Average coverage with deletions:        15.433
Standard deviation:                     4.545
Percent scaffolds with any coverage:    100.00
Percent of reference bases covered:     100.00
Total time:             150.129 seconds.
```

Results for the RRmetagenome_bin16 bin for its corresponding three metagenomes:
```
==> SRR5223441.log <==
Percent mapped:                         3.124
Percent proper pairs:                   3.005
Average coverage:                       23.674
Average coverage with deletions:        23.584
Standard deviation:                     8.962
Percent scaffolds with any coverage:    100.00
Percent of reference bases covered:     99.92
Total time:             35.155 seconds.

==> SRR5223442.log <==
Percent mapped:                         0.751
Percent proper pairs:                   0.722
Average coverage:                       5.782
Average coverage with deletions:        5.746
Standard deviation:                     3.084
Percent scaffolds with any coverage:    100.00
Percent of reference bases covered:     98.28
Total time:             30.410 seconds.

==> SRR5223443.log <==
Percent mapped:                         1.635
Percent proper pairs:                   1.569
Average coverage:                       12.764
Average coverage with deletions:        12.731
Standard deviation:                     5.548
Percent scaffolds with any coverage:    100.00
Percent of reference bases covered:     99.88
Total time:             29.782 seconds.
```

Both genome bins seem to have reasonable coverages. Thus, the lack of detection of PscA-like/FmoA sequences either means the HMMs do not cover 
the full diversity of this clade or that these genome bins do not encode Type I reaction center-associated genes. Given that PSI genes of 
possible chloroplasts were detected (which are highly divergent from the _Ca._ Chloroheliales RCI), the data imply that the HMMs (at least PscA) 
are able to detect remote PSI/RCI homologs among short read data. This makes the latter hypothesis (i.e., no Type I reaction center-associated genes 
in the bins) more likely.

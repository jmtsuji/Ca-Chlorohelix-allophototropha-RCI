# 16S ribosomal RNA gene amplicon analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

## Overview
As described in the paper, 16S rRNA gene amplicon sequencing data was analyzed in three different ways depending on the sample:
- Analysis of V4-V5 region amplicons via AXIOME3 (a wrapper of QIIME2)
- Analysis of V4 region amplicons via QIIME2 directly
- Analysis of V1-V9 region amplicons via NanoCLUST

The code for all three analysis methods is shown here.

## AXIOME3
Code will be run within the `AXIOME3` directory unless indicated otherwise.

### Install AXIOME3
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

### Get the data
See the `source_data` directory, specifically `culture_early_16S_rRNA_gene_data_accessions.tsv`, for download details (along with the README).

### Analyze the data
Analysis is done within the AXIOME3 folder cloned from Github above (a subfolder of the main `AXIOME3` directory).

For classification, you will need to download and train the Silva classifier, which is described 
[here](https://docs.qiime2.org/2019.10/tutorials/feature-classifier/) but not shown in this code summary.

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
Once finished, the ASV table will be available at `AXIOME3/AXIOME3/output/exported/ASV_table_combined.tsv`.

**In lieu of uploading all intermediate analysis files, a summary of the non-rarefied amplicon sequencing variant (ASV) table is provided at `ASV_table_non_rarefied.tsv`. **
**ASV table data was used to generate Supplementary Data 1 in the manuscript.**

## Direct use of QIIME2
Code will be run within the `direct_QIIME2` directory.

### Get the data
Data files will be made available online at BioProject PRJNA640240 before the paper is published.

TODO - add download instructions once publicly available.

After download, put files `Chx_S21.2c_R1.fastq.gz` and `Chx_S21.2c_R2.fastq.gz` here.

### Analyze with QIIME2
Assumes you have already installed QIIME2 2022.8 based on the instructions on the QIIME2 website (link [here](https://docs.qiime2.org/2022.8/install/native/#miniconda))

Make manifest
```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
CS21	input/Chx_S21.2c_R1.fastq.gz	input/Chx_S21.2c_R2.fastq.gz
```
Note that you might need absolute paths instead of relative paths.

Saved as `manifest.tsv` in `direct_QIIME2`.

Run
```bash
classifier_path="[ ADD PATH TO CLASSIFIER QZA FILE ]" # E.g., download 'Silva 138 99% OTUs from 515F/808R region of sequences' from the QIIME2 website: https://docs.qiime2.org/2023.2/data-resources/
bioproject_id="CS21"
primer_seq_f="GTGCCAGCMGCCGCGGTAA"
primer_seq_r="GGACTACHVGGGTWTCTAAT"
trunc_len_f=200
trunc_len_r=140
threads=20

conda activate qiime2-2022.8

mkdir -p "qiime2" && cd "$_"
mkdir -p "export" "analyze"

# Import FastQ
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "../manifest.tsv" \
  --output-path "01_demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

# Trim primers and throw away reads with no primer
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "01_demux.qza" \
  --o-trimmed-sequences "02_trim.qza" \
  --p-cores "${threads}" \
  --p-front-f "${primer_seq_f}" \
  --p-front-r "${primer_seq_r}" \
  --p-error-rate 0.1 \
  --p-discard-untrimmed \
  --verbose 2>&1 | \
  tee "02_trim.log"

# CHECKPOINT: evaluate primer trimming
grep -B 10 -A 3 "Read fate breakdown" "02_trim.log" | \
  grep -v "^Processing" | grep -v "^Finished" | grep -v "^$" | grep -v "^=" | grep -v "^Total" | \
  grep -v "^ " | sed 's\^Command.*/data/\\g' | sed 's/_[0-9]*_L[0-9]*_R[1-2]_[0-9]*.fastq.gz//g' \
  > "02_stats.txt"
# cat 02_stats.txt

# CHECKPOINT: evaluate read quality stats
qiime demux summarize \
  --i-data "02_trim.qza" \
  --o-visualization "02_qual.qzv"
qiime tools export --input-path 02_qual.qzv --output-path export/cutadapt_qual_stats
# Take a look at the resulting .qzv file at https://view.qiime2.org/ to identify good truncation lengths for preserving sequence regions before the sequence quality drops. Those lengths are set as trunc_len_f=200 and trunc_len_r=140 above (based on my manual inspection of the reads when I ran the analysis for the manuscript)

# Denoise (using the truncation lengths set above)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "02_trim.qza" \
  --o-table "03_table.qza" \
  --o-representative-sequences "03_seqs.qza" \
  --o-denoising-stats "03_stats.qza" \
  --p-trunc-len-f "${trunc_len_f}" \
  --p-trunc-len-r "${trunc_len_r}" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-min-overlap 12 \
  --p-pooling-method independent \
  --p-n-threads "${threads}" \
  --verbose 2>&1 | tee \
  "03_dada2.log"

# CHECKPOINT: evaluate DADA2 performance
qiime tools export --input-path 03_stats.qza --output-path export/dada2_stats
cat export/dada2_stats/stats.tsv | column -t -s $'\t'

# Classify
qiime feature-classifier classify-sklearn \
  --i-reads "03_seqs.qza" \
  --i-classifier "${classifier_path}" \
  --o-classification "04_taxonomy.qza" \
  --p-n-jobs "${threads}" \
  --verbose 2>&1 | \
  tee "04_taxonomy.log"

# Export
qiime tools export --input-path 03_table.qza --output-path export/table
qiime tools export --input-path 03_seqs.qza --output-path export/seqs
qiime tools export --input-path 04_taxonomy.qza --output-path export/taxonomy
biom convert --to-tsv -i export/table/feature-table.biom -o export/table/feature-table.tsv

# Output easy-to-understand table
git clone https://github.com/jmtsuji/qiime2-helpers.git
cd qiime2-helpers
git checkout 1275e4a
cd ..

conda deactivate
conda activate pandas # make this conda env with the following programs installed:
# python 3.9.2, pandas 1.2.3

qiime2-helpers/scripts/generate_combined_feature_table.py \
  -f export/table/feature-table.tsv \
  -s export/seqs/dna-sequences.fasta \
  -t export/taxonomy/taxonomy.tsv \
  -o export/combined_feature_table.tsv \
  -S 2>&1 | \
  tee export/combined_feature_table.log

cd ..
```

DADA2 stats - most of the reads passed the pipeline, which is good.
```
sample-id  input    filtered  percentage of input passed filter  denoised  merged   percentage of input merged  non-chimeric  percentage of input non-chimeric
#q2:types  numeric  numeric   numeric                            numeric   numeric  numeric                     numeric       numeric
CS21       57318    53649     93.6                               53629     53386    93.14                       49950         87.15
```

The `combined_feature_table.tsv` shown above is included in this repo and has the raw amplicon analysis results.

### Post analysis
There were a few low-count ASVs included in the data. Here, I check if those might be artefacts of high sequencing depth.

OTU clustering (99%) - as proof of concept that the low count ASVs are similar to the two high-count ones.

```bash
mkdir -p "otu_clustering" && cd "$_"

conda activate clustering # make this conda env with mmseqs 13.45111

mmseqs easy-cluster "../qiime2/export/seqs/dna-sequences.fasta" p99c99 tmp \
  -c 0.99 --cov-mode 5 --cluster-mode 2 --min-seq-id 0.99 --threads 10 2>&1 | \
  tee p99c99.log

ln -s p99c99_rep_seq.fasta repseqs.fasta

grep -c "^>" repseqs.fasta
# collapsed to 2 rep seqs
```

Repseq summary (`column -t p99c99_cluster.tsv `)
```
50997d1e6c398e7abaceb666004a2456  50997d1e6c398e7abaceb666004a2456
50997d1e6c398e7abaceb666004a2456  eed1df5d000ee2e59e32bd18e7f52365
07f4ce10c67f5c9585c783f63c14fdef  07f4ce10c67f5c9585c783f63c14fdef
07f4ce10c67f5c9585c783f63c14fdef  5bcefe19fdbf27d141ee3e85eb3d25fb
07f4ce10c67f5c9585c783f63c14fdef  b193ceb50cc19bb9ac1e7397aeb46537
07f4ce10c67f5c9585c783f63c14fdef  e04ee9556220165f6cfda8d2350708f1
```
In other words, the two high count ASVs became representative sequences, and the four low-count ASVs map to them at >99% identity.

Checked the ASV sequences against reference 16S rRNA gene sequences of "_Ca_. Chx. allophototropha" and _Geothrix_ sp. L227-G1 from their closed genomes.

Results: the `50997d1e6c398e7abaceb666004a2456` OTU is 100% identical across full length to both L227-G1 seqs.
         the `07f4ce10c67f5c9585c783f63c14fdef` OTU is 100% identical across full length to Ca. Chx. allophototropha's 16S seq.
         the other 4 low count sequences each has 1 bp mismatch to the high-count ASVs shown above.

I thus made a judgement call that, given the high sequencing depth of this sample, the 4 low-count sequences are likely just a collection of sequencing errors that couldn't be accurately corrected during ASV calling. I thus combined all the data into two ASVs based on >99% identity mapping to the "_Ca._ Chx. allophototropha" and _Geothrix_ L227-G1 reference sequences.

Made a simple pie chart in Excel using this collapsed data and exported it as a SVG (`pie-chart-CS21.svg`). This was added to Extended Data Fig. 1.

## NanoCLUST
Code will be run in the `NanoCLUST` subfolder.

Workflow overview:
- Get the data
- Run NanoCLUST
- Trim adapters and fix orientations using two rounds of cutadapt
- Cluster at 99% threshold
- Filter chimeric sequences
- Calculate relative abundances

### Get the data
The basecalled reads will be available at BioProject PRJNA640240 upon publication.

TODO: add info and download instructions once ready.

#### A note on basecalling
The basecalled reads will be downloaded directly from NCBI, but here I'll leave a few notes about how basecalling was performed that are outside of what is mentioned in the manuscript text.

Most samples were basecalled real-time using Guppy 5.1.12 with the SUP model.

One sample (2022.11.21 16S barcoding run for Chx S19.9) had to be basecalled after the run. Performed basecalling using Guppy 5.0.16 with the SUP model:
```bash
mkdir -p guppy5 && cd "$_"

# Make list of input fast5 files - assumes that the raw Fast5 files are in an adjacent folder called `raw`
find ../raw -name "*.fast5" | sort -h > fast5_input.list

# Assumes guppy is already installed at the location shown below... see instructions on Nanopore website
../ont-guppy/bin/guppy_basecaller \
  --input_file_list fast5_input.list \
  --save_path . \
  --config ../ont-guppy/data/dna_r9.4.1_450bps_sup.cfg \
  --num_callers 12 \
  --device 'auto' \
  --compress_fastq \
  --min_qscore 7 \
  --barcode_kits SQK-16S024 \
  --num_barcode_threads 2 \
  2>&1 | tee guppy5.log
```

#### Prepare metadata
Sample metadata - TODO - make sure this matches the NCBI download names once the files are available online.
```
sample-id	culture-id
20211112_barcode21	Chx-S19.9
20220216_barcode08	Chx-S21.2c
20220216_barcode11	Chx-S22.2a
20220216_barcode12	Chx-S22.2b
20220216_barcode13	Chx-S22.2c
20220216_barcode14	Chx-S22.2d
20221116_barcode08	G1-5.1b
```
Save as `sample-metadata.tsv`.

### Install NanoCLUST
#### Downloads
Download the Git repo:
```bash
git clone https://github.com/jmtsuji/NanoCLUST.git
cd NanoCLUST
git checkout expose_multithreading # commit a09991c
cd ..
```

Get BLAST 16S DB:
```bash
mkdir -p "db/taxdb"
cd db
wget -O - https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz | tar -xzvf -
cd taxdb
wget -O - https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar -xzvf -

cd ../..
```
Note: NanoCLUST may or may not work with never databases from NCBI without some manual modifications. (TODO: add link to issue) In the end, I do not use classifications from NanoCLUST for this work, so you might also be able to find a way to bypass this classification step entirely if you'd like.

#### Test
Run with docker as a test
```bash
# conda create -n nextflow -c bioconda nextflow
conda activate nextflow # nextflow version 20.10.0

nextflow run main.nf -profile test,docker
```
Confirm this finishes without errors.

#### Update
Update the version used for medaka so that it can use the SUP model. To do this, I'll need to initialize the docker container and then make a modified version with updated medaka. The docker container should already be there given that the test was run.
```bash
# Start the medaka docker container used by NanoCLUST
docker run -it hecrp/nanoclust-medaka_pass /bin/bash
```

Within the docker container
```bash
# Update medaka
conda update -n base conda -y
conda env remove -n medaka_pass
conda create -y -n medaka_pass -c bioconda -c conda-forge medaka=1.5.0
conda clean -a -y
```

OUTSIDE the docker repo, while the repo is running, i.e., in a separate terminal window, commit a new docker image
```bash
docker ps # see the ID of the currently running container and add it below.

container_id=7f90a01c879e # set manually after checking via docker ps

docker commit "${container_id}" microeco/nanoclust-medaka_pass
```

Then exit the container from the original window
```bash
exit
```

Finally, update the medaka model and docker image target in the NanoCLUST code. In the NanoCLUST Github folder cloned above, run:
```bash
# Add the SUP model to NanoCLUST
mv main.nf main.nf.backup
sed 's/r941_min_high_g303/r941_min_sup_g507/g' main.nf.backup > main.nf

# Add the docker container to NanoCLUST
mv nextflow.config nextflow.config.backup
sed 's\hecrp/nanoclust-medaka_pass\microeco/nanoclust-medaka_pass\g' nextflow.config.backup > nextflow.config
```
You're now ready to go!

### Run NanoCLUST
All samples used the 16S barcoding kit on a R9.4.1 Flongle flow cell with q-score limit of 7 and super high accuracy basecalling.

```bash
# Make sure you run these commands from inside the main `16S_rRNA_gene_amplicon/NanoCLUST` subfolder (not `16S_rRNA_gene_amplicon/NanoCLUST/NanoCLUST` where the Github repo is, for example)
input_dir="$(realpath input)" # put the downloaded FastQ files for the samples here
output_dir="$(realpath nc_output)"
nanoclust_dir="NanoCLUST" # Github repo dir
min_cluster_size=10

conda activate nextflow
# nextflow 20.10.0

mkdir -p "${output_dir}"
cd "${nanoclust_dir}"

input_fastqs=($(find "${input_dir}" -name "*.fastq" | sort -h))

# Clear the NanoCLUST history
rm -rf .nextflow work

for input_fastq in "${input_fastqs[@]}"; do

  sample_name="${input_fastq%.fastq}"
  sample_name="${sample_name##*/}"
  run_name="sample" # Doesn't really matter, because I will clear the NanoCLUST history after each run

  echo "[ $(date -u) ]: Analyzing sample '${sample_name}'"

  nextflow run main.nf \
    -profile docker \
    --reads "${input_fastq}" \
    --db "db/16S_ribosomal_RNA" \
    --tax "db/taxdb" \
    -name "${run_name}" \
    --max_cpus 40 \
    --min_cluster_size "${min_cluster_size}" \
    > "${sample_name}.log" 2>&1

  mv results "${output_dir}/${sample_name}"
  mv "${sample_name}.log" "${output_dir}"

  # Clear the NanoCLUST history
  rm -rf .nextflow work

done

# Organize output
cd "${output_dir}"

for input_fastq in "${input_fastqs[@]}"; do

  sample_name="${input_fastq%.fastq}"
  sample_name="${sample_name##*/}"

  echo "[ $(date -u) ]: Organizing sample '${sample_name}'"

  run_dir="${output_dir}/${sample_name}"

  # First, summarize the cluster sequences
  printf "" > "${run_dir}/cluster_sequences.fasta"

  cluster_dirs=($(find "${run_dir}/${sample_name}" -maxdepth 1 -type d -name "cluster*"))

  for cluster_dir in "${cluster_dirs[@]}"; do

    cluster_name="${cluster_dir##*/}"
    echo "[ $(date -u) ]: '${sample_name}': ${cluster_name}"

    cluster_fasta="${run_dir}/${sample_name}/${cluster_name}/consensus_medaka.fasta/consensus.fasta"

    printf ">${sample_name}_${cluster_name}\n" >> "${run_dir}/cluster_sequences.fasta"
    tail -n +2 "${cluster_fasta}" >> "${run_dir}/cluster_sequences.fasta"
  
  done

  # Next, append run name to the output table and change to a TSV
  head -n 1 "${run_dir}/${sample_name}/${sample_name}.nanoclust_out.txt" | sed 's/^/sample;/g' | \
    tr ";" $'\t' > "${run_dir}/cluster_counts.tsv"
  tail -n +2 "${run_dir}/${sample_name}/${sample_name}.nanoclust_out.txt" | sed "s/^/${sample_name};/g" | \
    tr ";" $'\t' >> "${run_dir}/cluster_counts.tsv"

done
```
Done. Check everything finished successfully based on NanoCLUST logs.

### Additional installations before subsequent analyses steps
Will use `uchime2_ref` as described here: https://drive5.com/usearch/manual/cmd_uchime2_ref.html (accessed Dec 28th, 2022)

#### Chimera DB
Prepare chimera DB
```bash
mkdir -p "chimera_db" && cd "$_"

conda activate general # install the following software in this conda env:
# blast 2.12.0
# seqtk 1.3-r106

# Downloaded NCBI 16S DB (on Dec 28 2022; DB last updated 2022-12-15)
mkdir -p 2022_12_28_NCBI_16S && cd "$_"
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz.md5
md5sum 16S_ribosomal_RNA.tar.gz > 16S_ribosomal_RNA.tar.gz.md5.local
cmp 16S_ribosomal_RNA.tar.gz.md5 16S_ribosomal_RNA.tar.gz.md5.local # OK
tar -xvzf 16S_ribosomal_RNA.tar.gz
cd ..

# Extract FastA
blastdbcmd -entry all -db 2022_12_28_NCBI_16S/16S_ribosomal_RNA -out 16S_ribosomal_RNA.fasta

grep -c "^>" 16S_ribosomal_RNA.fasta # 26809 entries
```

Pasted the Chlorohelix 16S rRNA gene sequence into this folder as `Chx_16S.fasta`:
```
>Chx_allophototropha 16S ribosomal RNA
TGAAGAGTTTGATCCTGGCTCAGGATAAACGCTGGCGGCGTGCCTAATGCATGCAAGTCG
AACGGGAGTAGCAATACTCAAGTGGCAAACGGGTGAGTAACACGTGGGAACCTGCCCTTT
AGTGGGGGACAACCTTTCGAAAGAGAGGCTAATACCGCATACGGTGGAGACACTAAAGCA
GCAATGCGCTGAAGGAGGGGCCTGCGTCTGATTAGCTAGTAGGTGGGGTAAAAGCCTACC
TAGGCGATGATCAGTAGCTGGTCTGAGAGGATGGTCAGCCACACTGGGACTGAGAAACGG
CCCAGACTCCTACGGGAGGCAGCAGCAAGGAATTTTCGGCAATGGGGGAAACCCTGACCG
AGCAACGCCGCGTGGGCGATGAAGTCTTTCGGGATGTAAAGCCCTTTTCTGTGGGAAGAG
CAAGGACGGTACCACAGGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATAC
GTAGGGTGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGCGCGCAGGCGGCTTGTTAA
GTTACTCGTGAAAGCCCCCGGCTAAACCGGGGAGGGTCGAGTGATACTGGCAGGCTAAGA
GAGCAGCAGAGGATAGTGGAATTCCCGGTGTAGTGGTGAAATGCGTAGATATCGGGAGGA
ACACCAGTGGCGAAAGCGACTATCTGGGCTGTGTCTGACGCTGAGGCGCGAAGGCTAGGG
GAGCAAACAGGATTAGATACCCTGGTAGTCCTAGCAGTAAACGATGAACACTAGGTGTGC
GGGGAAATTGACCCCCTGCGTGCCGTAGCTAACGCAATAAGTGTTCCGCCTGGGGAGTAC
GGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGACCCGCACAAGCAGCGGAGCGTGTG
GTTTAATTCGACGCAACGCGAAGAACCTTACCAGGGCTTGACATACTACGTTCAAGGGCG
GAAACGTTCTGGTCGCAAGACGAGTAGTACAGATGCTGCATGGCTGTCGTCAGCTCGTGT
CGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTAGTGTTAGTTGAAAGATCT
AGCGCGACTGCCGGTAGGAAACCGGAGGAAGGTGGGGATGACGTCAAGTCAGCATGGCTC
TTACGTCCTGGGCTACACACACGCTACAATGGCTGGTACAAAGCGAAGCGAGACGGCGAC
GTGGAGCGGATCGCAAAAAGCCAGTCTCAGTTCGGATAGCAGGCTGCAACTCGCCTGCTT
GAAGTCGGAGTTGCTAGTAATCGCGCATCAGCACGGTGCGATGAATATGTTCCCGGGTCT
TGTACACACCGCCCGTCACGTCATGGAAGCTGGTAACACCTGAAGCCGGTGTGCTAACCG
CAAGGAGGCAGCTGTCGAGGGTGGGATTAGTGACTGGGACGAAGTCGTAACAAGGTAGCC
GTAGCGGAAGCTGCGGCTGGATCACCTCCTTT
```

Then added to the 16S DB:
```bash
cp 16S_ribosomal_RNA.fasta 16S_ribosomal_RNA_with_Chx.fasta

seqtk seq -l 80 Chx_16S.fasta >> 16S_ribosomal_RNA_with_Chx.fasta

cd ..
```
The `16S_ribosomal_RNA_with_Chx.fasta` file will be used as a reference for chimera searching. 26810 entries including Chlorohelix.

#### UCHIME
Download UCHIME2 (free version)
```bash
# usearch11.0.667_i86linux32
wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz
chmod 755 usearch11.0.667_i86linux32
```

### Adapter trimming
Trim primers and adapters off the ends of both sequences. Flip sequences to the forward orientation at the same time.
```bash
mkdir -p "trim" && cd "$_"

source_dir="../nc_output"
fwd_primer_27f="AGAGTTTGATCMTGGCTCAG" # Lane, 1991
rev_primer_1492r="TACGGYTACCTTGTTACGACTT" # Lane, 1991

# Gather inputs
mkdir -p "01_inputs" && cd "$_"
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_sequences.fasta" | sort -h | \
  xargs cat > cluster_sequences_raw.fasta
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_counts.tsv" | sort -h | \
  head -n 1 | xargs head -n 1 > cluster_counts_raw.tsv
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_counts.tsv" | sort -h | \
  xargs tail -q -n +2 >> cluster_counts_raw.tsv
cd ..
# Note: saved a copy of these in the `intermediate` folder for reference

# Get 1492r reverse complement sequence
conda activate general # seqtk 1.3-r106
rev_primer_1492r_revcomp=$(printf ">1492r\n${rev_primer_1492r}\n" | seqtk seq -r | tail -n 1)
conda deactivate

# Trim
mkdir -p "cutadapt" && cd "$_"
conda activate cutadapt_3.4 # install cutadapt 3.4 in this env
# This code should work in theory to make the env, although the env was pre-existing: conda create -n cutadapt -c bioconda cutadapt=3.4
# cutadapt 3.4 with Python 3.9.6
# Setting O (minlength) to 15 instead of default 3 because I know the full primers should be there
# Using --revcomp to flip the sequences if the 27F is on the reverse complement strand
# Untrimmed sequences will be saved in a separate file and not used downstream (these sequences might be lacking a primer, for example, meaning that they are not full amplicons, or they might be garbage sequences)
cutadapt -o cluster_sequences_trimmed_fwd.fasta \
  -g "${fwd_primer_27f}" \
  --untrimmed-output untrimmed_fwd.fasta \
  --info-file cutadapt_trim_info_fwd.txt \
  --revcomp \
  ../01_inputs/cluster_sequences_raw.fasta 2>&1 | \
  tee cutadapt_fwd.log
# 57/61 sequences retained

# No need for --revcomp this time, because sequences have already been reverse complemented where needed
cutadapt -o cluster_sequences_trimmed_rev.fasta \
  -a "${rev_primer_1492r_revcomp}" \
  --untrimmed-output untrimmed_rev.fasta \
  --info-file cutadapt_trim_info_rev.txt \
  cluster_sequences_trimmed_fwd.fasta 2>&1 | \
  tee cutadapt_rev.log
# 55/57 sequences retained

# Link the final trimmed reads
ln -s cluster_sequences_trimmed_rev.fasta cluster_sequences_trimmed_complete.fasta

conda deactivate
cd ..
```
Put a copy of `cluster_sequences_trimmed_complete.fasta` into the `intermediate` folder.

### OTU clustering
```bash
mkdir -p "cluster" && cd "$_"

conda activate general # install:
# mmseqs 13.45111

mmseqs easy-cluster ../trim/cluster_sequences_trimmed_complete.fasta p99c99 tmp \
  -c 0.99 --cov-mode 5 --cluster-mode 2 --min-seq-id 0.99 --threads 2 2>&1 | \
  tee p99c99.log

# collapsed to 14 representative sequences

ln -s p99c99_rep_seq.fasta repseqs.fasta

conda deactivate
cd ..
```
Saved a copy of `p99c99_cluster.tsv` in the `intermediate` folder.

### Chimera filtration
Using UCHIME2 with the NCBI+Chx 16S DB downloaded above.

```bash
mkdir -p "chimera_filter" && cd "$_"

db_filepath="../chimera_db/16S_ribosomal_RNA_with_Chx.fasta"
usearch_filepath="../usearch11.0.667_i86linux32"

"${usearch_filepath}" -uchime2_ref "../cluster/repseqs.fasta" -db "${db_filepath}" \
  -mode high_confidence -strand plus -uchimeout chimera_report.tsv 2>&1 | \
  tee uchime2.log
```

Results: 9/14 sequences were identified as chimeras. 5 remain.

Summarized the chimera-filtered sequence set (i.e., the sequences not identified as chimeras):
```bash
# This was done manually. TODO - automate if doing this in high-throughput in the future.
printf "20211112_barcode21_cluster0\n20220216_barcode12_cluster5\n20220216_barcode08_cluster1\n20220216_barcode12_cluster0\n20220216_barcode14_cluster0\n" > chimera_filtered.list

conda activate general
# seqtk 1.3-r106

seqtk subseq "../cluster/repseqs.fasta" chimera_filtered.list > chimera_filtered.fasta
```
Saved a copy of this `chimera_filtered.fasta` file and the `chimera_report` in the `intermediate` folder for reference.

### Post processing
#### Relative abundances
Calculated relative abundances: did so using `summary/relative-abundance-calculations.ipynb`. 

Saved the output files in the same `summary` folder: `otu_table_counts.tsv` and `cluster_stats.tsv`. Added Excel-analyzed versions of both.

#### Sequence identification
In my analysis, the 5 representative sequences corresponded to the folllowing:
- `20211112_barcode21_cluster0`: 100% identical to Chx
- `20220216_barcode08_cluster1`: 100% identical to G1 rRNA gene `GEOCFX_001163`
- `20220216_barcode12_cluster5`: based on web BLASTN: near 100% identity to NR_180005.1	Acinetobacter oryzae
- `20220216_barcode12_cluster0`: based on web BLASTN: near 100% identity to NR_026163.1	Microbacterium testaceum
- `20220216_barcode14_cluster0 rc`: based on web BLASTN: near 100% identity to NR_116570.1	Sphingomonas hankookensis

The sequences that are now Chx or G1 have low relative abundances, looking at the `cluster_counts_raw.tsv` file:
```
sample	id	reads_in_cluster	used_for_consensus	reads_after_corr	draft_id	sciname	taxid	length	per_ident
20220216_barcode12	5	18	100	17	576ab204-ef44-4968-b7e4-9b39ee2af8de id=9	Acinetobacter johnsonii	40214	1501	99.467
20220216_barcode12	0	32	100	31	852e4339-43a8-497a-8bea-2f6dd9d17e01 id=19	Microbacterium testaceum	2033	1465	99.249
20220216_barcode14	0	43	100	41	d2ba43c4-3c27-421b-8f06-eed41d7e0224 id=35	Sphingomonas hankookensis	563996	1430	99.441
```
Reads in cluster ranged from 18-43 out of 40-50k total reads. This is less than 0.1% relative abundance.

These data were added into Extended Data Table 1.

Done!

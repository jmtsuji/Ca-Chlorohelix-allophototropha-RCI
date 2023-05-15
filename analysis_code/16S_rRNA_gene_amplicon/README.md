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
TODO

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
Data files will be made available online at BioProject ___TODO___ before the paper is published.

TODO - add download instructions once publicly available.

Got FastQ files
```bash
mdkir -p "input"
```
Put files `Chx_S21.2c_R1.fastq.gz` and `Chx_S21.2c_R2.fastq.gz` here.

### Analyze with QIIME2
Code is adapted from that used for the MGnify amplicon analyses in around Oct. 2022.

Make manifest
```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
CS20	/Analysis/jmtsuji/projects/culture/2023_01_11_Chx_S20.2c_MiSeq_16S_analysis/01_raw_fastq/CS20_R1.fastq.gz	/Analysis/jmtsuji/projects/culture/2023_01_11_Chx_S20.2c_MiSeq_16S_analysis/01_raw_fastq/CS20_R2.fastq.gz
```

Saved as `manifest.tsv` in `/Analysis/jmtsuji/projects/culture/2023_01_11_Chx_S20.2c_MiSeq_16S_analysis` as shown below.

Run
```bash
work_dir="/Analysis/jmtsuji/projects/culture/2023_01_11_Chx_S20.2c_MiSeq_16S_analysis"
manifest_dir="${work_dir}"
classifier_path="/Data/databases/qiime2/classifier/SILVA_138.1_SSURef_NR99/silva-138-ssu-nr99-classifier-uniq_515_806.qza"
bioproject_id="CS20"
primer_seq_f="GTGCCAGCMGCCGCGGTAA"
primer_seq_r="GGACTACHVGGGTWTCTAAT"
trunc_len_f=200
trunc_len_r=140
threads=20

conda activate qiime2-2022.8

mkdir -p "${work_dir}/${bioproject_id}" && cd "$_"
mkdir -p "export" "analyze"

# Import FastQ
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${manifest_dir}/manifest.tsv" \
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

# Denoise
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
conda activate jupyterlab
# python 3.9.2, pandas 1.2.3

qiime2-helpers/scripts/generate_combined_feature_table.py \
  -f export/table/feature-table.tsv \
  -s export/seqs/dna-sequences.fasta \
  -t export/taxonomy/taxonomy.tsv \
  -o export/combined_feature_table.tsv \
  -S 2>&1 | \
  tee export/combined_feature_table.log
```

DADA2 stats
```
sample-id  input    filtered  percentage of input passed filter  denoised  merged   percentage of input merged  non-chimeric  percentage of input non-chimeric
#q2:types  numeric  numeric   numeric                            numeric   numeric  numeric                     numeric       numeric
CS20       57318    53649     93.6                               53629     53386    93.14                       49950         87.15
```

### Post analysis
OTU clustering (99%)

```bash
mkdir -p "/Analysis/jmtsuji/projects/culture/2023_01_11_Chx_S20.2c_MiSeq_16S_analysis/03_post_analysis/otu_clustering" && cd "$_"

conda activate general2
# mmseqs 13.45111

mmseqs easy-cluster "../../CS20/export/seqs/dna-sequences.fasta" p99c99 tmp \
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
Combined this with the ASV table to make a 99% clustered file.

Checked against refs. Put `ref_sequences.fasta` in there with the 2 L227-G1 seqs (from `Geothrix_16S_rRNA_genes_ren_rc.fasta`) and the 1 Chx seq, full-length.
```bash
conda activate general2
# blast v2.12.0

blastn -query repseqs.fasta -subject ref_sequences.fasta -perc_identity 99 -qcov_hsp_perc 90 > repseqs_vs_refs.blastn.txt
```
Results: the `50997d1e6c398e7abaceb666004a2456` OTU is 100% identical across full length to both L227-G1 seqs.
         the `07f4ce10c67f5c9585c783f63c14fdef` OTU is 100% identical across full length to Ca. Chx. allophototropha's 16S seq.

Exported pie chart from Excel as SVG for relative abundances. About 50:50, two species.

## NanoCLUST
Workflow:
- Basecalling
- Run NanoCLUST
- Trim adapters and fix orientations using two rounds of cutadapt
- Cluster at 99% threshold
- Filter chimeric sequences

### Basecalling
Most samples were basecalled real-time using Guppy ____ with the SUP model.

One sample (2022/11/21 16S barcoding run) has to be basecalled after the run.

Performed basecalling using Guppy 5.0.16 with the SUP model:
```bash
mkdir -p /home/microeco2021/nanoclust/211112_DNA_cultures/guppy5 && cd "$_"

# Make list of input fast5 files
find ../raw -name "*.fast5" | sort -h > fast5_input.list

# Assumes guppy is already installed via instructions on Nanopore website
/home/microeco2021/guppy_5.0.16/ont-guppy/bin/guppy_basecaller \
  --input_file_list fast5_input.list \
  --save_path . \
  --config /home/microeco2021/guppy_5.0.16/ont-guppy/data/dna_r9.4.1_450bps_sup.cfg \
  --num_callers 12 \
  --device 'auto' \
  --compress_fastq \
  --min_qscore 7 \
  --barcode_kits SQK-16S024 \
  --num_barcode_threads 2 \
  2>&1 | tee guppy5.log
```

### Prepare NanoCLUST inputs
#### Prepare input reads
Download samples and rename
```bash
TODO
```

#### Prepare metadata
Sample metadata
```
sample-id	culture-id
20221112_barcode21  Chx-S18.9
20220216_barcode08	Chx-S20.2c
20220216_barcode11	Chx-S21.2a
20220216_barcode12	Chx-S21.2b
20220216_barcode13	Chx-S21.2c
20220216_barcode14	Chx-S21.2d
20221116_barcode08	G1-S5.1b
```
Saved as `sample-metadata.tsv` in `/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/02_combined_fastq`

TODO - make the reduced version shown above as a file in the repo directly.

### Install NanoCLUST
#### Downloads
Download the Git repo:
```bash
cd "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro"

git clone https://github.com/jmtsuji/NanoCLUST.git
cd NanoCLUST
git checkout expose_multithreading # commit a09991c
cd ..
```

Get BLAST 16S DB (downloaded on Dec. 28th 2022):
```bash
cd "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/NanoCLUST"

# mkdir -p "db/taxdb"
# cd db
# wget -O - https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz | tar -xzvf -
# cd taxdb
# wget -O - https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar -xzvf -

# Re-using DB downloaded previously (need to look up old README for this)
rsync -a ../../2022_02_17_Chx3_NatMicro_16S_barcoding/NanoCLUST/db_220216 .
mv db_220216 db

cd ../..
```
Note: NanoCLUST may or may not work with never databases from NCBI without some manual modifications. (TODO: add link to issue)

#### Test
Run with docker as a test
```bash
# conda create -n nextflow -c bioconda nextflow
conda activate nextflow # nextflow version 20.10.0

nextflow run main.nf -profile test,docker
```
Finished without error.

#### Update
Update the version used for medaka so that it can use the SUP model. To do this, I need to initialize the docker container and then make a modified version with updated medaka. The docker container should already be there given that the test was run.
```bash
# TODO - see previous notes
# hecrp/nanoclust-medaka_pass
```

Inside the docker container
```bash
# TODO - see previous notes
```

Then, outside the docker container, before shutting down the container, commited the modified Docker container as `microeco/nanoclust-medaka_pass`:
```bash
# TODO - see previous notes
```
The container can now be shut down from the inside using `exit`. # <- TODO: confirm command

Lastly, modify the NanoCLUST nextflow files to use the updated container and SUP model
```bash
cd "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/NanoCLUST"

# Add the SUP model to NanoCLUST
sed -i 's/r941_min_high_g303/r941_min_sup_g507/g' main.nf

# Add the docker container to NanoCLUST
sed -i 's\hecrp/nanoclust-medaka_pass\microeco/nanoclust-medaka_pass\g' nextflow.config
```

Now ready to go!

### Run NanoCLUST
All samples used the 16S barcoding kit on a R9.4.1 Flongle flow cell with q-score limit of 7 and super high accuracy basecalling.

```bash
work_dir="/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro"
input_dir="${work_dir}/02_combined_fastq"
output_dir="${work_dir}/03b_NanoCLUST_c10"
nanoclust_dir="${work_dir}/NanoCLUST"
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
  mv "${sample_name}.log" "${output_dir}" #/${sample_name}

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
Done. All finished successfully based on NanoCLUST logs.

### Installations
Will use `uchime2_ref` as described here: https://drive5.com/usearch/manual/cmd_uchime2_ref.html (accessed Dec 28th, 2022)

#### Chimera DB
Prepare chimera DB
```bash
mkdir -p "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/chimera_ref_db" && cd "$_"

conda activate general2
# blast 2.12.0
# seqtk 1.3-r106

# Download NCBI 16S DB on Dec 28 2022; DB last updated 2022-12-15
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
```
The `16S_ribosomal_RNA_with_Chx.fasta` file will be used as a reference for chimera searching. 26810 entries including Chlorohelix.

#### UCHIME
Download UCHIME2 (free version)
```bash
cd "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro"

# usearch11.0.667_i86linux32
wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

gunzip usearch11.0.667_i86linux32.gz
chmod 755 usearch11.0.667_i86linux32
```

### Adapter trimming
Trim primers and adapters off the ends of both sequences. Flip sequences to the forward orientation at the same time.
```bash
source_dir="${work_dir}/03b_NanoCLUST_c10"
fwd_primer_27f="AGAGTTTGATCMTGGCTCAG" # Lane, 1991
rev_primer_1492r="TACGGYTACCTTGTTACGACTT" # Lane, 1991

mkdir -p "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/04b_c10_analyzed" && cd "$_"

# Gather inputs
mkdir -p "01_inputs" && cd "$_"
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_sequences.fasta" | sort -h | \
  xargs cat > cluster_sequences_raw.fasta
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_counts.tsv" | sort -h | \
  head -n 1 | xargs head -n 1 > cluster_counts_raw.tsv
find "${source_dir}" -mindepth 2 -maxdepth 2 -name "cluster_counts.tsv" | sort -h | \
  xargs tail -q -n +2 >> cluster_counts_raw.tsv
cd ..

# Get 1492r reverse complement sequence
conda activate general2 # seqtk 1.3-r106
rev_primer_1492r_revcomp=$(printf ">1492r\n${rev_primer_1492r}\n" | seqtk seq -r | tail -n 1)
conda deactivate

# Trim
mkdir -p "02_cutadapt" && cd "$_"
conda activate cutadapt_3.4
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

### OTU clustering
```bash
cd "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/04b_c10_analyzed"
mkdir -p "03_cluster" && cd "$_"

conda activate general2
# mmseqs 13.45111

mmseqs easy-cluster ../02_cutadapt/cluster_sequences_trimmed_complete.fasta p99c99 tmp \
  -c 0.99 --cov-mode 5 --cluster-mode 2 --min-seq-id 0.99 --threads 2 2>&1 | \
  tee p99c99.log

# collapsed to 14 representative sequences

ln -s p99c99_rep_seq.fasta repseqs.fasta

conda deactivate
cd ..
```

### Chimera filtration
Using UCHIME2 with the NCBI+Chx 16S DB made in the `Pre-prep` section of the `Processed outputs` region of this README file.

```bash
db_filepath="../../chimera_ref_db/16S_ribosomal_RNA_with_Chx.fasta"
usearch_filepath="../../usearch11.0.667_i86linux32"

mkdir -p "/Analysis/jmtsuji/projects/culture/2022_12_28_Chx_Nanopore_16S_summary_NatMicro/04b_c10_analyzed/04_chimera_removal" && cd "$_"

"${usearch_filepath}" -uchime2_ref "../03_cluster/repseqs.fasta" -db "${db_filepath}" \
  -mode high_confidence -strand plus -uchimeout chimera_report.tsv 2>&1 | \
  tee uchime2.log
```

Results (cleaned version of UCHIME2 output) (`cat chimera_report.tsv | sed 's/ Geothrix fermentans strain H5 16S ribosomal RNA, partial sequence//g' | sed 's/ Holophaga foetida strain TMBS4 16S ribosomal RNA, partial sequence//g' | sed 's/ 16S ribosomal RNA//g' | column -s $'\t' -t`):
```
20211112_barcode21_cluster0       0.0000   N  *                    *                    Chx_allophototropha                                                      dqt=0;
20220216_barcode12_cluster5       0.0000   ?  *                    *                    NR_180005.1 Acinetobacter oryzae strain B23, partial sequence            dqt=5;why=nodiv;
20220216_barcode08_cluster1       0.5267   ?  NR_036779.1          NR_036891.1          (L)                                                                      dqt=25;dqm=28;L=64,13,15(1436);R=12,0,0(12);div=-1.0%;
20220216_barcode08_cluster2       30.9145  Y  Chx_allophototropha  NR_036779.1          (L)                                                                      dqt=154;dqm=12;L=216,0,0(749);R=151,10,3(689);div=10.1%;
20220216_barcode11_cluster6       30.9949  Y  Chx_allophototropha  NR_036779.1          (R)                                                                      dqt=161;dqm=13;L=162,0,0(510);R=204,10,4(954);div=10.2%;
20220216_barcode08_cluster6 rc    16.7630  Y  Chx_allophototropha  NR_036779.1          (R)                                                                      dqt=103;dqm=28;L=89,0,0(252);R=262,12,17(1206);div=5.3%;
20220216_barcode11_cluster8       25.9251  Y  Chx_allophototropha  NR_036779.1          (L)                                                                      dqt=107;dqm=12;L=263,0,0(1022);R=104,10,3(380);div=6.8%;
20220216_barcode08_cluster7 rc    17.4133  Y  Chx_allophototropha  NR_036779.1          (L)                                                                      dqt=69;dqm=13;L=302,0,0(1189);R=66,11,3(275);div=4.1%;
20220216_barcode08_cluster12 rc   34.9517  Y  NR_036779.1          Chx_allophototropha  (L)                                                                      dqt=173;dqm=27;L=195,6,21(784);R=161,0,0(668);div=10.2%;
20220216_barcode08_cluster3       33.9782  Y  NR_036779.1          Chx_allophototropha  (R)                                                                      dqt=157;dqm=26;L=137,6,20(486);R=220,0,0(960);div=9.1%;
20220216_barcode08_cluster4       41.9234  Y  NR_036779.1          Chx_allophototropha  (R)                                                                      dqt=88;dqm=11;L=81,4,7(252);R=291,0,0(1206);div=5.5%;
20220216_barcode13_cluster4       35.2690  Y  NR_036779.1          Chx_allophototropha  (L)                                                                      dqt=192;dqm=27;L=176,6,21(662);R=180,0,0(791);div=11.5%;
20220216_barcode12_cluster0       0.0000   ?  *                    *                    NR_026163.1 Microbacterium testaceum strain DSM 20166, partial sequence  dqt=10;why=nodiv;
20220216_barcode14_cluster0 rc    0.0000   ?  *                    *                    NR_116570.1 Sphingomonas hankookensis strain ODN7, partial sequence      dqt=7;why=nodiv;
```
9/14 sequences were identified as chimeras. 5 remain.

Summarized the chimera-filtered sequence set (i.e., the sequences not identified as chimeras):
```bash
# This was done manually. TODO - automate if doing this in high-throughput in the future.
printf "20211112_barcode21_cluster0\n20220216_barcode12_cluster5\n20220216_barcode08_cluster1\n20220216_barcode12_cluster0\n20220216_barcode14_cluster0\n" > chimera_filtered.list

conda activate general2
# seqtk 1.3-r106

seqtk subseq "../03_cluster/repseqs.fasta" chimera_filtered.list > chimera_filtered.fasta
```

Validate the remaining two sequences against the known Chx and G1 reference seqs, prepared in this folder as `culture_refs.fasta`
```bash
# blast 2.12.0
blastn -query culture_refs.fasta -subject chimera_filtered.fasta
```
Result: 
- `20211112_barcode21_cluster0`: 100% identical to Chx
one of the two is 100% identical to Chx. The other is 100% identical to one of the two Geothrix 16S seqs (and just has 3 mismatches to the other)
- `20220216_barcode08_cluster1`: 100% identical to G1 rRNA gene `GEOCFX_001163`
- `20220216_barcode12_cluster5`: unclear
- `20220216_barcode12_cluster0`: unclear
- `20220216_barcode14_cluster0 rc`: unclear

#### Interpretation of unclear sequences
Then ran the 5 chimera-filtered seqs against NCBI 16S DB via web BLAST (Dec 28 2022). Results for the 3 latter genes with "unclear" classification when compared to Chx and G1 via BLAST above:
- `20220216_barcode12_cluster5`: near 100% identity to NR_180005.1	Acinetobacter oryzae
- `20220216_barcode12_cluster0`: near 100% identity to NR_026163.1	Microbacterium testaceum
- `20220216_barcode14_cluster0 rc`: near 100% identity to NR_116570.1	Sphingomonas hankookensis

Low relative abundances, looking at the `cluster_counts_raw.tsv` file:
```
sample	id	reads_in_cluster	used_for_consensus	reads_after_corr	draft_id	sciname	taxid	length	per_ident
20220216_barcode12	5	18	100	17	576ab204-ef44-4968-b7e4-9b39ee2af8de id=9	Acinetobacter johnsonii	40214	1501	99.467
20220216_barcode12	0	32	100	31	852e4339-43a8-497a-8bea-2f6dd9d17e01 id=19	Microbacterium testaceum	2033	1465	99.249
20220216_barcode14	0	43	100	41	d2ba43c4-3c27-421b-8f06-eed41d7e0224 id=35	Sphingomonas hankookensis	563996	1430	99.441
```
Reads in cluster ranged from 18-43 out of 40-50k total reads. This is less than 0.1% relative abundance.

### Relative abundance calculations
TODO

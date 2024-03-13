# Genomics analyses
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper

Copyright Jackson M. Tsuji, Neufeld Research Group, 2024

**NOTE: for each code section provided below, the code ought to be run from within this `genomics` directory.**

This document section is related to the `culture_metagenomics` section, which focused on generating genome bins from metagenome data of non-pure bacterial cultures. Here, I generate closed genomes of "_Ca_. Chx. allophototropha" L227-S17 and _Geothrix_ sp. L227-G1. I also combine that data with genome bin data to determine the relative abundances of microbial populations in the enrichment cultures over time. Lastly, I perform some analyses of specific genes in the "_Ca_. Chx. allophototropha" L227-S17 and strain L227-5C genomes.

## Closed genome of _Ca_ Chx. allophototropha L227-S17
The genome of "_Ca_ Chx. allophototropha" L227-S17 was closed using short and long read data from a colony picked from the S19.9 culture.

Note: all code for this section should be run in the `Chlorohelix_genome_final` folder.

### Get the read data
Raw read data can be downloaded as described in the README in the `source_data` section of this repo, specifically the 
part of the README that mentions the `Chx_closed_genome_raw_data_accessions.tsv` file.

If you just want to see the final genome, you can download it from GenBank accession `GCA_030389965.1` on NCBI.

Save downloaded raw reads to the `input` folder. Make sure the downloaded files are named as follows:
- `ChxS19_R1.fastq.gz`: R1 short read file
- `ChxS19_R2.fastq.gz`: R2 short read file
- `ChxS19_Nanopore.fastq.gz`: long read file

### QC of short read data
#### Adapter trimming
Firstly, because Nextera adapters are bound to both sides of the reads, Nextera adapter removal was performed using CutAdapt.

```bash
conda activate cutadapt_3.4 # install cutadapt version 3.4 here

mkdir -p "input/trimmed" && cd "$_"

cutadapt -o ChxS19_R1.fastq.gz -p ChxS19_R2.fastq.gz \
  -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --discard-untrimmed -j 30 \
  ../ChxS19_R1.fastq.gz ../ChxS19_R2.fastq.gz \
  2>&1 | tee cutadapt.log

cd ../..
```

#### Prepare for further adapter trimming by ATLAS
Adapter: CTGTCTCTTATACACATCT  -- this is a partial sequence in several adapters in the ATLAS adapters.fa file. However, to be safe, I made a custom version of the adapters.fa file including it specifically as an entry.

After running ATLAS once, you will find the `adapters.fa` file in the database folder. Simply add the following entry at the end of the file and save a copy called `adapters_with_custom_nextera.fa`:
```
>custom_nextera
CTGTCTCTTATACACATCT
```
TODO: it might also be advisable to add the BGI adapters in the future, although in theory these should occur after the Nextera adapters and be trimmed off.

#### Read QC by ATLAS
Initialize ATLAS
```bash
conda create -n atlas_2.8.2 -c bioconda -c conda-forge mamba python=3.8
conda activate atlas_2.8.2
mamba install -c bioconda -c conda-forge metagenome-atlas=2.8.2

mkdir -p short_read_QC && cd "$_"
# -d should be the path to the database folder on your machine.
atlas init -d /Data/databases/atlas/2.8.x -w . --threads 30 ../input/trimmed
```

Edited config:
- Added `simplejob_threads: 30`
- `mem: 200`
- `large_threads: 30`
- `assembly_threads: 30`
- `preprocess_adapters: /Data/atlas/2.8.x/adapters_with_custom_nextera.fa` # put the custom path to your adapter file here

Ran:
```bash
# Continuing in the same folder and conda env as the above code block

atlas run -w . -c config.yaml -j 30 -n qc --printshellcmds --keep-going 2>&1 | tee atlas_steps.log
atlas run -w . -c config.yaml -j 30 qc --printshellcmds --keep-going 2>&1 | tee atlas.log

cd ..
```
QC'ed reads are located in `short_read_QC/ChxS19/sequence_quality_control` in the following files: `ChxS19_QC_R1.fastq.gz`, `ChxS19_QC_R2.fastq.gz`, and `ChxS19_QC_se.fastq.gz`.

### Hybrid assembly pipeline (rotary)
Rotary is a custom pipeline designed for circularization and polishing of long read genomes (with optional short read data).

Install and configure:
```bash
cd .. # i.e., go back to the main `genomics` folder

git clone https://github.com/jmtsuji/rotary.git
cd rotary
git checkout e636236 # commit e636236

conda env create -n rotary_e636236 --file=envs/toyako.yaml

cp config.yaml ../Chlorohelix_genome_final/rotary.yaml

# Return to the analysis dir
cd ../Chlorohelix_genome_final
```

Edited config `rotary.yaml`:
* `sample_id: "ChxS19"`
* `longreads: "input/ChxS19_Nanopore.fastq.gz"`
* `qc_short_r1: "short_read_QC/ChxS19/sequence_quality_control/ChxS19_QC_R1.fastq.gz"`
* `qc_short_r2: "short_read_QC/ChxS19/sequence_quality_control/ChxS19_QC_R2.fastq.gz"`
* `db_dir: "/Data/databases/rotary"` # The custom path to where you want to save rotary database files
* `flye_meta_mode: "True"`

NOTE: you might need to use absolute paths (instead of relative paths) above for rotary to work.

Run the pipeline
```bash
conda activate rotary_e636236

run_directory="rotary" # Wherever you want to store the run files
config="rotary.yaml"
conda_prefix="/Data/databases/rotary/conda_envs" # Wherever you want to store the conda envs
snakefile="../rotary/rules/toyako.smk"
jobs=40

mkdir -p "${run_directory}"

snakemake --dryrun --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread_steps.log"

snakemake --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread.log"
```
Done! You should now have a closed and polished genome of Chx in the `rotary/circularize/circularize.fasta`.

### Post analysis

#### Genome annotation via PGAP
Using NCBI PGAP pipeline instead of DFAST included in rotary:

##### Install PGAP
Install CWL language
```bash
conda create -n pgap -c bioconda -c conda-forge cwl_runner=1 cwltool=3.1.20220124184855
# cwl_runner version 1.0
# cwltool version 3.1.20220124184855
```

Install PGAP (note that Docker is also needed on your computer):
```bash
mkdir ../pgap && cd "$_" # i.e., install it in the `genomics` folder

# Note: I think there might be an easier way to do this now - see https://github.com/ncbi/pgap/wiki/Quick-Start
pgap_version="2022-02-10.build5872"
wget -O /dev/stdout "https://github.com/ncbi/pgap/archive/refs/tags/${pgap_version}.tar.gz" | \
  tar -xzf -
# saves as directory "pgap-2022-02-10.build5872"

# Move executable to an easy-to-use place
ln -s "pgap-${pgap_version}/scripts/pgap.py" .

# Set up Docker container
docker pull "ncbi/pgap-utils:${pgap_version}"
```
I then recommend to test the pipeline as described at https://github.com/ncbi/pgap/wiki/Quick-Start. It should auto download the annotation DB during testing.

You might also (or alternatively) need to run `./pgap.py --update` to download the DBs etc.


##### Prepare input genome
Before starting annotation: Here is the desired naming and order of the contigs for the final assembly (decided on manually):
```
#seq_name       length  NEW-NAME	location	name (chromosome or plasmid-name)
contig_5        2961143 contig1	chromosome	1
contig_2        2453369 contig2	chromosome	2
contig_3        375235  contig3	plasmid	unnamed1
contig_4        241397  contig4	plasmid	unnamed2
contig_6        54947   contig5	plasmid	unnamed3
```
Here, the 1st column is the original contig name, and the 3rd column is the new name I will give it. The 4th and 5th columns show the metadata labels I will give to each FastA entry according to the PGAP wiki: https://github.com/ncbi/pgap/wiki/Input-Files

Note: the original contig names and lengths might differ slightly in your rotary output compared to mine.

Manually modify contig names to match the above table
```bash
mkdir -p "post-analysis/pgap/input" && cd "$_"

# Semi-manually re-order and add labels
# Using seqtk 1.3-r106
seqtk seq ../../../rotary/circularize/circularize.fasta > tmp.fasta
grep -A 1 contig_5 tmp.fasta > ordered.fasta
grep -A 1 contig_2 tmp.fasta >> ordered.fasta
grep -A 1 contig_3 tmp.fasta >> ordered.fasta
grep -A 1 contig_4 tmp.fasta >> ordered.fasta
grep -A 1 contig_6 tmp.fasta >> ordered.fasta
rm tmp.fasta
cat ordered.fasta | sed 's/^>contig_5$/>contig1 [location=chromosome] [chromosome=1]/g' > renamed.fasta
sed -i 's/^>contig_2$/>contig2 [location=chromosome] [chromosome=2]/g' renamed.fasta
sed -i 's/^>contig_3$/>contig3 [location=plasmid] [plasmid-name=unnamed1]/g' renamed.fasta
sed -i 's/^>contig_4$/>contig4 [location=plasmid] [plasmid-name=unnamed2]/g' renamed.fasta
sed -i 's/^>contig_6$/>contig5 [location=plasmid] [plasmid-name=unnamed3]/g' renamed.fasta
rm ordered.fasta
seqtk seq -l 60 renamed.fasta > ChxS19_labelled.fasta
rm renamed.fasta

# Confirm everything looks OK
grep -n "^>" ChxS19_labelled.fasta
#1:>contig1 [location=chromosome] [chromosome=1]
#49355:>contig2 [location=chromosome] [chromosome=2]
#90247:>contig3 [location=plasmid] [plasmid-name=unnamed1]
#96503:>contig4 [location=plasmid] [plasmid-name=unnamed2]
#100528:>contig5 [location=plasmid] [plasmid-name=unnamed3]

cd ..
ln -s input/ChxS19_labelled.fasta .

cd ../..
```

##### Prepare PGAP config files
As described at: https://github.com/ncbi/pgap/wiki/Input-Files

You can find customized input files in the `Chlorohelix_genome_final` directory: 
- input.yaml # generic metadata file - made via manual input
- submol.yaml # annotation yaml - made via manual input. **For now, had to label the genus_species as 'Chloroflexus' because we don't yet have a valid name and 'Chloroflexi bacterium' is not accepted.
  - See https://github.com/ncbi/pgap/issues/164 -- this is probably an OK course of action.
  - Note that I streamlined the info in the `contact_info` section because this is a public repo.

TODO: confirm the `submol.yaml` version I put in this folder is up to date.

##### Perform the annotation
Note: PGAP should auto install/update the database. When I ran PGAP, it used DB version `2022-02-10.build5872`.

```bash
mkdir -p "post-analysis/pgap" && cd "$_"

conda activate pgap
pgap_dir="../../../../pgap" # Where PGAP was installed

"${pgap_dir}/pgap.py" -n -o ChxS19 -c 40 -m 200g input.yaml 2>&1 | tee ChxS19_pgap.log
# Takes several hours

conda deactivate
cd ../..
```
Done! You should now have genome annotation files at `post-analysis/pgap/ChxS19`.

Go ahead and copy `annot.fna` to the main `Chlorohelix_genome_final` folder as `Ca_Chx_allophototropha_L227-S17.fna` for reference.

#### Run CheckM
```bash
conda activate /Data/reference_databases/atlas/2.2.x/conda_envs/91946fa0
# checkm-genome=1.0.18   same as used in ATLAS 2.2.0; but you can make a custom conda env with the same checkM version and set/download the corresponding DB

mkdir -p "post-analysis/checkm/input" && cd "$_"

ln -s ../../pgap/ChxS19/annot.faa Ca_Chx_allophototropha_L227-S17.faa
cd ..

# --genes flag makes it expect a FAA file
checkm lineage_wf --genes --tab_table -x faa -t 40 --pplacer_threads 40 input output \
  > completeness.tsv 2> checkm_manual.log

# Oops - the above command still saves log info into completeness.tsv. Had to clean up:
grep -v "^\[" completeness.tsv > completeness_clean.tsv

conda deactivate
cd ..
```

Result:
```
Bin Id                           Marker lineage         # genomes  # markers  # marker sets  0  1    2  3  4  5+  Completeness  Contamination  Strain heterogeneity
Ca_Chx_allophototropha_L227-S17  k__Bacteria (UID1452)  924        160        109            8  150  2  0  0  0   95.66         1.38           0.00
```
This genome was then used for other downstream analyses (e.g., in the `phylogenomics` folder) for the paper.

## Closed genome of _Geothrix_ L227-G1
This is a slightly more complicated version of the above analysis. Unlike the above analysis, the L227-G1 culture was not pure, so genome binning is required.

Note: all code in this section should be run from the `Geothrix_genome_bin_final` folder.

### Get the read data
Download instructions are available in the README file in the `source_data` section of this Github repo.

Short read metagenome read files, from the `Capt_S15` metagenome, need to be downloaded using the subsection of the 
README that mentions the `culture_early_metagenome_data_accessions.tsv` file. For this section, we will assume that you 
have already downloaded these files and performed short read QC on them as described in the `culture_metagenomics` 
section of this repo - see below.

Meanwhile, the long read (Nanopore) data from the `ChxS15c_nanopore` metagenome, need to be downloaded using the 
subsection of the README that mentions the `Geothrix_closed_genome_bin_raw_data_accessions.tsv` file.

Save the raw Nanopore read data to the `input` folder.

(Note: if you just want to get the final assembled L227-G1 genome, you can download it from GenBank accession 
`GCA_030219325.1` on NCBI.)

### Short read QC
These reads are the same ones that were QC processed in the `culture_metagenomics` section, specifically the `Capt_S15` 
sample, as mentioned above. See the README in `culture_metagenomics` for more. Here, we will use the QC'ed `Capt_S15` 
read files that should have been generated in the `culture_metagenomics` section of the code via ATLAS 2.2.0.

Copy the two main QC'ed files (R1 and R2) to `short_read_qc` as `ChxS15_QC_R1.fastq.gz` and `ChxS15_QC_R2.fastq.gz`.

### Hybrid genome assembly
Assumes rotary was already installed during the Chx genome analysis above.

Update rotary version
```bash
# Update the rotary commit
cd ../rotary # see the Chx genome section above for the inital download of this repo

git checkout fd5acee # commit fd5acee

conda env create -n rotary_fd5acee --file=envs/toyako.yaml

cd ../Geothrix_genome_bin_final
```

Get config file
```bash
cp ../rotary/config.yaml rotary.yaml
```

Edited config:
* `sample_id: "L227G1"`
* `longreads: "input/ChxS15c_nanopore.fastq.gz"`
* `qc_short_r1: "short_read_QC/ChxS15_QC_R1.fastq.gz"`
* `qc_short_r2: "short_read_QC/ChxS15_QC_R2.fastq.gz"`
* `db_dir: "/Data/databases/rotary"` # where the rotary DB files are stored
* `flye_meta_mode: "True"`

Ran
```bash
conda activate rotary_fd5acee

run_directory="rotary" # Wherever you want to store the run files
config="rotary.yaml"
conda_prefix="/Data/databases/toyako/conda_envs" # Wherever you want to store the conda envs
snakefile="../rotary/rules/toyako.smk"
jobs=40

mkdir -p "${run_directory}"

# Ran until the end of the circularize module (no annotation)
snakemake --until circularize --dryrun --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread_steps.log"

snakemake --until circularize --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread.log"
```

### Genome binning
Using the assembled contig output of rotary at the end of the circularization module.

Map short reads to get relative abundances
```bash
cd "rotary"
mkdir -p "binning" && cd "$_"

# Custom params
contigs=../circularize/circularize.fasta 
short_r1=../../short_read_QC/ChxS15_QC_R1.fastq.gz
short=r2=../../short_read_QC/ChxS15_QC_R2.fastq.gz

mkdir -p "mapping"

conda activate general
# bbmap 37.99, samtools 1.13 - install these tools in the conda env

# map short reads
bbmap.sh -Xmx50G \
  ref="${contigs}" \
  in1="${short_r1}" \
  in2="${short_r2}" \
  outm="mapping/short_reads.sam" \
  minid=0.95 \
  pairedonly=f \
  ambiguous=best \
  sam=1.4 \
  nodisk=t \
  threads=30 2>&1 | \
  tee "mapping/short_reads.log"
samtools view -b -@ 4 "mapping/short_reads.sam" | \
  samtools sort -@ 4 -m 10G > "mapping/short_reads.bam"
samtools index "mapping/short_reads.bam"
rm "mapping/short_reads.sam"

cd ../..

conda deactivate
```

Perform binning step
```bash
mkdir -p "rotary/binning" && cd "$_"

mkdir -p "bins"
conda activate /Data/databases/atlas/2.8.x/conda_envs/0e8761dffcf86d0af63a275926c80d10 # only specified install is for metabat2=2.15, with conda-forge and bioconda and defaults channels; you can install this in your own conda env instead of using a conda env from ATLAS if you'd like

# metabat2 2.15
jgi_summarize_bam_contig_depths --outputDepth bins/contig_depths.txt mapping/short_reads.bam 2>&1 | tee bins/jgi_bam.log
metabat2 -a bins/contig_depths.txt -i "${contigs}" -o bins/MAG -t 30 2>&1 | tee metabat2_manual.log

conda deactivate

# Check out the composition of the bins
grep "^>" bins/*.fa # In my analysis, MAG3 seems to be Geothrix (contig_10 alone, 3.74 Mb)
cd ../..
```

Here is the contents of the contig depths file:
```
contigName  contigLen    totalAvgDepth  short_reads.bam  short_reads.bam-var
contig_10   3.73387e+06  45.6685        45.6685          294.912
contig_33   2.95871e+06  39.0052        39.0052          305.192
contig_7    2.45341e+06  32.8918        32.8918          172.73
contig_81   241391       14.1521        14.1521          62.9304
contig_82   375276       20.397         20.397           109.92
contig_87   54957        11.1464        11.1464          40.8853
contig_85   8142         4.85636        4.85636          10.6192
```

Here is the bin allocation info (from the grep command above):
```
bins/MAG.1.fa:>contig_81
bins/MAG.1.fa:>contig_82
bins/MAG.1.fa:>contig_87
bins/MAG.2.fa:>contig_33
bins/MAG.2.fa:>contig_7
bins/MAG.3.fa:>contig_10
```
Will continue with the Geothrix bin, `MAG.3.fa`. (I did several checks in the end and know this is the Geothrix genome and that it is a closed circular contig. MAG 2 seems to be the Chlorohelix main genome, and MAG 1 seems to be chromids/plasmids from Chlorohelix.)

```bash
cd "rotary/binning"
ln -s bins/MAG.3.fa L227G1.fasta
cd ../..
```

Swapped this FastA into the circularize output folder so that I can use it for annotation downstream
```bash
cd "rotary/circularize"
mv circularize.fasta circularize_all.fasta
ln -s ../binning/L227G1.fasta circularize.fasta
cd ../..
```

### Finish genome assembly/annotation pipeline
Continued pipeline for annotation steps using the genome bin:
```bash
conda activate rotary_fd5acee

run_directory="rotary" # Wherever you want to store the run files
config="rotary.yaml"
conda_prefix="/Data/databases/rotary/conda_envs" # Wherever you want to store the conda envs
snakefile="../rotary/rules/toyako.smk"
jobs=40

mkdir -p "${run_directory}"

## Commit fd5acee
snakemake --dryrun --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread_steps2.log"

snakemake --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee "${run_directory}/genome_longread2.log"
```
Finished without errors.

### Post analysis
#### PGAP annotation
Using NCBI PGAP pipeline - see the Chlorohelix section above for installation.

Annotate:
```bash
conda activate pgap # See Chlorohelix section for creation of this env
pgap_dir="../pgap" # See Chlorohelix section for download

mkdir -p "post-analysis/pgap" && cd "$_"

cp ../../rotary/binning/L227G1.fasta .
# MANUALLY adjust sequencer header to >contig1

# See https://github.com/ncbi/pgap/wiki/Input-Files - used template there
# See input.yaml and submol.yaml files added to the Geothrix_genome_bin_final folder.

# Now at PGAP version 2022-04-14.build6021 (was previously at 2022-02-10.build5872)

# Ran
"${pgap_dir}/pgap.py" -n -o L227G1 -c 30 -m 200g input.yaml 2>&1 | tee L227G1_pgap.log
# Takes several hours to run

conda deactivate
cd ../..
```
Done! You can go ahead and copy `post-analysis/pgap/L227G1/annot.fna` to `Geothrix_L227-G1.fna` in the main `Geothrix_genome_bin_final` folder.

## Relative abundance calculations of culture members over time
Unless otherwise indicated, please run all code in this section in the `culture_MAGs_intermediate` folder.

### Download the metagenome and genome bin data
To calculate relative abundances, you will need to download the sequences of the genome bins obtained from culture 
metagenomes that were analyzed in this study. You'll also need to download raw read files for the metagenomes. To do 
so, please see the `README` file in the `source_data` folder, along with the 
`culture_early_genome_bin_data_accessions.tsv` and `culture_early_metagenome_data_accessions.tsv` files in that folder.
Please also download the L227-S17 and L227-G1 final genome bins as described in the `source_data` README file, and save 
the fna.gz files with meaningful names: `L227-S17-closed.fna.gz` and `L227-G1-closed.fna.gz`, respectively.

After finishing download:
- Genome bins files: save to `input/genomes`. 
- Source metagenome read files: save to `input/reads`.

### Cluster and organize the MAGs
#### Calculate ANI
```bash
mkdir -p "input/genomes_clustered/ani" && cd "$_"

conda activate general2
# fastani version 1.33

find "../../genomes" -name "*.fna.gz" | sort -h > inputs.list

fastANI --ql inputs.list --rl inputs.list -t 10 -o fastani_raw.txt 2>&1 | tee fastani.log

cat fastani_raw.txt | \
  sed 's\../../genomes\\g' | \
  sed 's/\.fna\.gz//g' \
  > fastani_ren.tsv

cd ../../..
```

Based on the output (`fastani_ren.tsv`), only two genome pairs had >99% ANI:
- `L227-S17-closed` and `L227-S17` (99.9607% ANI) # L227-S17 genome
- `L227-G1-closed` and `Cfx3-01` (99.9677% ANI) # L227-G1 genome

For both of these, I will take the closed genome as the representative.

#### Summarize clustered MAGs
```bash
mkdir -p "input/genomes_clustered/clustered" && cd "$_"

# Manually weed out the two old MAGs of Chx and Geo
ln -s ../../genomes/*.fna.gz .
rm L227-S17.fna.gz Cfx3-01.fna.gz

cd ../../..
```
14 total

#### Rename MAGs
Just clarify the `-closed` MAGs by removing this suffix

```bash
mkdir -p "input/genomes_clustered/final" && cd "$_"

ln -s ../clustered/*.fna.gz .

mv L227-S17-closed.fna.gz L227-S17.fna.gz
mv L227-G1-closed.fna.gz L227-G1.fna.gz

cd ../../..
```

#### Get contig IDs
Got contig IDs associated with each MAG

```bash
mkdir -p "input/genomes_clustered/final" && cd "$_"

genomes=($(find -L . -name "*.fna.gz" | sort -h))

printf "genome\tcontig\n" > contig_mapping.tsv

for genome in "${genomes[@]}"; do

  genome_basename=${genome##*/}
  genome_basename=${genome_basename%.fna.gz}

  contigs=($(zgrep "^>" "${genome}" | cut -d " " -f 1 | cut -d ">" -f 2))

  echo "[ $(date -u) ]: ${genome_basename}: ${#contigs[@]} contigs"

  for contig in "${contigs[@]}"; do
    printf "${genome_basename}\t${contig}\n" >> contig_mapping.tsv
  done

done

# Semi-manual check for duplicate contig names
tail -n +2 contig_mapping.tsv | cut -f 2 | wc -l # 2866
tail -n +2 contig_mapping.tsv | cut -f 2 | sort -u | wc -l # 2866
# Row numbers match before and after duplicate filtration, so there are no duplicates contig IDs

cd ../../..
```

#### Make mapping reference
```bash
cd "input/genomes_clustered"

ln -s final/contig_mapping.tsv .

conda activate general
# seqtk 1.3-r106 should be installed in this conda env

# Standardize lines per row for FastA when combining
# Also drop FastA comments
find -L final -maxdepth 1 -name "*.fna.gz" | sort -h | xargs -I {} seqtk seq -C -l 80 {} | gzip > refs_combined.fna.gz

cd ../..
```

### Prepare reads for mapping
(Note: Source metagenome reads should be saved to `input/reads` - see the download section above)

#### Redo QC on the downloaded reads
```bash
mkdir -p "qc/short_batch1" && cd "$_"

conda activate atlas_2.8.2 # Same ATLAS version as used for Chx S19.9 short reads (above)

# Set the -d path to where your database files are stored
atlas init -w . -d /Data/databases/atlas/2.8.x --threads 40 ../../input/reads
```

Edit config.yaml:
- `large_mem`: 60
- `assembly_memory`: 60
- Added: `simplejob_threads: 40`

samples.tsv: no changes.

Ran ATLAS
```bash
# continuing in the same folder as above
conda activate atlas_2.8.2

atlas run -w . -c config.yaml -j 40 -n qc --printshellcmds --reason 2>&1 | tee atlas_230208_1635_steps.log
atlas run -w . -c config.yaml -j 40 qc --printshellcmds --reason 2>&1 | tee atlas_230208_1635.log

cd ../..
```
Finished without errors.

#### Get QC'ed reads from projects already analyzed in this section of the GitHub repo
```bash
mkdir -p "qc/short_batch2" && cd "$_"

ln -s ../../../Chlorohelix_genome_final/short_read_QC/ChxS19/sequence_quality_control/ChxS19_QC_R1.fastq.gz .
ln -s ../../../Chlorohelix_genome_final/short_read_QC/ChxS19/sequence_quality_control/ChxS19_QC_R2.fastq.gz .
ln -s ../../../Chlorohelix_genome_final/short_read_QC/stats/read_counts.tsv ChxS19_read_counts.tsv

cd ..

mkdir -p "qc/long" && cd "$_"

ln -s ../../../Chlorohelix_genome_final/rotary/qc_long/nanopore_qc.fastq.gz ChxS19_Nanopore.fastq.gz
grep "^Input:" ../../../Chlorohelix_genome_final/rotary/logs/qc_long.log | tr -s " " $'\t' | cut -f 2-3 # 509806 reads
grep "^Output:" ../../../Chlorohelix_genome_final/rotary/logs/qc_long.log | tr -s " " $'\t' | cut -f 2-3 # 131256 reads

ln -s ../../../Geothrix_genome_bin_final/rotary/qc_long/nanopore_qc.fastq.gz ChxS15c_Nanopore.fastq.gz
grep "^Input:" ../../../Geothrix_genome_bin_final/rotary/logs/qc_long.log | tr -s " " $'\t' | cut -f 2-3 # 230577 reads
grep "^Output:" ../../../Geothrix_genome_bin_final/rotary/logs/qc_long.log | tr -s " " $'\t' | cut -f 2-3 # 75646 reads

cd ../..
```

Made total read table for long reads, called `read_counts_Nanopore.tsv`, in the `long` folder:
```
Sample	Step	Reads_se
ChxS15c	raw	230577
ChxS15c	QC	75646
ChxS19	raw	509806
ChxS19	QC	131256
```

#### Make summary dir
Just use R1 and R2 files for short reads.

I will also standardize the file names here with the correct subculture codes

Note - please manually confirm that the names shown below are correct for your analysis before performing these 
commands. For example, `Capt-S01` might be named `Capt_S01` in your analysis.
```bash
mkdir -p "qc/short_summary" && cd "$_"

## Short reads
ln -s ../short_batch1/Capt-S01/sequence_quality_control/Capt-S01_QC_R1.fastq.gz L227-S17_S01_R1.fastq.gz
ln -s ../short_batch1/Capt-S01/sequence_quality_control/Capt-S01_QC_R2.fastq.gz L227-S17_S01_R2.fastq.gz

ln -s ../short_batch1/Capt-S15/sequence_quality_control/Capt-S15_QC_R1.fastq.gz L227-S17_S15_R1.fastq.gz
ln -s ../short_batch1/Capt-S15/sequence_quality_control/Capt-S15_QC_R2.fastq.gz L227-S17_S15_R2.fastq.gz

ln -s ../short_batch1/L227-5C/sequence_quality_control/L227-5C_QC_R1.fastq.gz L227-5C_S00_R1.fastq.gz
ln -s ../short_batch1/L227-5C/sequence_quality_control/L227-5C_QC_R2.fastq.gz L227-5C_S00_R2.fastq.gz

ln -s ../short_batch2/ChxS19_QC_R1.fastq.gz L227-S17_S19_R1.fastq.gz
ln -s ../short_batch2/ChxS19_QC_R2.fastq.gz L227-S17_S19_R2.fastq.gz

# Read stats
head -n 1 ../short_batch1/stats/read_counts.tsv > header.1.tmp
cat ../short_batch1/stats/read_counts.tsv | \
  sed 's/Capt-/L227-S17_/g' | \
  sed 's/L227-5C/L227-5C_S00/g' \
  > read_counts.tsv
head -n 1 ../short_batch2/ChxS19_read_counts.tsv > header.2.tmp
# MANUAL confirmation required here
cmp header.1.tmp header.2.tmp # identical, so can merge the tables
rm header.1.tmp header.2.tmp
cat ../short_batch2/ChxS19_read_counts.tsv | \
  sed 's/ChxS19/L227-S17_S19/g' | \
  tail -n +2 \
  >> read_counts.tsv

cd ../..

## Long reads
mkdir -p "qc/long_summary" && cd "$_"
ln -s ../long/ChxS15c_Nanopore.fastq.gz L227-S17_S15c_long.fastq.gz
ln -s ../long/ChxS19_Nanopore.fastq.gz L227-S17_S19_long.fastq.gz

# Read stats
cat ../long/read_counts_Nanopore.tsv | \
  sed 's/ChxS15c/L227-S17_S15c_long/g' | \
  sed 's/ChxS19/L227-S17_S19_long/g' \
  > read_counts.tsv

cd ../..
```
Now ready for read mapping

### Map the reads
#### Short reads

Map using bbmap version 37.99
```bash
mkdir -p "mapping/short" && cd "$_"

# Re-used a conda env from rotary
conda activate /Data/databases/rotary/conda_envs/2b03926f2933a42eda0b015b96779e22
# bbmap 37.99, samtools 1.15.1
# You can make a similar conda env yourself with the same tool versions

bbmap.sh ref=../../input/genomes_clustered/refs_combined.fna.gz 2>&1 | tee build_ref.log

metagenomes_R1=($(find -L ../../qc/short_summary -name "*_R1.fastq.gz" | sort -h))

printf "sample\tcontig\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > short_read_coverages.tsv
for metagenome_R1 in "${metagenomes_R1[@]}"; do

  metagenome_R2="${metagenome_R1%_R1.fastq.gz}_R2.fastq.gz"
  metagenome_basename="${metagenome_R1%_R1.fastq.gz}"
  metagenome_basename="${metagenome_basename##*/}"

  echo "[ $(date -u) ]: ${metagenome_basename}"

  # Map the reads
  bbmap.sh -Xmx50g in="${metagenome_R1}" in2="${metagenome_R2}" out="${metagenome_basename}.sam" \
    unpigz=t minid=0.90 threads=40 ambiguous=best pairedonly=f secondary=f > "${metagenome_basename}_bbmap.log" 2>&1

  printf "" > "${metagenome_basename}_samtools.log"
  # Reformat output
  samtools view -b -@ 10 "${metagenome_basename}.sam" 2>>"${metagenome_basename}_samtools.log" | \
    samtools sort -m 10G -@ 4 2>>"${metagenome_basename}_samtools.log" \
    > "${metagenome_basename}.bam"
  samtools index -@ 10 "${metagenome_basename}.bam" 2>>"${metagenome_basename}_samtools.log"
  rm "${metagenome_basename}.sam"

  # Calculate coverages
  samtools coverage "${metagenome_basename}.bam" > "${metagenome_basename}_coverage.tsv"
  cat "${metagenome_basename}_coverage.tsv" | tail -n +2 | \
    sed "s/^/${metagenome_basename}\t/g" >> short_read_coverages.tsv

done

cd ../..
```

#### Notes about the purity of the Chx S19 culture
After finishing the above read mapping, I checked the S19 short read culture sample more carefully to see how pure the culture was. 
I noticed that a few contigs (I guess 40 or something) from other genome bins (especially bin Chx3-03) had mapped reads from the Chx S19 metagenome. Could this be a sign that the S19 culture was not "pure", or could it mean that those genome bins had contamination in them from the true Chx genome? 

I checked the ~40 contigs in question, and they are almost all <5 kb and have >99% identity across >90% of the sequence to the closed Chx genome. Thus, they are likely cross-binning contamination of Chx genome fragments into the old MAG set.

I searched the fragments against the Chx genome. Almost all have a top hit of >=96.845% identity and 100% query coverage. Many have multiple copies. So maybe they were repetitive elements that didn't assemble well in the Chx genome?

Also, note that only 28 reads mapped to the L227-G1 genome. This is 0.00086% of the mapped reads... very low.

See the Excel file in this folder, `Chx-S19-short-contamination-check`, where I show that >99.9% of reads mapping to Cfx3-03 in the end are coming from these possible Chlorohelix contigs. Within that read set, I also show that contigs with >99% identity across 100% of their sequence to Chlorohelix allophototropha make up 96.86% of reads mapping to the Cfx-03 genome bin from the Chx S19 sample.

Conclusion: I can have high confidence that the Chx S19 metagenome was almost entirely Chx.

### Long reads
```bash
mkdir -p "mapping/long" && cd "$_"

# Env from rotary
conda activate /Data/databases/rotary/conda_envs/ec1a2684af35daa973f4f5bb6d6ab4cb
# minimap2 version 2.23, samtools version 1.15
# You could make a similar conda env having tools with the same versions as above

metagenomes=($(find -L ../../qc/long_summary -name "*.fastq.gz" | sort -h))

printf "sample\tcontig\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > long_read_coverages.tsv
for metagenome in "${metagenomes[@]}"; do

  metagenome_basename="${metagenome%.fastq.gz}"
  metagenome_basename="${metagenome_basename##*/}"

  echo "[ $(date -u) ]: ${metagenome_basename}"

  printf "" > "${metagenome_basename}_samtools.log"
  # Note: I exclude secondary and supplementary alignments so that I only get a max of 1 mapping per read. This is to get more accurate relative abundance measures downstream, but it is not good for other purposes like looking at genome quality etc. It's like using only the best hit as I've done in bbmap.
  # Help from https://github.com/lh3/minimap2/issues/416#issuecomment-499095442, accessed Feb. 8th, 2023
  minimap2 -t 40 -ax map-ont ../01f_mapping_reference_cat/refs_combined.fna.gz \
    "${metagenome}" 2> "${metagenome_basename}_minimap2.log" | \
    samtools view -b -@ 10 -F 0x900 2>> "${metagenome_basename}_samtools.log" | \
    samtools sort -@ 4 -m 10G 2>> "${metagenome_basename}_samtools.log" \
    > "${metagenome_basename}.bam"
  
  samtools index -@ 10 "${metagenome_basename}.bam" 2>> "${metagenome_basename}_samtools.log"

  samtools coverage "${metagenome_basename}.bam" > "${metagenome_basename}_coverage.tsv"
  cat "${metagenome_basename}_coverage.tsv" | tail -n +2 | \
    sed "s/^/${metagenome_basename}\t/g" >> long_read_coverages.tsv

done

cd ../..
```

#### Examine the purity of the Chx S19 culture of Chx via long read data
In general, all contigs from non-Chx genomes that had mapped reads overlapped with the set of contigs analyzed in the short read check above. There wre two exceptions:

```
#rname             startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
JACAAX010000047.1  1         53500   1         762       1.4243    0.014243   24         8
JACAUB010000123.1  1         1524    5         1524      100       3.3084     21         60
```
The 53 kb contig is potentially a spurious mapping due to only having a single mapped read. The 1.5 kb one is actually a 100% match with 100% qcovhsp against viral genomes based on online BLASTN! It is included with 100% match, 100% qcov in the genomes of many refseq entries as well (E. coli etc.). It match perfectly against a portion of the lambda phage genome (NC_001416.1), which was included during library prep as a spike-in. Thus, I conclude that the second contig is some kind of viral DNA fragment from the lambda phage DNA added to the library. Not real contamination.

Overall, then, the long read data also support that the Chx S19 colony was highly pure.

### Calculate relative abundances
```bash
mkdir -p "mapping/calculations" && cd "$_"

mkdir input && cd "$_"

ln -s ../../input/genomes_clustered/contig_mapping.tsv genome_contig_names.tsv
ln -s ../../qc/short_summary/read_counts.tsv read_counts_short.tsv
ln -s ../../qc/long_summary/read_counts.tsv read_counts_long.tsv
ln -s ../short/short_read_coverages.tsv coverage_short.tsv
ln -s ../long/long_read_coverages.tsv coverage_long.tsv

cd ../../..
```
Note that `GTDB_and_checkM_stats.tsv` is included in the main directory for `culture_MAGs_intermediate` and will also be used to generate summary files. Please copy this to the `input` folder above.

Then analyzed the data using the Jupyter notebook in this folder, `analyze_relabund.ipynb`.

Produces `MAG_rel_abundances.pdf` and `MAG_rel_abundances.tsv` in the `output` folder. The PDF was cleaned up to make a panel in Extended Data Fig. 1.

## Misc genome analyses
### Custom hidden Markov models
Custom hidden Markov models (HMMs) were made for the PscA-like and FmoA genes of '_Ca_. Chloroheliales' members:
```bash
# Dependencies you need to install:
# clustalo v1.2.3
# hmmbuild v3.1b2

input_dir="sequences_of_interest"
output_dir="hidden_markov_models"

gene=pscA
clustalo -i "${input_dir}/${gene}.faa" -o "${output_dir}/${gene}_aln.faa" -v -v \
  2>&1 | tee "${output_dir}/${gene}_aln.log"
hmmbuild -O "${output_dir}/${gene}_aln.stk.annotated" "${output_dir}/${gene}_Chloroheliales.hmm" \
  "${output_dir}/${gene}_aln.faa" 2>&1 | tee "${output_dir}/${gene}_hmmbuild.log"

gene=fmoA
clustalo -i "${input_dir}/${gene}.faa" -o "${output_dir}/${gene}_aln.faa" -v -v \
  2>&1 | tee "${output_dir}/${gene}_aln.log"
hmmbuild -O "${output_dir}/${gene}_aln.stk.annotated" "${output_dir}/${gene}_Chloroheliales.hmm" \
  "${output_dir}/${gene}_aln.faa" 2>&1 | tee "${output_dir}/${gene}_hmmbuild.log"
```

These HMMs are crude because they are only made with two sequences as input. 
They can be expanded further as the diversity of this clade becomes better understood.

The HMM files are pre-made and provided in the `hidden_markov_models` directory for reference.

### Summary of key photosynthesis-associated genes in '_Ca_. Chloroheliales' genomes
FastA files (non-aligned) are provided with the sequences of key photosynthes-associated genes and their predicted ORFs for the two genome bins. 
The genes summarized here are:
- _pscA_ encoding one of the core Type I reaction center components
- _fmoA_ encoding a the FMO energy-transfer protein between Type I reaction centers and chlorosomes
- _csmA_ encoding a key structural protein for chlorosomes
- _cbbL_ encoding the RuBisCO large subunit gene

FastA files can be found in the `sequences_of_interest` directory. `.ffn` files contain gene nucletide sequences. 
`.faa` files contain predicted amino acid sequences.

### Homology modelling (protein structural prediction)
Homology modelling was performed using I-TASSER for predicted amino acid sequences of the key genes above, with default settings.

The top-scoring predicted tertiary structure for each gene (PDB format) are available in the `homology_models` directory. 
These PDB files can be opened in a protein structural viewer such as [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/).

The full I-TASSER output for each gene can be found at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110).

One of these homology models ('_Ca_. Chx. allophototropha' PscA-like protein) is shown an Extended Data figure.

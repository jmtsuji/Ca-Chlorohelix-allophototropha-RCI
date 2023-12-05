# Culture metagenome analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

**NOTE: for each code section provided below, the code ought to be run from within this `culture_metagenomics` directory.**

## Overview
This analysis folder describes:
- Metagenome analysis of early to mid-phase enrichment cultures of "_Ca_. Chloroheliales" members
  - Note that the Ca. Chx. allophototropha genome was ultimately derived from a different dataset than this. See `genomics/Chlorohelix_genome_final`
  - However, these analyses resulted in the final genome bin for strain L227-5C that was further analyzed in the manuscript.
- Semi-related, minor analysis: searching for evidence of the "_Ca_. Chx. allophototropha"-like pscA-like gene in environmental metagenomes from other studies

## Short read or read cloud metagenomes for early- to mid-stage enrichment cultures

### Data download
See the `source_data` folder for details, specifically `culture_early_metagenome_data_accessions.tsv` and the `README` there.

### Install ATLAS
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


### Process the metagenome data
If it is your first time using ATLAS, note that a large amount of database files will be auto-downloaded during the run.

#### Ca. Chloroheliaceae bin L227-5C
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

#### Ca. Chx. allophototropha
Config files are also provided.

This case is slightly more complicated, though. 
We will use the Tell-Read QC processor and Tell-Link assembler to assemble the data ahead of time and then put these data into ATLAS partway through.

##### Assemble read cloud data separately
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
# No underscores in output_dir name!
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

##### Run ATLAS
Notes:
- make sure you modify the config.yaml file in the `atlas_L227_5C` folder so that `database_dir` is a real directory on your machine.
- you might have to modify the filepaths for the samples in the `samples.tsv` file so that they are correct on your system.

```bash
# Activate the environment by running: 
# conda activate atlas_2.2.0

cd "atlas_Chx_allophototropha"

# Finish the QC step - note this should use the downloaded short reads from NCBI rather than those processed via Tell-Link
# Note that the short read in NCBI are the ones demultiplexed via Picard
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


### Curation of the genome bins
I identified the two '_Ca_. Chloroheliales' genome bins from the above ATLAS run folders by looking for the bins classified to the _Chloroflexota_ phylum based on the `genomes/taxonomy/gtdb/gtdbtk.bac120.summary.tsv` file in each folder.

Then, I manually removed suspect scaffolds as described in the paper.

Curated genomes are available from NCBI, annotated using PGAP during submission. 
See download instructions in the README file in `source_data` section of this repo, specifically in the part of the 
README file that mentions the `culture_early_genome_bin_data_accessions.tsv` TSV file.
**Note**: the Ca. Chx. allophototropha genome analyzed in the manuscript was ultimately derived from a different dataset than this. See `genomics/Chlorohelix_genome_final`.

Done! These data will be used again later in the `genomics` folder under `culture_MAGs_intermediate`.


## Scanning metagenomes containing potential '_Ca._ Chloroheliales' members for photosynthesis genes
This was only briefly mentioned in the paper, but because two genomes from the GTDB were classified to the same order as '_Ca_. Chx. allophototropha', I examined the corresponding raw metagenome files that those two genome bins were derived from to search for signs of Type I reaction center-associated genes. I used the custom HMMs developed as decribed in the `genome_bin_analysis` folder.

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

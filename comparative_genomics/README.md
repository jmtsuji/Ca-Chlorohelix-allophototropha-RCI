# Comparative genomics
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020

**NOTE: for each code section provided below, the code ought to be run from within this `comparative_genomics` directory.**

## Genome datatset download
Download all genomes used for the comparative genomics analyses below.

### Downloaded NCBI genomes

Downloaded the genomes using the simple bash script `download_genome_collection.sh` (could take time):
```bash
./download_genome_collection.sh "genome_collection_info.tsv" "downloads/genome_collection" \
  | tee "downloads/genome_collection/download.log"
```

A few genomes of interest have no ORF predictions on NCBI. Thus, I downloaded the nucleotide files only and then predicted ORFs:
```bash
output_dir="downloads/genome_collection/no_faa"

# Run the downloader
./download_genome_collection_no_faa.sh "genome_collection_info_no_faa.tsv" "${output_dir}" \
  | tee "${output_dir}/download.log"

# Predict ORFs
# **Requires prodigal to be installed - v2.6.3 was used
gunzip "${output_dir}/*.fna.gz"

fna_files=($(find "${output_dir}" -type f -name "*.fna" | sort -h))

for fna_file in ${fna_files[@]}; do
  basename="${fna_file%.fna}"
  echo "[ $(date -u) ]: ${basename##*/}"

  # v2.6.3
  # N.B., using default translation table
  prodigal -i "${fna_file}" -a "${basename}.faa" -d "${basename}.ffn" -f gff \
      -o "${basename}.gff" > "${basename}_prodigal.log" 2>&1

done

gzip "${output_dir}/*.fna" "${output_dir}/*.faa" "${output_dir}/*.ffn" "${output_dir}/*.gff"
```

### Summarized genomes of interest not available on NCBI
Got Ca. Chloroheliales bins

```bash
# TODO
```

The other genomes were collected manually from various databases with the advice of Marcus and Vera Thiel (I think mostly Vera):
- Chloranaerofilum_corporosum_YNP-MS-B-OTU-15
- Roseilinea_gracile_YNP-MS-B-OTU-6
- Chlorothrix_halophila
**Note these three are missing GFF files.**

### Made full summary directory with ORF predictions only
Then combined into the same directory:
```bash
work_dir="downloads/genome_collection/all_faa"
cd "${work_dir}"

# Made subfolders for organization

# Made hard links of genome files from both directories into the subfolders
find ../02_ncbi_downloads -type f -name "*.faa.gz" | xargs -I {} ln {} faa
find ../02b_non_ncbi_genomes -type f -name "*.faa.gz" | xargs -I {} ln {} faa
```

### Made subdirectories based on usage of the genomes
```bash
work_dir="downloads/genomes_by_category"
mkdir -p "${work_dir}"

# Put guide file into this folder: `genome_subsets_info.tsv`
# Put simple linking script into this folder: `link_genomes_vs2.sh`. Expects this folder to be $PWD when run
./link_genomes.sh 2>&1 | tee link_genomes.log
```

## _Chloroflexota_ concatenated core protein phylogeny

```bash
work_dir="Chloroflexota_phylogeny"
output_dir="Chloroflexota_phylogeny/Cfx_taxonomy"
source_dir="downloads/genomes_by_category/Cfx_taxonomy"

mkdir -p "${output_dir}"
cd "${output_dir}"

# Summarize input FAA files
find "${source_dir}" -type f -name "*.faa.gz" > input_faa_filepaths.list

# Made concatenated core protein tree
# via gtotree 1.4.11
conda activate gtotree_1.4.11

# TODO - two omitted?
# TODO - check command
GToTree -A input_faa_filepaths.list -H Bacteria.hmm -o bacteria -T IQ-TREE -n 12 -j 12 -G 0.3
# Used model LG+F+R6
```

A copy of this phylogeny is already provided in the folder: `Chloroflexota_phylogeny/summary/Chloroflexota_core_gene_phylogeny.treefile`

## Bidirectional BLASTP among _Chloroflexota_ members
Using the same genomes downloaded for the section above.

Query files for bidirectional BLASTP are already provided in `Chloroflexota_bidirectional_BLASTP/queries`

Install BackBLAST version 2.0.0-alpha3:
```bash

```

Then run BackBLAST for each of the four query genomes against the entire reference set
```bash

```

The combine the output files and make the rough gene heatmap (with phylogenetic tree made above attached)
```bash

```
The output raw PDF file is then used to generate Figure 3 in the manuscript -- see more details in the Figure 3 folder in this repo.

## Alignments and phylogenies of photosynthesis-associated genes


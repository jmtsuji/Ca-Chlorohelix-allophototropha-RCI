# Sequencing data availability
Part of the larger '*Ca.* Chlorohelix allophototropha' Type I reaction center paper

Copyright Jackson M. Tsuji, Neufeld Research Group, 2024

## 16S rRNA gene amplicon data
### V4-V5 region amplicon data for early enrichment cultures
You can download the data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="culture_early_16S_rRNA_gene_data_accessions_IlluminaV4V5p1.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
r1_urls=($(cat "${accession_data_filepath}" | cut -f 6 | tail -n +2))
r2_urls=($(cat "${accession_data_filepath}" | cut -f 7 | tail -n +2))

# Iteratively download the FastQ files
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

### V4-V5 region amplicon data for later enrichment cultures
These data can be downloaded just like in the code block above, but 
`culture_early_16S_rRNA_gene_data_accessions_IlluminaV4V5p1.tsv` should be replaced with
`culture_early_16S_rRNA_gene_data_accessions_IlluminaV4V5p2.tsv`.

### V4 region amplicon data for later enrichment cultures

Available in BioProject PRJNA640240 on NCBI. To download, run the following Bash code in your working folder:

```bash
accession_data_filepath="culture_early_16S_rRNA_gene_data_accessions_IlluminaV4.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
r1_urls=($(cat "${accession_data_filepath}" | cut -f 6 | tail -n +2))
r2_urls=($(cat "${accession_data_filepath}" | cut -f 7 | tail -n +2))

# Iteratively download the FastQ files
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
As one limitation of these data, note that the reads will lack quality scores due to data compression measures being 
taken by NCBI. You can go to the NCBI website directly for this sample to get the version that still has quality scores 
if you would like.

### Long-read (Nanopore) amplicon data for later enrichment cultures, V1-V9 region
Available in BioProject PRJNA640240 on NCBI. To download, run the following Bash code in your working folder:

Install the NCBI SRA Toolkit - pre-compiled binaries are available at [this URL](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit#sra-toolkit)

```bash
accession_data_filepath="culture_early_16S_rRNA_gene_data_accessions_nanopore.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
run_ids=($(cat "${accession_data_filepath}" | cut -f 5 | tail -n +2))

# Iteratively download the FastQ files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  run_id=${run_ids[${j}]}

  # Download - note that fasterq-dump is part of the SRA Toolkit
  echo "[$(date -u)]: Downloading '${sample_id}'"
  fasterq-dump -o "${output_dir}/${sample_id}.fastq" "${run_id}" # makes a FastQ version with suffix .fastq
  gzip "${output_dir}/${sample_id}.fastq"
  
done

echo "[$(date -u)]: Finished."
```

## Metagenome and genome bin data for early enrichment cultures
### Short read metagenomes for early enrichment cultures
You can download the data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="culture_early_metagenome_data_accessions.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
r1_urls=($(cat "${accession_data_filepath}" | cut -f 6 | tail -n +2))
r2_urls=($(cat "${accession_data_filepath}" | cut -f 7 | tail -n +2))

# Iteratively download the FastQ files
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

You now have FastQ files for the three datasets. Note that these files might have their quality scores masked out 
as part of the NCBI/ENA's effort to reduce file size.

### Genome bin data from early enrichment cultures
#### Curated '_Ca_. Chloroheliales' genome bins
You can download the genome files using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="culture_early_genome_bin_data_accessions.tsv"
output_dir="downloads/Ca_Chloroheliales"
search_key="Ca_Chloroheliales_bin"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 2))
fna_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 7))
faa_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 8))
gff_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 9))

# Iteratively download the genome files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  fna_url=${fna_urls[${j}]}
  faa_url=${faa_urls[${j}]}
  gff_url=${gff_urls[${j}]}

  # Download
  echo "[$(date -u)]: Downloading '${sample_id}'"
  wget -O "${output_dir}/${sample_id}.fna.gz" "${fna_url}"
  wget -O "${output_dir}/${sample_id}.faa.gz" "${faa_url}"
  wget -O "${output_dir}/${sample_id}.gff.gz" "${gff_url}"
done

echo "[$(date -u)]: Finished."
```

#### Other uncurated bins recovered from the early enrichment culture metagenomes
```bash
accession_data_filepath="culture_early_genome_bin_data_accessions.tsv"
output_dir="downloads/other_recovered_bins"
search_key="other_recovered_bin"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 2))
fna_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 7))
faa_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 8))
gff_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 9))

# Iteratively download the genome files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  fna_url=${fna_urls[${j}]}
  faa_url=${faa_urls[${j}]}
  gff_url=${gff_urls[${j}]}

  # Download
  echo "[$(date -u)]: Downloading '${sample_id}'"
  wget -O "${output_dir}/${sample_id}.fna.gz" "${fna_url}"
  wget -O "${output_dir}/${sample_id}.faa.gz" "${faa_url}"
  wget -O "${output_dir}/${sample_id}.gff.gz" "${gff_url}"
done

echo "[$(date -u)]: Finished."
```

### Supplementary genome bin files for early enrichment cultures (if of interest)

#### Original prokka annotations of the curated '_Ca_. Chloroheliales' genome bins
These annotations were actually used when analyzing data for the manuscript.

Download the ORFs from Zenodo:
```bash
output_dir="downloads/extras/prokka"
zenodo_url_base="https://zenodo.org/record/3930588/files" # v1.0.1

mkdir -p "${output_dir}"

wget -O "${output_dir}/Ca_Chx_allophototropha_L227-S17_prokka_ORFs.faa.gz" "${zenodo_url_base}/Ca_Chx_allophototropha_L227-S17_prokka_ORFs.faa.gz"
wget -O "${output_dir}/Ca_Chloroheliaceae_bin_L227_5C_prokka_ORFs.faa.gz" "${zenodo_url_base}/Ca_Chloroheliaceae_bin_L227_5C_prokka_ORFs.faa.gz"
```

#### Non-curated '_Ca_. Chloroheliales' genome bins
In case you want to see what was lost during curation.

Download the nucleotide files from Zenodo:
```bash
output_dir="downloads/extras/prokka"
zenodo_url_base="https://zenodo.org/record/3930588/files" # v1.0.1

mkdir -p "${output_dir}"

wget -O "${output_dir}/Ca_Chx_allophototropha_L227-S17_uncurated.fna.gz" "${zenodo_url_base}/Ca_Chx_allophototropha_L227-S17_uncurated.fna.gz"
wget -O "${output_dir}/Ca_Chloroheliaceae_bin_L227_5C_uncurated.fna.gz" "${zenodo_url_base}/Ca_Chloroheliaceae_bin_L227_5C_uncurated.fna.gz"
```

## Genomics data for closed genomes of "_Ca_. Chx. allophototropha" L227-S17 and _Geothrix_ sp. L227-G1
### Raw reads for obtaining a closed genome of "_Ca_. Chx. allophototropha" L227-S17
Available in BioProject PRJNA909349 on NCBI. To download, run the following Bash code in your working folder:

Install the NCBI SRA Toolkit - pre-compiled binaries are available at [this URL](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit#sra-toolkit)

```bash
accession_data_filepath="Chx_closed_genome_raw_data_accessions.tsv"
output_dir="downloads/Chx_closed_genome"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
run_ids=($(cat "${accession_data_filepath}" | cut -f 5 | tail -n +2))

# Iteratively download the FastQ files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  run_id=${run_ids[${j}]}

  # Download - note that fasterq-dump is part of the SRA Toolkit
  echo "[$(date -u)]: Downloading '${sample_id}'"
  fasterq-dump -o "${output_dir}/${sample_id}.fastq" "${run_id}" # makes a FastQ version with suffix .fastq
  gzip "${output_dir}/${sample_id}.fastq"
  
done

echo "[$(date -u)]: Finished."
```

### Raw reads for obtaining a closed genome of _Geothrix_ sp. L227-G1
Available in BioProject PRJNA975665 on NCBI. To download, run the following Bash code in your working folder:

Install the NCBI SRA Toolkit - pre-compiled binaries are available at [this URL](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit#sra-toolkit)

```bash
accession_data_filepath="Geothrix_closed_genome_bin_raw_data_accessions.tsv"
output_dir="downloads/Geothrix_closed_genome"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
run_ids=($(cat "${accession_data_filepath}" | cut -f 5 | tail -n +2))

# Iteratively download the FastQ files
for i in $(seq 1 ${#sample_ids[@]}); do
  # Make a zero-ordered counter
  j=$((${i}-1))

  # Get variables
  sample_id=${sample_ids[${j}]}
  run_id=${run_ids[${j}]}

  # Download - note that fasterq-dump is part of the SRA Toolkit
  echo "[$(date -u)]: Downloading '${sample_id}'"
  fasterq-dump -o "${output_dir}/${sample_id}.fastq" "${run_id}" # makes a FastQ version with suffix .fastq
  gzip "${output_dir}/${sample_id}.fastq"
  
done

echo "[$(date -u)]: Finished."
```

### Closed genomes for both strains
These can be obtained on NCBI via GenBank accessions `GCA_030389965.1` (for L22-S17) 
and `GCA_030219325.1` (for L227-G1).

## Lake survey metagenomes and metatranscriptomes
You can download the metagenome data using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="lake_metagenome_data_accessions_ncbi.tsv"
output_dir="downloads"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | cut -f 2 | tail -n +2))
r1_urls=($(cat "${accession_data_filepath}" | cut -f 6 | tail -n +2))
r2_urls=($(cat "${accession_data_filepath}" | cut -f 7 | tail -n +2))

# Iteratively download the FastQ files
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

In addition, some of the metagenomes are stored in the JGI database. 
See `lake_metagenome_data_accessions_jgi.tsv` for the IMG Genome IDs of each of these metagenomes; you can manually 
download the metagenome reads from JGI using these IDs if you would like. 
In my case, when I downloaded the files, they were in interleaved FastQ format, and I had to de-interleave the FastQ 
files before starting ATLAS. 
TODO - there is also a way to download the raw files from NCBI now - post instructions on how to do this.

To download the metatranscriptome data from NCBI, repeat the above code block using the 
`lake_metatranscriptome_data_accessions.tsv` file in place of `lake_metagenome_data_accessions_ncbi.tsv`.

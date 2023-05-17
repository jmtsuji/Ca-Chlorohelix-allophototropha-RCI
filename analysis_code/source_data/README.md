# Sequencing data availability
Part of the larger '*Ca.* Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

## 16S rRNA gene amplicon data for early enrichment cultures
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


## Short read metagenomes for early enrichment cultures
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

## Genome bin data from early enrichment cultures
### Curated '_Ca_. Chloroheliales' genome bins
You can download the genome files using the provided TSV file by running the following code in Bash in your working folder:

```bash
accession_data_filepath="genome_bin_data_accessions.tsv"
output_dir="downloads/Ca_Chloroheliales"
search_key="Ca_Chloroheliales_bin"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 2))
fna_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 7))
faa_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 8))
gff_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 9))

# Iternatively download the genome files
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

### Other uncurated bins recovered from the metagenomes
```bash
accession_data_filepath="genome_bin_data_accessions.tsv"
output_dir="downloads/other_recovered_bins"
search_key="other_recovered_bin"

mkdir -p "${output_dir}"

# Read selected table columns
sample_ids=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 2))
fna_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 7))
faa_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 8))
gff_urls=($(cat "${accession_data_filepath}" | grep "${search_key}$" | cut -f 9))

# Iternatively download the genome files
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

### Additional files (if of interest)

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


## Lake survey metagenomes and metatranscriptomes
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
In my case, when I downloaded the files, they were in interleaved FastQ format, and I had to de-interleave the FastQ files before starting ATLAS.

To download the metatranscriptome data from NCBI, repeat the above code block using the `metatranscriptome_data_accessions.tsv` file 
in place of `metagenome_data_accessions_ncbi.tsv`.


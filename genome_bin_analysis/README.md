# Genome bin analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020

**NOTE: for each code section provided below, the code ought to be run from within this `genome_bin_analysis` directory.**

In this methods segment, the '_Ca_. Chloroheliales' genome bins that were generated as described in the `metagenome_analysis` folder are analyzed. 
In addition, download instructions are provided for the other genome bins recovered in this study.

## Data download
### Curated '_Ca_. Chloroheliales' genome bins
You can download the genome files using the provided TSV file by running the following code in Bash in your working folder:

**TODO**: Finish the `genome_bin_data_accessions.tsv` file once NCBI URLs are available.

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

## Key photosynthesis-associated genes in '_Ca_. Chloroheliales' genomes
FastA files (non-aligned) are provided with the sequences of key photosynthes-associated genes and their predicted ORFs for the two genome bins. 
The genes summarized here are:
- _pscA_ encoding one of the core Type I reaction center components
- _fmoA_ encoding a the FMO energy-transfer protein between Type I reaction centers and chlorosomes
- _csmA_ encoding a key structural protein for chlorosomes
- _cbbL_ encoding the RuBisCO large subunit gene

FastA files can be found in the `sequences_of_interest` directory. `.ffn` files contain gene nucletide sequences. 
`.faa` files contain predicted amino acid sequences.


## Homology modelling
Homology modelling was performed using I-TASSER for predicted amino acid sequences of the key genes above, with default settings.

The top-scoring predicted tertiary structure for each gene (PDB format) are available in the `homology_models` directory. 
These PDB files can be opened in a protein structural viewer such as [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/).

The full I-TASSER output for each gene can be found at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110).


## Custom hidden Markov models
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

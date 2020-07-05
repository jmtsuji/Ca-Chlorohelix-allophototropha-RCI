#!/usr/bin/env bash
set -euo pipefail
# Download genomes from the NCBI database - FNA files only
# Copyright Jackson M. Tsuji, 2020

# Get user variables
genome_info_filepath=$1 # "genome_collection_info_no_faa.tsv"
output_dir=$2 # "downloads/genome_collection"

# Get info from the genome info file
genome_basenames=($(cat "${genome_info_filepath}" | cut -d $'\t' -f 3 | tail -n +2))
ftp_dir_urls=($(cat "${genome_info_filepath}" | cut -d $'\t' -f 5 | tail -n +2))

for i in $(seq 1 ${#genome_basenames[@]}); do
  # Convert from 1 to 0 ordered
  j=$((${i}-1))

  # Get variable names
  genome_basename=${genome_basenames[${j}]}
  ftp_dir_url=${ftp_dir_urls[${j}]}

  ## Figure out download links
  # If directory url is: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1
  # Then looking for the following files:
  # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1/GCA_001898225.1_ASM189822v1_genomic.fna.gz - genome nucleotide
  # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1/GCA_001898225.1_ASM189822v1_protein.faa.gz - genome predicted protein
  # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1/GCA_001898225.1_ASM189822v1_genomic.gff.gz - genome flat file
  # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/898/225/GCA_001898225.1_ASM189822v1/GCA_001898225.1_ASM189822v1_rna_from_genomic.fna.gz - RNA

  ftp_basename=${ftp_dir_url##*/}
  ftp_url_no_extension="${ftp_dir_url}/${ftp_basename}"
  ftp_url_fna="${ftp_url_no_extension}_genomic.fna.gz"
  ftp_url_faa="${ftp_url_no_extension}_protein.faa.gz"
  ftp_url_gff="${ftp_url_no_extension}_genomic.gff.gz"
  #ftp_url_rna="${ftp_url_no_extension}_rna_from_genomic.fna.gz"

  echo "[ $(date -u) ]: Downloading '${ftp_basename}' as '${genome_basename}'"
  wget -nv -O - "${ftp_url_fna}" > "${genome_basename}.fna.gz"
  #wget -nv -O - "${ftp_url_faa}" > "${genome_basename}.faa.gz"
  #wget -nv -O - "${ftp_url_gff}" > "${genome_basename}.gff.gz"
  #wget -nv -O - "${ftp_url_rna}" > "${genome_basename}.rna.fna.gz" # Often isn't there

done

echo "[ $(date -u) ]: Finished."

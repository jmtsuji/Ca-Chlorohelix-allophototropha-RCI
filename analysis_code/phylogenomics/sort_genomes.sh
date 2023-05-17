#!/usr/bin/env bash
set -euo pipefail

# Variables
source_dir=$1
output_dir=$2
guide_filepath=$3

echo "[ $(date -u) ]: Running ORF file organizer"

mkdir -p "${output_dir}"
cd "${output_dir}"

mkdir -p Cfx_taxonomy RCI_PSI FmoA CsmA CbbL chlorophyll

# Guide file columns (3,13-18):
# genome_name	Cfx_taxonomy	RCI_PSI	FmoA	CsmA	CbbL	chlorophyll
# 'y' means to put the genome in that folder

# Check columns
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 3)" != "genome_name" ]; then
  echo "ERROR: column 3 does not equal 'genome_name'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 13)" != "Cfx_taxonomy" ]; then
  echo "ERROR: column 13 does not equal 'Cfx_taxonomy'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 14)" != "RCI_PSI" ]; then
  echo "ERROR: column 14 does not equal 'RCI_PSI'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 15)" != "FmoA" ]; then
  echo "ERROR: column 15 does not equal 'FmoA'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 16)" != "CsmA" ]; then
  echo "ERROR: column 16 does not equal 'CsmA'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 17)" != "CbbL" ]; then
  echo "ERROR: column 17 does not equal 'CbbL'. Exiting..."
  exit 1
fi
if [ "$(cat "${guide_filepath}" | head -n 1 | cut -f 18)" != "chlorophyll" ]; then
  echo "ERROR: column 18 does not equal 'chlorophyll'. Exiting..."
  exit 1
fi

# Get column contents
genome_names=($(cat "${guide_filepath}" | cut -f 3 | tail -n +2))
Cfx_taxonomys=($(cat "${guide_filepath}" | cut -f 13 | tail -n +2))
photosystem_Is=($(cat "${guide_filepath}" | cut -f 14 | tail -n +2))
fmos=($(cat "${guide_filepath}" | cut -f 15 | tail -n +2))
chlorosomes=($(cat "${guide_filepath}" | cut -f 16 | tail -n +2))
rubiscos=($(cat "${guide_filepath}" | cut -f 17 | tail -n +2))
chlorophylls=($(cat "${guide_filepath}" | cut -f 18 | tail -n +2))

for i in $(seq 1 ${#genome_names[@]}); do
  # Convert from 1 ordered to 0 ordered counter
  j=$((${i}-1))

  # Get variable names
  genome_name=${genome_names[${j}]}
  Cfx_taxonomy=${Cfx_taxonomys[${j}]}
  photosystem_I=${photosystem_Is[${j}]}
  fmo=${fmos[${j}]}
  chlorosome=${chlorosomes[${j}]}
  rubisco=${rubiscos[${j}]}
  chlorophyll=${chlorophylls[${j}]}

  printf "[ $(date -u) ]: ${genome_name}: "
  printf "pattern is '${Cfx_taxonomy}-${photosystem_I}-${fmo}-${chlorosome}-${rubisco}-${chlorophyll}': "

  copy_counter=0
  
  if [ "${Cfx_taxonomy}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "Cfx_taxonomy; "
    ln "${source_dir}/${genome_name}.faa.gz" Cfx_taxonomy
  fi
  
  if [ "${photosystem_I}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "photosystem_I; "
    ln "${source_dir}/${genome_name}.faa.gz" RCI_PSI
  fi

  if [ "${fmo}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "fmo; "
    ln "${source_dir}/${genome_name}.faa.gz" FmoA
  fi

  if [ "${chlorosome}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "chlorosome; "
    ln "${source_dir}/${genome_name}.faa.gz" CsmA
  fi
  
  if [ "${rubisco}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "rubisco; "
    ln "${source_dir}/${genome_name}.faa.gz" CbbL
  fi

  if [ "${chlorophyll}" = "y" ]; then
    copy_counter=$((${copy_counter}+1))
    printf "chlorophyll; "
    ln "${source_dir}/${genome_name}.faa.gz" chlorophyll
  fi

  if [ ${copy_counter} -eq 0 ]; then
    printf "No need to link.\n"
  elif [ ${copy_counter} -gt 0 ]; then
    printf "Linked ${copy_counter} times\n"
  else
    echo "[ $(date -u) ]: ERROR with copy_counter. Exiting..."
    exit 1
  fi

done
echo "[ $(date -u) ]: Done."

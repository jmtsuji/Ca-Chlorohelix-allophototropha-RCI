# Comparative genomics
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

**NOTE: for each code section provided below, the code ought to be run from within this `phylogenomics` directory.**

## Genome datatset download
Download all genomes used for the comparative genomics analyses below.

### Downloaded NCBI genomes

Downloaded the genomes using the simple bash script `download_genome_collection.sh` (could take time):
```bash
output_dir="downloads/genome_collection/with_faa"
mkdir -p "${output_dir}"

./download_genome_collection.sh "genome_collection_info.tsv" "${output_dir}" "ncbi" "all" \
  | tee "downloads/genome_collection/download.log"

gunzip "downloads/genome_collection/*.gz"
```

A few genomes of interest have no ORF predictions on NCBI. Thus, I downloaded the nucleotide files only and then predicted ORFs:
```bash
output_dir="downloads/genome_collection/no_faa"
mkdir -p "${output_dir}"

# Run the downloader
./download_genome_collection.sh "genome_collection_info_no_faa.tsv" "${output_dir}" "ncbi_no_faa" "only_fna" \
  | tee "${output_dir}/download.log"

# Predict ORFs using prodigal v2.6.3
# Install via: conda create -n prodigal_2.6.3 -c bioconda prodigal=2.6.3
# conda activate prodigal_2.6.3
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
```

### Summarized genomes of interest not available on NCBI
A few other genomes were also collected manually from various databases:
- `Chloranaerofilum_corporosum_YNP-MS-B-OTU-15` -- from RAST
- `Roseilinea_gracile_YNP-MS-B-OTU-6` -- from RAST
- `Chlorothrix_halophila` -- from the JGI GOLD database
- `Vulcanimicrobium_alpinus_WC8-2` -- from NCBI Genbank (no assemby entry): see https://www.ncbi.nlm.nih.gov/nuccore/AP025523.1
  - This was later made available under the NCBI assembly accession GCA_027923555.1 (so it can be downloaded like standard NCBI entries in future)

For _Chlorothrix halophila_, I had to export the contigs only and then predict the genes manually, using prokka 1.13.3 with default 
settings with respect to gene prediction.

Save unzipped files in `downloads/genome_collection/manual_downloads` with the extension `.faa` for the code below to work.

### Summarized Ca. Chloroheliales genomes
These genomes were prepared as part of this specific study. Get the predicted amino acid files:

- `Ca_Chx_allophototropha_L227-S17`: using the PGAP annotated version of the complete genome

- `Ca_Chloroheliaceae_bin_L227_5C`: using the PGAP annotated version of the manually curated genome bin
  - Note: this genome can now also be downloaded from NCBI assembly accession `GCA_013390945.1`

TODO: in future, point to how to download the L227-S17 genome from here.

As above, save unzipped files in `downloads/genome_collection/manual_downloads` with the extension `.faa` for the code below to work.

### Made full summary directory with ORF predictions only
Then combined into the same directory:
```bash
input_dir="downloads/genome_collection"
output_dir="downloads/genome_collection_faa_summary"

mkdir -p "${output_dir}"

# Made hard links of genome files from both directories into the subfolders
find "${input_dir}" -type f -name "*.faa" | xargs -I {} ln {} "${output_dir}"
```

### Made subdirectories based on usage of the genomes
```bash
input_dir="downloads/genome_collection_faa_summary"
output_dir="downloads/genomes_by_category"
guide_filepath="genome_collection_info.tsv"

mkdir -p "${output_dir}"

# Run code to sort genomes into relevant folders (via hard linking)
./sort_genomes.sh "${input_dir}" "${output_dir}" "${guide_filepath}" 2>&1 | tee sort_genomes.log
```

## _Chloroflexota_ concatenated core protein phylogeny
```bash
work_dir="Chloroflexota_phylogeny"
output_dir="${work_dir}/Cfx_taxonomy"
source_dir="downloads/genomes_by_category/Cfx_taxonomy"
threads=12

mkdir -p "${work_dir}"
cd "${work_dir}"

# Summarize input FAA files
find "${source_dir}" -type f -name "*.faa" > "${work_dir}/input_faa_filepaths.list"

# Made concatenated core protein tree via gtotree 1.4.11 
# First, make the conda env: conda create -n gtotree_1.4.11 -c conda-forge -c bioconda -c astrobiomike gtotree=1.4.11
# Then activate: conda activate gtotree_1.4.11

GToTree -A "${work_dir}/input_faa_filepaths.list" -H Bacteria.hmm -o "${output_dir}" \
  -T IQ-TREE -n "${threads}" -j "${threads}" -G 0.3
# Note that this was originally run with gzipped .faa files as input, but for simplicity in the overall workflow here, I show the command run with
#   non-gzipped .faa files. This should not affect the run.
```

A copy of this phylogeny is already provided in the folder: `Chloroflexota_phylogeny/Chloroflexota_core_gene_phylogeny.treefile`

## Bidirectional BLASTP among _Chloroflexota_ members
Using the same genomes downloaded for the section above.

Query files for bidirectional BLASTP are already provided in `Chloroflexota_bidirectional_BLASTP/queries`. 
In addition, the metadata for plotting are provided as `gene_metadata.tsv` and `genome_metadata.tsv`.

### Install BackBLAST version 2.0.0-alpha3:
```bash
# Get the repo
wget -O - https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST/archive/v2.0.0-alpha3.tar.gz | tar -xzf -
cd BackBLAST_Reciprocal_BLAST-2.0.0-alpha3

# Install requirements
conda env create -n backblast_2.0.0-alpha3 --file=envs/conda_requirements.yaml
conda activate backblast_2.0.0-alpha3

# Add to conda env PATH
mkdir -p ${CONDA_PREFIX}/share/backblast
cp -r * ${CONDA_PREFIX}/share/backblast
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
if [[ ! -f ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh ]]; then
  echo '#!/bin/sh' > ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
fi
echo "export PATH=\${PATH}:${CONDA_PREFIX}/share/backblast" \
  >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
chmod 755 ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# Re-activate the repo to apply the changes
conda deactivate
conda activate backblast_2.0.0-alpha3

# Delete the downloaded repo
cd ..
rm -rf BackBLAST_Reciprocal_BLAST-2.0.0-alpha3
```

### Run BackBLAST
Then run BackBLAST for each of the four query genomes against the entire reference set
```bash
work_dir="Chloroflexota_bidirectional_BLASTP"
run_dir="${work_dir}/backblast"
subject_dir="downloads/genomes_by_category/Cfx_taxonomy"
query_dir="${work_dir}/queries" # already provided in this repo
tree_filepath="Chloroflexota_phylogeny/Chloroflexota_core_gene_phylogeny.treefile"
jobs=50
evalue=1e-3
pident=20
qcov=50

# FYI, Make sure you are in the conda env. To activate:
# conda activate backblast_2.0.0-alpha3

# Cfl. aurantiacus
output_dir="${run_dir}/01a_Cfl_aurantiacus"
genome_prefix="Chloroflexus_aurantiacus_J-10-fl"
mkdir -p "${output_dir}"
backblast auto -j ${jobs} -e ${evalue} -p ${pident} -c ${qcov} \
  ${query_dir}/${genome_prefix}_gene_targets.faa \
  ${subject_dir}/${genome_prefix}.faa ${subject_dir} \
  ${output_dir} --until finished_blast 2>&1 | \
  tee ${output_dir}/backblast_vs2a.log

# Osc. trichoides
output_dir="01b_Osc_trichoides"
genome_prefix="Oscillochloris_trichoides_DG-6"
mkdir -p "${output_dir}"
backblast auto -j ${jobs} -e ${evalue} -p ${pident} -c ${qcov} \
  ${query_dir}/${genome_prefix}_gene_targets.faa \
  ${subject_dir}/${genome_prefix}.faa ${subject_dir} \
  ${output_dir} --until finished_blast 2>&1 | \
  tee ${output_dir}/backblast_vs2b.log

# Ca. Chx. allophototropha
output_dir="01c_Chx_allophototropha"
genome_prefix="Ca_Chx_allophototropha_L227-S17"
mkdir -p "${output_dir}"
backblast auto -j ${jobs} -e ${evalue} -p ${pident} -c ${qcov} \
  ${query_dir}/${genome_prefix}_gene_targets.faa \
  ${subject_dir}/${genome_prefix}.faa ${subject_dir} \
  ${output_dir} --until finished_blast 2>&1 | \
  tee ${output_dir}/backblast_vs2c.log

# Roseiflexus_castenholzii_DSM_13941
output_dir="01d_Rfl_castenholzii"
genome_prefix="Roseiflexus_castenholzii_DSM_13941"
mkdir -p ${output_dir}
backblast auto -j ${jobs} -e ${evalue} -p ${pident} -c ${qcov} \
  ${query_dir}/${genome_prefix}_gene_targets.faa \
  ${subject_dir}/${genome_prefix}.faa ${subject_dir} \
  ${output_dir} --until finished_blast 2>&1 | \
  tee ${output_dir}/backblast_vs2d.log

### Combine blast tables
output_dir="${run_dir}/02_combined"
mkdir -p "${output_dir}"
cp ${run_dir}/01a_Cfl_aurantiacus/blast/blast_tables_combined.csv ${output_dir}/blast_tables_combined.csv
tail -n +2 ${run_dir}/01b_Osc_trichoides/blast/blast_tables_combined.csv >> ${output_dir}/blast_tables_combined.csv
tail -n +2 ${run_dir}/01c_Chx_allophototropha/blast/blast_tables_combined.csv >> ${output_dir}/blast_tables_combined.csv
tail -n +2 ${run_dir}/01d_Rfl_castenholzii/blast/blast_tables_combined.csv >> ${output_dir}/blast_tables_combined.csv

### Make final viz
# Genome and gene metadata are already provided in the base folder
backblast generate_heatmap \
  -m genome_metadata.tsv \
  -g gene_metadata.tsv \
  -r "Methylobacterium_extorquens_AM1_outgroup" \
  -w 400 -z 250 \
  ${tree_filepath} \
  ${output_dir}/blast_tables_combined.csv \
  ${output_dir}/Figure_02_raw_Blues.pdf 2>&1 | \
  tee ${output_dir}/Figure_02_raw_Blues.log
```

### Custom plots
However, because I wanted each segment of the heatmap to use a different colour scale, I then re-generated different versions of the heatmap using
different colour scales to combine them manually later:
```bash
# First make sure to run: conda activate backblast_2.0.0-alpha3
heatmap_script="${CONDA_PREFIX}/share/backblast/generate_heatmap.R"

work_dir="Chloroflexota_bidirectional_BLASTP"
run_dir="${work_dir}/backblast"
output_dir="${run_dir}/02_combined"

# Make a copy of the heatmap script
cp "${heatmap_script}" "${heatmap_script}.backup"

palettes=(Purples Greens Oranges)
for pal in ${palettes[@]}; do
  # Change the colour palette option directly in the script
  sed -e "s/name = \"Blues\"/name = \"${pal}\"/" "${heatmap_script}.backup" > "${heatmap_script}"

  # Make the heatmap
  backblast generate_heatmap \
    -m Chloroflexota_bidirectional_BLASTP/genome_metadata.tsv \
    -g Chloroflexota_bidirectional_BLASTP/gene_metadata.tsv \
    -r "Methylobacterium_extorquens_AM1_outgroup" \
    -w 400 -z 250 \
    ${tree_filepath} \
    ${output_dir}/blast_tables_combined.csv \
    ${output_dir}/Figure_02_raw_${pal}.pdf 2>&1 | \
    tee ${output_dir}/Figure_02_raw_${pal}.log
done

# Restore the original script
mv "${heatmap_script}.backup" "${heatmap_script}"
```

The four raw plots (Blues, Purples, Greens, Oranges) are then manually combined to produce Fig. 2. 
The `blast_tables_combined.csv` file is modified to produce Supplementary Data 2.

## Alignments and phylogenies of photosynthesis-associated genes
In reality, I scanned the downloaded genomes above (in `genomes_by_category`) via BackBLAST to search for the genes of interest. 
I then manually checked the sequence alignments and prototype phylogenies to curate the list of BackBLAST hits for each gene 
to the final lists used below.

For simplicity, for this code, I give the final summary accessions that were decided upon for each gene and then download those directly 
from NCBI. However, you could cross-check these results with the downloaded genomes above. Note that a few sequences in each phylogeny 
do not have a corresponding NCBI full genome.

### Fe-S type reaction center (RCI/PSI)
A set of 68 RCI/PSI sequences were used for the phylogeny; accessions are summarized in `RCI_PSI_accessions.tsv`.

Download the sequences:
```bash
run_name="PSI_RCI"
cd "alignments_and_phylogenies/${run_name}"

# Summarize gene hits
filenames=($(cut -f 1 "${run_name}_accessions.tsv" | tail -n +2))
accessions=($(cut -f 2 "${run_name}_accessions.tsv" | tail -n +2))

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

output_filepath_1="${run_name}_raw.faa"
output_filepath_2="${run_name}_renamed.faa"
printf "" > "${output_filepath_1}"
printf "" > "${output_filepath_2}"

# You'll need to instal seqtk; I used version 1.3-r106
# conda create -n seqtk -c bioconda seqtk=1.3
# conda activate seqtk

for i in $(seq 1 ${#filenames[@]}); do
  j=$((${i}-1))
  file=${filenames[${j}]}
  acc=${accessions[${j}]}
  echo "[ $(date -u) ]: ${file}: ${acc}"
  
  # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
  echo "[ $(date -u) ]: ${file}: ${url}"
  
  wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
done

cd ../..
```

Align, mask, and make the PSI/RCI tree:
```bash
# User variables
threads=10
run_name="PSI_RCI"

cd "alignments_and_phylogenies/${run_name}"

# Align using Clustal Omega v.1.2.3
# Install via: conda create -n clustalo_1.2.3 -c bioconda clustalo=1.2.3
# conda activate clustalo_1.2.3
clustalo -i "${run_name}_renamed.faa" -o "${run_name}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${run_name}_aligned.log"

# Mask using Gblocks v0.91b
# Install via: conda create -n gblocks_0.91b -c bioconda gblocks=0.91b
# conda activate gblocks_0.91b
Gblocks "${run_name}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee "${run_name}_aligned_masked.log"
mv "${run_name}_aligned.faa_GB01" "${run_name}_aligned_masked.faa"

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
# Install via: conda create -n iqtree_1.6.11 -c bioconda iqtree=1.6.11
# conda activate iqtree_1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

# Supplement: try tree again with no mask for comparison
mkdir -p "phylogeny_unmasked"
iqtree -s "${run_name}_aligned.faa" -nt ${threads} -pre "phylogeny_unmasked/${run_name}_aligned" -bb 1000 -m TEST

cd ../..
```
The resuting tree file `PSI_RCI_aligned_masked.treefile` is plotted in Fig. 1b and Extended Data Fig. 2.

### FMO protein
A set of 30 FMO sequences were used for the phylogeny; accessions are summarized in `FmoA_accessions.tsv`.

Download the sequences:
```bash
run_name="FmoA"
cd "alignments_and_phylogenies/${run_name}"

# Summarize gene hits
filenames=($(cut -f 1 "${run_name}_accessions.tsv" | tail -n +2))
accessions=($(cut -f 2 "${run_name}_accessions.tsv" | tail -n +2))

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

output_filepath_1="${run_name}_raw.faa"
output_filepath_2="${run_name}_renamed.faa"
printf "" > "${output_filepath_1}"
printf "" > "${output_filepath_2}"

# You'll need to instal seqtk; I used version 1.3-r106
# conda create -n seqtk -c bioconda seqtk=1.3
# conda activate seqtk

for i in $(seq 1 ${#filenames[@]}); do
  j=$((${i}-1))
  file=${filenames[${j}]}
  acc=${accessions[${j}]}
  echo "[ $(date -u) ]: ${file}: ${acc}"
  
  # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
  echo "[ $(date -u) ]: ${file}: ${url}"
  
  wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
done

cd ../..
```

Align, mask, and make the FmoA phylogeny:
```bash
# User variables
threads=10
run_name="FmoA"

cd "alignments_and_phylogenies/${run_name}"

# Align using Clustal Omega v.1.2.3
# Install via: conda create -n clustalo_1.2.3 -c bioconda clustalo=1.2.3
# conda activate clustalo_1.2.3
clustalo -i "${run_name}_renamed.faa" -o "${run_name}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${run_name}_aligned.log"

# Mask using Gblocks v0.91b
# Install via: conda create -n gblocks_0.91b -c bioconda gblocks=0.91b
# conda activate gblocks_0.91b
Gblocks "${run_name}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee "${run_name}_aligned_masked.log"
mv "${run_name}_aligned.faa_GB01" "${run_name}_aligned_masked.faa"

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
# Install via: conda create -n iqtree_1.6.11 -c bioconda iqtree=1.6.11
# conda activate iqtree_1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

# Supplement: try tree again with no mask for comparison
mkdir -p "phylogeny_unmasked"
iqtree -s "${run_name}_aligned.faa" -nt ${threads} -pre "phylogeny_unmasked/${run_name}_aligned" -bb 1000 -m TEST

cd ../..
```
The resuting tree file `FmoA_aligned_masked.treefile` is plotted in Extended Data Fig. 4.

### Chlorosom protein (CsmA)
A set of 41 CsmA sequences were used for the phylogeny; accessions are summarized in `CsmA_accessions.tsv`.

Download the sequences:
```bash
run_name="CsmA"
cd "alignments_and_phylogenies/${run_name}"

# Summarize gene hits
filenames=($(cut -f 1 "${run_name}_accessions.tsv" | tail -n +2))
accessions=($(cut -f 2 "${run_name}_accessions.tsv" | tail -n +2))

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

output_filepath_1="${run_name}_raw.faa"
output_filepath_2="${run_name}_renamed.faa"
printf "" > "${output_filepath_1}"
printf "" > "${output_filepath_2}"

# You'll need to instal seqtk; I used version 1.3-r106
# conda create -n seqtk -c bioconda seqtk=1.3
# conda activate seqtk

for i in $(seq 1 ${#filenames[@]}); do
  j=$((${i}-1))
  file=${filenames[${j}]}
  acc=${accessions[${j}]}
  echo "[ $(date -u) ]: ${file}: ${acc}"
  
  # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
  echo "[ $(date -u) ]: ${file}: ${url}"
  
  wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
done

cd ../..
```

Align, mask, and make the CsmA phylogeny:
```bash
# User variables
threads=10
run_name="CsmA"

cd "alignments_and_phylogenies/${run_name}"

# Align using Clustal Omega v.1.2.3
# Install via: conda create -n clustalo_1.2.3 -c bioconda clustalo=1.2.3
# conda activate clustalo_1.2.3
clustalo -i "${run_name}_renamed.faa" -o "${run_name}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${run_name}_aligned.log"

# Mask using Gblocks v0.91b
# Install via: conda create -n gblocks_0.91b -c bioconda gblocks=0.91b
# conda activate gblocks_0.91b
Gblocks "${run_name}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee "${run_name}_aligned_masked.log"
mv "${run_name}_aligned.faa_GB01" "${run_name}_aligned_masked.faa"

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
# Install via: conda create -n iqtree_1.6.11 -c bioconda iqtree=1.6.11
# conda activate iqtree_1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

# Supplement: try tree again with no mask for comparison
mkdir -p "phylogeny_unmasked"
iqtree -s "${run_name}_aligned.faa" -nt ${threads} -pre "phylogeny_unmasked/${run_name}_aligned" -bb 1000 -m TEST

cd ../..
```
The resuting tree file `CsmA_aligned_masked.treefile` is plotted in Extended Data Fig. 8, and the multiple sequence alignment is shown 
in Extended Data Fig. 5.

I also summarized the phylogenetic distances between the taxa in this phylogeny using the provided R script `CsmA_distances.R`.

### RuBisCO large subunit (CbbL)
A set of 106 CbbL sequences were used for the phylogeny; accessions are summarized in `CbbL_accessions.tsv`.

Download the sequences:
```bash
run_name="CbbL"
cd "alignments_and_phylogenies/${run_name}"

# Summarize gene hits
filenames=($(cut -f 1 "${run_name}_accessions.tsv" | tail -n +2))
accessions=($(cut -f 2 "${run_name}_accessions.tsv" | tail -n +2))

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

output_filepath_1="${run_name}_raw.faa"
output_filepath_2="${run_name}_renamed.faa"
printf "" > "${output_filepath_1}"
printf "" > "${output_filepath_2}"

# You'll need to instal seqtk; I used version 1.3-r106
# conda create -n seqtk -c bioconda seqtk=1.3
# conda activate seqtk

for i in $(seq 1 ${#filenames[@]}); do
  j=$((${i}-1))
  file=${filenames[${j}]}
  acc=${accessions[${j}]}
  echo "[ $(date -u) ]: ${file}: ${acc}"
  
  # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
  echo "[ $(date -u) ]: ${file}: ${url}"
  
  wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
done

cd ../..
```

Align, mask, and make the CbbL phylogeny:
```bash
# User variables
threads=10
run_name="CbbL"

cd "alignments_and_phylogenies/${run_name}"

# Align using Clustal Omega v.1.2.3
clustalo -i "${run_name}_renamed.faa" -o "${run_name}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${run_name}_aligned.log"

# Mask using Gblocks v0.91b
Gblocks "${run_name}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee "${run_name}_aligned_masked.log"
mv "${run_name}_aligned.faa_GB01" "${run_name}_aligned_masked.faa"

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

# Supplement: try tree again with no mask for comparison
mkdir -p "phylogeny_unmasked"
iqtree -s "${run_name}_aligned.faa" -nt ${threads} -pre "phylogeny_unmasked/${run_name}_aligned" -bb 1000 -m TEST

cd ../..
```
The resulting tree file `CbbL_aligned_masked.treefile` is plotted in Extended Data Fig. 6. Note that I realized after producing this phylogeny that 
three pairs of sequences (between the reference genome set and the gene set derived from Tabita and colleagues, 2008) in the Group III CbbL set had 
been duplicated: the CbbL sequences of _Hyperthermus butylicus_, _Methanocaldococcus jannaschii_, and _Archaeoglobus fulgidus_. I collapsed those 
three branches into single nodes when producing the final figure for clarity.


### BchIDH concatenated primary sequence phylogeny
A set of 82 groups of three BchIDH sequences were used for the phylogeny; accessions are summarized in `BchIDH_accessions.tsv`.

Download the sequences:  
```bash
run_name="BchIDH"
gene_list=(BchI BchD BchH)

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

  # Make a file specific for those hits
  mkdir -p "${gene}"
  head -n 1 "${run_name}_accessions.tsv" > "${gene}/${gene}_accessions.tsv"
  grep "^${gene}" "${run_name}_accessions.tsv" >> "${gene}/${gene}_accessions.tsv"

  # Summarize gene hits
  filenames=($(cut -f 1 "${gene}/${gene}_accessions.tsv" | tail -n +2))
  accessions=($(cut -f 2 "${gene}/${gene}_accessions.tsv" | tail -n +2))

  output_filepath_1="${gene}/${gene}_raw.faa"
  output_filepath_2="${gene}/${gene}_renamed.faa"
  printf "" > "${output_filepath_1}"
  printf "" > "${output_filepath_2}"

  # You'll need to instal seqtk; I used version 1.3-r106
  # conda create -n seqtk -c bioconda seqtk=1.3
  # conda activate seqtk

  for i in $(seq 1 ${#filenames[@]}); do
    j=$((${i}-1))
    file=${filenames[${j}]}
    acc=${accessions[${j}]}
    echo "[ $(date -u) ]: ${file}: ${acc}"
  
    # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
    echo "[ $(date -u) ]: ${gene}: ${file}: ${url}"
  
    wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
  done
  
done

cd ../..
```

Align, mask, and make the phylogenies for each of BchI, BchD, and BchH:
```bash
# User variables
threads=10
run_name="BchIDH"
gene_list=(BchI BchD BchH)

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

# Tree generation below is not needed per se, but I generated them anyway and checked them manually against the final concatenated tree for reference.

echo "[ $(date -u) ]: Working on gene ${gene}"

  # Align using Clustal Omega v.1.2.3
  clustalo -i "${gene}/${gene}_renamed.faa" -o "${gene}/${gene}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${gene}/${gene}_aligned.log"

  # Mask using Gblocks v0.91b
  Gblocks "${gene}/${gene}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
    2>&1 | tee "${gene}/${gene}_aligned_masked.log"
  mv "${gene}/${gene}_aligned.faa_GB01" "${gene}/${gene}_aligned_masked.faa"

  # Made maximum likelihood phylogeny using IQ-TREE v1.6.11
  mkdir -p "${gene}/phylogeny"
  iqtree -s "${gene}/${gene}_aligned_masked.faa" -nt ${threads} -pre "${gene}/phylogeny/${gene}_aligned_masked" -bb 1000 -m TEST
  cp "${gene}/phylogeny/${gene}_aligned_masked.treefile" .

  # Supplement: try tree again with no mask for comparison
  mkdir -p "${gene}/phylogeny_unmasked"
  iqtree -s "${gene}/${gene}_aligned.faa" -nt ${threads} -pre "${gene}/phylogeny_unmasked/${gene}_aligned" -bb 1000 -m TEST

done

cd ../..
```

Concatenated the multiple sequence alignments into a single alignment via `concatenate_alignment_BchIDH.ipynb`, which was run in a Jupyter notebook. Ran this code to produce `${run_name}_aligned_masked.faa`.  

Then built a tree from this concatenated alignment:
```bash
run_name="BchIDH"
cd "alignments_and_phylogenies/${run_name}"
threads=26

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

cd ../..
```
The resuting tree file `BchIDH_aligned_masked.treefile` is plotted in Extended Data Fig. 9.

### BchLNB concatenated primary sequence phylogeny
A set of 90 groups of three BchLNB sequences were used for the phylogeny; accessions are summarized in `BchLNB_accessions.tsv`.

Download the sequences:  
```bash
run_name="BchLNB"
gene_list=(BchL BchN BchB)

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

  # Make a file specific for those hits
  mkdir -p "${gene}"
  head -n 1 "${run_name}_accessions.tsv" > "${gene}/${gene}_accessions.tsv"
  grep "^${gene}" "${run_name}_accessions.tsv" >> "${gene}/${gene}_accessions.tsv"

  # Summarize gene hits
  filenames=($(cut -f 1 "${gene}/${gene}_accessions.tsv" | tail -n +2))
  accessions=($(cut -f 2 "${gene}/${gene}_accessions.tsv" | tail -n +2))

  output_filepath_1="${gene}/${gene}_raw.faa"
  output_filepath_2="${gene}/${gene}_renamed.faa"
  printf "" > "${output_filepath_1}"
  printf "" > "${output_filepath_2}"

  # You'll need to instal seqtk; I used version 1.3-r106
  # conda create -n seqtk -c bioconda seqtk=1.3
  # conda activate seqtk

  for i in $(seq 1 ${#filenames[@]}); do
    j=$((${i}-1))
    file=${filenames[${j}]}
    acc=${accessions[${j}]}
    echo "[ $(date -u) ]: ${file}: ${acc}"
  
    # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
    echo "[ $(date -u) ]: ${gene}: ${file}: ${url}"
  
    wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
  done
  
done

cd ../..
```

Align, mask, and make the phylogenies for each of BchL, BchN, and BchB:
```bash
# User variables
threads=10
run_name="BchLNB"
gene_list=(BchL BchN BchB)

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

# Tree generation below is not needed per se, but I generated them anyway and checked them manually against the final concatenated tree for reference.

echo "[ $(date -u) ]: Working on gene ${gene}"

  # Align using Clustal Omega v.1.2.3
  clustalo -i "${gene}/${gene}_renamed.faa" -o "${gene}/${gene}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${gene}/${gene}_aligned.log"

  # Mask using Gblocks v0.91b
  Gblocks "${gene}/${gene}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
    2>&1 | tee "${gene}/${gene}_aligned_masked.log"
  mv "${gene}/${gene}_aligned.faa_GB01" "${gene}/${gene}_aligned_masked.faa"

  # Made maximum likelihood phylogeny using IQ-TREE v1.6.11
  mkdir -p "${gene}/phylogeny"
  iqtree -s "${gene}/${gene}_aligned_masked.faa" -nt ${threads} -pre "${gene}/phylogeny/${gene}_aligned_masked" -bb 1000 -m TEST
  cp "${gene}/phylogeny/${gene}_aligned_masked.treefile" .

  # Supplement: try tree again with no mask for comparison
  mkdir -p "${gene}/phylogeny_unmasked"
  iqtree -s "${gene}/${gene}_aligned.faa" -nt ${threads} -pre "${gene}/phylogeny_unmasked/${gene}_aligned" -bb 1000 -m TEST

done

cd ../..
```

Concatenated the multiple sequence alignments into a single alignment via `concatenate_alignment_BchLNB.ipynb`, which was run in a Jupyter notebook. Ran this code to produce `${run_name}_aligned_masked.faa`.  

Then built a tree from this concatenated alignment:
```bash
run_name="BchLNB"
cd "alignments_and_phylogenies/${run_name}"
threads=26

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

cd ../..
```
The resuting tree file `BchLNB_aligned_masked.treefile` is plotted in Extended Data Fig. 10.

### BchXYZ concatenated primary sequence phylogeny
A set of 67 groups of three BchXYZ sequences were used for the phylogeny; accessions are summarized in `BchXYZ_accessions.tsv`.

Download the sequences:  
```bash
run_name="BchXYZ"
gene_list=(BchX BchY BchZ)

# **MANUALLY enter your NCBI API key here to ensure your download is faster and less interrupted:
# API_KEY=[INSERT_HERE]

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

  # Make a file specific for those hits
  mkdir -p "${gene}"
  head -n 1 "${run_name}_accessions.tsv" > "${gene}/${gene}_accessions.tsv"
  grep "^${gene}" "${run_name}_accessions.tsv" >> "${gene}/${gene}_accessions.tsv"

  # Summarize gene hits
  filenames=($(cut -f 1 "${gene}/${gene}_accessions.tsv" | tail -n +2))
  accessions=($(cut -f 2 "${gene}/${gene}_accessions.tsv" | tail -n +2))

  output_filepath_1="${gene}/${gene}_raw.faa"
  output_filepath_2="${gene}/${gene}_renamed.faa"
  printf "" > "${output_filepath_1}"
  printf "" > "${output_filepath_2}"

  # You'll need to instal seqtk; I used version 1.3-r106
  # conda create -n seqtk -c bioconda seqtk=1.3
  # conda activate seqtk

  for i in $(seq 1 ${#filenames[@]}); do
    j=$((${i}-1))
    file=${filenames[${j}]}
    acc=${accessions[${j}]}
    echo "[ $(date -u) ]: ${file}: ${acc}"
  
    # If you don't want to use an API key, then just remove the '&api_key=${API_KEY}' part of the URL below.
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text&api_key=${API_KEY}"
    echo "[ $(date -u) ]: ${gene}: ${file}: ${url}"
  
    wget -nv -O - "${url}" | seqtk seq | tee -a "${output_filepath_1}" | sed -e "1s/^>.*$/>${file}__${acc}/g" >> "${output_filepath_2}"
  done
  
done

cd ../..
```

Align, mask, and make the phylogenies for each of BchX, BchY, and BchZ:
```bash
# User variables
threads=10
run_name="BchXYZ"
gene_list=(BchX BchY BchZ)

cd "alignments_and_phylogenies/${run_name}"

for gene in ${gene_list[@]}; do

# Tree generation below is not needed per se, but I generated them anyway and checked them manually against the final concatenated tree for reference.

echo "[ $(date -u) ]: Working on gene ${gene}"

  # Align using Clustal Omega v.1.2.3
  clustalo -i "${gene}/${gene}_renamed.faa" -o "${gene}/${gene}_aligned.faa" --threads=${threads} -v -v 2>&1 | tee "${gene}/${gene}_aligned.log"

  # Mask using Gblocks v0.91b
  Gblocks "${gene}/${gene}_aligned.faa" -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
    2>&1 | tee "${gene}/${gene}_aligned_masked.log"
  mv "${gene}/${gene}_aligned.faa_GB01" "${gene}/${gene}_aligned_masked.faa"

  # Made maximum likelihood phylogeny using IQ-TREE v1.6.11
  mkdir -p "${gene}/phylogeny"
  iqtree -s "${gene}/${gene}_aligned_masked.faa" -nt ${threads} -pre "${gene}/phylogeny/${gene}_aligned_masked" -bb 1000 -m TEST
  cp "${gene}/phylogeny/${gene}_aligned_masked.treefile" .

  # Supplement: try tree again with no mask for comparison
  mkdir -p "${gene}/phylogeny_unmasked"
  iqtree -s "${gene}/${gene}_aligned.faa" -nt ${threads} -pre "${gene}/phylogeny_unmasked/${gene}_aligned" -bb 1000 -m TEST

done

cd ../..
```

Concatenated the multiple sequence alignments into a single alignment via `concatenate_alignment_BchXYZ.ipynb`, which was run in a Jupyter notebook. Ran this code to produce `${run_name}_aligned_masked.faa`.  

Then built a tree from this concatenated alignment:
```bash
run_name="BchXYZ"
cd "alignments_and_phylogenies/${run_name}"
threads=26

# Made maximum likelihood phylogeny using IQ-TREE v1.6.11
mkdir -p "phylogeny"
iqtree -s "${run_name}_aligned_masked.faa" -nt ${threads} -pre "phylogeny/${run_name}_aligned_masked" -bb 1000 -m TEST
cp "phylogeny/${run_name}_aligned_masked.treefile" .

cd ../..
```
The resuting tree file `BchXYZ_aligned_masked.treefile` is plotted in Extended Data Fig. 10.

### (Bacterio)chlorophyll synthase (ChlG/BchG)
This analysis forms part of the supplementary materials and was done a bit differently/independently of the analyses above.

Here, genes were independently collected/downloaded and were then analyzed. See the unaligned input sequences in the `input` sub-folder.

#### Setup
Install needed software

```bash
conda create -n phylogenies -c bioconda -c conda-forge seqtk=1.3 clustalo=1.2.4 gblocks=0.91b iqtree=2.2.0.3 blast=2.12.0
```

#### Candidate genes from the L227-S17 and L227-5C genomes
Based on a BLASTP search of the Chlorobaculum tepidum ChlG sequence (AAM72500.1) against the strain L227-S17 genome:
```bash
conda activate phylogenies

# Save AAM72500.1 as AAM72500.1_ChlG.faa

# Assumes you have downloaded the closed genome of strain L227-S17 as described in the README in the source_data section. Make a copy here.
blastp -query AAM72500.1_ChlG.faa -subject L227-S17.faa -outfmt "6 qseqid sseqid pident evalue bitscore qcovhsp qlen qstart qend slen sstart send" -evalue 1e-1
```

Results
```
qseqid      sseqid                  pident  evalue   bits qcov qlen qst qend slen sst send notes
AAM72500.1  gnl|extdb|OZ401_002374  34.251  4.15e-57  186  87  367  41  359  333  15  329 # ChlG?
AAM72500.1  gnl|extdb|OZ401_001677  36.301  6.15e-56  183  78  367  75  362  323  30  320 # BchG based on BackBLAST previously
AAM72500.1  gnl|extdb|OZ401_002285  34.722  9.26e-45  153  76  367  75  353  306  7   287 # BchK based on BackBLAST previously
AAM72500.1  gnl|extdb|OZ401_003162  34.014  6.53e-40  141  78  367  75  359  321  23  308 # Extra gene? UbiA?
```
(More evidence than just this, e.g., previous BackBLAST searches, were also used to identify these 4 gene candidates)

Equivalents in L227-5C were also identified in Extended Data Table 2 and by BLASTP.

Summary of candidates by locus tag:
```
gene-candidate  L227-S17      L227-5C
chlG?           OZ401_002374  HXX20_02935
bchG            OZ401_001677  HXX20_04355
bchK            OZ401_002285  HXX20_10895
Second bchK?    OZ401_003162  HXX20_17010
```

These candidates will be added into the unaligned input files in the `input` folder below.

#### ChlG/BchG phylogeny that also includes BchK
This first phylogeny was done to confirm whether one of the chlorophyll synthases in the strain L227-S17 genome (and L227-5C genome) was 
likely BchK.

Align
```bash
conda activate phylogenies

# strip off comments
seqtk seq -C input/ChlG_BchG_BchK.faa > ChlG_BchG_BchK_ren.faa

# clustalo 1.2.4
clustalo -v -v -i ChlG_BchG_BchK_ren.faa -o ChlG_BchG_BchK_aln.faa 2>&1 | tee ChlG_BchG_BchK_aln.log

# Gblocks v0.91b
Gblocks ChlG_BchG_BchK_aln.faa -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee ChlG_BchG_BchK_mask.log
mv ChlG_BchG_BchK_aln.faa_GB01 ChlG_BchG_BchK_mask.faa
```

Gblocks output summary
```
Original alignment: 451 positions
Gblocks alignment:  215 positions (47 %) in 9 selected block(s)
```

Make tree
```bash
mkdir -p ChlG_BchG_BchK
iqtree -s "ChlG_BchG_BchK_mask.faa" -T 4 --prefix "ChlG_BchG_BchK/masked" -B 1000 -m TEST --boot-trees
# Best-fit model: LG+F+I+G4 chosen according to BIC

# Supplement: try tree again with no mask for comparison
iqtree -s "ChlG_BchG_BchK_ren.faa" -T 4 --prefix "ChlG_BchG_BchK/unmasked" -B 1000 -m TEST --boot-trees
# Best-fit model: LG+F+I+G4 chosen according to BIC
```

Results:
- The topology is stable between masked and unmasked. Bootstraps are reasonable.
- Confirmed as BchK based on phylogeny: OZ401_002285, HXX20_10895
 - Also seem to be BchK-like: OZ401_003162, HXX20_17010
- ChlG seems paraphyletic. Separate clades for oxygenic, anoxygenic (Ctepi+CAB), and anoxygenic (Chx). I wonder if BchK arose from within the ChlG clade, disrupting the phylogeny? Removing the BchK clade might help.

#### ChlG/BchG phylogeny (BchK-free)
Also omitted the 2 L227-S17 and 2 L227-5C sequences associated with the BchK group.

Align
```bash
conda activate phylogenies

# strip off comments
seqtk seq -C input/ChlG_BchG.faa > ChlG_BchG_ren.faa

# clustalo 1.2.4
clustalo -v -v -i ChlG_BchG_ren.faa -o ChlG_BchG_aln.faa 2>&1 | tee ChlG_BchG_aln.log

# Gblocks v0.91b
Gblocks ChlG_BchG_aln.faa -t=p -b3=40 -b4=4 -b5=h -e=_GB01 \
  2>&1 | tee ChlG_BchG_mask.log
mv ChlG_BchG_aln.faa_GB01 ChlG_BchG_mask.faa
```

Gblocks output summary
```
Original alignment: 416 positions
Gblocks alignment:  255 positions (61 %) in 7 selected block(s)
```

Make tree
```bash
mkdir -p ChlG_BchG
iqtree -s "ChlG_BchG_mask.faa" -T 4 --prefix "ChlG_BchG/masked" -B 1000 -m TEST --boot-trees
# Best-fit model: LG+F+G4 chosen according to BIC

# Supplement: try tree again with no mask for comparison
iqtree -s "ChlG_BchG_ren.faa" -T 4 --prefix "ChlG_BchG/unmasked" -B 1000 -m TEST --boot-trees
# Best-fit model: LG+F+I+G4 chosen according to BIC
```

This time:
- ChlG and BchG are well separated
- Topology seems stable between phylogenies made with masked and unmasked alignments

Loaded in Dendroscope to make Supplementary Figure 2.

## Annotree GTDB visualization
From the Annotree website, I exported a phylogeny of the Genome Tree Database, release 89, summarized to the class level, for an overview of 
the tree of life. The raw SVG is available in `annotree/annotree_GTDB_r89_class_raw.svg`.

This raw SVG was edited to produce Fig. 1a.

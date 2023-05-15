# Genome bin analysis
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

**NOTE: for each code section provided below, the code ought to be run from within this `genome_bin_analysis` directory.**

In this methods segment, the '_Ca_. Chloroheliales' genome bins that were generated as described in the `metagenome_analysis` folder are analyzed. 
In addition, download instructions are provided for the other genome bins recovered in this study.

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

One of these homology models ('_Ca_. Chx. allophototropha' PscA-like protein) is shown in Fig. 1c.
Another ('_Ca_. Chx. allophototropha' putative CsmA protein) is shown at the bottom of Extended Data Fig. 4.


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

## Gene context around the detected _pscA_-like genes
Checked on the closest BLASTP hits against RefSeq for the 20 genes up/downstream of the detected _pscA_-like genes in the 
two '_Ca._ Chloroheliales' genome bins.

First, installed a convenient BLASTP wrapper:
```bash
# Create the conda env
conda create -n blast -c bioconda blast=2.9.0 entrez-direct=11.0
conda activate blast

# Install the wrapper
git clone https://github.com/jmtsuji/basic-sequence-analysis.git
cd basic-sequence-analysis
git checkout 447ccc3
PATH="${PATH}:${PWD}/scripts"
cd ..

# Install a wrapper dependency
git clone https://github.com/dib-lab/2018-ncbi-lineages.git
cd 2018-ncbi-lineages
git checkout 0d41546
cd ..

# Download the RefSeq database if you don't already have it
#   For this paper, RefSeq was downloaded on Oct. 10th, 2019
#   You might want to download this in a shared area on your server given that it is a large file and might be used by others
#   See https://www.ncbi.nlm.nih.gov/books/NBK279680/
output_dir="refseq_protein"
mkdir -p "${output_dir}"
cd "${output_dir}"

# update_blastdb.pl version 581818
#   Part of the entrez-direct conda install
update_blastdb.pl \
  --decompress \
  --verbose \
  refseq_protein 2>&1 | \
  tee update_blastdb.log

# Get taxonomy files
update_blastdb.pl \
  --decompress \
  --verbose \
  taxdb 2>&1 | \
  tee update_blastdb_taxdb.log

# Also get taxdump files; these are not used here but can be handy in some circumstances
mkdir -p taxdump
cd taxdump
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xvzf -
cd ..

# Lastly, set the hard-coded variables in the wrapper script to match the location of your BLAST DBs and 2018-ncbi-lineages repo
which run_blastp_parallel.sh # get the file location
# Then edit the global variables in the first few lines of the code using your favourite text editor.
```

Then run BLASTP on the provided ORF prediction subsets (that represent predicted proteins of the 20 up/downstream genes of pscA):
```bash
cd gene_context

run_blastp_parallel.sh Capt_00887_context_20.faa Capt_00887_context_20_blastp.tsv 1e-10 5 10 \
  2>&1 | tee run_blastp_parallel_Capt.log

run_blastp_parallel.sh Chx_L227_5C_00166_context_20.faa Chx_L227_5C_00166_context_20_blastp.tsv 1e-10 5 10 \
  2>&1 | tee run_blastp_parallel_L227_5C.log
```
The output TSV files (included in this folder) are then combined to make Supplementary Data 1.
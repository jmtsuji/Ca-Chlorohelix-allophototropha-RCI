# Supplementary files
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

This folder contains (or points to) supplementary data files referenced in the manuscript text that could be particularly helpful to understand the analyses performed.

### Supplementary data tables
- `Data_01_RCI_gene_context.xlsx` - Gene context of the the novel _pscA_-like genes detected among '_Ca_. Chloroheliales' members. BLASTP results are shown of the predicted primary sequences of 20 genes up/downstream of the novel _pscA_-like genes.
- `Data_02_bidirectional_BLASTP_results.xlsx` - Raw results of the data presented in Fig. 2 in the manuscript
- `Data_03_MAG_relative_abundance_table.xlsx` - Relative abundances of environmental MAGs based on read recruitment to the set of 44 lake metagenomes.
- `Data_04_MAG_expression_tables.xlsx` - Relative expression of environmental MAGs based on read recruitment to the three lake metatranscriptomes (each with three replicates)
- `Data_05_Chloroheliales_MAG_gene_expression.xlsx` - Gene up/down regulation relative to the single-copy marker gene _dnaK_ for the two environmental MAGs classified to the "_Ca._ Chloroheliales" order, based on read recruitment from Lake 221/304 metatranscriptomes
- `Data_06_ASV_table.xlsx` - results of 16S rRNA gene amplicon sequencing of the enrichment cultures presented in this study

### Important reference files
In addition, you might be interested in the following reference files:
- Custom **profile hidden Markov models (HMMs)** for the RCI protein and FMO protein of "_Ca_. Chx. allophototropha" are available in the `genome_bin_analysis/hidden_markov_models` folder in this repo.
- Nucleotide and amino acid sequences of **key photosynthesis genes** in "_Ca_. Chx. allophototropha" can be found in the `genome_bin_analysis/sequences_of_interest` folder in this repo.
- **Homology models** of key photosynthesis genes can be found under the `genome_bin_analysis/homology_models` folder in this repo. 
- **Amino acid sequence alignments** of key photosynthesis genes can be found in the `comparative_genomics/alignments_and_phylogenies` folder in this repo.
- The set of **756 MAGs** with original gene accession numbers and annotations are available at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110) under `lake_survey_MAGs.tar.gz`.

### Detailed additional files
These files might also be of interest if you are looking for fine details:
- The full output of I-TASSER for homology modelling (for those who really want to dive deep into the data) is available at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110).
- ATLAS config files containing the settings used for enrichment culture metagenome analysis are available under `metagenome_analysis/atlas_Chx_allophototropha` and `metagenome_analysis/atlas_L227_5C`
- ATLAS config files for lake survey metagenome/metatranscriptome analysis are available under `lake_survey_analysis/lake_metagenomes` and `lake_survey_analysis/lake_metatranscriptomes`

Enjoy!

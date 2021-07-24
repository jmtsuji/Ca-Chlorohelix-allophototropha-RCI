# Supplementary files
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper  
Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

This folder contains supplementary data files referenced in the manuscript text that could be particularly helpful to understand the analyses performed:
- `Data_01_RCI_gene_context.xlsx` - Gene context of the the novel _pscA_-like genes detected among '_Ca_. Chloroheliales' members. BLASTP results are shown of the predicted primary sequences of 20 genes up/downstream of the novel _pscA_-like genes.
- `Data_02_bidirectional_BLASTP_results.xlsx` - Raw results of the data presented in Fig. 2 in the manuscript
- `Data_03_MAG_relative_abundance_table.xlsx` - Relative abundances of environmental MAGs based on read recruitment to the set of 44 lake metagenomes.
- `Data_04_MAG_expression_tables.xlsx` - Relative expression of environmental MAGs based on read recruitment to the three lake metatranscriptomes (each with three replicates)
- `Data_05_Chloroheliales_MAG_gene_expression.xlsx` - Gene up/down regulation relative to the single-copy marker gene _dnaK_ for the two environmental MAGs classified to the "_Ca._ Chloroheliales" order, based on read recruitment from Lake 221/304 metatranscriptomes
- `Data_06_ASV_table.xlsx` - results of 16S rRNA gene amplicon sequencing of the enrichment cultures presented in this study

In addition:
- Custom profile hidden Markov models and homology models (mentioned in the manuscript) can be found under the `genome_bin_analysis` folder in this repo. 
- The full output of I-TASSER for homology modelling (for those who really want to dive deep into the data) is available at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110).
- ATLAS config files containing the settings used for enrichment culture metagenome analysis are available under `metagenome_analysis/atlas_Chx_allophototropha` and `metagenome_analysis/atlas_L227_5C`
- The set of 756 MAGs with original gene accession numbers and annotations are available at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110)
- ATLAS config files for lake survey metagenome/metatranscriptome analysis are available under `lake_survey_analysis/atlas_metagenomes` and `lake_survey_analysis/atlas_metatranscriptomes`

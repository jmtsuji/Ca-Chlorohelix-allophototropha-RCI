# Code repository for the cultivation of '_Ca._ Chlorohelix allophototropha'

[![DOI](https://zenodo.org/badge/273164101.svg)](https://zenodo.org/badge/latestdoi/273164101)
[![biorxiv](https://img.shields.io/badge/biorxiv-10.1101%2F2020.07.07.190934-green)](https://doi.org/10.1101/2020.07.07.190934)

Copyright Jackson M. Tsuji, Neufeld Research Group, 2021

This Github repository provides additional intermediate data files and code related to the publication by Tsuji and colleagues, 2021, 
"Type I photosynthetic reaction center in an anoxygenic phototrophic member of the _Chloroflexota_," which is currently deposited in [biorxiv](https://doi.org/10.1101/2020.07.07.190934).

This repo also relies on a supplementary Zenodo data repository containing larger files, available [here](https://doi.org/10.5281/zenodo.3930110).

The main focus of this repo is the genomic analysis of '_Candidatus_ Chlorohelix allophototropha', the first cultivated representative of a novel 
order (the '_Ca_. Chloroheliales') within the _Chloroflexota_ phylum. '_Ca_. Chx. allophototropha' -- and one of its relatives that we also managed 
to cultivate for a short time -- uses Type I photosynthetic reaction centers instead of Type II RC's like other phototrophs in the _Chloroflexota_ 
phylum, and this has important implications for the evolution of photosynthesis. (Take a look at the manuscript mentioned above for more details.)

## Repo contents
The repo is organized into four main folders representing the analyses performed in this study:
- `16S_rRNA_gene_amplicon_analysis`: download and analysis of 16S rRNA gene amplicon data from '_Ca_. Chloroheliales' enrichment cultures
  - Used to produce Supplementary Data 6
- `metagenome_analysis`: download, assembly, and binning of metagenome data from '_Ca_. Chloroheliales' enrichment cultures. Unassembled 
read-based analysis are also performed. Lastly, additional environmental metagenomes are also scanned for '_Ca_. Chloroheliales' sequences.
  - Used to produce Extended Data Fig. 1
- `genome_bin_analysis`: download of all genome bins produced in this study, including the two bins classified to the novel '_Ca_. Chloroheliales' order. 
Analysis of the '_Ca_. Chloroheliales' genome bins, including detection of the novel Type I reaction center-associated genes, homology modelling, and 
HMM development, are also included.
  - Used to produce Supplementary Data 1, Fig. 1c and the bottom portion of Extended Data Fig. 5.
- `comparative_genomics`: phylogenetic and functional gene comparisons between the novel '_Ca_. Chloroheliales' bins and other phototrophs or members 
of the _Chloroflexota_ phlylum.
  - Used to produce the bulk of the data presentations in the paper: Fig. 1a-b; Fig. 2, Supplementary Data 2, and Extended Data Table 1; Extended Data Figs. 2,4-6,8-10, and Extended Data Table 2.
- `lake_survey_analysis`: analyze Boreal Shield lake metagenome and metatranscriptome data
  - Used to produce Fig. 3, Extended Data Fig. 7, and Supplementary Data 3-5

In addition, summary folders contain data of particular relevance to readers:
- The `data_presentations` folder contains the individual tables and figures (vector format where possible) that are 
included in the manuscript, and
- The `supplementary_files` folder contains Suppementary Data summaries that are referenced in the manuscript (such as ASV tables and raw results of bidirectional BLASTP).

## How to use this repo
Each of the four subfolders above contains a README file with the code used to analyze the data. Key intermediate files are included in folders 
for your reference. Some key usage cases for the repo include:
- Quickly download the datasets used in this paper (instructions are provided)
- Check how a certain analysis was performed
- Analyze a specific intermediate data file included in the repo

In general, it helpful to have the manuscript itself (see link at top of repo) on hand for comparison while reading this repo.

## Issues
To make this repo, I streamlined several aspects of the raw code to make it more presentable and easy to follow. As such, it is possible that 
syntax errors were introduced along the way. If you notice that any code block does not run, please post a Github issue, and I'll get back to 
you as soon as possible. You can also post an issue if you have questions/concerns about an aspect of the data analysis.

Enjoy!

Jackson M. Tsuji, PhD

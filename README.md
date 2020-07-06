# Code repository for the cultivation of '_Ca._ Chlorohelix allophototropha'
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020

TODO: Add Zenodo button

This Github repository provides additional intermediate data files and code related to the publication by Tsuji and colleagues, 2020, 
"Anoxygenic phototrophic _Chloroflexota_ member uses a Type I reaction center", which is currently deposited in 
[biorxiv]().

This repo also relies on a supplementary Zenodo data repository containing larger files, available [here](https://doi.org/10.5281/zenodo.3930110).

## Repo contents
The repo is organized into four main folders representing the analyses performed in this study:
- `16S_rRNA_gene_amplicon_analysis`: download and analysis of 16S rRNA gene amplicon data from '_Ca_. Chloroheliales' enrichment cultures
  - Used to produce Supplementary Data 1
- `metagenome_analysis`: download, assembly, and binning of metagenome data from '_Ca_. Chloroheliales' enrichment cultures. Unassembled 
read-based analysis are also performed. Lastly, additional environmental metagenomes are also scanned for '_Ca_. Chloroheliales' sequences.
  - Used to produce Extended Data Fig. 1
- `genome_bin_analysis`: download of all genome bins produced in this study, including the two bins classified to the novel '_Ca_. Chloroheliales' order. 
Analysis of the '_Ca_. Chloroheliales' genome bins, including detection of the novel Type I reaction center-associated genes, homology modelling, and 
HMM development, are also included.
  - Used to produce Supplementary Data 2, Figure 2c, and the bottom portion of Extended Data Fig. 4.
- `comparative_genomics`: phylogenetic and functional gene comparisons between the novel '_Ca_. Chloroheliales' bins and other phototrophs or members 
of the _Chloroflexota_ phlylum.
  - Used to produce the bulk of the data presentations in the paper: Figure 2a-b; Figure 3 and Supplementary Data 3; Extended Data Figs. 2-8.

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
Enjoy!

Jackson M. Tsuji  
PhD Candidate, University of Waterloo, Canada

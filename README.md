# Code repository for the cultivation of '_Ca._ Chlorohelix allophototropha'

[![DOI](https://zenodo.org/badge/273164101.svg)](https://zenodo.org/badge/latestdoi/273164101)
[![biorxiv](https://img.shields.io/badge/biorxiv-10.1101%2F2020.07.07.190934-green)](https://doi.org/10.1101/2020.07.07.190934)

Copyright Jackson M. Tsuji, Neufeld Research Group, 2023

This Github repository provides supplementary data files and code related to the publication by Tsuji and colleagues, 
2023, "Anoxygenic phototrophic _Chloroflexota_ member uses a Type I reaction center," which is currently deposited in 
[bioRxiv](https://doi.org/10.1101/2020.07.07.190934).

This repo also relies on a supplementary Zenodo data repository containing larger files, available 
[here](https://doi.org/10.5281/zenodo.3930110).

The main focus of this repo is the genomic analysis of "_Candidatus_ Chlorohelix allophototropha", the first cultivated 
representative of a novel order (the "_Ca_. Chloroheliales") within the _Chloroflexota_ phylum. "_Ca_. Chx. 
allophototropha" uses Type I photosynthetic reaction centers instead of Type II RC's like other phototrophs in the 
_Chloroflexota_ phylum, and this has important implications for the evolution of photosynthesis. (Take a look at the 
manuscript mentioned above for more details.)

## Repo contents
Three key folders contain the content you might have come to this repo for:
- `tables_figures_etc` - raw versions of the tables, figures, and supplementary data files shown in the paper. These 
  are helpful if you want to directly view the vector-format versions of our graphics or want to directly copy and 
  paste from tables that aren't in a PDF.
- `supporting_files` - contains some key files used in the paper that aren't included in the supplementary data, like 
  custom HMMs built for photosynthesis proteins associated with "_Ca_. Chx. allophototropha".
  - Some supporting files are too large for Github and are instead stored on Zenodo 
    [here](https://doi.org/10.5281/zenodo.3930110)
- `analysis_code` - shows code for analyses that we performed that weren't out-of-the-box. We also include 
  miscellaneous code for data viz etc. performed in the paper, for those who are interested.

### Analysis code: contents
The `analysis_code` folder is further divided into the following sub-folders:

- `source_data` - contains links to the different publicly available sequences datasets we used
- `physiology` - contains code and some small data files for analysis of spectroscopy data
- `16S_rRNA_gene_amplicon` - shows the three different methods used to process 16S rRNA gene amplicon sequencing data 
  in this work, and includes count table data.
- `culture_metagenomics` - includes code used to obtain genome bins from metagenomes of early enrichment cultures
- `genomics` - includes code for: hybrid assembly approaches used to obtain closed genome (bins), and read mapping to 
  determine relative abundances across all genome bins from `culture_metagenomics`
- `lake_omics` - shows the analysis methods used to process lake metagenome/metatranscriptome data
- `phylogenomics` - code used to build the _Chloroflexota_ species tree, search for photosynthesis genes among 
  _Chloroflexota_ members, and build phylogenies of photosynthesis-related genes

## How to use this repo
Each of the major subfolders above contains a README file details about the folder contents and/or code. Some key 
usage cases for the repo include:
- Find supplementary files not included in the supplemental data in the paper
- Obtain slightly more raw (e.g., vector format or Word/Excel format) data presentations
- Find how to download the datasets used in this paper
- Check how a certain analysis was performed
- Analyze a specific intermediate data file included in the repo

In general, it helpful to have the manuscript itself (see link at top of repo) on hand for comparison while reading 
this repo.

## Versions of this repo
The Chlorohelix manuscript has been updated several times since the initial discovery of this exciting strain as we 
have added more experimental data. The following Github repo versions correspond to major bioRxiv versions of the 
manuscript:
- repo [v3.0.0](https://github.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/tree/v3.0.0): corresponds to a second 
  major manuscript revision [here](https://www.biorxiv.org/content/10.1101/2020.07.07.190934v4) (2023)
- repo [v2.0.1](https://github.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/tree/v2.0.1): corresponds to a major 
  manuscript revision [here](https://www.biorxiv.org/content/10.1101/2020.07.07.190934v3) (2021)
- repo [v1.0.0](https://github.com/jmtsuji/Ca-Chlorohelix-allophototropha-RCI/tree/v1.0.0): corresponds to the initial 
  bioRxiv submission [here](https://www.biorxiv.org/content/10.1101/2020.07.07.190934v1) (2020)

## Issues
To make this repo, I streamlined several aspects of the raw code to make it more presentable and easy to follow. As 
such, it is possible that syntax errors were introduced along the way. If you notice that any code block does not run, 
please post a Github issue, and I'll get back to you as soon as possible. You can also post an issue if you have 
questions/concerns about an aspect of the data analysis.

Enjoy!

Jackson M. Tsuji, PhD  
On behalf of the manuscript's coauthor team

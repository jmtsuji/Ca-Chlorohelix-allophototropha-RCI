# Supplementary files/data
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper

Copyright Jackson M. Tsuji, Neufeld Research Group, 2024

This folder contains (or points to) supplementary data files referenced in the manuscript (or supplementary) text that
aren't included in the Extended Data files of the manuscript.

## Folder contents
- `hidden_markov_models`: contains custom **profile hidden Markov models (HMMs)** for the RCI protein and FMO protein of
  "_Ca_. Chx. allophototropha" that were used to perform sequence searches for the manuscript.
- `homology_models`: contains PDB-format **predicted strutures** for the PscA-like proteins of "_Ca_. Chx.
  allophototropha" strain L227-S17 and strain L227-5C; the predicted structure of the PscA-like protein for "_Ca_. Chx.
  allophototropha" was visualized in an Extended Data figure in the manuscript.
  - As a bonus, although not described in the manuscript, predicted structures for FmoA, CsmA, and CbbL (i.e., RbcL) are
    also included. These were predicted using the I-TASSER webserver in May 2019 (FmoA, strain L227-S17), March 2020
    (FmoA, strain L227-5C; also CsmA and CbbL, strain L227-S17), and April 2020 (CsmA and CbbL, strain L227-5C).
- `BchXYZ_phylogenies`: individual phylogenies (Newick format) and sequence alignments, with and without masking of
  non-informative regions, for BchX, BchY, and BchZ genes as described in the manuscript text.
  - As indicated below, detailed phylogeny files for all gene phylogenies shown in the manuscript (e.g., Fig. 5) can be
    found in the `analysis_code/phylogenomics/alignments_and_phylogenies` folder in this repo.
- `sequences_of_interest`: sequences of key photosynthesis-associated gene homologs from "_Ca_. Chloroheliales" members.
  Nucleotide gene sequences and predicted protein sequences are included.

## Other files of interest
- **Amino acid sequence alignments and phylogenies** of key photosynthesis genes can be found in the
  `analysis_code/phylogenomics/alignments_and_phylogenies` folder in this repo.
- The set of **756 MAGs** with original gene accession numbers and annotations are available at the
  [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110) under `lake_survey_MAGs.tar.gz`.
- ATLAS config files containing the settings used for enrichment culture metagenome analysis are available under
  `culture_metagenomics/atlas_Chx_allophototropha` and `culture_metagenomics/atlas_L227_5C`
- ATLAS config files for lake survey metagenome/metatranscriptome analysis are available under
  `lake_omics/lake_metagenomes` and `lake_omics/lake_metatranscriptomes`
- The full output of I-TASSER for homology modelling (for those who really want to dive deep into the data) is available
  at the [Zenodo data repository corresponding to this code repo](https://doi.org/10.5281/zenodo.3930110).

Enjoy!

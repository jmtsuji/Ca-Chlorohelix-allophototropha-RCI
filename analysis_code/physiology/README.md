# Phototrophic physiology
Part of the larger '_Ca._ Chlorohelix allophototropha' Type I reaction center paper

Copyright Jackson M. Tsuji, Neufeld Research Group, 2024

This folder contains analysis code and some small raw data files for spectroscopic analyses of the 
"_Ca_. Chx. allophototropha" L227-S17 culture compared to reference strains. Three major spectroscopic 
analyses were performed:

- `in-vivo`: _in vivo_ absorption spectra of cultures, presented in Fig. 1a and Extended Data Fig. 2a
- `hplc`: HPLC-based separation of bacteriochlorophyll _c_ and _a_, presented in Fig. 1b-c and Supplementary Fig. 1.
- `light-dark-test`: partial _in vivo_ spectra of cultures incubated in the light vs. dark.

For each folder, you will find the following contents:
- `input`: folder containing small input data files. For `hplc`, the data are larger and so are stored on the 
[Zenodo repository](https://zenodo.org/doi/10.5281/zenodo.3930110) associated with this work, rather than in a
`input` folder.
- iPython notebooks containing analysis code; these are named after their associated data presentation (e.g., Fig1a)
- `output_raw`: raw outputs of the iPython notebook analyses. These were cleaned up and used for figures etc. in the 
manuscript.

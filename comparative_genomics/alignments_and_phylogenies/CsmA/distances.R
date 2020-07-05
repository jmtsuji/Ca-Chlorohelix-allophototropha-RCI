setwd("D:/share/manuscripts/2020_04_14_Ca_Chx_allophototropha_paper/02_edits/05_github/Ca_Chlorohelix_allophototropha_RCI/alignments_and_phylogenies/CsmA/")
library(ape)
library(dplyr)
library(tibble)

tree <- ape::read.tree("CsmA_aligned_masked.treefile")
tbl <- ape::cophenetic.phylo(tree) %>%
  tibble::as_tibble(rownames = "genome")
write.table(tbl, "distances.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

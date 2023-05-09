# Set the working directory as the same directory that this script is in.
library(ape)
library(dplyr)
library(tibble)

tree <- ape::read.tree("CsmA_aligned_masked.treefile")
tbl <- ape::cophenetic.phylo(tree) %>%
  tibble::as_tibble(rownames = "genome")
write.table(tbl, "distances.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Figure S1 plotter

# Load libraries
library(metannoviz) # Used v1.0.2
library(egg)
# N.B., install metannoviz via: 
#   library(devtools)
#   devtools::install_github("metannotate/metannoviz@v1.0.2")
#   You might then need to restart R to complete installation.

# User variables
# N.B., Your working directory should be the current directory of the script.
metannotate_annotations_filepath <- "rpoB_all_annotations.tsv"
bin_table_allophototropha_filepath <- "Capt_MAG_table_to_assembled.tsv"
bin_table_L227_5C_filepath <- "L227_5C_MAG_table_to_assembled.tsv"
bin_table_sample_naming_info_filepath <- "guides/dataset_info_MAGs.tsv"
output_metannotate_info_filepath <- "output/metannotate_data_processed.tsv"
output_plot_filepath <- "output/Figure_S1_raw.pdf"

# N.B., the two MetAnnotate guide files here, your `raw_name` column might have slight spelling differences depending on the details of your MetAnnotate run. You can just look at the first few columns of the "rpoB_all_annotations.tsv" to check.
hmm_info_filepath <- "guides/hmm_info_metannotate.tsv"
dataset_info_filepath <- "guides/dataset_info_metannotate.tsv"

### Process MetAnnotate data
metannotate_data <- metannoviz::read_metannotate_data(metannotate_annotations_filepath)
metannotate_data_mapped <- metannoviz::map_naming_information(metannotate_data, 
                                                              hmm_naming_info_filename = hmm_info_filepath,
                                                              dataset_naming_info_filename = dataset_info_filepath)

# Customize the MetAnnotate bubble plot
fill_values = RColorBrewer::brewer.pal(8, "Set1")[c(1,8,5,3)]
metannotate_plot_raw <- metannoviz::explore(metannotate_data_mapped, plot_type = "bubble")

metannotate_plot <- ggplot(metannotate_plot_raw$data, aes(x = Dataset, y = Closest.Homolog.Family)) +
  geom_point(aes(size = percent_abundance, fill = Closest.Homolog.Phylum), shape = 21,
             colour = "black", alpha = 0.7) +
  geom_text(aes(label = label), size = 2.5) +
  scale_size_continuous(range = c(3,20), guide = "none") +
  scale_fill_manual(values = fill_values) +
  guides(fill = guide_legend(override.aes = list(size = 8), title = "Phylum")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, colour = "black", face = "italic"), 
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_blank()) +
  xlab(NULL) +
  ylab("Family")

### Process bin data
# Read and combine the data
bin_table_allophototropha <- read.table(bin_table_allophototropha_filepath, header = TRUE,
                                        sep = "\t", stringsAsFactors = FALSE) %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(cols = c(Cfx3_2018, Cfx3_2019), names_to = "Sample_ID",
                      values_to = "Relative_abundance")
bin_table_allophototropha$MAG_ID <- paste0("Cfx3_", bin_table_allophototropha$MAG_ID)
bin_table_allophototropha$Culture <- "Cfx3"

bin_table_L227_5C <- read.table(bin_table_L227_5C_filepath, header = TRUE,
                                        sep = "\t", stringsAsFactors = FALSE) %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(cols = c(Cfx_L227_5C), names_to = "Sample_ID",
                      values_to = "Relative_abundance")
bin_table_L227_5C$MAG_ID <- paste0("Cfx2_", bin_table_L227_5C$MAG_ID)
bin_table_L227_5C$Culture <- "Cfx2"

# Combine and filter out low abundance taxa
bin_table <- dplyr::bind_rows(bin_table_allophototropha, bin_table_L227_5C) %>%
  dplyr::filter(Relative_abundance >= 0.05)

# Make plotting labels
bin_table$Relative_abundance_label <- round(bin_table$Relative_abundance, digits = 1)
bin_table$bin_label <- paste0(bin_table$MAG_ID, " (f__", bin_table$Family,
                              "; g__", bin_table$Genus, "; ", bin_table$Completeness,
                              "/", bin_table$Contamination, ")")

# Sort genomes
bin_table <- dplyr::arrange(bin_table, Sample_ID, Phylum, Class, Order, Family, Genus, Species)
bin_table$bin_label <- factor(bin_table$bin_label, levels = rev(unique(bin_table$bin_label)),
                              ordered = TRUE)

# Samples
bin_table_sample_naming_info <- read.table(bin_table_sample_naming_info_filepath, sep = "\t",
                                           header = TRUE, stringsAsFactors = FALSE) %>%
  tibble::as_tibble()
bin_table$Sample_ID <- plyr::mapvalues(bin_table$Sample_ID, from = bin_table_sample_naming_info$raw_name,
                                       to = bin_table_sample_naming_info$Dataset) %>%
  factor(levels = bin_table_sample_naming_info$Dataset, ordered = TRUE)


# Plot
fill_values = RColorBrewer::brewer.pal(8, "Set1")[c(1,8,5,2,6,3)]
bin_plot <- ggplot(bin_table, aes(x = Sample_ID, y = bin_label)) +
  geom_point(aes(size = Relative_abundance, fill = Phylum), 
             shape = 21, alpha = 0.7, colour = "black") +
  geom_text(aes(label = Relative_abundance_label), size = 2.5) +
  scale_fill_manual(values = fill_values) +
  scale_size_continuous(range = c(3, 20), guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 6))) + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(face = "italic")) +
  ylab(NULL) +
  xlab(NULL)

combined_plots <- egg::ggarrange(plots = list(metannotate_plot, bin_plot), nrow = 2, 
                                 newpage = TRUE, heights = c(5,3.5))

# Save output
pdf(file = output_plot_filepath, width = 200 / 25.4, height = 240 / 25.4)
print(combined_plots)
dev.off()

write.table(metannotate_plot$data, file = output_metannotate_info_filepath, quote = FALSE, 
            sep = "\t", col.names = TRUE, row.names = FALSE)

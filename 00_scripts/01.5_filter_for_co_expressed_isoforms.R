setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

# Load necessary libraries
library(dplyr)
library(readr)

# Read the CSV file
data <- read_csv("./01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv")

# Filter the data for isoform-specific peptides and at least one '1' in the sample columns
filtered_data <- data %>%
  filter(category == "isoform-specific") %>% # Select isoform-specific peptides
  filter(sn_pep31811_pep_07022024 == 1 | sn_pep31812_pep_07022024 == 1 | sn_pep31813_pep_07022024 == 1 | sn_pep31814_pep__07022024 == 1) %>% # At least one '1' in the last four columns
  group_by(gene_name) %>% # Group by gene
  filter(n_distinct(transcript_name) >= 2) %>% # Keep genes with two or more unique isoforms
  ungroup()

# Check if the directory exists
output_dir <- "./01_filter_data/01.5_filtered_co_expressed_isoforms"
if (!dir.exists(output_dir)) {
  # Create the directory
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory already exists:", output_dir, "\n")
}
# Write the filtered data to a CSV file
write_csv(filtered_data, paste0(output_dir, "/filtered_co_expressed_isoforms.csv"))

# Get the number of non-redundant entries in each column
non_redundant_counts <- filtered_data %>%
  summarize(across(everything(), n_distinct))

# Write the non-redundant counts to a CSV file
write_csv(non_redundant_counts, paste0(output_dir, "/non_redundant_counts.csv"))

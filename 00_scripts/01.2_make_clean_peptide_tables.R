setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

# Make clean peptide tables for input for peptide-protein mapping
# Removes flanking AAs as well as non-alpha characters


# Load necessary library
library(dplyr)

# Function to extract base peptide sequence and remove non-alphabet characters
extract_base_sequence <- function(peptide) {
  peptide <- gsub("^[A-Z]\\.|\\.[A-Z]$", "", peptide)
  gsub("[^A-Za-z]", "", peptide) # Remove non-alphabet characters
}

# Function to clean a single file
clean_file <- function(input_path, output_path) {
  data <- read.csv(input_path)
  cleaned_data <- data %>%
    select(ProteinId, GeneSymbol, PeptideSequence) %>%
    mutate(PeptideSequence = sapply(PeptideSequence, extract_base_sequence))
  write.csv(cleaned_data, output_path, row.names = FALSE)
}

# Directory containing the input files
input_directory <- "./00_peptide_data_20240702"
# Directory to save the cleaned files
output_directory <- "./01_filter_data/01.2_make_clean_peptide_tables"

# Ensure the output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# List all CSV files in the input directory
file_list <- list.files(input_directory, pattern = "*.csv", full.names = TRUE)

# Process each file
for (input_path in file_list) {
  # Create the output file path
  output_path <- file.path(output_directory, basename(input_path))
  # Clean the file
  clean_file(input_path, output_path)
}



# Summary stats of the peptide data


# ---- Import libraries ----
library(dplyr)
library(readr)

# ---- Input and Output directories ----
input_dir <- '../mouse_tissue_peptide_hits'
output_dir <- '../output/03_view_peptide_data'

# ---- Create output directory ----
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ---- Get all the files in the input directory ----
# This corresponds to the mouse tissue peptide results from Tian
file_list <- list.files(input_dir, full.names = TRUE)


# ---- Loop through each file and get the summary stats ----

process_file <- function(file_path) {
  # Read the data
  data <- read_csv(file_path)
  
  # Clean the peptide sequences
  data <- data %>%
    mutate(PeptideSequenceClean = sub("^.\\.(.*)\\..$", "\\1", PeptideSequence))
  
  # Get statistics
  distinct_proteins <- data %>% distinct(ProteinId) %>% nrow()
  distinct_proteins_tr <- data %>% filter(grepl("^tr\\|", ProteinId)) %>% distinct(ProteinId) %>% nrow()
  distinct_proteins_sp <- data %>% filter(grepl("^sp\\|", ProteinId)) %>% distinct(ProteinId) %>% nrow()
  total_genes <- data %>% distinct(GeneSymbol) %>% nrow()
  total_genes_tr <- data %>% filter(grepl("^tr\\|", ProteinId)) %>% distinct(GeneSymbol) %>% nrow()
  total_genes_sp <- data %>% filter(grepl("^sp\\|", ProteinId)) %>% distinct(GeneSymbol) %>% nrow()
  unique_group_ids <- data %>% distinct(GroupId) %>% nrow()
  cleaned_peptides <- data %>% select(PeptideSequenceClean)
  
  # extract filename without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Define the output file path
  output_file <- file.path(output_dir, paste0(file_name, '_peptide_summary.txt'))
  
  
  # Open a connection to a file
  sink(output_file)
  
  # Your analysis and printing results
  cat("Number of distinct ProteinId:", distinct_proteins, "\n")
  cat("Number of distinct ProteinId (tr):", distinct_proteins_tr, "\n")
  cat("Number of distinct ProteinId (sp):", distinct_proteins_sp, "\n")
  cat("Total number of genes:", total_genes, "\n")
  cat("Number of genes (tr):", total_genes_tr, "\n")
  cat("Number of genes (sp):", total_genes_sp, "\n")
  cat("Unique number of GroupId:", unique_group_ids, "\n")
  
  # If you want to print the cleaned peptide sequences as well
  cat("Cleaned Peptide Sequences:\n")
  print(cleaned_peptides)
  
  # Close the connection to the file
  sink()
}

# Loop through each file
lapply(file_list, process_file)
#----
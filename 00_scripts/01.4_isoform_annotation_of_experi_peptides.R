setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

# annoate isoform specificity of experimentally identified peptides

# ------------------ load required libraries ------------------

library(tidyverse)

# read in file for the 

# read in in silico tryptic peptide annotations
pep_map <- read_csv('./01_filter_data/01.5_peptide_to_isoform_mapping/tryptic_peptide_to_mouse_protein_mapping.csv')

# read in experimentally identified peptides
experiment_files <- list.files(path='./01_filter_data/01.2_make_clean_peptide_tables/', pattern='p*csv', full.names=TRUE)

# initialize list to hold the data frames
experiment_dfs <- list()

# Read each experimental file into a data frame and store in the list
for (file in experiment_files) {
  experiment_name <- tools::file_path_sans_ext(basename(file))
  experiment_df <- read_csv(file)
  experiment_df <- experiment_df[1:3] %>%
    select(PeptideSequence) %>%
    rename(!!experiment_name := PeptideSequence)
  experiment_dfs[[experiment_name]] <- experiment_df
}

# Function to check if peptides are present in the experimental data frames
check_presence <- function(peptide, experiment_df) {
  ifelse(peptide %in% experiment_df[[1]], 1, 0)
}


# Add columns to pep_map for each experiment
for (experiment_name in names(experiment_dfs)) {
  print(paste("Adding column for", experiment_name))
  
  pep_map[[experiment_name]] <- sapply(pep_map$pep_seq, check_presence, experiment_df = experiment_dfs[[experiment_name]])
  
  print(paste("Column added for", experiment_name))
  print(head(pep_map[[experiment_name]]))
}


# write out peptide map with the peptides found in Tian's experiments marked
write_csv(pep_map, './01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv')


#%%

# Read each experimental file into a data frame and store in the list
for (file in experiment_files) {
  print(paste("Processing file:", file))
  
  experiment_name <- tools::file_path_sans_ext(basename(file))
  print(paste("Experiment name:", experiment_name))
  
  experiment_df <- read_csv(file)
  print("First few rows of experiment_df before selection:")
  print(head(experiment_df))
  
  experiment_df <- experiment_df %>%
    select(PeptideSequence) %>%
    rename(!!experiment_name := PeptideSequence)
  print("First few rows of experiment_df after selection and renaming:")
  print(head(experiment_df))
  
  experiment_dfs[[experiment_name]] <- experiment_df
  print(paste("Stored data frame for", experiment_name))
  
}


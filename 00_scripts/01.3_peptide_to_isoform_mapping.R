setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

# Peptide to isoform mapping using the cleaver package
# Modified from: https://github.com/sheynkman-lab/protein-isoform-detection/tree/e3a53bbec5bdf933d57745e0b87c2de4897fa03e/Depreciated-code


# ---- Load necessary libraries ----
library(cleaver)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(progress)
library(stringr)
library(tidyr) 
library(data.table)
library(janitor)

# read in fasta
fasta <- readAAStringSet("./00_gencode_mouse_models/gencode.vM35.pc_translations.fa")

# tryptic digest
peptides <- cleave(fasta, enzym = "trypsin", missedCleavages= 0, unique=TRUE)
peptides_df <- as.data.frame(peptides)
pep.df.raw <- as.data.table(peptides) # data table for faster processing

# tidy digested peptides
pep.df <- pep.df.raw %>%
  select(-group) %>%
  separate(col = group_name, into = paste0("col", 1:8), sep = "\\|") %>%
  select(gene_name = col7, transcript_name = col6, pep_seq = value) %>%
  mutate("pep_length" = str_length(pep_seq)) %>%
  filter(pep_length >= 5 & pep_length <= 75)

# get number of isoform mappings
pep.df <- pep.df %>%
  group_by(gene_name) %>%
  mutate(iso_total = n_distinct(transcript_name)) %>%
  ungroup()

# categorize peptides by number of isoform mappings 
collapsed.pep.df <- pep.df %>% 
  group_by(gene_name, iso_total, pep_seq) %>% 
  summarize(transcript_name= toString(transcript_name), count = n()) %>% 
  mutate(category = case_when(
    count == iso_total ~ "constitutive", 
    count == 1 ~ "isoform-specific",
    TRUE ~ "isoform-informative"
  ))

# get gene specificity of peptide mappings
df.gene.nonspecific <- collapsed.pep.df %>% 
  get_dupes(pep_seq) %>% 
  group_by(pep_seq) %>% 
  summarize(multigene.IDs = toString(transcript_name), genes = toString(gene_name)) %>% 
  ungroup()

# join isoform mapping and gene specificity results
collapsed.pep.df.joined <- collapsed.pep.df %>% 
  left_join(df.gene.nonspecific, by = "pep_seq")

collapsed.pep.df.joined <- collapsed.pep.df.joined %>% 
  mutate(gene.specificity = ifelse(is.na(multigene.IDs), "gene-specific", "not gene-specific"))

# write out results
output_dir <- "./01_filter_data/01.5_peptide_to_isoform_mapping/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
output_file <- file.path(output_dir, 'tryptic_peptide_to_mouse_protein_mapping.csv')
write.csv(collapsed.pep.df.joined, file = output_file, row.names = FALSE)







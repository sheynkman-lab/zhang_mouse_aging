# Format an output file for input into Pogo to make Genome Browser tracks
# Restructure the dataframe for input into Pogo
# Outputs are named first by TMT or TOM, their ID number, WTC11, then unfrac or 8frac
setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

library(readr)
library(dplyr)
library(stringr)

# 31811
pq31811_df <- read_csv("01_filter_data/01.2_make_clean_peptide_tables/sn_pep31811_pep_07022024.csv")

restructured_pogo_df <- pq31811_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31811", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/sn_pq31811_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 31812
pq31812_df <- read_csv("01_filter_data/01.2_make_clean_peptide_tables/sn_pep31812_pep_07022024.csv")

restructured_pogo_df <- pq31812_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31812", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/sn_pq31812_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 31813
pq31813_df <- read_csv("01_filter_data/01.2_make_clean_peptide_tables/sn_pep31813_pep_07022024.csv")

restructured_pogo_df <- pq31813_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31813", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/sn_pq31813_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 31814
pq31814_df <- read_csv("01_filter_data/01.2_make_clean_peptide_tables/sn_pep31814_pep__07022024.csv")

restructured_pogo_df <- pq31814_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31814", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/sn_pq31814_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# annotated isoforms
experi_df <- read_csv("01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv")

restructured_pogo_df <- experi_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "experimental_pep", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = pep_seq, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/experimental_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# coexpressed
coexpressed_df <- read_csv("01_filter_data/01.5_filtered_co_expressed_isoforms/filtered_co_expressed_isoforms.csv")

restructured_pogo_df <- coexpressed_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "coexpressed", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = pep_seq, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "02_Peptides2Pogo/coexpressed_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)













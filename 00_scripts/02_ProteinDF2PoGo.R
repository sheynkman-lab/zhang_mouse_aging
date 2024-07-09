# Format an output file for input into Pogo to make Genome Browser tracks
# Restructure the dataframe for input into Pogo
# Outputs are named first by TMT or TOM, their ID number, WTC11, then unfrac or 8frac
setwd("/Volumes/sheynkman/projects/zhang_mouse_aging")

# pq31153
library(readr)
pq31153_df <- read_csv("00_2024-05-29_Tian_Aging_Mouse_SpliceProt/mouse_tissue_peptide_hits/pq31153.csv")

library(dplyr)
library(stringr)

restructured_pogo_df <- pq31153_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31153", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "pq31153_peptide_POGO.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
#PoGo prefers tables to CSVs
#write.csv(restructured_pogo_df, "pq31153_peptide_output_POGO.csv", quote= FALSE)



# pq31193
pq31193_df <- read_csv("00_2024-05-29_Tian_Aging_Mouse_SpliceProt/output/04_make_clean_peptide_tables/pq31193.csv")

restructured_pogo_df <- pq31193_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31193", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "pq31193_peptide_POGO.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# pq31194
pq31194_df <- read_csv("00_2024-05-29_Tian_Aging_Mouse_SpliceProt/output/04_make_clean_peptide_tables/pq31194.csv")

restructured_pogo_df <- pq31194_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31194", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "pq31194_peptide_POGO.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)


# pq31190
pq31190_df <- read_csv("2024-05-29_Tian_Aging_Mouse_SpliceProt/output/04_make_clean_peptide_tables/pq31190.csv")

restructured_pogo_df <- pq31190_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "pq31190", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = PeptideSequence, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "pq31190_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# all cleaved
all_df <- read_csv("00_2024-05-29_Tian_Aging_Mouse_SpliceProt/output/07_filtered_co_expressed_isoforms/filtered_co_expressed_isoforms.csv")

restructured_pogo_df <- all_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "all", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = pep_seq, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "all_peptide_POGO.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# all aggregate
all_df <- read_csv("2024-05-29_Tian_Aging_Mouse_SpliceProt/output/06_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv")

restructured_pogo_df <- all_df %>%
  mutate( #mutate just allows us to modify and add columns!
    Experiment = "all", #create column called "Experiment" that is filled with the experiment type
    PSMs = 1, #create column called "PSMs" and fill with 1
    Quant = 1 #create column called "Quant" and fill with 1
  ) %>%
  select(Experiment, Distinct_Peptide = pep_seq, PSMs, Quant) #the dataframe we're creating will only have these columns

file_path <- "aggregate_peptides.txt"
write.table(restructured_pogo_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)













import pandas as pd

# Load the CSV files
aggregate = '/project/sheynkman/projects/zhang_mouse_aging/2024-05-29_Tian_Aging_Mouse_SpliceProt/output/06_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv'
cleaved = '/project/sheynkman/projects/zhang_mouse_aging/2024-05-29_Tian_Aging_Mouse_SpliceProt/output/07_filtered_co_expressed_isoforms/filtered_co_expressed_isoforms.csv'

df1 = pd.read_csv(aggregate)
df2 = pd.read_csv(cleaved)

# Function to count unique values
def count_unique(df):
    unique_gene_names = df['gene_name'].nunique()
    unique_pep_seqs = df['pep_seq'].nunique()
    unique_transcript_names = df['transcript_name'].apply(lambda x: x.split(', ')).explode().nunique()
    return unique_gene_names, unique_pep_seqs, unique_transcript_names

# Count unique values for each file
unique_counts_aggregate = count_unique(df1)
unique_counts_cleaved = count_unique(df2)

# Print the results
print(f"Aggregate - Unique gene_names: {unique_counts_aggregate[0]}, Unique pep_seqs: {unique_counts_aggregate[1]}, Unique transcript_names: {unique_counts_aggregate[2]}")
print(f"Cleaved - Unique gene_names: {unique_counts_cleaved[0]}, Unique pep_seqs: {unique_counts_cleaved[1]}, Unique transcript_names: {unique_counts_cleaved[2]}")

 import pandas as pd

# Load the Excel file with multiple sheets (peptide data)
peptide_file_path = '00_peptide_data_20240702/unique_peptide_effects_20240702.xlsx'
xls = pd.ExcelFile(peptide_file_path)

# Load the gene expression data
gene_expression_file_path = '01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv'
gene_expression_df = pd.read_csv(gene_expression_file_path)

# Check the columns of gene_expression_df to ensure 'gene_name' exists or find the relevant column
print("Gene Expression Data Columns:", gene_expression_df.columns)

# Define the threshold for significance
padj_threshold = 0.05

# Initialize an empty list to store significant peptides from all sheets
all_significant_peptides = []

# Process each sheet in the Excel file
for sheet_name in xls.sheet_names:
    # Load the data from the current sheet
    df = pd.read_excel(xls, sheet_name=sheet_name)
    
    # Filter for significant changes
    significant_peptides = df[
        (df['sex_padj'] < padj_threshold) |
        (df['age_padj'] < padj_threshold) |
        (df['sex_by_age_padj'] < padj_threshold)
    ]
    
    # Add a column to indicate the tissue type using .loc
    significant_peptides['tissue'] = sheet_name
    
    # Append the significant peptides to the list
    all_significant_peptides.append(significant_peptides)

# Concatenate all significant peptides into a single DataFrame
all_significant_peptides_df = pd.concat(all_significant_peptides, ignore_index=True)

# Check the columns of all_significant_peptides_df
print("Peptide Data Columns:", all_significant_peptides_df.columns)

# Define function to determine change in expression
def determine_change(effect):
    if effect > 0:
        return '+'
    elif effect < 0:
        return '-'
    else:
        return 'same'

# Identify the correct column for gene expression change if 'effect' is not present
# Assuming 'count' is a proxy for expression level in this example
if 'count' in gene_expression_df.columns:
    gene_expression_df['effect'] = gene_expression_df.groupby('gene_name')['count'].diff().fillna(0)
    gene_expression_df['gene_change'] = gene_expression_df['effect'].apply(determine_change)
else:
    print("Error: Relevant column for gene expression change not found in gene_expression_df")
    raise KeyError("Relevant column for gene expression change not found in gene_expression_df")

# Determine change in isoform (peptide) expression
if 'age_effect' in all_significant_peptides_df.columns:
    all_significant_peptides_df['isoform_change'] = all_significant_peptides_df['age_effect'].apply(determine_change)
else:
    print("Error: 'age_effect' column not found in all_significant_peptides_df")
    raise KeyError("'age_effect' column not found in all_significant_peptides_df")

# Assuming 'symbol' should be replaced with the correct column name from all_significant_peptides_df
# Inspect the columns and replace 'symbol' if needed
merge_left_on = 'symbol'  # Adjust if necessary
merge_right_on = 'gene_name'  # Adjust if necessary

if merge_left_on not in all_significant_peptides_df.columns:
    print(f"Error: Column '{merge_left_on}' not found in peptide data")
    raise KeyError(f"Column '{merge_left_on}' not found in peptide data")

if merge_right_on not in gene_expression_df.columns:
    print(f"Error: Column '{merge_right_on}' not found in gene expression data")
    raise KeyError(f"Column '{merge_right_on}' not found in gene expression data")

# Merge the peptide data with gene expression data
merged_df = pd.merge(all_significant_peptides_df, gene_expression_df, left_on=merge_left_on, right_on=merge_right_on)

# Extract relevant columns and rename them
summary_df = merged_df[[merge_left_on, 'gene_change', 'isoform_change', 'pep_seq']]

# Rename columns to 'gene' and 'peptide'
summary_df = summary_df.rename(columns={merge_left_on: 'gene', 'pep_seq': 'peptide'})

# Optionally drop duplicates if necessary
summary_df = summary_df.drop_duplicates()

# Output the results
output_file_path = '04_candidate_peps/gene_peptide_changes_significant.csv'
summary_df.to_csv(output_file_path, index=False)
print(f"Summary of gene and peptide changes has been saved to '{output_file_path}'")
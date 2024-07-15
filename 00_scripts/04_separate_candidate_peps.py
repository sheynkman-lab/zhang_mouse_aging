import pandas as pd
import numpy as np

# Load the Excel file with multiple sheets (peptide data)
peptide_file_path = '00_peptide_data_20240702/unique_peptide_effects_20240702.xlsx'
xls = pd.ExcelFile(peptide_file_path)

# Load the gene expression data
gene_expression_file_path = '01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv'
gene_expression_df = pd.read_csv(gene_expression_file_path)

# Initialize a dictionary to store tables for each tissue type
tissue_tables = {}

# Process each sheet in the Excel file (each sheet corresponds to a tissue type)
for sheet_name in xls.sheet_names:
    # Load the data from the current sheet
    df = pd.read_excel(xls, sheet_name=sheet_name)
    
    # Filter columns that contain the word "effect"
    effect_columns = [col for col in df.columns if 'effect' in col.lower()]
    padj_columns = [col for col in df.columns if 'padj' in col.lower()]
    
    # Skip processing if no 'effect' columns are found
    if not effect_columns:
        print(f"No columns with 'effect' found in '{sheet_name}' sheet. Skipping...")
        continue
    
    # Merge with gene expression data to add gene expression change
    merge_left_on = 'symbol'  # Assuming 'symbol' is the gene identifier column in df
    merge_right_on = 'gene_name'  # Assuming 'gene_name' is the gene identifier column in gene_expression_df
    
    merged_df = pd.merge(df, gene_expression_df, left_on=merge_left_on, right_on=merge_right_on, how='left')
    
    # Determine change in gene expression if 'count' column is present
    if 'count' in merged_df.columns:
        # Calculate gene change based on all isoforms of the same gene
        def calculate_gene_change(x):
            if (x > 0).all():
                return '+'
            elif (x < 0).all():
                return '-'
            else:
                return 'same'
        
        gene_change = merged_df.groupby('gene_name')['count'].apply(calculate_gene_change).reset_index()
        gene_change.columns = ['gene_name', 'gene_change']
        
        # Merge gene_change back to merged_df
        merged_df = pd.merge(merged_df, gene_change, on='gene_name', how='left')
    else:
        print(f"Error: 'count' column not found in merged_df for '{sheet_name}'")
        continue
    
    # Calculate isoform change based on peptide-level data
    def calculate_isoform_change(row):
        # Extract peptide-level data
        peptides = row[['sn_pep31811_pep_07022024', 'sn_pep31812_pep_07022024', 'sn_pep31813_pep_07022024', 'sn_pep31814_pep__07022024']]
        
        # Count unique isoforms
        unique_isoforms = set()
        for peptide_list in peptides:
            if isinstance(peptide_list, str):  # Check if peptide_list is a string
                isoforms = peptide_list.split(', ')
                unique_isoforms.update(isoforms)
        
        # Determine isoform change
        if len(unique_isoforms) == 1:
            return 'same'
        else:
            return 'altered'  # Change to 'altered' if different isoforms are detected
    
    merged_df['isoform_change'] = merged_df.apply(calculate_isoform_change, axis=1)
    
    # Prepare summary tables for each effect column
    sheet_tables = {}
    for effect_column in effect_columns:
        if effect_column in merged_df.columns:
            effect_df = merged_df.copy()
            
            # Filter based on absolute value of the current effect column
            effect_df = effect_df[(effect_df[effect_column].abs() >= 0.5) & (effect_df[effect_column].abs() <= 1)]
            
            # Check if necessary columns are present
            required_columns = ['symbol', 'pep_seq', 'gene_change', 'isoform_change', effect_column]
            # Include padj column if it exists
            padj_column = [col for col in padj_columns if effect_column in col]
            if padj_column:
                required_columns += padj_column
            
            if all(col in effect_df.columns for col in required_columns):
                effect_df = effect_df[required_columns].rename(columns={'symbol': 'gene', 'pep_seq': 'peptide'}).drop_duplicates()
                sheet_tables[effect_column] = effect_df
            else:
                print(f"Error: Required columns {required_columns} not found in df for '{sheet_name}'")
        else:
            print(f"Error: '{effect_column}' column not found in df for '{sheet_name}'")
    
    # Store tables for the current tissue type
    tissue_tables[sheet_name] = sheet_tables

# Output the results for each tissue type as a single Excel file with each effect type as a separate sheet/tab
for tissue_name, sheet_tables in tissue_tables.items():
    output_file_path = f'{tissue_name}_summary.xlsx'
    with pd.ExcelWriter(output_file_path) as writer:
        for effect_name, table_df in sheet_tables.items():
            table_df.to_excel(writer, sheet_name=effect_name, index=False)
            print(f"Sheet '{effect_name}' for '{tissue_name}' has been saved.")
    print(f"All sheets for '{tissue_name}' have been saved to '{output_file_path}'")

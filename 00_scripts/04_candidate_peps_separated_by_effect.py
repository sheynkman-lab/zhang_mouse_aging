import pandas as pd

# Load the data from both CSV files
peptide_isoform_df = pd.read_csv('01_filter_data/01.4_isoform_annotation_of_experi_peptides/peptide_to_isoform_mapping_with_experimental_peptides.csv')
peptide_effects_df = pd.read_excel('00_peptide_data_20240702/unique_peptide_effects_20240702.xlsx', sheet_name=None)

# Initialize an empty DataFrame to store combined effects from all tissues
combined_effects_df = pd.DataFrame()

# Process each tissue tab
for tissue, df in peptide_effects_df.items():
    df = df.copy()
    df['tissue'] = tissue  # Add a tissue column
    combined_effects_df = pd.concat([combined_effects_df, df], ignore_index=True)

# Ensure the correct column name is used for merging
peptide_seq_column = 'peptide_seq' if 'peptide_seq' in combined_effects_df.columns else 'peptide_id'

# Standardize the peptide sequences in combined_effects_df
combined_effects_df['standardized_peptide_seq'] = combined_effects_df[peptide_seq_column].apply(lambda x: x.split('.')[1] if '.' in x else x)

# Merge the datasets on the standardized peptide sequence
merged_df = pd.merge(peptide_isoform_df, combined_effects_df, left_on='pep_seq', right_on='standardized_peptide_seq', how='inner')

# Define a function to determine if the change is significant based on p-values and effect sizes
def map_change(p_value, effect_size):
    try:
        p_value = float(p_value)
        effect_size = float(effect_size)
    except ValueError:
        return '0'
    
    if p_value < 0.05:
        if effect_size > 0:
            return '1'
        elif effect_size < 0:
            return '-1'
        else:
            return '0'
    else:
        return '0'

# Apply the change mapping to the dataset
merged_df['sex_change'] = merged_df.apply(lambda row: map_change(row['sex_pval'], row['sex_effect']), axis=1)
merged_df['age_change'] = merged_df.apply(lambda row: map_change(row['age_pval'], row['age_effect']), axis=1)
merged_df['age_intercept_change'] = merged_df.apply(lambda row: map_change(row['age_cat_pval'], row['age_intercept_effect']), axis=1)
merged_df['age_m12_change'] = merged_df.apply(lambda row: map_change(row['age_cat_pval'], row['age_m12_effect']), axis=1)
merged_df['age_m20_change'] = merged_df.apply(lambda row: map_change(row['age_cat_pval'], row['age_m20_effect']), axis=1)

# Aggregate peptides for each gene
def aggregate_data(group):
    result = []

    # Group isoform-specific peptides together
    isoform_specific = group[group['category'] == 'isoform-specific']
    if not isoform_specific.empty:
        peptides = isoform_specific['pep_seq'].tolist()
        changes = isoform_specific[['sex_change', 'age_change', 'age_intercept_change', 'age_m12_change', 'age_m20_change']].astype(str).apply(lambda x: '; '.join(x), axis=1).tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate average fraction of significant changes
            fraction_gene_change = sum(1 for v in '; '.join(changes).split('; ') if v != '0') / num_peptides

            result.append({
                'gene': isoform_specific['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'isoform-specific',
                'mapped_isoform': '; '.join(isoform_specific['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(isoform_specific.apply(lambda x: '; '.join(map(str, [x['sex_effect'], x['age_effect'], x['age_intercept_effect'], x['age_m12_effect'], x['age_m20_effect']])), axis=1).tolist()),
                'fraction_gene_change': fraction_gene_change
            })

    # Group isoform-informative peptides together
    isoform_informative = group[group['category'] == 'isoform-informative']
    if not isoform_informative.empty:
        peptides = isoform_informative['pep_seq'].tolist()
        changes = isoform_informative[['sex_change', 'age_change', 'age_intercept_change', 'age_m12_change', 'age_m20_change']].astype(str).apply(lambda x: '; '.join(x), axis=1).tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate average fraction of significant changes
            fraction_gene_change = sum(1 for v in '; '.join(changes).split('; ') if v != '0') / num_peptides

            result.append({
                'gene': isoform_informative['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'isoform-informative',
                'mapped_isoform': '; '.join(isoform_informative['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(isoform_informative.apply(lambda x: '; '.join(map(str, [x['sex_effect'], x['age_effect'], x['age_intercept_effect'], x['age_m12_effect'], x['age_m20_effect']])), axis=1).tolist()),
                'fraction_gene_change': fraction_gene_change
            })

    # Group constitutive peptides together
    constitutive = group[group['category'] == 'constitutive']
    if not constitutive.empty:
        peptides = constitutive['pep_seq'].tolist()
        changes = constitutive[['sex_change', 'age_change', 'age_intercept_change', 'age_m12_change', 'age_m20_change']].astype(str).apply(lambda x: '; '.join(x), axis=1).tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate average fraction of significant changes
            fraction_gene_change = sum(1 for v in '; '.join(changes).split('; ') if v != '0') / num_peptides

            result.append({
                'gene': constitutive['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'constitutive',
                'mapped_isoform': '; '.join(constitutive['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(constitutive.apply(lambda x: '; '.join(map(str, [x['sex_effect'], x['age_effect'], x['age_intercept_effect'], x['age_m12_effect'], x['age_m20_effect']])), axis=1).tolist()),
                'fraction_gene_change': fraction_gene_change
            })

    return pd.DataFrame(result)

# Aggregate data by gene
aggregated_data = merged_df.groupby('gene_name', as_index=False, group_keys=False).apply(aggregate_data).reset_index(drop=True)

# Create Excel writer
writer = pd.ExcelWriter('aggregated_peptide_data_final.xlsx', engine='xlsxwriter')

# Create separate sheets for each effect
for effect in ['sex', 'age', 'age_intercept', 'age_m12', 'age_m20']:
    # Filter data for the specific effect
    effect_df = aggregated_data[aggregated_data['effect_sizes'].str.contains(effect)]

    # Write to Excel
    effect_df.to_excel(writer, sheet_name=f'{effect}_effect', index=False)

# Save and close the Excel file
writer._save()
writer.close()

print("Final aggregated peptide data saved to 'aggregated_peptide_data_final.xlsx'.")

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

# Print the merged dataframe shape and a sample
print("Merged DataFrame shape:", merged_df.shape)
print("Merged DataFrame sample:")
print(merged_df.head())

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
merged_df['change'] = merged_df.apply(
    lambda row: '; '.join(map_change(row[p], row[e]) for p, e in [
        ('sex_pval', 'sex_effect'),
        ('age_pval', 'age_effect'),
        ('age_cat_pval', 'age_intercept_effect'),
        ('age_cat_pval', 'age_m12_effect'),
        ('age_cat_pval', 'age_m20_effect')
    ]), axis=1
)

# Print the changes column for debugging
print("Changes column sample:")
print(merged_df[['pep_seq', 'change']].head())

# Aggregate peptides for each gene
# Define a function to aggregate data for each gene
def aggregate_data(group):
    result = []

    # Group isoform-specific peptides together
    isoform_specific = group[group['category'] == 'isoform-specific']
    if not isoform_specific.empty:
        peptides = isoform_specific['pep_seq'].tolist()
        changes = isoform_specific['change'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in '; '.join(changes).split('; ') if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides

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
        changes = isoform_informative['change'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in '; '.join(changes).split('; ') if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides

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
        changes = constitutive['change'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in '; '.join(changes).split('; ') if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides

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
aggregated_data = merged_df.groupby('gene_name').apply(aggregate_data).reset_index(drop=True)

# Print the aggregated data sample
print("Aggregated data sample:")
print(aggregated_data.head())

# Define the filter condition function
def filter_condition(group):
    shared_peptides = group[group['peptide_cat'].str.contains('constitutive|semishared')]
    isoform_specific_peptides = group[group['peptide_cat'].str.contains('isoform-specific')]

    if not shared_peptides.empty and not isoform_specific_peptides.empty:
        shared_no_change_fraction = (shared_peptides['change'].str.count('0') / shared_peptides['change'].str.split('; ').str.len()).mean()
        isoform_specific_change_present = isoform_specific_peptides['change'].str.contains('1|\\-1').any()
        
        print(f"Gene: {group['gene'].iloc[0]}")
        print(f"Shared no change fraction: {shared_no_change_fraction}")
        print(f"Isoform-specific change present: {isoform_specific_change_present}")

        return shared_no_change_fraction > 0.8 and isoform_specific_change_present
    return False

# Filter the aggregated data
filtered_data = aggregated_data.groupby('gene').filter(filter_condition)

# Print some debug information
print("Filtered data sample:")
print(filtered_data.head())

# Save the result to a new CSV file
filtered_data.to_csv('filtered_peptide_data.csv', index=False)

print("Filtered peptide data saved to 'filtered_peptide_data.csv'.")

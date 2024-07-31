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

# Apply the change mapping to the dataset for all effects
merged_df['change'] = merged_df.apply(
    lambda row: map_change(row['sex_pval'], row['sex_effect']) if 'sex_pval' in row else '0', axis=1
)

# Define a function to aggregate data by peptide category
def aggregate_data(group):
    result = []

    # Group isoform-specific peptides together
    isoform_specific = group[group['category'] == 'isoform-specific']
    if not isoform_specific.empty:
        peptides = isoform_specific['pep_seq'].tolist()
        changes = isoform_specific['change'].tolist()
        effect_sizes = isoform_specific['sex_effect'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in changes if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides
            avg_effect_size = sum(float(e) for e in effect_sizes) / num_peptides

            result.append({
                'gene': isoform_specific['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'isoform-specific',
                'mapped_isoform': '; '.join(isoform_specific['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(map(str, effect_sizes)),
                'fraction_gene_change': fraction_gene_change,
                'avg_effect_size': avg_effect_size
            })

    # Group isoform-informative peptides together
    isoform_informative = group[group['category'] == 'isoform-informative']
    if not isoform_informative.empty:
        peptides = isoform_informative['pep_seq'].tolist()
        changes = isoform_informative['change'].tolist()
        effect_sizes = isoform_informative['sex_effect'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in changes if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides
            avg_effect_size = sum(float(e) for e in effect_sizes) / num_peptides

            result.append({
                'gene': isoform_informative['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'isoform-informative',
                'mapped_isoform': '; '.join(isoform_informative['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(map(str, effect_sizes)),
                'fraction_gene_change': fraction_gene_change,
                'avg_effect_size': avg_effect_size
            })

    # Group constitutive peptides together
    constitutive = group[group['category'] == 'constitutive']
    if not constitutive.empty:
        peptides = constitutive['pep_seq'].tolist()
        changes = constitutive['change'].tolist()
        effect_sizes = constitutive['sex_effect'].tolist()
        num_peptides = len(peptides)

        if num_peptides > 0:
            # Calculate fraction of significant changes
            non_zero_changes = sum(1 for v in changes if v != '0')
            fraction_gene_change = non_zero_changes / num_peptides
            avg_effect_size = sum(float(e) for e in effect_sizes) / num_peptides

            result.append({
                'gene': constitutive['gene_name'].iloc[0],
                'peptides': '; '.join(peptides),
                'peptide_cat': 'constitutive',
                'mapped_isoform': '; '.join(constitutive['transcript_name'].tolist()),
                'change': '; '.join(changes),
                'effect_sizes': '; '.join(map(str, effect_sizes)),
                'fraction_gene_change': fraction_gene_change,
                'avg_effect_size': avg_effect_size
            })

    return pd.DataFrame(result)

# Aggregate data by gene
aggregated_data = merged_df.groupby('gene_name').apply(aggregate_data).reset_index(drop=True)

# Filter the aggregated data to keep only genes with both isoform-specific and constitutive peptides
isoform_specific_genes = aggregated_data[aggregated_data['peptide_cat'] == 'isoform-specific']['gene'].unique()
constitutive_genes = aggregated_data[aggregated_data['peptide_cat'] == 'constitutive']['gene'].unique()
genes_with_both = set(isoform_specific_genes) & set(constitutive_genes)

# Define the filter condition function
def filter_condition(group):
    isoform_specific_peptides = group[group['peptide_cat'] == 'isoform-specific']
    constitutive_peptides = group[group['peptide_cat'] == 'constitutive']

    if not isoform_specific_peptides.empty and not constitutive_peptides.empty:
        isoform_specific_changes = isoform_specific_peptides['change'].str.split('; ').apply(lambda x: set(x))
        if all(len(changes) == 1 for changes in isoform_specific_changes):
            constitutive_no_change_fraction = constitutive_peptides['change'].str.split('; ').apply(lambda x: x.count('0') / len(x)).mean()
            return constitutive_no_change_fraction >= 0.95

    return False

# Filter the aggregated data to include only genes with both isoform-specific and constitutive peptides, and meet the condition
filtered_data = aggregated_data[
    (aggregated_data['gene'].isin(genes_with_both)) & (
        (aggregated_data['peptide_cat'] == 'isoform-specific') & 
        (aggregated_data['fraction_gene_change'] != 0) |
        (aggregated_data['peptide_cat'] == 'constitutive') &
        (aggregated_data['fraction_gene_change'] == 0)
    )
]

# Apply additional filter condition
final_filtered_data = filtered_data.groupby('gene').filter(filter_condition)

# Include isoform-informative peptides in the final data
final_filtered_data_with_informative = pd.concat([
    final_filtered_data,
    aggregated_data[(aggregated_data['gene'].isin(final_filtered_data['gene'].unique())) & (aggregated_data['peptide_cat'] == 'isoform-informative')]
])

# Save the result to a new CSV file
final_filtered_data_with_informative.to_csv('./04_candidate_peps/subset_filtered_peptides_sex_effects.csv', index=False)

print("Filtered peptide data saved to 'subset_filtered_peptides_sex_effects.csv'.")

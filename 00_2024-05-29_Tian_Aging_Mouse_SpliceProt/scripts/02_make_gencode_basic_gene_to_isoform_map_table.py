# Tian shared 4 MS runs representing peptides identified from mouse tissue samples.
# The goal is to analyze the peptides that are specific to each isoform.


# First, need to obtain the mouse basic GENCODE protein fasta file.
# A translations fasta for comprehensive GENCODE is available only.
# Extract basic GENCODE gene/isoforms from GTF.

# NOTE - this script takes ~1 hour to run - could be optimized

#%%

import gffutils
import pandas as pd
import os

#%%


# gencode mouse gtf file - full
gtf_file = '../gencode_mouse_models/gencode.vM35.basic.annotation.gtf'

# Create a gffutils database
db = gffutils.create_db(gtf_file, dbfn=':memory:', force=True, keep_order=True, 
                        merge_strategy='merge', sort_attribute_values=True)

data = [] # will hold accession groupings

for transcript in db.features_of_type('transcript'):
    if 'gene_type' in transcript.attributes and 'protein_coding' in transcript.attributes['gene_type']:
        gene_name = transcript.attributes.get('gene_name', [''])[0]
        transcript_name = transcript.attributes.get('transcript_name', [''])[0]
        transcript_id = transcript.attributes.get('transcript_id', [''])[0]
        gene_id = transcript.attributes.get('gene_id', [''])[0]
        protein_id = transcript.attributes.get('protein_id', [''])[0]
        gencode_basic = 'basic' in transcript.attributes.get('tag', [])
        data.append((gene_name, transcript_name, transcript_id, gene_id, protein_id, gencode_basic))
        

# Convert to DataFrame
df = pd.DataFrame(data, columns=['gene_name', 'transcript_name', 'transcript_id', 'gene_id', 'protein_id', 'gencode_basic'])

# Write out to tsv file in output directory
def write_dataframe(directory, filename, df):
    os.makedirs(directory, exist_ok=True)
    file_path = os.path.join(directory, filename)
    df.to_csv(file_path, sep='\t', index=False)

write_dataframe('../output/02_make_gencode_basic_gene_to_isoform_map_table', 'gencode_basic_gene_to_isoform_map_table.tsv', df)
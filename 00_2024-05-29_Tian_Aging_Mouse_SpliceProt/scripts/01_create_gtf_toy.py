#%%
import os

# Create a small mouse gtf to test read in of gtf with gffutils

def create_small_gtf(input_gtf, output_gtf, num_genes=10):
    gene_count = 0
    gene_ids = set()
    
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if fields[2] == 'transcript':
                attributes = {key: value.strip('"') for key, value in 
                              (attr.split() for attr in fields[8].split('; ') if attr)}
                
                if 'gene_type' in attributes and attributes['gene_type'] == 'protein_coding':
                    gene_id = attributes.get('gene_id')
                    if gene_id not in gene_ids:
                        gene_ids.add(gene_id)
                        gene_count += 1
                        if gene_count > num_genes:
                            break
                    outfile.write(line)
                else:
                    if attributes.get('gene_id') in gene_ids:
                        outfile.write(line)

# Replace 'path_to_your_gtf_file.gtf' with the actual path to your original GTF file
input_gtf = '../gencode_mouse_models/gencode.vM35.basic.annotation.gtf'

os.makedirs('../output/01_create_gtf_toy', exist_ok=True)

output_gtf = '../output/01_create_gtf_toy/gencode.vM35.basic.annotation.toy.gtf'
create_small_gtf(input_gtf, output_gtf, num_genes=10)


import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Extract ENSG to gene name map from a GTF file.')
parser.add_argument('--gtf_file', type=str, required=True, help='Path to the Ensembl GTF file')
parser.add_argument('--output_file', type=str, default='ensg_to_gene_name_map.txt', help='Output file for the ENSG to gene name map')

args = parser.parse_args()

# Assign arguments to variables
gtf_file_path = args.gtf_file
output_file_path = args.output_file

# Initialize a dictionary to store the mapping of ENSG to gene names
ensg_to_genename = {}

# Open and read the GTF file
with open(gtf_file_path, "r") as gtf_file:
    for line in gtf_file:
        # Skip header lines (lines that start with "#")
        if line.startswith("#"):
            continue
        
        # Split the GTF line into fields
        fields = line.strip().split("\t")
        
        # Check if this line represents a gene (feature type is "gene")
        if fields[2] == "gene":
            # Extract the gene ID (ENSG) and gene name from the attributes field
            attributes = fields[8].split(";")
            for attr in attributes:
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    ensg = attr.split('"')[1]
                elif attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
            
            # Store the mapping in the dictionary
            ensg_to_genename[ensg] = gene_name

# Write out the ENSG to gene name map
with open(output_file_path, 'w') as ofile:
    for ensg, gene_name in ensg_to_genename.items():
        ofile.write(f"{ensg}\t{gene_name}\n")

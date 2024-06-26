#!/bin/bash

#SBATCH --job-name=six_frame
#SBATCH --cpus-per-task=10 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --account=sheynkman_lab
#SBATCH --output=%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yqy3cu@virginia.edu

# Load necessary modules (if needed)
module purge
module load gcc/11.4.0
module load mamba/22.11.1-4
module load bioconda/py3.10
module load anaconda/2023.07-py3.11
module load openmpi/4.1.4
module load python/3.11.4

# It is assumed that your GTF file and FASTA file, along with ensg_to_gene_name_map.txt are in a folder called "data"
  # If you need to generate the ensg_to_gene_name_map.txt file, see the script extract_ensg_to_gene_map_gencode_gtf.py
  # Update the name of your GTF and FASTA files in the script map_peptides_to_pacbio_database.py 

# First, make PoGo compatible input GTF/FASTA files

conda activate PoGo

python map_peptides_to_pacbio_database.py

conda deactivate
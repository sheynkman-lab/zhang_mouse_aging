#!/bin/bash

#SBATCH --job-name=r_scripts
#SBATCH --cpus-per-task=20 # number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00 # amount of time for the whole job
#SBATCH --partition=standard # the queue/partition to run on
#SBATCH --account=sheynkman_lab
#SBATCH --output=%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yqy3cu@virginia.edu
#SBATCH --mem=800G # memory per node 

# Load necessary modules (if needed)
module purge 

module load gcc/11.4.0
module load mamba/22.11.1-4
module load bioconda/py3.10
module load anaconda/2023.07-py3.11
module load openmpi/4.1.4
module load python/3.11.4
module load git-lfs/2.10.0
module load apptainer/1.2.2
module load R/4.3.1

cd /project/sheynkman/projects/zhang_mouse_aging

Rscript ./2024-05-29_Tian_Aging_Mouse_SpliceProt/scripts/06_isoform_annotation_of_experi_peptides.R
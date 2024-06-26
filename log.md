# Step 1 - create directory on Rivanna & load modules
I created a directory on Rivanna called `zhang_mouse_aging` that corresponds to [this GitHub repository](https://github.com/sheynkman-lab/zhang_mouse_aging). <br/>
This protocol follows the ones set up in [this GitHub repository](https://github.com/efwatts/PoGo2GenomeBrowser) with some minor changes. The environment I loaded was created using that module <br/>
```
cd /project/sheynkman/projects
git clone https://github.com/sheynkman-lab/zhang_mouse_aging.git

cd zhang_mouse_aging

module load gcc/11.4.0
module load mamba/22.11.1-4
module load bioconda/py3.10
module load anaconda/2023.07-py3.11
module load openmpi/4.1.4
module load python/3.11.4
module load git-lfs/2.10.0

conda env create -f PoGo.yml
conda activate PoGo
```

# Step 2 - convert GTF and FASTA files to be PoGo compatible 
First, we need an ensg to gene name map from the Gencode database, then we can convert the GTF and FASTA files to be PoGo compatible. <br/>
We need a GTF file with CDS and a FASTA file with ORF information to create the PoGo database. <br/>
```
python ./GTF_FASTA_2PoGo/extract_ensg_to_gene_map_gencode_gtf.py --gtf_file ./gencode_mouse_models/gencode.vM35.annotation.gtf --output_file ./GTF_FASTA_2PoGo/ensg_to_gene_name_map.txt

python ./GTF_FASTA_2PoGo/map_peptides_to_pacbio_database.py --pacbio_fasta ./gencode_mouse_models/gencode.vM35.pc_translations.fa --pacbio_gtf ./gencode_mouse_models/gencode.vM35.annotation.gtf --ensg_to_gene ./GTF_FASTA_2PoGo/ensg_to_gene_name_map.txt --output_dir ./GTF_FASTA_2PoGo

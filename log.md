# Step 0 - create directory on Rivanna & load modules
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
module load apptainer/1.2.2

conda env create -f PoGo.yml
conda activate PoGo
```

When coming back to this script, move to the directory and activate the environment
```
cd /project/sheynkman/projects/zhang_mouse_aging
conda activate PoGo
```


# Step 1 - Peptides 2 PoGo
This is an R script that takes the peptides and converts them to a PoGo database. <br/>
I just ran it through RStudio on my local machine. <br/>
./01_Peptides2PoGo/ProteinDF2PoGo.R

# Step 1 - LRP 
Reference Tables
```
conda activate reference_tab

apptainer pull docker://gsheynkmanlab/generate-reference-tables:latest

apptainer exec generate-reference-tables_latest.sif /bin/bash -c "\
    python ./LRP/00_scripts/01_prepare_reference_tables.py \
        --gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
        --fa ./00_gencode_mouse_models/gencode.vM35.pc_transcripts.fa \
        --ensg_gene ./LRP/01_reference_tables/ensg_gene.tsv \
        --enst_isoname ./LRP/01_reference_tables/enst_isoname.tsv \
        --gene_ensp ./LRP/01_reference_tables/gene_ensp.tsv \
        --gene_isoname ./LRP/01_reference_tables/gene_isoname.tsv \
        --isoname_lens ./LRP/01_reference_tables/isoname_lens.tsv \
        --gene_lens ./LRP/01_reference_tables/gene_lens.tsv \
        --protein_coding_genes ./LRP/01_reference_tables/protein_coding_genes.txt
"
conda deactivate
```

SQANTI
```
conda activate SQANTI3.env

chmod +x /project/sheynkman/programs-needs_attentionEFW/SQANTI3-5.2/utilities/gtfToGenePred
export PYTHONPATH=$PYTHONPATH:/project/sheynkman/programs-needs_attentionEFW/SQANTI3-5.2/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/project/sheynkman/programs-needs_attentionEFW/SQANTI3-5.2/cDNA_Cupcake/

python /project/sheynkman/programs-needs_attentionEFW/SQANTI3-5.2/sqanti3_qc.py \
./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
./00_gencode_mouse_models/GRCm39.primary_assembly.genome.fa \
--skipORF \
-o mouse \
-d ./LRP/02_sqanti/ 

conda deactivate
```

CPAT
```
export PATH="$HOME/.local/bin:$PATH"

cpat \
   -x ./LRP/00_input_data/Mouse_Hexamer.tsv \
   -d ./LRP/00_input_data/Mouse_logitModel.RData \
   -g ./00_gencode_mouse_models/gencode.vM35.pc_transcripts.fa \
   --min-orf=50 \
   --top-orf=50 \
   -o mouse_CPAT \
   2> cpat.error
```

Transcriptome Summary 
```
conda activate transcriptome_sum

python ./LRP/00_scripts/04_transcriptome_summary_gene_table_only.py \
--sq_out ./LRP/02_sqanti/mouse_classification.txt \
--ensg_to_gene ./LRP/01_reference_tables/ensg_gene.tsv \
--enst_to_isoname ./LRP/01_reference_tables/enst_isoname.tsv \
--odir ./LRP/04_transcriptome_summary/

conda deactivate
```

ORF-Calling
```
apptainer pull docker://gsheynkmanlab/orf_calling:latest

conda activate orf-calling

apptainer exec orf_calling_latest.sif /bin/bash -c "\
   python ./LRP/00_scripts/05_orf_calling.py \
      --orf_coord ./LRP/04_CPAT/mouse_CPAT.ORF_prob.tsv \
      --orf_fasta ./LRP/04_CPAT/mouse_CPAT.ORF_seqs.fa \
      --gencode ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
      --sample_gtf ./LRP/02_sqanti/mouse_corrected.gtf \
      --pb_gene ./LRP/04_transcriptome_summary/pb_gene.tsv \
      --classification ./LRP/02_sqanti/mouse_classification.txt \
      --sample_fasta ./LRP/02_sqanti/mouse_corrected.fasta \
      --output ./LRP/05_orf_calling/mouse_best_ORF.tsv
"
conda deactivate
```

Refine ORF database
```
conda activate refined-database-generation

python ./LRP/00_scripts/06_refine_orf_database.py \
--name ./LRP/06_refine_orf_database/mouse \
--orfs ./LRP/05_orf_calling/mouse_best_ORF.tsv \
--pb_fasta ./LRP/02_sqanti/mouse_corrected.fasta \
--coding_score_cutoff 0.3 

conda deactivate
```

Make CDS GTF
```
conda activate reference_tab

apptainer exec pb-cds-gtf_latest.sif /bin/bash -c "\
      python ./LRP/00_scripts/07_make_pacbio_cds_gtf.py \
      --sample_gtf ./LRP/02_sqanti/mouse_corrected.gtf \
      --agg_orfs ./LRP/06_refine_orf_database/mouse_orf_refined.tsv \
      --refined_orfs ./LRP/05_orf_calling/mouse_best_ORF.tsv \
      --pb_gene ./LRP/04_transcriptome_summary/pb_gene.tsv \
      --output_cds ./LRP/07_make_cds_gtf/mouse_cds.gtf
"
conda deactivate
```
# Step 2 - convert GTF and FASTA files to be PoGo compatible 
First, we need an ensg to gene name map from the Gencode database, then we can convert the GTF and FASTA files to be PoGo compatible. <br/>
We need a GTF file with CDS and a FASTA file with ORF information to create the PoGo database. <br/>
Ok, it actually looks like the gencode files may work with Pogo.
```
conda activate PoGo 

python ./02_GTF_FASTA_2PoGo/extract_ensg_to_gene_map_gencode_gtf.py --gtf_file ./00_gencode_mouse_models/gencode.vM35.annotation.gtf --output_file ./02_GTF_FASTA_2PoGo/ensg_to_gene_name_map.txt

python ./02_GTF_FASTA_2PoGo/map_peptides_to_pacbio_database.py --pacbio_fasta ./LRP/06_refine_orf_database/mouse_orf_refined.fasta --pacbio_gtf ./LRP/07_make_cds_gtf/mouse_cds.gtf --ensg_to_gene ./02_GTF_FASTA_2PoGo/ensg_to_gene_name_map.txt --output_dir ./02_GTF_FASTA_2PoGo
```


# Step 3 - Run PoGo
Be sure PoGo is in your PATH. <br/>
So, this works, but I need to run it through some parts of the LRP to get the right files types. for fasta and gtf pogo.
```
export PATH=$PATH:/project/sheynkman/programs-needs_attentionEFW/PoGo 

PoGo -fasta ./00_gencode_mouse_models/human/gencode.v46.pc_translations.fa -gtf ./00_gencode_mouse_models/human/gencode.v46.annotation.gtf -in all_peptide_POGO.txt -format BED -species MOUSE 

PoGo -fasta ./00_gencode_mouse_models/gencode.vM34.pc_translations.fa -gtf ./00_gencode_mouse_models/gencode.vM34.annotation.gtf -in all_peptide_POGO.txt -format BED

PoGo -fasta ./00_gencode_mouse_models/gencode.vM32.pc_translations.fa -gtf ./00_gencode_mouse_models/gencode.vM32.basic.annotation.gtf -in all_peptide_POGO.txt -format BED


```

Trying to use the LRP scripts to make the bed files for genome browser
Multiregion BED generation
```
conda activate visualization

python ./00_scripts/18_make_region_bed_for_ucsc.py \
--name mouse_multi \
--sample_gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
--reference_gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
--output_dir ./LRP_tracks/multiregion_bed
```
Peptide Track Visualization
```
python ./00_scripts/18_make_peptide_gtf_file.py \
--name ./LRP_tracks/mouse_peptides \
--sample_gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
--reference_gtf ./00_input_data/gencode.v35.annotation.canonical.gtf \
--peptides ./01_Peptides2Pogo/all_peptide_POGO.txt \
--pb_gene ./LRP/04_transcriptome_summary/pb_gene.tsv \
--gene_isoname ./LRP/01_reference_tables/gene_isoname.tsv \
--refined_fasta ./00_gencode_mouse_models/gencode.vM35.pc_transcripts.fa 

gtfToGenePred ./18_track_visualization/peptide/jurkat_hybrid_peptides.gtf ./18_track_visualization/peptide/jurkat_hybrid_peptides.genePred
genePredToBed ./18_track_visualization/peptide/jurkat_hybrid_peptides.genePred ./18_track_visualization/peptide/jurkat_hybrid_peptides.bed12
# add rgb to colorize specific peptides
python ./00_scripts/18_finalize_peptide_bed.py \
--bed ./18_track_visualization/peptide/jurkat_hybrid_peptides.bed12 \
--name ./18_track_visualization/jurkat_hybrid

conda deactivate
```

This is actually a new version of the script...
```
export PATH=$PATH:/project/sheynkman/programs-needs_attentionEFW/PoGo
python ./LRP_tracks/peptide_gtf.py --name mouse --gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf --peptides all_peptide_POGO.txt --fasta ./00_gencode_mouse_models/gencode.vM35.pc_translations.fa
```

# Hyrbid
python ./LRP_tracks/peptide_to_bed.py \
--name ./LRP_tracks/mouse \
--sample_gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
--reference_gtf ./00_gencode_mouse_models/gencode.vM35.annotation.gtf \
--peptides all_peptide_POGO.txt \
--pb_gene ./LRP/04_transcriptome_summary/pb_gene.tsv \
--gene_isoname ./LRP/01_reference_tables/gene_isoname.tsv \
--refined_fasta ./00_gencode_mouse_models/gencode.vM35.pc_translations.fa 

gtfToGenePred ./18_track_visualization/peptide/jurkat_hybrid_peptides.gtf ./18_track_visualization/peptide/jurkat_hybrid_peptides.genePred
genePredToBed ./18_track_visualization/peptide/jurkat_hybrid_peptides.genePred ./18_track_visualization/peptide/jurkat_hybrid_peptides.bed12
# add rgb to colorize specific peptides
python ./00_scripts/18_finalize_peptide_bed.py \
--bed ./18_track_visualization/peptide/jurkat_hybrid_peptides.bed12 \
--name ./18_track_visualization/jurkat_hybrid

conda deactivate

# Trying a new approach
## Step 1 - Download newest version of PoGo & data
```
wget https://github.com/cschlaffner/PoGo/releases/download/v1.2.3/PoGo_v1.2.3.zip
unzip PoGo_v1.2.3.zip
export PATH=$PATH:/project/sheynkman/projects/zhang_mouse_aging/PoGo_v1.2.3/Linux
```
This version of PoGo will work with mouse data, but it does not support Gencode. It only supports Ensambl. So I'm downloading the new dataset.

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in all_peptide_POGO.txt -format BED
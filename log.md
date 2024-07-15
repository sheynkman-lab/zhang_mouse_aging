# Step 0 - load modules and set working data
```
cd /project/sheynkman/projects/zhang_mouse_aging

module load gcc/11.4.0
module load mamba/22.11.1-4
module load bioconda/py3.10
module load anaconda/2023.07-py3.11
module load openmpi/4.1.4
module load python/3.11.4
module load git-lfs/2.10.0
module load apptainer/1.2.2
module load R/4.3.1
```

# Step 1 - sort through peptide data
I used the following scripts in R to create these:
`01.1_view_peptide_data.R`
`01.2_make_clean_peptide_tables.R`
`01.3_peptide_to_isoform_mapping.R`
`01.4_isoform_annotation_of_experi_peptides.R`
`01.5_filter_for_co_expressed_isoforms.R`

# Step 2 - prepare peptide data for PoGo
I used this script:
`02_ProteinDF2PoGo.R`

# Step 3 - Run PoGo
```
export PATH=$PATH:/project/sheynkman/programs-needs_attentionEFW/PoGo_v1.2.3/Linux

PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/coexpressed_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/experimental_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31811_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31812_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31813_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31814_peptides.txt -format BED
```

# Step 3 - Convert to bigBED
```
module load perl/5.36.0

perl 00_scripts/TrackHubGenerator.pl /project/sheynkman/projects/zhang_mouse_aging/ mm39 /project/sheynkman/projects/zhang_mouse_aging/03_pogo_out/ /project/sheynkman/programs-needs_attentionEFW/ watts.emily.f@virginia.edu
```

# Step 4 - Find candidate peptides
```
conda env create -f ./00_scripts/candidate_peps.yml

conda activate candidate_peps

python ./00_scripts/04_candidate_peps_sig.py
python ./00_scripts/04_candidate_peps_all.py
python ./00_scripts/04_separate_candidate_peps.py
python ./00_scripts/04_candidate_peps_summary.py
python ./00_scripts/04_candidate_peps_summary_reduced.py
python ./00_scripts/04_candidate_peps_separated_by_effect.py
```
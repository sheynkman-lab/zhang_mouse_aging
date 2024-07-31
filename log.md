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
export PATH=$PATH:/project/sheynkman/programs/PoGo_v1.2.3/Linux

PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/coexpressed_isoform_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/all_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31811_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31812_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31813_peptides.txt -format BED
PoGo -fasta ./00_ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./00_gencode_mouse_models/gencode.vM35.basic.annotation.gtf -in ./02_Peptides2Pogo/sn_pq31814_peptides.txt -format BED
```

# Step 3 - Convert to bigBED
```
module load perl/5.36.0

perl 00_scripts/TrackHubGenerator.pl /project/sheynkman/projects/zhang_mouse_aging/ mm39 /project/sheynkman/projects/zhang_mouse_aging/03_pogo_out/ /project/sheynkman/programs/ watts.emily.f@virginia.edu
```

# Step 4 - Find candidate peptides
All 3 scripts apply the following steps:
1. Filter for peptides with a significant (>0.05) change
2. Determine "fraction gene change." 0 if no [significant] change, -1 if decreases with age, sex, etc., and 1 if increases
3. Calculate average effect size if multiple peptides map to the same gene and are in the same category (constitutive & isoform-infomative; isoform-specific is always separate)

```
conda env create -f ./00_scripts/candidate_peps.yml

conda activate candidate_peps
```
This script filters for genes (across all tissue types and effect types) where more than 90% (and 80%) of the shared peptides (constitutive and isoform-informative) have no change (fraction gene change = 0) and at least one isoform-specific peptide have a change (faction gene change = -1 or 1). The file here is very large.
```
python ./00_scripts/04_candidate_peps_summary_reduced.py
```
This script filters for genes (across all tissue types but only sex effect) where the constitutive paptides have 0 fraction gene change, and the isofom-specific peptides have a change. This file is still very large.
```
python ./00_scripts/04_sex_subset_candidate_peps.py
```
This script filters for genes (across all tissue types but only age effect) where the constitutive paptides have 0 fraction gene change, and the isofom-specific peptides have a change. This file is still very large. 
```
python ./00_scripts/04_age_subset_candidate_peps.py
```
Those didn't include any coexpressed isoform, so I modified the first script to only include genes with coexpressed isoforms.
```
python ./00_scripts/04_0.9_coexpressed.py
```
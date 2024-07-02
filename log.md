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

## Step 2 - Download newest version of PoGo & data
```
cd /project/sheynkman/programs-needs_attentionEFW
wget https://github.com/cschlaffner/PoGo/releases/download/v1.2.3/PoGo_v1.2.3.zip
unzip PoGo_v1.2.3.zip
export PATH=$PATH:/project/sheynkman/programs-needs_attentionEFW/PoGo_v1.2.3/Linux
```
This version of PoGo will work with mouse data, but it does not support Gencode. It only supports Ensambl. So I'm downloading the new dataset.

# Step 3 - Run PoGo
```
cd /project/sheynkman/projects/zhang_mouse_aging

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/all_peptides.txt -format BED

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/pq31153_peptides.txt -format BED

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/pq31190_peptides.txt -format BED

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/pq31193_peptides.txt -format BED

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/pq31194_peptides.txt -format BED

PoGo -fasta ./ensambl_mouse/Mus_musculus.GRCm39.pep.all.fa -gtf ./ensambl_mouse/Mus_musculus.GRCm39.112.gtf -in ./01_Peptides2Pogo/aggregate_peptides.txt -format BED
```

# Step 3 - Convert to bigBED
First, I'll download it then do the conversion
```
cd /project/sheynkman/programs-needs_attentionEFW
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x bedToBigBed

wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x fetchChromSizes

export PATH=$PATH:/project/sheynkman/programs-needs_attentionEFW/

cd /project/sheynkman/projects/zhang_mouse_aging/00_scripts
wget https://github.com/cschlaffner/TrackHubGenerator/blob/master/TrackHubGenerator.pl

cd ..
module load perl/5.36.0

perl 00_scripts/TrackHubGenerator.pl /project/sheynkman/projects/zhang_mouse_aging/ mm39 /project/sheynkman/projects/zhang_mouse_aging/PoGo_output/ /project/sheynkman/programs-needs_attentionEFW/ watts.emily.f@virginia.edu

```

# Count the number of peptides, genes, etc. in each dataset 
```
python ./00_scripts/count_data.py
```
Aggregate - Unique gene_names: 22024, Unique pep_seqs: 717490, Unique transcript_names: 66178
Cleaved - Unique gene_names: 271, Unique pep_seqs: 1977, Unique transcript_names: 581
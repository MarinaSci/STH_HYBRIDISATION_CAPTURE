# Housekeeping of data and reference files 
Author: Marina Papaiakovou, mpapaiakovou[at]gmail.com 

## Contents: 
- Trimming of raw data(using trim-galore)
- Counting reads

```bash 

#conda activate adapter-removal on HPC
conda install -c bioconda adapterremoval
#version
#AdapterRemoval ver. 2.3.2

#running it
#########################################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_adapterremoval_non_collapsed_%j_%a.out
#SBATCH --error=job_adapterremoval_non_collapsed_%j_%a.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-24

source activate adapter-removal

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt | cut -f 1)
#FWD=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt | cut -f 2)
#REV=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt | cut -f 3)

mkdir -p adapterremoval_non_collapsed

AdapterRemoval \
--file1 ${SAMPLE}_1.fq.gz \
--file2 ${SAMPLE}_2.fq.gz \
--gzip \
--qualitymax 64 \
--trimns \
--trimqualities \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
--minlength 30 \
--minquality 20 \
--minadapteroverlap 1 \
--basename adapterremoval_non_collapsed/${SAMPLE}

echo Complete!
  03_adapter_removal.sh (END)
#########################################################

#Installing trim-galore as well 
conda install bioconda::trim-galore
#version
#v0.6.2

#I only had one sample in that list because I put them on a separate folder as a test
#########################################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_trim_galore_%j_%a.out
#SBATCH --error=job_trim_galore_%j_%a.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-24 

source activate trim-galore

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt)

trim_galore --paired ${SAMPLE}_1.fq.gz  ${SAMPLE}_2.fq.gz

echo "I am done trimming"
#########################################################
#I ran trim-galore for both datasets, WGS and hybrid capture 

#counting raw reads
/home/marip3/mbl_genome_skimming/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/02.TRIMMED_DATA
/home/marip3/mbl_genome_skimming/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/02.TRIMMED_DATA

#there is a folder with multiqc stats as output you can get the raw reads as well 
#or can do the following: 
for file in *.fq.gz; do zcat "$file" | echo $(($(wc -l)/4)); done > CAP_read_counts.txt
for file in *.fq.gz; do zcat "$file" | echo $(($(wc -l)/4)); done > WGS_read_counts.txt

#the above will generate a list of read numbers for both R1 and R2. Because they are identical, 
#manually remove the duplicates and import only the one

########################
###### FASTQC ---- #
########################

#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --output=job_fastqc_%j_%a.out
#SBATCH --error=job_fastqc_%j_%a.err
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-48

#SCRIPT TO RUN FASTQC NOT ON A LOOP, BUT IN AN ARRAY

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

#make a list of samples present before you run the script!
#ls -1 *_trimmed_?.fq.gz > allfiles.txt
#because each R1, R2 are separate files for fast qc

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p allfiles.txt)

echo "I will now run FASTQC"

fastqc -o ./FASTQC -t 8  ${SAMPLE}
#remember to have already created that directory

echo "I AM DONE QC-ING!"

#04_fastqc_hpc.sh (END)

#generate names of samples 
ls -1 *_CIGAR_final_nosambamba_nosamclip.bam | cut -d'_' -f1,2 > CAPTURE_FILENAMES.txt
ls -1 *_CIGAR_final_nosambamba_nosamclip.bam | cut -d'_' -f1,2 > WGS_FILENAMES.txt

```

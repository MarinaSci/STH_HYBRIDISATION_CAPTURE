# Indexing and mapping of hybrid capture & WGS datasets
Author: Marina Papaiakovou, mpapaiakovou[at]gmail.com 

## Contents: 
- Indexing of two types of refs: mito and the one with all the bait targets (to calculate percentage of on target mapping)
- Mapping of two datasets, hybrid capture AND WGS to mito (human_mito_ref.fasta)  and the bait reference (PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta)
- PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta, containing all probe/bait sequences obtained by Arbor - was used to generate Figure 1A (might change), to show what percentage of the reads was mapped to on-target 
- human_mito_ref.bed, containing coordinates (start/end) of mitogenomes, was used to map both HYBRID CAPTURE AND WGS DATA TO TRIMMED READS 
- Gettting basic mapping stats with bedtools multicov

- Indexing mitochondrial genome reference
```bash
########################################
## UPDATED SCRIPT FOR INDEXING MITO REF ----
#######################################
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --output=job_index_file_%j.out
#SBATCH --error=job_index_file_%j.err
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

bwa index human_mito_ref.fasta
#05_bwam_index_mito.sh (END)
#######################################
```
- Index reference with all baits/targets 

```bash
##########################################################
## UPDATED SCRIPT FOR INDEXING REF WITH ALL BAITS ----
#########################################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_index_file_%j.out
#SBATCH --error=job_index_file_%j.err
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

bwa index PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta

echo "Done indexing"

#05_bwa_index_bait_targets.sh (END)
#######################################
```

- Script to map hybrid capture probe data to mitogenome reference 

```bash
######################################################################
## UPDATED SCRIPT FOR MITO MAPPING FOR PROBE DATA TO MITO  REFs ----
######################################################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bwa_simplified_%j_%a.out
#SBATCH --error=job_bwa_simplified_%j_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

#TESTING NOW TO RUN THE SCRIPT FROM THE SAME FOLDER ΑS THE TRIMMED FILES - IT WORKS !!!

#make a list of samples present
#ls -1 *_1.fq.gz | sed 's/_1.fq.gz//' > sample_list.txt # sample list needs to only have sample_name_trimmed

#REF_DIR=/mbl/share/workspaces/groups/genome-skimming/03.GLOBAL_SKIM/04.ANALYSIS
TRIM_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/02.TRIMMED_DATA
MAPPING_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/04.ANALYSIS/01.MTDNA_MAPPING


SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt) #I got this working !!!!
echo "I will now run the bwa mem"

bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" human_mito_ref.fasta ${TRIM_DIR}/${SAMPLE}_1.fq.gz ${TRIM_DIR}/${SAMPLE}_2.fq.gz > ${MAPPING_DIR}/${SAMPLE}.sam
#added the -R flag with RG and SM tag so I can easily remove duplicates from the samples down the line

#doublechck with samtools view -H your.bam #and look for the sample ID
echo "get rid of unmapped reads"
samtools view -q 30 -F 4 -S -h ${MAPPING_DIR}/${SAMPLE}.sam > ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam #keep the header -h

echo "filter them by length"
samtools view -h ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam | awk 'length($10) > 80 || $1 ~ /^@/' > ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam

samtools view -S -b ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam > ${MAPPING_DIR}/${SAMPLE}.bam

#sort
samtools sort -o ${MAPPING_DIR}/${SAMPLE}_sorted.bam ${MAPPING_DIR}/${SAMPLE}.bam
#remove duplicates

echo "tag the duplicate reads"
#tag the duplicates here with picard. It will give a warning/error but it does not stpo picard from running
picard -Xmx4G  MarkDuplicates I=${MAPPING_DIR}/${SAMPLE}_sorted.bam REMOVE_DUPLICATES=TRUE O=${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam M=${MAPPING_DIR}/${SAMPLE}_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT

samtools index  ${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam  #you need to sort before you index and you need the sorting before calling duplicates too

samtools view -h ${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam | awk '$6 !~ /H|S/ || $1 ~ /@/' | samtools view -bS - > ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam

samtools index ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam

samtools idxstats ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam > ${MAPPING_DIR}/${SAMPLE}_mapped_reads_no_samclip_no_sambamba_yes_CIGAR_filter.txt

mv *.err  ${MAPPING_DIR}/ERR_FILES
mv *.out ${MAPPING_DIR}/OUT_FILES

echo "I am done"
#06.bwa_mem_mito_simplified.sh
##################
```
- Script to map hybrid capture probe data to probe/bait reference 

```bash
################################################################
## UPDATED SCRIPT FOR MAPPING FOR PROBE DATA TO BAIT REFs ----
#################################################################

#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bwa_targets_only_%j_%a.out
#SBATCH --error=job_bwa_targets_only_%j_%a.err
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

TRIM_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/02.TRIMMED_DATA
MAPPING_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/04.ANALYSIS/05.TARGET_MAPPING

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt) #I got this working !!!!
echo "I will now run the bwa mem"

bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta  ${TRIM_DIR}/${SAMPLE}_1.fq.gz ${TRIM_DIR}/${SAMPLE}_2.fq.gz > ${MAPPING_DIR}/${SAMPLE}_mapped_to_targets_only.sam
#added the -R flag with RG and SM tag so I can easily remove duplicates from the samples down the line

mv *targets_only*.err  ${MAPPING_DIR}/ERR_FILES
mv *targets_only*.out ${MAPPING_DIR}/OUT_FILES

echo "Done mapping!"
#06.bwa_mem_bait_mapping.sh (END)

###################################################################################
```
- Now map the WGS data to the same refs are above
- Indexing of the same references would not change here, so copy the fasta and indexed files to the WGS data 
- Script to map WGS data to mito reference 

```bash

######################################################################
## UPDATED SCRIPT FOR MITO MAPPING FOR PROBE DATA TO MITO  REFs ----
######################################################################

#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bwa_simplified_%j_%a.out
#SBATCH --error=job_bwa_simplified_%j_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

#TESTING NOW TO RUN THE SCRIPT FROM THE SAME FOLDER ΑS THE TRIMMED FILES - IT WORKS !!!

#make a list of samples present
#ls -1 *_1.fq.gz | sed 's/_1.fq.gz//' > sample_list.txt # sample list needs to only have sample_name_trimmed

#REF_DIR=/mbl/share/workspaces/groups/genome-skimming/03.GLOBAL_SKIM/04.ANALYSIS
TRIM_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/02.TRIMMED_DATA
MAPPING_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/04.ANALYSIS/01.MTDNA_MAPPING


SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt) #I got this working !!!!
echo "I will now run the bwa mem"


bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" human_mito_ref.fasta ${TRIM_DIR}/${SAMPLE}_1.fq.gz ${TRIM_DIR}/${SAMPLE}_2.fq.gz > ${MAPPING_DIR}/${SAMPLE}.sam
#added the -R flag with RG and SM tag so I can easily remove duplicates from the samples down the line

#doublechck with samtools view -H your.bam #and look for the sample ID
echo "get rid of unmapped reads"
samtools view -q 30 -F 4 -S -h ${MAPPING_DIR}/${SAMPLE}.sam > ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam #keep the header -h

echo "filter them by length"
samtools view -h ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam | awk 'length($10) > 80 || $1 ~ /^@/' > ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam

samtools view -S -b ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam > ${MAPPING_DIR}/${SAMPLE}.bam

#sort
samtools sort -o ${MAPPING_DIR}/${SAMPLE}_sorted.bam ${MAPPING_DIR}/${SAMPLE}.bam
#remove duplicates

echo "tag the duplicate reads"
#tag the duplicates here with picard. It will give a warning/error but it does not stpo picard from running
picard -Xmx4G  MarkDuplicates I=${MAPPING_DIR}/${SAMPLE}_sorted.bam REMOVE_DUPLICATES=TRUE O=${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam M=${MAPPING_DIR}/${SAMPLE}_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT

samtools index  ${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam  #you need to sort before you index and you need the sorting before calling duplicates too

#remove funny CIGARS
samtools view -h ${MAPPING_DIR}/${SAMPLE}_sorted_dup.bam | awk '$6 !~ /H|S/ || $1 ~ /@/' | samtools view -bS - > ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam

samtools index ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam

#calculate stats
samtools idxstats ${MAPPING_DIR}/${SAMPLE}_sorted_dup_filtered_CIGAR_final_nosambamba_nosamclip.bam > ${MAPPING_DIR}/${SAMPLE}_mapped_reads_no_samclip_no_sambamba_yes_CIGAR_filter.txt

mv *.err  ${MAPPING_DIR}/ERR_FILES
mv *.out ${MAPPING_DIR}/OUT_FILES

echo "I am done"
#06.bwa_mem_mito_simplified.sh
########################################
```
- Script to map WGS data to probe/bait reference 

```bash
################################################################
## UPDATED SCRIPT FOR MAPPING FOR PROBE DATA TO BAIT REFs ----
#################################################################

#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bwa_targets_only_%j_%a.out
#SBATCH --error=job_bwa_targets_only_%j_%a.err
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

TRIM_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/02.TRIMMED_DATA
MAPPING_DIR=/mbl/share/workspaces/groups/marip3-phd/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/04.ANALYSIS/05.TARGET_MAPPING

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt) #I got this working !!!!
echo "I will now run the bwa mem"

bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta ${TRIM_DIR}/${SAMPLE}_1.fq.gz ${TRIM_DIR}/${SAMPLE}_2.fq.gz > ${MAPPING_DIR}/${SAMPLE}_mapped_to_targets_only.sam
#added the -R flag with RG and SM tag so I can easily remove duplicates from the samples down the line

mv *targets_only*.err  ${MAPPING_DIR}/ERR_FILES
mv *targets_only*.out ${MAPPING_DIR}/OUT_FILES

echo "Done mapping!"
#06.bwa_mem_bait_mapping.sh (END)
##################################################

```
- Getting stats now from mapping of both datasets (hybrid capture and WGS) to the mitogenome references

```bash
#get stats
/home/marip3/mbl_genome_skimming/04.PROBE_BAIT/01.PROBE_CAPTURE_DATA/04.ANALYSIS/01.MTDNA_MAPPING
/home/marip3/mbl_genome_skimming/04.PROBE_BAIT/02.SHOTGUN_PROBE_SAMPLES_DATA/04.ANALYSIS/01.MTDNA_MAPPING

bedtools multicov -bams *filtered_CIGAR_final_nosambamba_nosamclip.bam -bed human_mito_ref.bed > ALL_CAPTURE_SAMPLES_BEDTOOLS_MULTICOV.txt
bedtools multicov -bams *filtered_CIGAR_final_nosambamba_nosamclip.bam -bed human_mito_ref.bed > ALL_WGS_SAMPLES_BEDTOOLS_MULTICOV.txt

```
- Ready for some plotting!!
- Note that the mapping script here is simplified compared to the one in the global skim script. The sambamba and samclip filtering was not needed, filtering by duplicates and H/S was taking care of the funny CIGARs in the data
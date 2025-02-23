# Variant calling for hybrid capture and WGS data 
- Variants were called using Bcftools mpileup, adding flag for orphan reads
- Filtering varians for min/max alleles 
- Filtering VCFs for depth 

# Contents 

- Bash script for variant calling using bcftools for CAPTURE and WGS data
- Variant filtering with vcftools for *Ascaris lumbricoides* and *Trichuris trichiura* 
- VCF file filtering for depth (DP>10)

- Variant calling for hybrid capture data 

```bash 
#################################
#### VARIANT CALLING ---- #######

#new environment on the HPC
conda activate bcftools-env
conda install -c bioconda bcftools openssl=1.0
#version that it installed 
#Version: 1.9 (using htslib 1.9)
#I have an old version of conda but I am not updating anything until I am done with the phd
#answer: https://www.biostars.org/p/9480029/

#install vcftools
conda install bioconda::vcftools
#VCFtools (0.1.16)

#CAPTURE DATA
##########################################
####### SCRIPT FOR VARIANT_CALLING.sh ----
##########################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bcftools_%j.out
#SBATCH --error=job_bcftools_%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate bcftools-env

#UPDATE 11 Nov 2024
#I needed to add the '-A' flag to allow for variants that exist on orphaned reads. During my mapping, I remove funny CIGARS, which means I am getting reads of reads that map improperly 
#Which then means that I am left with 'orphaned' reads. That's why bcftools was giving me no variants on the WGS data. 
#WHen I tested, there was no difference in the capture data. But there was a difference on the WGS data. Having paired reads increases the confidence overall, but will check those variants and 
#see what they are. Because I know get Trichuris variants .... 

bcftools mpileup -A --annotate FORMAT/AD -Ov -f human_mito_ref.fasta -d 100000 -b bamlist_for_variant_call | bcftools call --ploidy 1 -mv --skip-variants indels -o CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf

#split by species
#VCFtools (0.1.16)

#take the species separately 
#isolate ALUM
vcftools --vcf CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf --chr NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome --recode --recode-INFO-all --out ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads
#After filtering, kept 510 out of a possible 1947 Sites

#min max alleles 
vcftools --vcf ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER
#After filtering, kept 510 out of a possible 510 Sites


#Isolate TT
vcftools --vcf CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf --chr NC_017750_Trichuris_trichiura_mitochondrion_complete_genome --recode --recode-INFO-all --out TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads
#After filtering, kept 1434 out of a possible 1947 Sites

#min max alleles 
vcftools --vcf TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER
#After filtering, kept 1406 out of a possible 1434 Sites


```

- Variant calling for WGS data 

```bash
#SHOTGUN DATA
##########################################
####### SCRIPT FOR VARIANT_CALLING.sh ----
##########################################
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_bcftools_%j.out
#SBATCH --error=job_bcftools_%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate bcftools-env

bcftools mpileup  -A --annotate FORMAT/AD -Ov -f human_mito_ref.fasta -d 100000 -b bamlist_for_variant_call | bcftools call --ploidy 1 -mv --skip-variants indels -o SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf

#split by species
#VCFtools (0.1.16)

#take the species separately 
#isolate ALUM
vcftools --vcf SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf --chr NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome --recode --recode-INFO-all --out ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads
#After filtering, kept 491 out of a possible 1186 Sites
#min max alleles 
vcftools --vcf ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER
#After filtering, kept 490 out of a possible 491 Sites

#isolate TT
vcftools --vcf SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.vcf --chr NC_017750_Trichuris_trichiura_mitochondrion_complete_genome --recode --recode-INFO-all --out TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads
#After filtering, kept 483 out of a possible 1186 Sites

#min max alleles 
vcftools --vcf TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER
#After filtering, kept 481 out of a possible 483 Sites

echo "I am  done calling variants"

#The above are TOTAL VARIANTS. FOR Ascaris and Trichuris. 
```

- Remove the duplicates (if any)

```bash
#Remove the duplicate variants from the VCf or check if they are all unique
#how to check, using the non depth filtered data as an example
cut -f1,2 TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf | sort | uniq -c | sort
cut -f1,2 ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf | sort | uniq -c | sort
#same for shotgun data (filtered/non filtered) and for capture data (filtered)

#checked and they are unique variants
cut -f1,2 CAPTURE_DATA_ALL_SAMPLES_bcftools_NOMINIMUMALLELEFREQ_mtDNA_coding_genes.recode.vcf | sort | uniq -c | sort
cut -f1,2 SHOTGUN_DATA_ALL_SAMPLES_bcftools_NOMINIMUMALLELEFREQ_mtDNA_coding_genes.recode.vcf | sort | uniq -c | sort

```
- Handling VCF files in R
- Extract allele depth to calculate allele frequencies for CAPTURE & WGS DATA that have not been filtered for depth 
- VCFs below without any filtering 

```{r, warning=FALSE, message=FALSE}
suppressMessages(library(vcfR))
suppressMessages(library(tidyverse))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(stringr))
suppressMessages(library(patchwork))
#NO DEPTH FILTER !!!!!!!! 
####################################
# CAPTURE DATA NO DEPTH FILTER ! ---- 
####################################
# based on this assessemtn: ~/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/07.CODE/Probe_capture_STHs_analysis/XX_DEPTH_OF_COVERAGE_THRESHOLDS.R
#decided to keep the 10X as the depth of choice so will compare variants as 'raw' without any filter on DEPTH and with filter of DP > 10

setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/05.CAPTURE_DATA/02_TRIMMED_DATA/06_VARIANT_CALLING/01_MITOGENOME_VARS/01_WITH_ORPHAN_FLAG/")

#Ascaris 
ALUM_CAPTURE_all_mito_SNPs_bcftools <- read.vcfR ("ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf")

ALUM_CAPTURE_all_mito_SNPs_bcftools <- vcfR2tidy(ALUM_CAPTURE_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

ALUM_CAPTURE_all_mito_SNPs_bcftools_2 <- ALUM_CAPTURE_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
ALUM_CAPTURE_all_mito_SNPs_bcftools_3 <- ALUM_CAPTURE_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

ALUM_CAPTURE_all_mito_SNPs_bcftools_3$CHROM[ALUM_CAPTURE_all_mito_SNPs_bcftools_3$CHROM == 'NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome'] <- 'Ascaris lumbricoides'

#Trichuris 
TT_CAPTURE_all_mito_SNPs_bcftools <- read.vcfR ("TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf")

TT_CAPTURE_all_mito_SNPs_bcftools <- vcfR2tidy(TT_CAPTURE_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

TT_CAPTURE_all_mito_SNPs_bcftools_2 <- TT_CAPTURE_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
TT_CAPTURE_all_mito_SNPs_bcftools_3 <- TT_CAPTURE_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

TT_CAPTURE_all_mito_SNPs_bcftools_3$CHROM[TT_CAPTURE_all_mito_SNPs_bcftools_3$CHROM == 'NC_017750_Trichuris_trichiura_mitochondrion_complete_genome'] <- 'Trichuris trichiura'

CAPTURE_ALL_NO_DEPTH_FILTER <- rbind(ALUM_CAPTURE_all_mito_SNPs_bcftools_3, TT_CAPTURE_all_mito_SNPs_bcftools_3)

####################################
# SHOTGUN DATA NO DEPTH FILTER ! ---- 
####################################
# based on this assessemtn: ~/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/07.CODE/Probe_capture_STHs_analysis/XX_DEPTH_OF_COVERAGE_THRESHOLDS.R
#decided to keep the 10X as the depth of choice so will compare variants as 'raw' without any filter on DEPTH and with filter of DP > 10

setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/06.SHOTGUN_DATA/02_TRIMMED_DATA/06_VARIANT_CALLING/01_MITOGENOME_VARS/01_WITH_ORPHAN_FLAG/")

#Ascaris
ALUM_SHOTGUN_all_mito_SNPs_bcftools <- read.vcfR ("ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf")

ALUM_SHOTGUN_all_mito_SNPs_bcftools <- vcfR2tidy(ALUM_SHOTGUN_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

ALUM_SHOTGUN_all_mito_SNPs_bcftools_2 <- ALUM_SHOTGUN_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
ALUM_SHOTGUN_all_mito_SNPs_bcftools_3 <- ALUM_SHOTGUN_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

ALUM_SHOTGUN_all_mito_SNPs_bcftools_3$CHROM[ALUM_SHOTGUN_all_mito_SNPs_bcftools_3$CHROM == 'NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome'] <- 'Ascaris lumbricoides'

#Trichuris
TT_SHOTGUN_all_mito_SNPs_bcftools <- read.vcfR ("TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf")

TT_SHOTGUN_all_mito_SNPs_bcftools <- vcfR2tidy(TT_SHOTGUN_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

TT_SHOTGUN_all_mito_SNPs_bcftools_2 <- TT_SHOTGUN_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
TT_SHOTGUN_all_mito_SNPs_bcftools_3 <- TT_SHOTGUN_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

TT_SHOTGUN_all_mito_SNPs_bcftools_3$CHROM[TT_SHOTGUN_all_mito_SNPs_bcftools_3$CHROM == 'NC_017750_Trichuris_trichiura_mitochondrion_complete_genome'] <- 'Trichuris trichiura'


#COMBINE THE DATASETS TO COMPARE THEIR SNPS
SHOTGUN_ALL_NO_DEPTH_FILTER <- rbind(ALUM_SHOTGUN_all_mito_SNPs_bcftools_3, TT_SHOTGUN_all_mito_SNPs_bcftools_3)

```

- Filtering VCFs for DP>10 
- BCFtools steps

```bash 

####################################
# CAPTURE DATA WITH DEPTH FILTER ! ---- 
# DP > 10
####################################
# based on this assessment : ~/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/07.CODE/Probe_capture_STHs_analysis/06_DEPTH_OF_COVERAGE_THRESHOLDS.R
#decided to keep the 10X as the depth of choice so will compare variants as 'raw' without any filter on DEPTH and with filter of DP > 10

#bcftools-env in the HPC
#Version: 1.9 (using htslib 1.9)

#CAPTURE DATA 
#filter both files for DP>10 and assess number of SNPs post filter in each 
bcftools view -i 'INFO/DP>10' ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf -o ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf
#499 NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome

bcftools view -i 'INFO/DP>10' TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf -o TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf
#1406 NC_017750_Trichuris_trichiura_mitochondrion_complete_genome

#SHOTGUN DATA
bcftools view -i 'INFO/DP>10' ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf -o ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf
#484 NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome

bcftools view -i 'INFO/DP>10' TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_NO_DEPTH_FILTER.recode.vcf -o TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf
#0


```
- Handling VCFs that have been filtered for DP >10 
- Hybrid capture data & WGS data


```{r, warning = FALSE, message=FALSE}

setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/05.CAPTURE_DATA/02_TRIMMED_DATA/06_VARIANT_CALLING/01_MITOGENOME_VARS/01_WITH_ORPHAN_FLAG/")

#Ascaris
ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools <- read.vcfR ("ALUM_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf")

ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools <- vcfR2tidy(ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 <- ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3 <- ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM[ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM == 'NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome'] <- 'Ascaris lumbricoides'
#CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM[CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM == 'NC_017750_Trichuris_trichiura_mitochondrion_complete_genome'] <- 'Trichuris trichiura'

#Trichuris
TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools <- read.vcfR ("TT_CAPTURE_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf")

TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools <- vcfR2tidy(TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 <- TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3 <- TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)

TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM[TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM == 'NC_017750_Trichuris_trichiura_mitochondrion_complete_genome'] <- 'Trichuris trichiura'


CAPTURE_ALL_DEPTH_FILTERED <- rbind(ALUM_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3, TT_CAPTURE_DEPTH_FILTERED_all_mito_SNPs_bcftools_3)

####################################
# SHOTGUN DATA WITB DEPTH FILTER ! ---- 
# DP > 10
####################################
# based on this assessment: ~/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/07.CODE/Probe_capture_STHs_analysis/06_DEPTH_OF_COVERAGE_THRESHOLDS.R
#decided to keep the 10X as the depth of choice so will compare variants as 'raw' without any filter on DEPTH and with filter of DP > 10

setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/06.SHOTGUN_DATA/02_TRIMMED_DATA/06_VARIANT_CALLING/01_MITOGENOME_VARS/01_WITH_ORPHAN_FLAG/")

#Ascaris
ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools <- read.vcfR ("ALUM_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf")

ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools <- vcfR2tidy(ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 <- ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools$dat

#will need to plot the QUALITY SCORES pper species
ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3 <- ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_2 %>%
  dplyr::select(CHROM, POS, QUAL, REF, ALT,  DP4, Indiv)


ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM[ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM == 'NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome'] <- 'Ascaris lumbricoides'
#SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM[SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3$CHROM == 'NC_017750_Trichuris_trichiura_mitochondrion_complete_genome'] <- 'Trichuris trichiura'

#Trichuris
TT_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools <- read.vcfR ("TT_SHOTGUN_DATA_ALL_SAMPLES_bcftools_with_orphan_reads_NOMINIMUMALLELEFREQ_DP10_FILTER.recode.vcf")

#not working because it's empty ....
#TT_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools <- vcfR2tidy(TT_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools, single_frame = TRUE, toss_INFO_column = TRUE, alleles =TRUE)

SHOTGUN_ALL_DEPTH_FILTERED <- ALUM_SHOTGUN_DEPTH_FILTERED_all_mito_SNPs_bcftools_3

```

- SNPs from either of VCFs (pre and post filter) will be compared between datasets (hybrid capture and WGS)
- Keep reading! :) 





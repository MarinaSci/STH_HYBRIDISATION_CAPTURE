---
title: "STH_hybridisation_capture_02_ON_TARGET_ENRICHMENT"
author: "Marina Papaiakovou"
date: "2025-02-22"
output: html_document
---

# On-target percentage WGS vs hybridisation capture, absolute and fold-enrichment

- 

## Contents: 

- Bash code to get basic stats with flagstats from both hybrid capture and WGS mapped to probe/baits
- 

- You need the sam files from mapping hybrid capture and WGS (separately) to the probe/bait fasta (PARASITE_BAIT_TARGETS_KRILL_FORMATTED.fasta) 

```bash
for i in *mapped_to_targets_only.sam; do   base_name=${i%%_trimmed*};   samtools flagstat "$i" > "${base_name}_flagstat.txt"; done

```
- Need to collate all the stats in a single txt file 

```bash
############################################################################################################
#COLLATE FLAGSTAT RESULTS
##############################################################
#!/bin/bash

# Directory containing samtools flagstat outputs
OUTPUT_DIR=.

# Create a file to store the results
RESULT_FILE="mapping_percentages.txt"
echo -e "Sample\tTotal_Reads_Mapped(%)" > $RESULT_FILE

# Loop through all flagstat output files
for file in "$OUTPUT_DIR"/*.txt; do
# Extract the sample name from the filename
sample_name=$(basename "$file" .txt)

# Use grep and awk to extract the percentage from the "mapped" line
mapped_percentage=$(grep "mapped (" "$file" | awk -F'[()%]' '{print $2}')

# Append the sample name and percentage to the results file
echo -e "${sample_name}\t${mapped_percentage}" >> $RESULT_FILE
done

# Print the results
cat $RESULT_FILE
#COLLATE_MAPPED_PERCENTAGE.sh (END)
#####################################################################################################################

```
- Import them in R for plotting 

```{r, warning = FALSE,message = FALSE}
library(reshape2)
library(tidyverse)

#import the files with the results from samtools flagstat
CAPTURE_MAPPED_READS <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/05.CAPTURE_DATA/02_TRIMMED_DATA/07_CONTAMINATION_ASSESSMENT/CAPTURE_PARASITE_KRILL_mapping_percentages.txt", sep = "\t", header =T )
colnames(CAPTURE_MAPPED_READS) <- c('sample_id', 'reads_mapped_percentage')
CAPTURE_MAPPED_READS$dataset <- 'capture'

#remove last row
CAPTURE_MAPPED_READS <- CAPTURE_MAPPED_READS %>% filter(row_number() <= n()-1)


#isolate the numbers from the sample_ids
CAPTURE_MAPPED_READS_2 <- as.data.frame(sapply(strsplit(CAPTURE_MAPPED_READS$sample_id, "_"), function(x) x[2]))
colnames(CAPTURE_MAPPED_READS_2) <- 'sample_id'

#collate the datasets
CAPTURE_MAPPED_READS <- CAPTURE_MAPPED_READS %>%
  dplyr::select(2) #remove the previous sample_id
CAPTURE_MAPPED_READS_3 <- cbind(CAPTURE_MAPPED_READS_2,CAPTURE_MAPPED_READS )
CAPTURE_MAPPED_READS_3$dataset <- 'capture'


#SHOTGUN
SHOTGUN_MAPPED_READS <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/05.Hybridization_probe/06.SHOTGUN_DATA/02_TRIMMED_DATA/07_CONTAMINATION_ASSESSMENT/SHOTGUN_PARASITE_KRILL_mapping_percentages.txt", sep = "\t", header =T)
colnames(SHOTGUN_MAPPED_READS) <- c('sample_id', 'reads_mapped_percentage')
SHOTGUN_MAPPED_READS$dataset <- 'wgs'

#remove first row
SHOTGUN_MAPPED_READS = SHOTGUN_MAPPED_READS[-1,]


#isolate the numbers from the sample_ids
SHOTGUN_MAPPED_READS$sample_id <- as.character(SHOTGUN_MAPPED_READS$sample_id)

SHOTGUN_MAPPED_READS_2 <- as.data.frame(sapply(strsplit(SHOTGUN_MAPPED_READS$sample_id, "_"), function(x) x[2]))
colnames(SHOTGUN_MAPPED_READS_2) <- 'sample_id'

#collate the datasets
SHOTGUN_MAPPED_READS <- SHOTGUN_MAPPED_READS %>%
  dplyr::select(2) #remove the previous sample_id
SHOTGUN_MAPPED_READS_3 <- cbind(SHOTGUN_MAPPED_READS_2,SHOTGUN_MAPPED_READS )
SHOTGUN_MAPPED_READS_3$dataset <- 'wgs'


#COLLATE BOTH DATASETS 
ALL_MAPPED_READS <- rbind(CAPTURE_MAPPED_READS_3, SHOTGUN_MAPPED_READS_3)

#Create a species vector
species <- c("ascaris lumbricoides", "trichuris trichiura")

# Expand the data frame by duplicating the rows and adding species
ALL_MAPPED_READS_EXPANDED <- ALL_MAPPED_READS[rep(1:nrow(ALL_MAPPED_READS), each = length(species)), ]

# Add the 'species' column
ALL_MAPPED_READS_EXPANDED$species <- rep(species, times = nrow(ALL_MAPPED_READS))


#make the factors of dataset so you can print the wgs first 
ALL_MAPPED_READS$dataset <- factor(ALL_MAPPED_READS$dataset, levels = c("wgs", "capture"))


# Create the contamination plot to call it later
CONTAMINATION_PLOT <- ggplot(ALL_MAPPED_READS, aes(x = dataset, y = reads_mapped_percentage, fill = dataset)) +
  geom_boxplot(outlier.shape = NA) + #THIS WILL TURN OFF THE OUTLIARS THAT WERE COMING UP AS A POINT!!! 
 # geom_point()+
 # geom_jitter(aes(shape = ifelse(sample_id == "Ascaris", "worm", "faecal")), 
    #          width = 0.2, height = 0) + 
  geom_point(position = position_jitter(seed = 1), aes(shape = ifelse(sample_id == "Ascaris", "worm", "faecal")), 
             width = 0.2, height = 0, size = 2)+ #setting a seed = 1, or seed=42 makes the code reproducible. So you avoid generating plots with jitters that are in different spots on the x axis
  theme_bw() +
  labs(
    title = "Percentage of reads mapped with BWA MEM - no filters",
    x = "Dataset: wgs vs capture",
    y = "Mapped reads percentage", shape = "shape"
  ) +
  scale_fill_manual(values = c("wgs" = "#1f77b4", "capture" = "#ff7f0e")) +
  theme(
    axis.text.x = element_text(color = "black"),  # Explicitly set x-axis text color to black
    axis.text.y = element_text(color = "black")   # Set y-axis text color to black (optional)
  )
# Custom colors for points
CONTAMINATION_PLOT
#for the paper, to report the mean percentage etcf
summary_stats <- ALL_MAPPED_READS %>%
  filter(sample_id != "Ascaris") %>%  # Exclude Ascaris
  group_by(dataset) %>%
  summarise(
    mean_percentage = round(mean(reads_mapped_percentage, na.rm = TRUE), 2),
    sd_percentage = round(sd(reads_mapped_percentage, na.rm = TRUE), 2),
    min_value = round(min(reads_mapped_percentage, na.rm = TRUE), 2),
    max_value = round(max(reads_mapped_percentage, na.rm = TRUE), 2)
  ) %>%
  mutate(range = round(max_value - min_value, 2))

#just for Ascaris 
summary_stats_Ascaris <- ALL_MAPPED_READS %>%
  filter(sample_id == "Ascaris") %>%  # Exclude Ascaris
  group_by(dataset) %>%
  summarise(
    mean_percentage = mean(reads_mapped_percentage, na.rm = TRUE),
    sd_percentage = sd(reads_mapped_percentage, na.rm = TRUE),
    min_value = min(reads_mapped_percentage, na.rm = TRUE),
    max_value = max(reads_mapped_percentage, na.rm = TRUE)
  ) %>%
  mutate(range = max_value - min_value)


```

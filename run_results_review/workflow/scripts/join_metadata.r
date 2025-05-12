library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
outlier_status_fp <- args[1]
counts_fp <- args[2]
#metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
metadata_fp <- args[3]

outlier_status <- read_csv(outlier_status_fp)
counts <- read_csv(counts_fp)
#metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
metadata <- read_csv(metadata_fp)

output_file <- args[4]
RIN_file <- args[5]

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch, sampleID)

#####################-----Low RIN-----########################
low_RIN <- metadata %>% select(sampleID, RIN) %>% filter(RIN < 7)
write_csv(low_RIN, RIN_file)

#####################-----Counts-----########################
joined_counts <- left_join(counts, metadata, by="sampleID")
joined_counts <- joined_counts %>% filter(RIN >= 7)

#####################-----Outlier Status-----########################
joined_outlier_status <- left_join(outlier_status, metadata, by="sampleID")
joined_status <- joined_outlier_status %>% filter(RIN >= 7)

#####################-----combined-----########################
joined <- left_join(joined_outlier_status, counts)
write_csv(joined, output_file)

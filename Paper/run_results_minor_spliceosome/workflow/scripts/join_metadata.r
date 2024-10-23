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

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch, sampleID)

#####################-----Counts-----########################
joined_counts <- left_join(counts, metadata, by="sampleID")

write_csv(joined_counts, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

#####################-----Outlier Status-----########################
joined_outlier_status <- left_join(outlier_status, metadata, by="sampleID")

write_csv(joined_outlier_status, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/outlier_status_metadata_filtered_with_sibs.csv")

#####################-----combined-----########################
joined <- left_join(joined_outlier_status, joined_counts)
write_csv(joined, output_file)

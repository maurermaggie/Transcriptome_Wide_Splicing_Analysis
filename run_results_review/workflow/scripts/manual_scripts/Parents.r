library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_Jan_6/DataFrames/metadata_counts_outlier_joined.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")
family_id_RD380 <- metadata %>% filter(RDID == "RD380") %>% pull(family_identifier)
family_RD380 <- metadata %>% filter(family_identifier == family_id_RD380) %>% select(RDID)

joined %>% filter(sampleID %in% c("RD380", "RD381", "RD382")) %>% select(no_MIGs_with_theta_juncs)

joined %>% filter(sampleID %in% c("RD380", "RD381", "RD382")) 

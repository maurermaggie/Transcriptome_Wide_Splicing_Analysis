library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

FRASER_metadata_important <- read_csv("/home/maurertm/smontgom/shared/UDN/FRASER_analysis/scripts/FRASER_snakemake/output/FRASER_all/FRASER_metadata.csv")

######################################################
################----Subset_Data----###################
######################################################
FRASER_metadata_subset <- FRASER_metadata_important %>% filter(!is.na(outlier_type)) %>% select(outlier_type, )
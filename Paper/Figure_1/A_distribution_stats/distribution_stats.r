library(readxl)
library(magrittr)
require(data.table)
library(tidyverse)

metadata_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

#################----Affected Status----#################
metadata_joined %>% group_by(affected_status) %>% tally
metadata_joined %>% filter(affected_status == "Unknown") %>% select(RDID, UDN_ID, GSS_ID)

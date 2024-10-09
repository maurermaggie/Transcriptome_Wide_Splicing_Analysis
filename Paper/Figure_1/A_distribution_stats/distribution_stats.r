library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")

#################----Affected Status----#################
metadata %>% group_by(affected_status) %>% tally
metadata %>% filter(affected_status == "Unknown") %>% select(RDID, UDN_ID, GSS_ID)

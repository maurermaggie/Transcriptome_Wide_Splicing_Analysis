library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv") %>% select(GSS_ID, UDN_ID, RDID, batch)
batchz <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchz.txt")
batchy <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchy.txt")
batchx <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchx.txt")
new_batch_info <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_batch_new.csv")
batch_info <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/Kevins_Batch_plot.csv") %>% select(PPID, `RD - ID`, `SEQ RUN #`)
metadata_UDN <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/2017to2022_metadata.tsv") %>% select(indv_id, wetlab_id, sample_id, seq_batch...8)

no_batch <- metadata %>% filter(is.na(batch)) %>% select(RDID, GSS_ID, UDN_ID)
intersect(no_batch$UDN_ID, metadata_UDN$indv_id)
intersect(no_batch$UDN_ID, metadata_UDN$wetlab_id)
intersect(no_batch$UDN_ID, metadata_UDN$sample_id)
intersect(no_batch$UDN_ID, batch_info$PPID)
intersect(no_batch$UDN_ID, batch_info$`RD - ID`)
intersect(no_batch$UDN_ID, new_batch_info$UDN_ID)
intersect(no_batch$UDN_ID, new_batch_info$GSS_ID)
intersect(no_batch$UDN_ID, new_batch_info$RDID)
intersect(no_batch$UDN_ID, batchx$Batch)
intersect(no_batch$UDN_ID, batchy$Batch)
intersect(no_batch$UDN_ID, batchz$Batch)

batch_info_1_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 1) %>% pull(PPID)
one <- metadata %>% filter(UDN_ID %in% batch_info_1_UDN_ID)
wrong_UDN <- one %>% filter(batch != 1) %>% pull(UDN_ID)
metadata_UDN %>% filter(indv_id %in% wrong_UDN) %>% select(batch)

batch_info_2_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 2) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_2_UDN_ID)

batch_info_3_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 3) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_3_UDN_ID)

batch_info_4_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 4) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_4_UDN_ID)

batch_info_5_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 5) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_5_UDN_ID)

batch_info_6_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 6) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_6_UDN_ID)

batch_info_7_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 7) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_7_UDN_ID)

batch_info_8_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 8) %>% pull(PPID)
metadata %>% filter(UDN_ID %in% batch_info_8_UDN_ID)

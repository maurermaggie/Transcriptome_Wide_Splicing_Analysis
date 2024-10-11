library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

outlier_status <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/outlier_status_filtered_with_sibs.csv")
counts <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_filtered_with_sibs.csv")
#metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch)

#####################-----Counts-----########################
names(counts)[names(counts) == 'sampleID'] <- 'RDID'
all_filter_RD <- filter(counts, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)
names(all_RDID)[names(all_RDID) == 'RDID'] <- 'sampleID'

names(counts)[names(counts) == 'RDID'] <- 'GSS_ID'
all_filter_GSS <- filter(counts, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)
names(all_GSS)[names(all_GSS) == 'GSS_ID'] <- 'sampleID'

names(counts)[names(counts) == 'GSS_ID'] <- 'UDN_ID'
all_filter_UDN <- filter(counts, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)
names(all_UDN)[names(all_UDN) == 'UDN_ID'] <- 'sampleID'

names(counts)[names(counts) == 'UDN_ID'] <- 'sampleID'

joined_counts <- bind_rows(all_RDID, all_GSS, all_UDN)

write_csv(joined_counts, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

#####################-----Outlier Status-----########################
names(outlier_status)[names(outlier_status) == 'sampleID'] <- 'RDID'
all_filter_RD <- filter(outlier_status, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)
names(all_RDID)[names(all_RDID) == 'RDID'] <- 'sampleID'

names(outlier_status)[names(outlier_status) == 'RDID'] <- 'GSS_ID'
all_filter_GSS <- filter(outlier_status, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)
names(all_GSS)[names(all_GSS) == 'GSS_ID'] <- 'sampleID'

names(outlier_status)[names(outlier_status) == 'GSS_ID'] <- 'UDN_ID'
all_filter_UDN <- filter(outlier_status, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)
names(all_UDN)[names(all_UDN) == 'UDN_ID'] <- 'sampleID'

names(outlier_status)[names(outlier_status) == 'UDN_ID'] <- 'sampleID'

joined_outlier_status <- bind_rows(all_RDID, all_GSS, all_UDN)
names(joined_outlier_status)[names(joined_outlier_status) == 'RDID'] <- 'sampleID'

write_csv(joined_outlier_status, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/outlier_status_metadata_filtered_with_sibs.csv")

#####################-----combined-----########################
joined <- left_join(joined_outlier_status, joined_counts)
write_csv(joined, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

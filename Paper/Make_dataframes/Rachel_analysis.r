library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

Rachel_metadata <- fread("/home/maurertm/smontgom/shared/UDN/Data/SampleTables/sample_table_UDN_sep2023.txt") %>% select(INDIVIDUAL, SAMPLE, TISSUE, INSTITUTION)
my_metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")

colnames(Rachel_metadata) <- c("UDN_ID", "RDID", "Tissue_Type", "Source")
Rachel_metadata <- filter(Rachel_metadata, Tissue_Type == "Blood")
Rachel_metadata$RDID <- gsub("run18", "", Rachel_metadata$RDID)

find_old_UDN <- intersect(my_metadata$UDN_ID, Rachel_metadata$UDN_ID)
find_old_RDID <- intersect(my_metadata$RDID, Rachel_metadata$RDID)

find_old_UDN_ID <- filter(my_metadata, UDN_ID %in% find_old_UDN) %>% pull(ID)
find_old_RDID_ID <- filter(my_metadata, RDID %in% find_old_RDID) %>% pull(ID)

find_old <- c(find_old_RDID_ID, find_old_UDN_ID) %>% unique

my_metadata$old <- ifelse(my_metadata$ID %in% find_old, "Old", "Novel")
my_metadata %>% group_by(old) %>% tally()

old_in_Rachel <- filter(Rachel_metadata, UDN_ID %in% find_old_UDN) 
old_in_Rachel_utah <- old_in_Rachel %>% filter(Source == "UDN Utah")

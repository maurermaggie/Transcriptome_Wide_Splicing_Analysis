library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_minor_spliceosome/output/output_Nov_18_try5/DataFrames/metadata_counts_outlier_joined.csv")
leafcutter <- fread("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/LeafCutterMD/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz")
leafcutter$V5 %>% unique

intersect(metadata$RDID, leafcutter$V5) %>% unique %>% length #270 in length
setdiff(metadata$RDID, leafcutter$V5)
#NA

novel <- metadata %>% filter(is.na(RDID)) %>% select(GSS_ID, UDN_ID, batch)
old <- metadata %>% filter(!is.na(RDID)) %>% select(GSS_ID, UDN_ID, RDID, batch)
write_csv(novel, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/novel_samples.csv")
write_csv(old, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/old_samples.csv")

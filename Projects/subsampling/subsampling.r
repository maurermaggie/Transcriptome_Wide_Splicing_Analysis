library(biomaRt)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(GO.db)
library(readxl)
require(data.table)
library(tidyverse)

all_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/config/symlink_blood.csv", col_names = FALSE)
colnames(all_files) <- "filepath"
basenames <- lapply(all_files, basename) %>% as.data.frame
colnames(basenames) <- "filename"
base <-  basenames %>% mutate(base = str_replace_all(filename, ".Aligned.sortedByCoord.out.bam.bai|.Aligned.sortedByCoord.out.bam|.Aligned.sortedByCoord.out.sorted_opticalDup2500_minMQ30.bam.bai|.Aligned.sortedByCoord.out.sorted_opticalDup2500_minMQ30.bam|.Aligned.sortedByCoord.sorted_opticalDup2500_minMQ30.bam.bai|.Aligned.sortedByCoord.sorted_opticalDup2500_minMQ30.bam|_star_hg38", ""))

base <- bind_cols(all_files, base)

MIGs <- base %>% filter(base %in% c("RD380", "RD268", "GSS225379"))
no_MIGs <- base %>% filter(! base %in% c("RD380", "RD268", "GSS225379")) %>% group_by(base) %>% arrange(desc(base))
RD268 <- base %>% filter(base == "RD268")

twenty4 <- no_MIGs %>% head(n=48)
fourty9 <- no_MIGs %>% head(n=98)
ninety9 <- no_MIGs %>% head(n=198)
one49 <- no_MIGs %>% head(n=298)
one99 <- no_MIGs %>% head(n=398)
two99 <- no_MIGs %>% head(n=598)
three99 <- no_MIGs %>% head(n=798)

twenty5 <- bind_rows(twenty4, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
fifty <- bind_rows(fourty9, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
one00 <- bind_rows(ninety9, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
one50 <- bind_rows(one49, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
two00 <- bind_rows(one99, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
three00 <- bind_rows(two99, RD268) %>% dplyr::pull(filepath) %>% as.data.frame
four00 <- bind_rows(three99, RD268) %>% dplyr::pull(filepath) %>% as.data.frame

write_csv(twenty5, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_25.csv", col_names=FALSE)
write_csv(fifty, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_50.csv", col_names=FALSE)
write_csv(one00, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_100.csv", col_names=FALSE)
write_csv(one50, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_150.csv", col_names=FALSE)
write_csv(two00, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_200.csv", col_names=FALSE)
write_csv(three00, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_300.csv", col_names=FALSE)
write_csv(four00, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/subsampling/subset_400.csv", col_names=FALSE)

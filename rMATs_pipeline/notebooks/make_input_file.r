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

file <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/config/symlink_blood.csv", col_names = FALSE)
basenames <- lapply(file, basename) %>% as.data.frame
output_file <- bind_cols(file, basenames)
colnames(output_file) <- c("File_path", "ID")
output_file <- output_file[- grep("bai", output_file$File_path),]

write_csv(output_file, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/rMATS_pipeline/config/filepaths_all.csv")

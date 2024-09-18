library(tidyverse)
library(readxl)
require(data.table)
library(magrittr)
library(RColorBrewer) # for a colourful plot
library(ggrepel)

files <- list.files("/home/maurertm/smontgom/shared/UDN/Gateway/Downloads/FASTQ", pattern = "\\R1.fastq.gz$", ignore.case = TRUE)

blood <- files[grepl("blood", files)]
fibro <- files[grepl("fibroblast", files)]
not_fibro <- files[!grepl("fibroblast", files)]
not_type <- not_fibro[!grepl("blood", not_fibro)]

not_type_df <- data.frame(ID = gsub('_R1.fastq.gz', '', not_type), Type = "None Provided")

fibro_cleaned <- fibro %>% gsub('-M-fibroblast_R1.fastq.gz', '', .) %>%
                    gsub("-F-fibroblast_R1.fastq.gz", '', .) %>%
                    gsub("-P-fibroblast_R1.fastq.gz", '', .)
fibro_df <- data.frame(ID = fibro_cleaned, Type = "Fibroblast")

blood_cleaned <- blood %>% gsub('-F-blood_R1.fastq.gz', '', .) %>%
                    gsub("-P-blood_R1.fastq.gz", '', .) %>%
                    gsub("-M-blood_R1.fastq.gz", '', .) %>%
                    gsub("-O-blood_R1.fastq.gz", '', .)
blood_df <- data.frame(ID = blood_cleaned, Type = "Blood")

dataframe <- bind_rows(not_type_df, fibro_df, blood_df)
dataframe$issues <- c("None", "Low Uniquely mapped reads; high % of reads unmapped b/c too short","Low Uniquely mapped reads; high % of reads unmapped b/c too short",
"Low Uniquely mapped reads; high % of reads unmapped b/c too short", "Never finished; high multimapped; low mapped unique", "Low Uniquely mapped reads; high % of reads unmapped b/c too short",
"Low Uniquely mapped reads; high % of reads unmapped b/c too short", "Low Uniquely mapped reads; high % of reads unmapped b/c too short",
"Low Uniquely mapped reads; high % of reads unmapped b/c too short", "None", "Low Uniquely mapped reads; high % of reads unmapped b/c too short",
"None", "None", "Can't Find", "None", "None", "None", "None", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "0% uniquely mapped reads", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique",
"Never finished; high multimapped; low mapped unique", "Never finished; high multimapped; low mapped unique")

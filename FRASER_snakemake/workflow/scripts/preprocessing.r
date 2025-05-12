library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

args <- commandArgs(TRUE)
input_directory <- args[1]
output_directory <- args[2]

files <- list.files(input_directory)
bam_files <- str_subset(files, "bam$")
full_bam_files <- paste(input_directory, "/", bam_files, sep="")
file_names <- trimws(basename(bam_files))
samples <- gsub("\\_.*", "", file_names)
samples <- gsub(".Aligned.sortedByCoord.out.sorted", "", samples)
samples <- gsub(".Aligned.sortedByCoord.sorted", "", samples)

FRASER_input <- data.frame(sampleID = samples, bamFile = full_bam_files, pairedEnd=TRUE) %>% as.tibble()
filepath <- paste(output_directory, "/FRASER_input.csv", sep="")
write_csv(FRASER_input, filepath)

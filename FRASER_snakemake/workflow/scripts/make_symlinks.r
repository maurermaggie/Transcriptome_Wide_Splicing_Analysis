library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

args <- commandArgs(TRUE)
file_list <- args[1]
input_directory <- args[2]

print(file_list)
print(input_directory)
 
all_filepaths <- read_csv(file_list, col_names = FALSE) %>% as.data.frame
colnames(all_filepaths) <- c("filepaths")

print(all_filepaths)

all_filepaths_missing_beginning <- gsub('.*/ ?(\\w+)', '\\1', all_filepaths$filepaths)
all_filepaths_destination <- paste(input_directory, "/", all_filepaths_missing_beginning, sep = "") %>% as.data.frame
all_filepaths_command <- data.frame("ln -s", all_filepaths, " ", all_filepaths_destination)

print(all_filepaths_command)

output_file <- paste(input_directory, "/", "FRASER_symlinking_instructions_unclean.sh", sep="")
write.table(all_filepaths_command, output_file, col.names = FALSE, row.names = FALSE)

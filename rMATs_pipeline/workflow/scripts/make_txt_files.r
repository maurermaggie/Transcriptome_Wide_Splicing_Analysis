library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

args <- commandArgs(TRUE)
file_path <- args[1]
data_directory <- args[2]

base_name <- basename(file_path) %>% unlist()

file_name <- paste0(data_directory, "/", base_name, ".txt")

write.table(file_path, file_name, col.names=FALSE, row.names=FALSE, sep="", na ="", quote=FALSE)

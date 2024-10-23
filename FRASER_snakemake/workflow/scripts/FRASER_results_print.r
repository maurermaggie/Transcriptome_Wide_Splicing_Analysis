library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

args <- commandArgs(TRUE)
input_file <- args[1]
FRASER_type <- args[2]
output_file <- args[3]

inputs <- read_csv(input_file) %>% pull(sampleID) %>% list()

if (FRASER_type == "FRASER1") {
    message <- paste("You have successfully run FRASER1! on", inputs, sep = " ")
    write_lines(message, output_file)
} else if (FRASER_type == "FRASER2") {
    message <- paste("You have successfully run FRASER2! on", inputs, sep = " ")
    write_lines(message, output_file)
} else {
    message <- paste("You have successfully run FRASER1 and FRASER2! on", inputs, sep = " ")
    write_lines(message, output_file)
}
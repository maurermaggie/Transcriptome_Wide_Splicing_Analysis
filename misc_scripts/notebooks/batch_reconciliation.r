library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")

no_batch <- metadata %>% filter(is.na(batch))

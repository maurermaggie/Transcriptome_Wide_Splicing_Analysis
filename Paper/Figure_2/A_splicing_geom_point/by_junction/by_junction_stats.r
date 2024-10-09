library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

metadata_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

########################################################################
####################-----By Junction per Person-----####################
########################################################################
joined <- metadata_joined %>% select(sampleID, ends_with("Junctions"))

two_sd_psi3 <- mean(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions)
two_sd_psi5 <- mean(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions)
two_sd_theta <- mean(joined$Theta_Junctions) + sd(joined$Theta_Junctions) + sd(joined$Theta_Junctions)

filter(joined, Psi3_Junctions >= two_sd_psi3) %>% nrow
filter(joined, Psi5_Junctions >= two_sd_psi5) %>% nrow
filter(joined, Theta_Junctions >= two_sd_theta) %>% nrow

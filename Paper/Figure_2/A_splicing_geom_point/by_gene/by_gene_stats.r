library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

########################################################################
######################-----By Gene per Person-----######################
########################################################################
two_sd_psi3 <- mean(joined$Psi3_Genes) + sd(joined$Psi3_Genes) + sd(joined$Psi3_Genes)
two_sd_psi5 <- mean(joined$Psi5_Genes) + sd(joined$Psi5_Genes) + sd(joined$Psi5_Genes)
two_sd_theta <- mean(joined$Theta_Genes) + sd(joined$Theta_Genes) + sd(joined$Theta_Genes)

filter(joined, Psi3_Genes >= two_sd_psi3) %>% nrow
filter(joined, Psi5_Genes >= two_sd_psi5) %>% nrow
filter(joined, Theta_Genes >= two_sd_theta) %>% nrow

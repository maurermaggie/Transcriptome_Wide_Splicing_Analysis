library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
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

psi3_outs <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID)
psi5_outs <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID)
theta_outs <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID)

stats_list <- list(`Psi3 Outliers`= psi3_outs, `Psi5 Outliers`=psi5_outs, `Theta Outliers` =theta_outs )

stats <- ggVennDiagram(stats_list, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)

stats_venn <- ggVennDiagram(stats_list,  order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)

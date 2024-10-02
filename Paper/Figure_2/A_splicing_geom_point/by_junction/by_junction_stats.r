library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

########################################################################
######################-----By Gene per Person-----######################
########################################################################
filter_junctions <- function(dataframe, filter, column_name, significance){
      dataframe_filtered <- dataframe %>% filter(padjust <= significance) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
      if (filter == "psi3") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi3")
        print("filtering by psi3")
      } else if (filter == "psi5") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi5")
        print("filtering by psi5")
      } else if (filter == "theta") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "theta")
        print("filtering by theta")
      } else {
         dataframe_filtered < dataframe_filtered
      }
      dataframe_count <- dataframe_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

########################################################################
####################-----By Junction per Person-----####################
########################################################################

################-----Get Number Junctions per Person-----#################
#0.05
psi3_count <- filter_junctions(all_uncompiled, "psi3", "Psi3", 0.05)
psi5_count <- filter_junctions(all_uncompiled, "psi5", "Psi5", 0.05)
theta_count <- filter_junctions(all_uncompiled, "theta", "Theta", 0.05)

#join
joined <- left_join(psi3_count, psi5_count) %>% left_join(theta_count)

two_sd_psi3 <- mean(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3)
two_sd_psi5 <- mean(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5)
two_sd_theta <- mean(joined$Theta) + sd(joined$Theta) + sd(joined$Theta)

filter(psi3_count, Psi3 >= two_sd_psi3) %>% nrow
filter(psi5_count, Psi5 >= two_sd_psi5) %>% nrow
filter(theta_count, Theta >= two_sd_theta) %>% nrow

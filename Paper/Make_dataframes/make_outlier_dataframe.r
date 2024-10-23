library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

counts <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_filtered_with_sibs.csv")

########################################################################
######################-----By Gene per Person-----######################
########################################################################
################-----Get Outlier Status-----#################
get_two_sd_status <- function(dataframe, column_name){
        two_sd <- mean(dataframe[[column_name]]) + sd(dataframe[[column_name]]) + sd(dataframe[[column_name]])
        save_column_name <- paste(column_name, "Outlier_Status", sep="_")
        dataframe <- dataframe %>% select(sampleID, column_name)
        dataframe[[save_column_name]] <- ifelse(dataframe[[column_name]] >= two_sd, 1, 0)
        dataframe <- dataframe %>% select(-column_name)
}

#junctions
All_Junctions <- get_two_sd_status(counts, "All_Junctions")
All_Junctions_001 <- get_two_sd_status(counts, "All_Junctions_001")
All_Junctions_01 <- get_two_sd_status(counts, "All_Junctions_01")
All_Junctions_1e_6 <- get_two_sd_status(counts, "All_Junctions_1e_6")

Psi3_Junctions <- get_two_sd_status(counts, "Psi3_Junctions")
Psi3_Junctions_001 <- get_two_sd_status(counts, "Psi3_Junctions_001")
Psi3_Junctions_01 <- get_two_sd_status(counts, "Psi3_Junctions_01")
Psi3_Junctions_1e_6 <- get_two_sd_status(counts, "Psi3_Junctions_1e_6")

Psi5_Junctions <- get_two_sd_status(counts, "Psi5_Junctions")
Psi5_Junctions_001 <- get_two_sd_status(counts, "Psi5_Junctions_001")
Psi5_Junctions_01 <- get_two_sd_status(counts, "Psi5_Junctions_01")
Psi5_Junctions_1e_6 <- get_two_sd_status(counts, "Psi5_Junctions_1e_6")

Theta_Junctions <- get_two_sd_status(counts, "Theta_Junctions")
Theta_Junctions_001 <- get_two_sd_status(counts, "Theta_Junctions_001")
Theta_Junctions_01 <- get_two_sd_status(counts, "Theta_Junctions_01")
Theta_Junctions_1e_6 <- get_two_sd_status(counts, "Theta_Junctions_1e_6")

Jaccard_Junctions <- get_two_sd_status(counts, "Jaccard_Junctions")
Jaccard_Junctions_01 <- get_two_sd_status(counts, "Jaccard_Junctions_01")
Jaccard_Junctions_001 <- get_two_sd_status(counts, "Jaccard_Junctions_001")
Jaccard_Junctions_1e_6 <- get_two_sd_status(counts, "Jaccard_Junctions_1e_6")

#genes
All_Genes <- get_two_sd_status(counts, "All_Genes")
All_Genes_001 <- get_two_sd_status(counts, "All_Genes_001")
All_Genes_01 <- get_two_sd_status(counts, "All_Genes_01")
All_Genes_1e_6 <- get_two_sd_status(counts, "All_Genes_1e_6")

Psi3_Genes <- get_two_sd_status(counts, "Psi3_Genes")
Psi3_Genes_001 <- get_two_sd_status(counts, "Psi3_Genes_001")
Psi3_Genes_01 <- get_two_sd_status(counts, "Psi3_Genes_01")
Psi3_Genes_1e_6 <- get_two_sd_status(counts, "Psi3_Genes_1e_6")

Psi5_Genes <- get_two_sd_status(counts, "Psi5_Genes")
Psi5_Genes_001 <- get_two_sd_status(counts, "Psi5_Genes_001")
Psi5_Genes_01 <- get_two_sd_status(counts, "Psi5_Genes_01")
Psi5_Genes_1e_6 <- get_two_sd_status(counts, "Psi5_Genes_1e_6")

Theta_Genes <- get_two_sd_status(counts, "Theta_Genes")
Theta_Genes_001 <- get_two_sd_status(counts, "Theta_Genes_001")
Theta_Genes_01 <- get_two_sd_status(counts, "Theta_Genes_01")
Theta_Genes_1e_6 <- get_two_sd_status(counts, "Theta_Genes_1e_6")

Jaccard_Genes <- get_two_sd_status(counts, "Jaccard_Genes")
Jaccard_Genes_01 <- get_two_sd_status(counts, "Jaccard_Genes_01")
Jaccard_Genes_001 <- get_two_sd_status(counts, "Jaccard_Genes_001")
Jaccard_Genes_1e_6 <- get_two_sd_status(counts, "Jaccard_Genes_1e_6")

MIG_count <- get_two_sd_status(counts, "no_MIGs_theta_juncs")

joined_filtered <- left_join(All_Genes, Psi3_Genes) %>% left_join(Psi5_Genes) %>% left_join(Theta_Genes) %>% left_join(Jaccard_Genes) %>% 
                left_join(All_Genes_001) %>% left_join(Psi3_Genes_001) %>% left_join(Psi5_Genes_001) %>% left_join(Theta_Genes_001) %>% left_join(Jaccard_Genes_001) %>% 
                left_join(All_Genes_01) %>% left_join(Psi3_Genes_01) %>% left_join(Psi5_Genes_01) %>% left_join(Theta_Genes_01) %>% left_join(Jaccard_Genes_01) %>% 
                left_join(All_Genes_1e_6) %>% left_join(Psi3_Genes_1e_6) %>% left_join(Psi5_Genes_1e_6) %>% left_join(Theta_Genes_1e_6) %>% left_join(Jaccard_Genes_1e_6) %>%                 
                left_join(All_Junctions) %>% left_join(Psi3_Junctions) %>% left_join(Psi5_Junctions) %>% left_join(Theta_Junctions) %>% left_join(Jaccard_Junctions) %>% 
                left_join(All_Junctions_001) %>% left_join(Psi3_Junctions_001) %>% left_join(Psi5_Junctions_001) %>% left_join(Theta_Junctions_001) %>% left_join(Jaccard_Junctions_001) %>% 
                left_join(All_Junctions_01) %>% left_join(Psi3_Junctions_01) %>% left_join(Psi5_Junctions_01) %>% left_join(Theta_Junctions_01) %>% left_join(Jaccard_Junctions_01) %>% 
                left_join(All_Junctions_1e_6) %>% left_join(Psi3_Junctions_1e_6) %>% left_join(Psi5_Junctions_1e_6) %>% left_join(Theta_Junctions_1e_6) %>% left_join(Jaccard_Junctions_1e_6) %>%
                left_join(MIG_count)

write_csv(joined_filtered, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/outlier_status_filtered_with_sibs.csv")


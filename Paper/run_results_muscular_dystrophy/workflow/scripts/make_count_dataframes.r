library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

psi3 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi3_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
psi5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi5_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
theta <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_theta_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_muscle_jaccard_results_jaccard.csv")
output_file <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv"
number_samples <- 185

all <- bind_rows(psi3, psi5, theta, jaccard)

########################################################################
######################-----By Gene per Person-----######################
########################################################################
filter_genes <- function(dataframe, filter, column_name, significance){
      if (filter == "jaccard"){
          dataframe_filtered <- dataframe %>% filter(padjust < significance) %>% filter(abs(deltaPsi) >= 0.3)
      } else {
          dataframe_filtered <- dataframe %>% filter(padjust < significance) %>% filter(abs(deltaPsi) >= 0.3)
      }

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
      dataframe_by_gene <- dataframe_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
      dataframe_count <- dataframe_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      dataframe_number_samples <- dataframe_count %>% select(sampleID) %>% unique %>% nrow()
      print("dataframe number samples:")
      print(dataframe_number_samples)
      missing <- number_samples - as.numeric(dataframe_number_samples)
      print(missing[1])
      
      if (missing[1] == 0) {
            colnames(dataframe_count) <- c("sampleID", column_name)
      } else {
              missing_dataframe_df <- data.frame(sampleID = as.character(1:missing), n = 0)
              dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
              colnames(dataframe_count) <- c("sampleID", column_name)
      }
      dataframe_count
}

filter_junctions <- function(dataframe, filter, column_name, significance){
      if (filter == "jaccard"){
          dataframe_filtered <- dataframe %>% filter(padjust < significance) %>% filter(abs(deltaPsi) >= 0.3)
      } else {
          dataframe_filtered <- dataframe %>% filter(padjust < significance) %>% filter(abs(deltaPsi) >= 0.3)
      }     
      
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

      dataframe_number_samples <- dataframe_count$sampleID %>% unique %>% length()
      missing <- number_samples - as.numeric(dataframe_number_samples)

      if (missing[1] == 0) {
          colnames(dataframe_count) <- c("sampleID", column_name)
      } else {
              missing_dataframe_df <- data.frame(sampleID = as.character(1:missing), n = 0)
              dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
              colnames(dataframe_count) <- c("sampleID", column_name)
      }

      dataframe_count
}

################-----Get Number by Gene per Person-----#################
#0.05
psi3_count_genes <- filter_genes(psi3, "psi3", "Psi3_Genes", 0.05)
psi5_count_genes <- filter_genes(psi5, "psi5", "Psi5_Genes", 0.05)
theta_count_genes <- filter_genes(theta, "theta", "Theta_Genes", 0.05)
all_count_genes <- filter_genes(all, "all", "All_Genes", 0.05)
jaccard_count_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes", 0.05)

#0.01
psi3_count_01_genes <- filter_genes(psi3, "psi3", "Psi3_Genes_01", 0.01)
psi5_count_01_genes <- filter_genes(psi5, "psi5", "Psi5_Genes_01", 0.01)
theta_count_01_genes <- filter_genes(theta, "theta", "Theta_Genes_01", 0.01)
all_count_01_genes <- filter_genes(all, "all", "All_Genes_01", .01)
jaccard_count_01_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_01", 0.01)

#.001
psi3_count_001_genes <- filter_genes(psi3, "psi3", "Psi3_Genes_001", .001)
psi5_count_001_genes <- filter_genes(psi5, "psi5", "Psi5_Genes_001", .001)
theta_count_001_genes <- filter_genes(theta, "theta", "Theta_Genes_001", .001)
all_count_001_genes <- filter_genes(all, "all", "All_Genes_001", .001)
jaccard_count_001_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_001", 0.001)


#1e-6
psi3_count_1e_6_genes <- filter_genes(psi3, "psi3", "Psi3_Genes_1e_6", 1e-6)
psi5_count_1e_6_genes <- filter_genes(psi5, "psi5", "Psi5_Genes_1e_6", 1e-6)
theta_count_1e_6_genes <- filter_genes(theta, "theta", "Theta_Genes_1e_6", 1e-6)
all_count_1e_6_genes <- filter_genes(all, "all", "All_Genes_1e_6", 1e-6)
jaccard_count_1e_6_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_1e_6", 1e-6)

#join
all_count_genes <- left_join(psi3_count_genes, psi5_count_genes) %>% left_join(theta_count_genes) %>%  left_join(all_count_genes) %>% left_join(jaccard_count_genes) %>%
            left_join(psi3_count_01_genes) %>% left_join(psi5_count_01_genes) %>% left_join(theta_count_01_genes) %>% left_join(all_count_01_genes) %>% left_join(jaccard_count_01_genes) %>%
            left_join(psi3_count_001_genes) %>% left_join(psi5_count_001_genes) %>% left_join(theta_count_001_genes) %>% left_join(all_count_001_genes) %>% left_join(jaccard_count_001_genes) %>%
            left_join(psi3_count_1e_6_genes) %>% left_join(psi5_count_1e_6_genes) %>% left_join(theta_count_1e_6_genes) %>% left_join(all_count_1e_6_genes) %>% left_join(jaccard_count_1e_6_genes)

all_count_genes$Psi3_Psi5_Genes <- all_count_genes$Psi3_Genes + all_count_genes$Psi5_Genes
all_count_genes$Psi3_Psi5_Genes_01 <- all_count_genes$Psi3_Genes_01 + all_count_genes$Psi5_Genes_01
all_count_genes$Psi3_Psi5_Genes_001 <- all_count_genes$Psi3_Genes_001 + all_count_genes$Psi5_Genes_001
all_count_genes$Psi3_Psi5_Genes_1e_6 <- all_count_genes$Psi3_Genes_1e_6 + all_count_genes$Psi5_Genes_1e_6

#
################-----Get Number by Junction per Person-----#################
psi3_count <- filter_junctions(psi3, "psi3", "Psi3_Junctions", 0.05)
psi5_count <- filter_junctions(psi5, "psi5", "Psi5_Junctions", 0.05)
theta_count <- filter_junctions(theta, "theta", "Theta_Junctions", 0.05)
all_count <- filter_junctions(all, "all", "All_Junctions", 0.05)
jaccard_count <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions", 0.05)

psi3_count_001 <- filter_junctions(psi3, "psi3", "Psi3_Junctions_001", .001)
psi5_count_001 <- filter_junctions(psi5, "psi5", "Psi5_Junctions_001", .001)
theta_count_001 <- filter_junctions(theta, "theta", "Theta_Junctions_001", .001)
all_count_001 <- filter_junctions(all, "all", "All_Junctions_001", .001)
jaccard_count_001 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_001", 0.001)

psi3_count_01 <- filter_junctions(psi3, "psi3", "Psi3_Junctions_01", .01)
psi5_count_01 <- filter_junctions(psi5, "psi5", "Psi5_Junctions_01", .01)
theta_count_01 <- filter_junctions(theta, "theta", "Theta_Junctions_01", .01)
all_count_01 <- filter_junctions(all, "all", "All_Junctions_01", .01)
jaccard_count_01 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_01", 0.01)

psi3_count_1e_6 <- filter_junctions(psi3, "psi3", "Psi3_Junctions_1e_6", 1e-6)
psi5_count_1e_6 <- filter_junctions(psi5, "psi5", "Psi5_Junctions_1e_6", 1e-6)
theta_count_1e_6 <- filter_junctions(theta, "theta", "Theta_Junctions_1e_6", 1e-6)
all_count_1e_6 <- filter_junctions(all, "all", "All_Junctions_1e_6", 1e-6)
jaccard_count_1e_6 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_1e_6", 1e-6)


all_count_junctions <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(all_count) %>% left_join(jaccard_count) %>% 
                        left_join(psi3_count_01) %>% left_join(psi5_count_01) %>% left_join(theta_count_01) %>% left_join(all_count_01) %>% left_join(jaccard_count_01) %>% 
                        left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% left_join(theta_count_001) %>% left_join(all_count_001) %>% left_join(jaccard_count_001) %>% 
                        left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% left_join(all_count_1e_6) %>% left_join(jaccard_count_1e_6)

all_count_junctions$Psi3_Psi5_Junctions <- all_count_junctions$Psi3_Junctions + all_count_junctions$Psi5_Junctions
all_count_junctions$Psi3_Psi5_Junctions_01 <- all_count_junctions$Psi3_Junctions_01 + all_count_junctions$Psi5_Junctions_01
all_count_junctions$Psi3_Psi5_Junctions_001 <- all_count_junctions$Psi3_Junctions_001 + all_count_junctions$Psi5_Junctions_001
all_count_junctions$Psi3_Psi5_Junctions_1e_6 <- all_count_junctions$Psi3_Junctions_1e_6 + all_count_junctions$Psi5_Junctions_1e_6

################-----Combined Junctions to Genes-----#################
all_count <- left_join(all_count_junctions, all_count_genes)
all_count[is.na(all_count)] <- 0

write_csv(all_count, output_file)

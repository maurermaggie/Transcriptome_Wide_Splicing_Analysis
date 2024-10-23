library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled_fp <- args[1]
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files_fp <- args[2]
files_jaccard_fp <- args[3]
low_RIN_fp <- args[4]
MIGs_fp <- args[5]
jaccard_fp <- args[6]
output_file <- args[7]
missing_metadata_fp <- args[8]
missing_metadata <- read_csv(missing_metadata_fp)

all_uncompiled <- read_csv(all_uncompiled_fp)
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv(files_fp)
files_jaccard <- read_csv(files_jaccard_fp)
low_RIN <- read_csv(low_RIN_fp)
MIGs <- read_csv(MIGs_fp)
jaccard <- read_csv(jaccard_fp)

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
            if (filter == "jaccard"){
              missing_dataframe <- setdiff(files_jaccard$sampleID, dataframe_count$sampleID)
      } else {
              missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)

      }

      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
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

      if (filter == "jaccard"){
              missing_dataframe <- setdiff(files_jaccard$sampleID, dataframe_count$sampleID)
      } else {
              missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)

      }
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

################-----Get Number by Gene per Person-----#################
#0.05
psi3_count_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes", 0.05)
psi5_count_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes", 0.05)
theta_count_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes", 0.05)
all_count_genes <- filter_genes(all_uncompiled, "all", "All_Genes", 0.05)
jaccard_count_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes", 0.05)

#0.01
psi3_count_01_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_01", 0.01)
psi5_count_01_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_01", 0.01)
theta_count_01_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_01", 0.01)
all_count_01_genes <- filter_genes(all_uncompiled, "all", "All_Genes_01", .01)
jaccard_count_01_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_01", 0.01)

#.001
psi3_count_001_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_001", .001)
psi5_count_001_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_001", .001)
theta_count_001_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_001", .001)
all_count_001_genes <- filter_genes(all_uncompiled, "all", "All_Genes_001", .001)
jaccard_count_001_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_001", 0.001)


#1e-6
psi3_count_1e_6_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_1e_6", 1e-6)
psi5_count_1e_6_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_1e_6", 1e-6)
theta_count_1e_6_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_1e_6", 1e-6)
all_count_1e_6_genes <- filter_genes(all_uncompiled, "all", "All_Genes_1e_6", 1e-6)
jaccard_count_1e_6_genes <- filter_genes(jaccard, "jaccard", "Jaccard_Genes_1e_6", 1e-6)

#join
all_count_genes <- left_join(psi3_count_genes, psi5_count_genes) %>% left_join(theta_count_genes) %>%  left_join(all_count_genes) %>% left_join(jaccard_count_genes) %>%
            left_join(psi3_count_01_genes) %>% left_join(psi5_count_01_genes) %>% left_join(theta_count_01_genes) %>% left_join(all_count_01_genes) %>% left_join(jaccard_count_01_genes) %>%
            left_join(psi3_count_001_genes) %>% left_join(psi5_count_001_genes) %>% left_join(theta_count_001_genes) %>% left_join(all_count_001_genes) %>% left_join(jaccard_count_001_genes) %>%
            left_join(psi3_count_1e_6_genes) %>% left_join(psi5_count_1e_6_genes) %>% left_join(theta_count_1e_6_genes) %>% left_join(all_count_1e_6_genes) %>% left_join(jaccard_count_1e_6_genes)

#
################-----Get Number by Junction per Person-----#################
psi3_count <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions", 0.05)
psi5_count <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions", 0.05)
theta_count <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions", 0.05)
all_count <- filter_junctions(all_uncompiled, "all", "All_Junctions", 0.05)
jaccard_count <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions", 0.05)

psi3_count_001 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_001", .001)
psi5_count_001 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_001", .001)
theta_count_001 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_001", .001)
all_count_001 <- filter_junctions(all_uncompiled, "all", "All_Junctions_001", .001)
jaccard_count_001 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_001", 0.001)

psi3_count_01 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_01", .01)
psi5_count_01 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_01", .01)
theta_count_01 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_01", .01)
all_count_01 <- filter_junctions(all_uncompiled, "all", "All_Junctions_01", .01)
jaccard_count_01 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_01", 0.01)

psi3_count_1e_6 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_1e_6", 1e-6)
psi5_count_1e_6 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_1e_6", 1e-6)
theta_count_1e_6 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_1e_6", 1e-6)
all_count_1e_6 <- filter_junctions(all_uncompiled, "all", "All_Junctions_1e_6", 1e-6)
jaccard_count_1e_6 <- filter_junctions(jaccard, "jaccard", "Jaccard_Junctions_1e_6", 1e-6)


all_count_junctions <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(all_count) %>% left_join(jaccard_count) %>% 
                        left_join(psi3_count_01) %>% left_join(psi5_count_01) %>% left_join(theta_count_01) %>% left_join(all_count_01) %>% left_join(jaccard_count_01) %>% 
                        left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% left_join(theta_count_001) %>% left_join(all_count_001) %>% left_join(jaccard_count_001) %>% 
                        left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% left_join(all_count_1e_6) %>% left_join(jaccard_count_1e_6)

################-----Combined Junctions to Genes-----#################
all_count <- left_join(all_count_junctions, all_count_genes)

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG")

#get theta dataframe
results_filtered_theta <- all_uncompiled %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol")
colnames(results_filtered_theta) <- c("sampleID", "gene_symbol")

#combine dataframes
sig_genes_type_all_samples <- left_join(MIG_table_select, results_filtered_theta) %>% filter(!is.na(sampleID))

joined <- sig_genes_type_all_samples %>% select(sampleID, gene_symbol) %>% 
    group_by(sampleID) %>% tally() %>% arrange(desc(n))

missing_all <- setdiff(files$sampleID, joined$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
MIG_count <- bind_rows(joined, missing_all_df)
MIG_count <- MIG_count %>% filter(! sampleID %in% low_RIN$sampleID)

colnames(MIG_count) <- c("sampleID", "no_theta_juncs_in_MIGs")

all_count <- left_join(all_count, MIG_count)

########################################################################
##########-----Get Significant MIG Jaccard Index Events-----#########
########################################################################

#get jaccard dataframe
results_filtered_jaccard <- jaccard %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "jaccard")
results_filtered_jaccard <- results_filtered_jaccard %>% select("sampleID", "hgncSymbol")
colnames(results_filtered_jaccard) <- c("sampleID", "gene_symbol")

#combine dataframes
sig_genes_type_all_samples_jaccard <- left_join(MIG_table_select, results_filtered_jaccard) %>% filter(!is.na(sampleID))

joined_jaccard <- sig_genes_type_all_samples_jaccard %>% select(sampleID, gene_symbol) %>% 
    group_by(sampleID) %>% tally() %>% arrange(desc(n))

missing_all_jaccard <- setdiff(files_jaccard$sampleID, joined_jaccard$sampleID)
missing_all_jaccard_df <- data.frame(sampleID = missing_all_jaccard, n = 0)
MIG_count_jaccard <- bind_rows(joined_jaccard, missing_all_jaccard_df)
MIG_count_jaccard <- MIG_count_jaccard %>% filter(! sampleID %in% low_RIN$sampleID)

colnames(MIG_count_jaccard) <- c("sampleID", "no_jaccard_juncs_in_MIGs")
all_count <- left_join(all_count, MIG_count_jaccard)

########################################################################
################-----Get Significant MIGs with Theta-----###############
########################################################################
#get theta dataframe
results_filtered_theta <- all_uncompiled %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% 
                             filter(type == "theta")
results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol") %>% filter(hgncSymbol %in% MIG_table_select$gene_symbol) %>% unique()
colnames(results_filtered_theta) <- c("sampleID", "gene_symbol")
results_filtered_MIGs <- results_filtered_theta %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
colnames(results_filtered_MIGs) <- c("sampleID", "no_MIGs_with_theta_juncs")

all_count <- left_join(all_count, results_filtered_MIGs)

all_count <- all_count %>% filter(! sampleID %in% missing_metadata$sampleID)
all_count[is.na(all_count)] <- 0

write_csv(all_count, output_file)

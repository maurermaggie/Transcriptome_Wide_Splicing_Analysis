require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################
filter_genes <- function(dataframe, filter, column_name, significance){
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
      dataframe_by_gene <- dataframe_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
      dataframe_count <- dataframe_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

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

################-----Get Number by Junction per Person-----#################
psi3_count <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions", 0.05)
psi5_count <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions", 0.05)
theta_count <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions", 0.05)
all_count <- filter_junctions(all_uncompiled, "all", "All_Junctions", 0.05)

psi3_count_001 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_001", .001)
psi5_count_001 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_001", .001)
theta_count_001 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_001", .001)
all_count_001 <- filter_junctions(all_uncompiled, "all", "All_Junctions_001", .001)

psi3_count_01 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_01", .01)
psi5_count_01 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_01", .01)
theta_count_01 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_01", .01)
all_count_01 <- filter_junctions(all_uncompiled, "all", "All_Junctions_01", .01)

psi3_count_1e_6 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_1e_6", 1e-6)
psi5_count_1e_6 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_1e_6", 1e-6)
theta_count_1e_6 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_1e_6", 1e-6)
all_count_1e_6 <- filter_junctions(all_uncompiled, "all", "All_Junctions_1e_6", 1e-6)

all_count_junctions <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(all_count) %>% left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% 
              left_join(theta_count_001) %>% left_join(all_count_001) %>% left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% 
              left_join(all_count_1e_6) %>% left_join(psi3_count_01) %>% left_join(psi5_count_01) %>% left_join(theta_count_01) %>% left_join(all_count_01)

################-----Get Number by Genes per Person-----#################
psi3_count_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes", 0.05)
psi5_count_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes", 0.05)
theta_count_genes<- filter_genes(all_uncompiled, "theta", "Theta_Genes", 0.05)
all_count_genes <- filter_genes(all_uncompiled, "all", "All_Genes", 0.05)

psi3_count_001_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_001", .001)
psi5_count_001_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_001", .001)
theta_count_001_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_001", .001)
all_count_001_genes <- filter_genes(all_uncompiled, "all", "All_Genes_001", .001)

psi3_count_01_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_01", .01)
psi5_count_01_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_01", .01)
theta_count_01_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_01", .01)
all_count_01_genes <- filter_genes(all_uncompiled, "all", "All_Genes_01", .01)

psi3_count_1e_6_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_1e_6", 1e-6)
psi5_count_1e_6_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_1e_6", 1e-6)
theta_count_1e_6_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_1e_6", 1e-6)
all_count_1e_6_genes <- filter_genes(all_uncompiled, "all", "All_Genes_1e_6", 1e-6)

all_count_genes <- left_join(psi3_count_genes, psi5_count_genes) %>% left_join(theta_count_genes) %>% left_join(all_count_genes) %>% left_join(psi3_count_001_genes) %>% left_join(psi5_count_001_genes) %>% 
              left_join(theta_count_001_genes) %>% left_join(all_count_001_genes) %>% left_join(psi3_count_1e_6_genes) %>% left_join(psi5_count_1e_6_genes) %>% left_join(theta_count_1e_6_genes) %>% 
              left_join(all_count_1e_6_genes) %>% left_join(psi3_count_01_genes) %>% left_join(psi5_count_01_genes) %>% left_join(theta_count_01_genes) %>% left_join(all_count_01_genes)

################-----Combined Junctions to Genes-----#################
all_count <- left_join(all_count_junctions, all_count_genes)

################-----Join to Metadata-----#################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch)

names(all_count)[names(all_count) == 'sampleID'] <- 'RDID'
all_filter_RD <- filter(all_count, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)

names(all_count)[names(all_count) == 'RDID'] <- 'GSS_ID'
all_filter_GSS <- filter(all_count, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)

names(all_count)[names(all_count) == 'GSS_ID'] <- 'UDN_ID'
all_filter_UDN <- filter(all_count, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)

names(all_count)[names(all_count) == 'UDN_ID'] <- 'sampleID'

joined_filtered <- bind_rows(all_RDID, all_GSS, all_UDN)

########################################################################
#######################-----Linear Model-----###########################
########################################################################
make_estimate_df <- function(row, dataframe) {
    Esitmate_lm_df <- lm(dataframe[[row]] ~ RIN + age + sex+ batch, data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
    Value <- summary(Esitmate_lm_df)$coefficients[,0]
    df$Value <- c("Intercept", "RIN", "age", "sexM", "batch")
    df$Numbers <- ifelse(df$Value == "Intercept", 5,
                            ifelse(df$Value == "RIN", 4,
                            ifelse(df$Value == "age", 3,
                            ifelse(df$Value == "sexM", 2, 1))))
   df
}

################-----Linear Model Genes-----#################
#.005
All_05_Genes <- make_estimate_df("All_Genes", joined_filtered)
Psi3_05_Genes <- make_estimate_df("Psi3_Genes", joined_filtered)
Psi5_05_Genes <- make_estimate_df("Psi5_Genes", joined_filtered)
Theta_05_Genes <- make_estimate_df("Theta_Genes", joined_filtered)


################-----Linear Model Junctions-----#################
#.005
All_05_Junctions <- make_estimate_df("All_Junctions", joined_filtered)
Psi3_05_Junctions <- make_estimate_df("Psi3_Junctions", joined_filtered)
Psi5_05_Junctions <- make_estimate_df("Psi5_Junctions", joined_filtered)
Theta_05_Junctions <- make_estimate_df("Theta_Junctions", joined_filtered)

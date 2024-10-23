library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(broom)
library(FRASER)
library("org.Hs.eg.db")
library(tidyverse)

psi3 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi3_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
psi5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi5_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
theta <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_theta_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_muscle_jaccard_results_jaccard.csv")
output_file <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/GSEA_analysis.pdf"
all_count <- bind_rows(psi3, psi5, theta, jaccard)

count <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv")

FSHD_genes <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/FSHD_genes.xlsx")

FSHD2_genes <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/FSHD_2.csv") %>% filter(`DTU adjusted p-value` < 0.05)
FSHD3_genes <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/FSHD3.csv") %>% filter(`DTU adjusted p-value` < 0.05)
Expression <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/expression_pos_neg.csv")
FSHD4_genes <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/FSHD_4.xlsx") %>% filter(`adj. p-value` < 0.05)
count %>% arrange(desc(All_Junctions)) %>% select(sampleID, All_Junctions)
FSHD5_genes <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/splicing_schatzl.xlsx")
################-----Get ENSG IDs-----#################
Genes <- FSHD_genes %>% filter(event == "RI") %>% pull(geneSym)
ensg_psi <- FSHD_genes %>% filter(event %in% c("3ss", "5ss")) %>% pull(geneID)
ensg_RI <- FSHD_genes %>% filter(event %in% c("RI", "SE", "MXE")) %>% pull(geneID)

################-----Linear regression of # juncs in those genes-----#################
FSHD_RI <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("theta")) %>% filter(hgncSymbol %in% ensg_RI) %>% group_by(sampleID) %>% tally() 
colnames(FSHD_RI) <- c("sampleID", "RI")
FSHD_psi <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg_psi) %>% group_by(sampleID) %>% tally() 
colnames(FSHD_psi) <- c("sampleID", "PSI")
FSHD <- left_join(FSHD_RI, FSHD_psi)
FSHD$n <- FSHD$RI + FSHD$PSI

non_FSHD_RI <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("theta")) %>% filter(hgncSymbol %in% ensg_RI) %>% group_by(sampleID) %>% tally() 
colnames(non_FSHD_RI) <- c("sampleID", "RI")
non_FSHD_psi <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg_psi) %>% group_by(sampleID) %>% tally() 
colnames(non_FSHD_psi) <- c("sampleID", "PSI")
non_FSHD <- left_join(non_FSHD_RI, non_FSHD_psi)
non_FSHD$n <- non_FSHD$RI + non_FSHD$PSI

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg <- FSHD2_genes %>% pull(`Ensembl ID`)

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg <- FSHD3_genes %>% pull(`Ensembl ID`)

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg1 <- FSHD3_genes %>% pull(`Ensembl ID`)
ensg2 <- FSHD2_genes %>% pull(`Ensembl ID`)
ensg <- intersect(ensg1, ensg2)
################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df


################-----Get ENSG IDs-----#################
ensg1 <- FSHD3_genes %>% pull(`Ensembl ID`)
ensg2 <- FSHD2_genes %>% pull(`Ensembl ID`)
ensg <- intersect(ensg1, ensg2)
################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg <- Expression %>% pull(Ensemble_positive)

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg <- Expression %>% pull(Ensemble_negative)

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df
################-----Get ENSG IDs-----#################
ensg <- FSHD4_genes %>% pull(`ENSEMBL ID`) 

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

################-----Get ENSG IDs-----#################
ensg <- FSHD5_genes %>% pull(`ENSEMBL ID`) %>% str_replace_all("gene:", "")

################-----Linear regression of # juncs in those genes-----#################
FSHD <- all_count %>% filter(sampleID == "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_FSHD <- all_count %>% filter(sampleID != "1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5", "theta")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 

FSHD$FSHD <- 1
non_FSHD$FSHD <- 0

all <- rbind(FSHD, non_FSHD)

Esitmate_lm_df <- lm(n ~ FSHD, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "FSHD")
df

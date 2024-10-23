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

################-----Get ENSG IDs-----#################
Genes <- c("CLCN1", "CACNA1S", "RYR1", "ATP2A1", "INSR", "MBNL1", "MBNL2", "TNNT2", "DMD", "LDB3", "BIN1", "SOS1", "ANK2", "PHKA1", "KIF13A", "ARHGEF7")


ensg = mapIds(org.Hs.eg.db,
                    keys=Genes, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")

################-----Linear regression of # juncs in those genes-----#################
number_DMPK_junctions_psi3_psi5 <- all_count %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))

number_DMPK_junctions_psi3_psi5_reads <-  all_count %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally(meanTotalCounts) %>% arrange(desc(n))

DMPK <- all_count %>% filter(sampleID == "BEG_1230-1_T1227") %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
DMPK$DMPK <- 1
non_DMPK <- all_count %>% filter(sampleID != "BEG_1230-1_T1227") %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() 
non_DMPK$DMPK <- 0

all <- rbind(DMPK, non_DMPK)
missing <- setdiff(count$sampleID, all$sampleID)

missing_df <- data.frame(sampleID= missing, n=0)

all <- bind_rows(all, missing_df)

Esitmate_lm_df <- lm(n ~ DMPK, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "DMPK")
df



################-----Linear regression of # meanTotalCount in those genes-----#################
DMPK <- all_count %>% filter(sampleID == "BEG_1230-1_T1227") %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally(totalCounts) 
DMPK$DMPK <- 1
non_DMPK <- all_count %>% filter(sampleID != "BEG_1230-1_T1227") %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally(totalCounts) 
non_DMPK$DMPK <- 0

all <- rbind(DMPK, non_DMPK)

Esitmate_lm_df <- lm(n ~ DMPK, data = all)
df <- data.frame(matrix(nrow=2, ncol=0))
df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
total <- sum(df$Absolute_Effect_Estimate)
df$variance_explained <- df$Absolute_Effect_Estimate/total
df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
Value <- summary(Esitmate_lm_df)$coefficients[,0]
df$Value <- c("Intercept", "DMPK")
df

################-----Contingency Table of the DMPK vs FSHD cases-----#################
DMPKlist_DMPKcases <- filter(all_count, sampleID=="BEG_1230-1_T1227") %>% filter(padjust < 0.05) %>% filter(type %in% c("psi3", "psi5")) %>% filter(abs(deltaPsi) > 0.3) %>% filter(hgncSymbol %in% ensg) %>% nrow
DMPKlist_FSHDcases <- filter(all_count, sampleID=="1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% nrow
DMPKnonlist_DMPKcases <- filter(all_count, sampleID=="BEG_1230-1_T1227") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(!hgncSymbol %in% ensg) %>% nrow
DMPKnonlist_FSHDcases <- filter(all_count, sampleID=="1258-1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(!hgncSymbol %in% ensg) %>% nrow


DMPK_ct <- data.frame(
            "DMPK_case" = c(DMPKlist_DMPKcases, DMPKnonlist_DMPKcases),
            "FSHD_case" = c(DMPKlist_FSHDcases, DMPKnonlist_FSHDcases),
            row.names = c("DMPK_list", "DMPK_non_list"),
            stringsAsFactors = FALSE
        )

colnames(DMPK_ct) <- c("DMPK_case", "FSHD_case")
            
DMPK_f <- fisher.test(DMPK_ct)
DMPK_tidy <- tidy(DMPK_f)

################-----Contingency Table of the DMPK vs FSHD cases-----#################
DMPKlist_DMPKcases <- filter(all_count, sampleID=="BEG_1230-1_T1227") %>% filter(padjust < 0.05) %>% filter(type %in% c("psi3", "psi5")) %>% filter(abs(deltaPsi) > 0.3) %>% filter(hgncSymbol %in% ensg) %>% nrow
DMPKlist_MANcases <- filter(all_count, sampleID=="MAN_1438-01-M1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% nrow
DMPKnonlist_DMPKcases <- filter(all_count, sampleID=="BEG_1230-1_T1227") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(!hgncSymbol %in% ensg) %>% nrow
DMPKnonlist_MANcases <- filter(all_count, sampleID=="MAN_1438-01-M1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(type %in% c("psi3", "psi5")) %>% filter(!hgncSymbol %in% ensg) %>% nrow

DMPK_ct <- data.frame(
            "DMPK_case" = c(DMPKlist_DMPKcases, DMPKnonlist_DMPKcases),
            "Non_DMPK_case" = c(DMPKlist_MANcases, DMPKnonlist_MANcases),
            row.names = c("DMPK_list", "DMPK_non_list"),
            stringsAsFactors = FALSE
        )

colnames(DMPK_ct) <- c("DMPK_case", "FSHD_case")
            
DMPK_f <- fisher.test(DMPK_ct)
DMPK_tidy <- tidy(DMPK_f)














































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































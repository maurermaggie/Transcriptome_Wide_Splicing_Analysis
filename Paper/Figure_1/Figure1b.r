library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(GO.db)
library(tidyverse)
library(RColorBrewer)

display.brewer.all
pal <- brewer.pal(10, "Paired")

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

################-----Get Number by Gene per Person-----#################
all_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
all_by_gene <- all_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
all_count <- all_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_all <- setdiff(files$sampleID, all_count$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
all_count <- bind_rows(all_count, missing_all_df)
colnames(all_count) <- c("sampleID", "Combined")

psi3_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_by_gene <- psi3_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
psi3_count <- psi3_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi3 <- setdiff(files$sampleID, psi3_count$sampleID)
missing_psi3_df <- data.frame(sampleID = missing_psi3, n = 0)
psi3_count <- bind_rows(psi3_count, missing_psi3_df)
colnames(psi3_count) <- c("sampleID", "Psi3")

psi5_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_by_gene <- psi5_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
psi5_count <- psi5_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi5 <- setdiff(files$sampleID, psi5_count$sampleID)
missing_psi5_df <- data.frame(sampleID = missing_psi5, n = 0)
psi5_count <- bind_rows(psi5_count, missing_psi5_df)
colnames(psi5_count) <- c("sampleID", "Psi5")

theta_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_by_gene <- theta_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
theta_count <- theta_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_theta <- setdiff(files$sampleID, theta_count$sampleID)
missing_theta_df <- data.frame(sampleID = missing_theta, n = 0)
theta_count <- bind_rows(theta_count, missing_theta_df)
colnames(theta_count) <- c("sampleID", "Theta")

joined <- left_join(all_count, psi5_count) %>% left_join(psi3_count) %>% left_join(theta_count)
joined_long <- melt(setDT(joined), id.vars = c("sampleID"), variable.name = "genes_with_aberrant_junction_of_type")

################-----Violin Plot-----#################
p <- ggplot(joined_long, aes(x=genes_with_aberrant_junction_of_type, y=value)) + 
  geom_violin()

################-----Box Plot-----#################
title_lab <- paste("Number of Genes with Aberrant Splcie Junctions", "\n", "of each FRASER Type", sep="")

psi3 <- expression(psi ~ 3)
psi5 <- expression(psi ~ 5)
theta <- expression(theta)

FRASER_boxplot <- ggplot(joined_long, aes(fill = genes_with_aberrant_junction_of_type ,x=genes_with_aberrant_junction_of_type, y=value)) + 
  geom_boxplot(show.legend = FALSE) +
  scale_x_discrete(labels = c('Psi3' = psi3,
                                      'Theta' = theta,
                                          'Psi5' = psi5))+
  xlab("Splicing Aberration Type") +
  ylab("Number of Genes with Aberrant Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(10,5,2,4)])

FRASER_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Figure1b.pdf", plot=FRASER_boxplot,  limitsize = FALSE, units = "in")

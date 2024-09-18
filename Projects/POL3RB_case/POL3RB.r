library(tidyverse)
library(readxl)
require(data.table)
library(magrittr)
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

GREGoR_expression <- fread("/oak/stanford/groups/smontgom/gss_prospective_rnaseq/aggregate/expression/eoutliers_withUDN_blood.txt.gz")
GTEx_expression <- read_tsv("/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv")
FRASER_metadata_important <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_metadata_combined.csv")

################-----GREGoR Expression-----#################
mean_GREGoR <- mean(GREGoR_expression$zscore, na.rm=TRUE)
sd_GREGoR <- sd(GREGoR_expression$zscore, na.rm=TRUE)
GREGoR_expression$pvalue <- pnorm(GREGoR_expression$zscore, mean = mean_GREGoR, sd = sd_GREGoR, lower.tail=TRUE)
GREGoR_expression$diffexpressed <- ifelse(GREGoR_expression$zscore >=2, "UP",
                                        ifelse(GREGoR_expression$zscore <= -2, "DOWN", "NO"))                                   

filter(GREGoR_expression, gene == "ENSG00000013503") 
filter(GREGoR_expression, sample == "RD143")
filter(GREGoR_expression, sample == "RD143") %>% filter(gene == "ENSG00000013503")

filter(GREGoR_expression, gene == "ENSG00000013503") %>% arrange(zscore)

RD143 <- GREGoR_expression %>% filter(sample == "RD143")
POL3RB <- RD143 %>% filter(gene == "ENSG00000013503")
POL3RB$delabel <- POL3RB$gene
rest <- RD143 %>% filter(gene != "ENSG00000013503")
rest$delabel <- NA
RD143 <- rbind(rest, POL3RB)    

RD143_plot <- ggplot(data = RD143, aes(x = zscore, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)<br />
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('Outliers in RD143') + # Plot title
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  geom_label_repel()

RD143_plot


POL3RB <- GREGoR_expression %>% filter(gene == "ENSG00000013503")
g1 <- subset(POL3RB, sample == "RD143")

function_POL3RB_GREGoR <- ggplot(POL3RB, aes(x = reorder(sample, zscore), y = zscore)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab("Z-score") +
  labs(title="Zscore distribution of POL3RB expression")  +
  geom_point(data=g1, aes(colour= sample)) +
  labs(colour="Sample_ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_POL3RB_GREGoR


################-----GTEx Expression-----#################
mean_GTEx <- mean(GTEx_expression$normalized_resid, na.rm=TRUE)
sd_GTEx <- sd(GTEx_expression$normalized_resid, na.rm=TRUE)
GTEx_expression$pvalue <- pnorm(GTEx_expression$normalized_resid, mean = mean_GREGoR, sd = sd_GREGoR, lower.tail=TRUE)
GTEx_expression$diffexpressed <- ifelse(GTEx_expression$normalized_resid >=2, "UP",
                                        ifelse(GTEx_expression$normalized_resid <= -2, "DOWN", "NO"))                                  
GTEx_expression$affected_status <- ifelse(grepl("GTEX",GTEx_expression$Ind),"Control","Case")

FRASER_affected_status <- FRASER_metadata_important %>% select(UDN_ID, affected_status)
FRASER_case <- FRASER_affected_status %>% filter(affected_status == "Case") %>% pull(UDN_ID)
FRASER_control <- FRASER_affected_status %>% filter(affected_status == "Control") %>% pull(UDN_ID)

setdiff(GTEx_expression$Ind, FRASER_affected_status$UDN_ID)
GTEx_expression$affected_status <- ifelse(GTEx_expression$Ind %in% FRASER_case,"Case",
                                        ifelse(GTEx_expression$Ind %in% FRASER_control, "Control", GTEx_expression$affected_status))


filter(GTEx_expression, gene == "ENSG00000013503") 
filter(GTEx_expression, Ind == "UDN716672")

filter(GTEx_expression, Ind == "UDN716672") %>% filter(gene == "ENSG00000013503")
filter(GTEx_expression, gene == "ENSG00000013503") %>% arrange(normalized_resid)

RD143_GTEx <- GTEx_expression %>% filter(Ind == "UDN716672")
POL3RB_GTEx <- RD143_GTEx %>% filter(gene == "ENSG00000013503")
POL3RB_GTEx$delabel <- POL3RB_GTEx$gene
rest <- RD143_GTEx %>% filter(gene != "ENSG00000013503")
rest$delabel <- NA
RD143_GTEx <- rbind(rest, POL3RB_GTEx)    

RD143_GTEx_plot <- ggplot(data = RD143_GTEx, aes(x = normalized_resid, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)<br />
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('Outliers in RD143') + # Plot title
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  geom_label_repel()

RD143_GTEx_plot

POL3RB_GTEx <- GTEx_expression %>% filter(gene == "ENSG00000013503")
g1 <- subset(POL3RB_GTEx, Ind == "UDN716672")

function_POL3RB_GTEx <- ggplot(POL3RB_GTEx, aes(x = reorder(Ind, normalized_resid), y = normalized_resid)) +
  geom_point(aes(color = affected_status)) +
  xlab("Sample_ID") +
  ylab("Normalized residual") +
  labs(title="Normalized residual distribution of POL3RB expression")  +
  geom_point(data=g1, aes(colour= Ind)) +
  geom_label(data=g1, aes(label = "UDN716672"), nudge_x=100) +
  labs(colour="Sample_ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_POL3RB_GTEx
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/POL3RB_case/POL3RB_GTEx.pdf", plot=function_POL3RB_GTEx,  limitsize = FALSE, units = "in")

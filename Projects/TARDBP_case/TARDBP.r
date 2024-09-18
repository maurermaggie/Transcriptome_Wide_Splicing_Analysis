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

filter(GREGoR_expression, gene == "ENSG00000120948") 
filter(GREGoR_expression, sample == "RD465")
filter(GREGoR_expression, sample == "RD465") %>% filter(gene == "ENSG00000120948")

filter(GREGoR_expression, gene == "ENSG00000120948") %>% arrange(zscore)
filter(GREGoR_expression, sample == "RD465") %>% arrange(zscore)


RD465 <- GREGoR_expression %>% filter(sample == "RD465")
TARDBP <- RD465 %>% filter(gene == "ENSG00000120948")
TARDBP$delabel <- TARDBP$gene
rest <- RD465 %>% filter(gene != "ENSG00000120948")
rest$delabel <- NA
RD465 <- rbind(rest, TARDBP)    

RD465_plot <- ggplot(data = RD465, aes(x = zscore, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)<br />
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('Outliers in RD465') + # Plot title
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  geom_label_repel()

RD465_plot


TARDBP <- GREGoR_expression %>% filter(gene == "ENSG00000120948")
g1 <- subset(TARDBP, sample == "RD465")

function_TARDBP_GREGoR <- ggplot(TARDBP, aes(x = reorder(sample, zscore), y = zscore)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab("Z-score") +
  labs(title="Zscore distribution of TARDBP expression")  +
  geom_point(data=g1, aes(colour= sample)) +
  labs(colour="Sample_ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_TARDBP_GREGoR


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


filter(GTEx_expression, gene == "ENSG00000120948") 
filter(GTEx_expression, Ind == "UDN955248")

filter(GTEx_expression, Ind == "UDN955248") %>% filter(gene == "ENSG00000120948")
filter(GTEx_expression, gene == "ENSG00000120948") %>% arrange(normalized_resid)

RD143_GTEx <- GTEx_expression %>% filter(Ind == "UDN955248")
TARDBP_GTEx <- RD143_GTEx %>% filter(gene == "ENSG00000120948")
TARDBP_GTEx$delabel <- TARDBP_GTEx$gene
rest <- RD143_GTEx %>% filter(gene != "ENSG00000120948")
rest$delabel <- NA
RD143_GTEx <- rbind(rest, TARDBP_GTEx)    

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

TARDBP_GTEx <- GTEx_expression %>% filter(gene == "ENSG00000120948")
g1 <- subset(TARDBP_GTEx, Ind == "UDN955248")

function_TARDBP_GTEx <- ggplot(TARDBP_GTEx, aes(x = reorder(Ind, normalized_resid), y = normalized_resid)) +
  geom_point(aes(color = affected_status)) +
  xlab("Sample_ID") +
  ylab("Normalized residual") +
  labs(title="Normalized residual distribution of TARDBP expression")  +
  geom_point(data=g1, aes(colour= Ind)) +
  geom_label(data=g1, aes(label = "UDN955248"), nudge_x=100) +
  labs(colour="Sample_ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_TARDBP_GTEx
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/TARDBP_case/TARDBP_GTEx.pdf", plot=function_TARDBP_GTEx,  limitsize = FALSE, units = "in")

################-----GTEx Expression-----#################
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- GTEx_expression$gene
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
merge(GTEx_expression,G_list,by.x="gene",by.y="ensembl_gene_id")
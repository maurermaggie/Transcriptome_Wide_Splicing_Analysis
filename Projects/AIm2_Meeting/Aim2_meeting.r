#library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library (biomaRt)
library("org.Hs.eg.db")
library(R.utils)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
fibro <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/fibroblasts/outputs_6_26/compiled_counts/FRASER_results_raw.csv")

filepath_blood_MIGs <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/blood_MIGs.pdf"
filepath_blood <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/blood_intron_retention.pdf"
filepath_blood_sig <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/blood_intron_retention_highly_sig_pval_psi.pdf"
filepath_volcano_MIG <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/Volcano_MIG.pdf"
constraint_plot_filepath <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/constraint_plot.pdf"

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
transcriptomics <- read_tsv("/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
#gunzip("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/gnomad.v2.1.1.lof_metrics.by_gene.txt")
constraint <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/AIm2_Meeting/gnomad.v2.1.1.lof_metrics.by_gene.txt")

################-----Number Theta MIG Blood-----#################
y_lab <- "Number of Minor Intron Retention Events"
title_lab <- paste("Number of Minor Intron Retention Events per Sample in the", "\n", "UDN and GREGoR Stanford Sites", sep=" ")

top_10 <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(Sample_ID, no_MIGs_theta_juncs) %>% head(2) %>% pull(Sample_ID)
g1 <- subset(joined, Sample_ID %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood_MIGs <- ggplot(joined, aes(x = reorder(Sample_ID, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Sample ID Ordered by Number of Minor Intron Retention Events") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample_ID") +
  geom_text_repel(data=(joined %>% filter(Sample_ID %in% top_10)), size=3, aes(x = reorder(Sample_ID, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = Sample_ID), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)-sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_MIGs_theta_juncs)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) - sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 12*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_blood_MIGs

ggsave(filename=filepath_blood_MIGs, plot=function_plot_blood_MIGs,  limitsize = FALSE, units = "in")

################-----Number Theta Blood-----#################
y_lab <- "Number of Intron Retention Events"
title_lab <- paste("Number of Intron Retention Events per Sample in the", "\n", "UDN and GREGoR Stanford Sites", sep=" ")

top_10 <- joined %>% filter(Sample_ID %in% c("RD268", "GSS225379")) %>% select(Sample_ID, no_theta_juncs) %>% head(10) %>% pull(Sample_ID)
g1 <- subset(joined, Sample_ID %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood <- ggplot(joined, aes(x = reorder(Sample_ID, no_theta_juncs), y = no_theta_juncs)) +
  geom_point() +
  xlab("Sample ID Ordered by Number of Intron Retention Events") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample_ID") +
  geom_text_repel(data=(joined %>% filter(Sample_ID %in% top_10)), size=3, aes(x = reorder(Sample_ID, no_theta_juncs), y = no_theta_juncs, label = Sample_ID), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs)+sd(no_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs)+sd(no_theta_juncs)+sd(no_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs)-sd(no_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs)+sd(no_theta_juncs)+sd(no_theta_juncs)+sd(no_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs) + sd(joined$no_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs) + sd(joined$no_theta_juncs) + sd(joined$no_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs) + sd(joined$no_theta_juncs) + sd(joined$no_theta_juncs) + sd(joined$no_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs) - sd(joined$no_theta_juncs)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_blood

ggsave(filename=filepath_blood, plot=function_plot_blood,  limitsize = FALSE, units = "in")

################-----Number Theta Blood Adj P-val-----#################
y_lab <- "Number of Highly Significant Intron Retention Events"
title_lab <- paste("Number of Highly Significant Intron Retention Events per Sample in the", "\n", "UDN and GREGoR Stanford Sites", "\n", "(adjusted pval < 1E-6, abs(deltaPsi) > 0.3)", sep=" ")

top_10 <- joined %>% filter(Sample_ID %in% c("RD268", "GSS225379")) %>% select(Sample_ID, no_theta_juncs) %>% head(10) %>% pull(Sample_ID)
g1 <- subset(joined, Sample_ID %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood_sig <- ggplot(joined, aes(x = reorder(Sample_ID, no_theta_pval_1E_neg6_abs_deltaPsi), y = no_theta_pval_1E_neg6_abs_deltaPsi)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample_ID") +
  geom_text_repel(data=(joined %>% filter(Sample_ID %in% top_10)), size=3, aes(x = reorder(Sample_ID, no_theta_pval_1E_neg6_abs_deltaPsi), y = no_theta_pval_1E_neg6_abs_deltaPsi, label = Sample_ID), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_pval_1E_neg6_abs_deltaPsi), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_pval_1E_neg6_abs_deltaPsi)-sd(no_theta_pval_1E_neg6_abs_deltaPsi), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi)+sd(no_theta_pval_1E_neg6_abs_deltaPsi), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_pval_1E_neg6_abs_deltaPsi)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi) + sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_pval_1E_neg6_abs_deltaPsi) - sd(joined$no_theta_pval_1E_neg6_abs_deltaPsi)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_blood_sig

ggsave(filename=filepath_blood_sig, plot=function_plot_blood_sig,  limitsize = FALSE, units = "in")

################-----MIG and NMD RD268-----#################
RD268_tx <- transcriptomics %>% filter(Ind == "UDN889016")
RD268_splicing <- all_uncompiled %>% filter(sampleID == "RD268")

RD268_splicing_theta <- RD268_splicing %>% filter(type == "theta")

MIGs_hgnc <- MIGs %>% filter(gene_class == "MIG") %>% pull(gene_symbol)
MIGs_ens <- MIGs %>% filter(gene_class == "MIG") %>% pull(ensembl_gene_id)

RD268_underexp <- RD268_tx %>% filter(normalized_resid <= -2 )

RD268_MIG_underexp <- RD268_tx %>% filter(gene %in% MIGs_ens) %>% filter(normalized_resid <= -2 )

RD268_MIG <- RD268_tx %>% filter(gene %in% MIGs_ens)
mean_GTEx <- mean(transcriptomics$normalized_resid, na.rm=TRUE)
sd_GTEx <- sd(transcriptomics$normalized_resid, na.rm=TRUE)
RD268_MIG$diffexpressed <- ifelse(RD268_MIG$normalized_resid >=2, "UP",
                                        ifelse(RD268_MIG$normalized_resid <= -2, "DOWN", "NO"))
RD268_MIG$pvalue <- 2*pnorm(abs(RD268_MIG$normalized_resid), mean = mean_GTEx, sd = sd_GTEx, lower.tail=FALSE)
RD268_MIG$delabel <- ifelse(RD268_MIG$gene %in% head(RD268_MIG[order(RD268_MIG$pvalue), "gene"], 30), RD268_MIG$gene, NA)


RD268_splicing_theta_MIG <- RD268_splicing_theta %>% filter(hgncSymbol %in% MIGs_hgnc)

RD268_splicing_theta_MIG$ensembl_gene_id = mapIds(org.Hs.eg.db,
                    keys=RD268_splicing_theta_MIG$hgncSymbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")

#find those that are lowly expressed AND have intron retention events
intersect(RD268_MIG_underexp$gene, RD268_splicing_theta_MIG$ensembl_gene_id)

#find those that are lowly expressed AND do NOT have intron retention events
NMD <- setdiff(RD268_MIG_underexp$gene, RD268_splicing_theta_MIG$ensembl_gene_id)

myvolcanoplot <- ggplot(data = RD268_MIG, aes(x = normalized_resid, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 37.5), xlim = c(-9.5, 13.5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-15, 15, 2)) + # to customise the breaks in the x axis
  ggtitle('Gene Expression Levels of Minor Intron containing Genes (MIGs)') + # Plot title
  geom_text_repel(max.overlaps = Inf) # To show all labels 

myvolcanoplot
ggsave(filename=filepath_volcano_MIG, plot=myvolcanoplot,  limitsize = FALSE, units = "in")

RD268_MIG %>% filter(normalized_resid <= 2) %>% filter(pvalue >= 0.05) %>% nrow

################-----pLi Constraint Score-----#################
NMD_Gene_Name = mapIds(org.Hs.eg.db,
                    keys=NMD, 
                    keytype="ENSEMBL",
                    column="SYMBOL",
                    multiVals="first")

constraint_pruned <- constraint %>% filter(gene %in% NMD_Gene_Name) %>% dplyr::select(gene, pLI)
constraint_pruned_below_0.33 <- constraint_pruned %>% filter(pLI <= 0.33) %>% nrow()
constraint_pruned_0.33_0.66 <- constraint_pruned %>% filter(between(pLI, 0.34, 0.66)) %>% nrow()
constraint_pruned_0.67 <- constraint_pruned %>% filter(pLI > 0.66) %>% nrow()

constraint_counts <- c(constraint_pruned_below_0.33, constraint_pruned_0.33_0.66, constraint_pruned_0.67)
constraint_cats <- c("0-0.33", "0.34-0.66", "0.67-1")
constraint_df <- data.frame(constraint_counts, constraint_cats)

constraint_plot <- ggplot(constraint_df, aes(x=constraint_cats, y=constraint_counts)) + 
  geom_bar(stat = "identity")+
  xlab("Constraint pLI Score Category") +
  ylab("Number of Genes with Low Expression and no Minor Intron Retention Events") +
  labs(title="Number of Genes with Constraints (pLI) from 0-.33, .34-.66, and .67-1")

constraint_plot

ggsave(filename=constraint_plot_filepath, plot=constraint_plot,  limitsize = FALSE, units = "in")

constraint_high <- constraint_pruned %>% filter(pLI > 0.66)

################-----pLi Constraint Score Distribution-----#################
constraint_below_0.33 <- constraint %>% filter(pLI <= 0.33) %>% nrow()
constraint_0.33_0.66 <- constraint %>% filter(between(pLI, 0.34, 0.66)) %>% nrow()
constraint_0.67 <- constraint %>% filter(pLI > 0.66) %>% nrow()
constraint_total <- constraint %>% nrow()

constraint_counts <- c(constraint_below_0.33, constraint_0.33_0.66, constraint_0.67, constraint_total)
constraint_cats <- c("0-0.33", "0.34-0.66", "0.67-1", "All")
constraint_df <- data.frame(constraint_counts, constraint_cats)

constraint_df$proportion <- constraint_df$constraint_counts/constraint_total
7/45

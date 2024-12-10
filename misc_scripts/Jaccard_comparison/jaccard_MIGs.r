library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

all_outliers <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake_old_filters/output/FRASER2_RD112_try3_prpf8/FRASER2/FRASER_output_filtered.csv", col_names=FALSE)
colnames(all_outliers) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi",
                      "counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

MIGs_Jaccard <- all_outliers %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) 

MIGs_Jaccard_count <- MIGs_Jaccard %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
mean(MIGs_Jaccard_count$n) + sd(MIGs_Jaccard_count$n) + sd(MIGs_Jaccard_count$n)

RD268 <- (MIGs_Jaccard_count %>% filter(sampleID == "RD268") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
GSS225379 <- (MIGs_Jaccard_count %>% filter(sampleID == "GSS225379") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
UDN550488  <- (MIGs_Jaccard_count %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
UDN238929  <- (MIGs_Jaccard_count %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
RD380  <- (MIGs_Jaccard_count %>% filter(sampleID == "RD380") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
RD143 <- (MIGs_Jaccard_count %>% filter(sampleID == "RD143") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
RD140 <- (MIGs_Jaccard_count %>% filter(sampleID == "RD140") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)


RD268_genes <- all_outliers %>% filter(sampleID == "RD268") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
GSS225379_genes <- all_outliers %>% filter(sampleID == "GSS225379") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
UDN550488_genes  <- all_outliers %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
UDN238929_genes  <- all_outliers %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
RD380_genes  <- all_outliers %>% filter(sampleID == "RD380") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
RD143_genes  <- all_outliers %>% filter(sampleID == "RD143") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
RD140_genes  <- all_outliers %>% filter(sampleID == "RD140") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length

MIGs_Jaccard <- all_outliers %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) 

MIGs_Jaccard_count_genes <- MIGs_Jaccard %>% select(sampleID, hgncSymbol) %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
MIGs_Jaccard_count_genes$zscore <- (MIGs_Jaccard_count_genes$n - mean(MIGs_Jaccard_count_genes$n))/sd(MIGs_Jaccard_count_genes$n)

########################################################################
######################-----Plot-----##########################
########################################################################

######################-----Geom_Point-----##########################
x_lab <- expression(atop("Samples Ordered by Number of",
                       "Jaccard Index Outlier Junctions in MIGs"))
y_lab <- expression("Number of Jaccard Index Outlier Junctions in MIGs")
MIGs_Jaccard_count$numbers <- rownames(MIGs_Jaccard_count)
colnames(MIGs_Jaccard_count) <- c("sampleID", "no_theta_juncs_in_MIGs", "numbers")

RNU4ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
RNU6ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "RNU6ATAC-opathy"
non_RNU4ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC, RNU6ATAC)

g1 <- subset(joined, type %in% c("RNU4ATAC-opathy", "RNU6ATAC-opathy"))

#check which distribution this is considering
#give a z-score
ninet <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*19)

colours <- c("#D7DCD8", "#DA786C", "#D3B054")
names(colours) <- c(NA, "RNU4ATAC-opathy", "RNU6ATAC-opathy")

MIGs_RNU4ATAC <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=ninet, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 3, y=mean(joined$no_theta_juncs_in_MIGs)+6, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 3, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+6),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 3, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+6),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 3, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+6),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 19*SD", x = 1, y=ninet +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)

MIGs_RNU4ATAC

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/Jaccard_comparison/geom_point.pdf", plot=MIGs_RNU4ATAC,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----Boxplot-----##########################
colnames(MIGs_Jaccard_count_genes) <- c("sampleID", "no_MIGs_with_theta_juncs", "zscore")
title_lab <- expression(atop("Comparison of the Number of" ~ theta ~ "Outlier Junctions in MIGs",
                     "Per Sample between Samples with and without an Excess",
                     "Number of" ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of Jaccard Index Outliers Junctions in MIGs")
joined_RNU4ATAC <- MIGs_Jaccard_count_genes %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_RNU4ATAC$MS <- "RNU4atac-opathy"

joined_RNU6ATAC <- MIGs_Jaccard_count_genes %>% filter(sampleID == "RD380")
joined_RNU6ATAC$MS <- "RNU6atac-opathy"

joined_non_MS <- MIGs_Jaccard_count_genes %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_RNU4ATAC, joined_non_MS, joined_RNU6ATAC)
colours <- c("#D7DCD8", "#DA786C", "#D3B054")
names(colours) <- c(NA, "RNU4atac-opathy", "RNU6atac-opathy")

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_with_theta_juncs, fill=MS)) + 
  geom_boxplot(show.legend=FALSE, outlier.shape=21, outlier.alpha=.7, outlier.size=10) +
  ylab(y_lab) +
  xlab("Groups")+
  #labs(title=title_lab)  +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text( vjust = 1))+
  scale_fill_manual(values=colours)


MIG_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/Jaccard_comparison/boxplot.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in", height=10, width=10)


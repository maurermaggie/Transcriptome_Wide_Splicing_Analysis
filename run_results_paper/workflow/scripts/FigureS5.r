library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

args <- commandArgs(TRUE)
all_outliers_fp <- args[1]
all_outliers <- read_csv(all_outliers_fp, col_names=FALSE)
colnames(all_outliers) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi",
                      "counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")

MIGs_fp <- args[2]
MIGs <- read_csv(MIGs_fp)

jaccard_boxplot_geom_pt <- args[3]

#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

MIGs_Jaccard <- all_outliers %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) 

MIGs_Jaccard_count <- MIGs_Jaccard %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
mean(MIGs_Jaccard_count$n) + sd(MIGs_Jaccard_count$n) + sd(MIGs_Jaccard_count$n)

A1 <- (MIGs_Jaccard_count %>% filter(sampleID == "A1") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
B1 <- (MIGs_Jaccard_count %>% filter(sampleID == "B1") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
UDN550488  <- (MIGs_Jaccard_count %>% filter(sampleID == "C1") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
UDN238929  <- (MIGs_Jaccard_count %>% filter(sampleID == "C2") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
D1  <- (MIGs_Jaccard_count %>% filter(sampleID == "D1") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
RD143 <- (MIGs_Jaccard_count %>% filter(sampleID == "RD143") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)
RD140 <- (MIGs_Jaccard_count %>% filter(sampleID == "RD140") %>% pull(n) - mean(MIGs_Jaccard_count$n))/ sd(MIGs_Jaccard_count$n)


A1_genes <- all_outliers %>% filter(sampleID == "A1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
B1_genes <- all_outliers %>% filter(sampleID == "B1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
UDN550488_genes  <- all_outliers %>% filter(sampleID == "C1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
UDN238929_genes  <- all_outliers %>% filter(sampleID == "C2") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique %>% length
D1_genes  <- all_outliers %>% filter(sampleID == "D1") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% 
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

RNU4ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("A1", "B1", "C1", "C2"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
RNU6ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("D1"))
RNU6ATAC$type <- "RNU6ATAC-opathy"
non_RNU4ATAC <- MIGs_Jaccard_count %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("A1", "B1", "C1", "C2", "D1"))
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

ggsave(filename=jaccard_boxplot_geom_pt, plot=MIGs_RNU4ATAC,  limitsize = FALSE, units = "in", height=12, width=10)
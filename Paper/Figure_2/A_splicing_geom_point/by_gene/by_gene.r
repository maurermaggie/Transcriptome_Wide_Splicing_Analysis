library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

metadata_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

########################################################################
######################-----By Gene per Person-----######################
########################################################################
joined <- metadata_joined %>% select(sampleID, ends_with("Genes"))
psi3_count <- joined %>% select(sampleID, Psi3_Genes)
psi5_count <- joined %>% select(sampleID, Psi5_Genes)
theta_count <- joined %>% select(sampleID, Theta_Genes)
jaccard_count <- joined %>% select(sampleID, Jaccard_Genes)

colnames(psi3_count) <- c("sampleID", "n")
colnames(psi5_count) <- c("sampleID", "n")
colnames(theta_count) <- c("sampleID", "n")
colnames(jaccard_count) <- c("sampleID", "n")

########################################################################
#########################-----By Gene Plots-----########################
########################################################################

################-----Psi3 Plot-----#################
joined <- psi3_count

title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi3 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  labs(colour="Sample ID") +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi3

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/A_splicing_geom_point/by_gene/psi3_by_gene.pdf", plot=function_plot_psi3,  limitsize = FALSE, units = "in")

################-----Psi3 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi5 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi5

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/A_splicing_geom_point/by_gene/psi5_by_gene.pdf", plot=function_plot_psi5,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_theta <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001, vjust = -.01) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001, vjust =-.01) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_theta

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/A_splicing_geom_point/by_gene/theta_by_gene.pdf", plot=function_plot_theta,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- jaccard_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- "Jaccard Index"
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant Jaccard Index Events"))
y_lab <- expression("Number of Genes with Significant Jaccard Index Events")

#check which distribution this is considering
#give a z-score
function_plot_jaccard <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001, vjust = -.01) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001, vjust =-.01) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_jaccard

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/A_splicing_geom_point/by_gene/jaccard_by_gene.pdf", plot=function_plot_jaccard,  limitsize = FALSE, units = "in")




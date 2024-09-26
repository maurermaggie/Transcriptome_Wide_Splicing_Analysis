#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(reshape)
library(GO.db)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")

################-----Get Number by Gene per Person-----#################
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

################-----Psi3 Plot-----#################
joined <- psi3_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
x_lab <- paste("Samples ordered by nuber of", "\n", "Genes with Significant Psi3 Events")

#check which distribution this is considering
#give a z-score
function_plot_psi3 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab("Number of Genes with Significant Psi3 Events") +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/psi3.pdf", plot=function_plot_psi3,  limitsize = FALSE, units = "in")

################-----Psi3 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
x_lab <- paste("Samples ordered by nuber of", "\n", "Genes with Significant Psi5 Events")

#check which distribution this is considering
#give a z-score
function_plot_psi5 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab("Number of Genes with Significant Psi5 Events") +
  #labs(title=title_lab)  +
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/psi5.pdf", plot=function_plot_psi5,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
x_lab <- paste("Samples ordered by nuber of", "\n", "Genes with Significant Theta Events")

#check which distribution this is considering
#give a z-score
function_plot_theta <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab("Number of Genes with Significant Theta Events") +
  #labs(title=title_lab)  +
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/theta.pdf", plot=function_plot_theta,  limitsize = FALSE, units = "in")


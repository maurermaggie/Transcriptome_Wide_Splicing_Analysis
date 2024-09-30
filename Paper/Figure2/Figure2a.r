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

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_output.csv")

########################################################################
######################-----By Gene per Person-----######################
########################################################################
get_count <- function(filtered_df, colname) {
    by_gene <- filtered_df %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
    count <- by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
    missing <- setdiff(files$sampleID, count$sampleID)
    missing_df <- data.frame(sampleID = missing, n = 0)
    count <- bind_rows(count, missing_df)
    colnames(count) <- c("sampleID", colname)
    count <- count %>% filter(! sampleID %in% low_RIN$sampleID)
    count
}

################-----Get Number by Gene per Person-----#################
#0.05
psi3_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_count <- get_count(psi3_filtered, "Psi3")

psi5_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_count <- get_count(psi5_filtered, "Psi5")

theta_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_count <- get_count(theta_filtered, "Theta")

#1e-6
psi3_filtered_1e_6 <- all_uncompiled %>% filter(padjust <= 1e-6) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_count_1e_6 <- get_count(psi3_filtered_1e_6, "Psi3_1e_6")

psi5_filtered_1e_6 <- all_uncompiled %>% filter(padjust <= 1e-6) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_count_1e_6 <- get_count(psi5_filtered_1e_6, "Psi5_1e_6")

theta_filtered_1e_6 <- all_uncompiled %>% filter(padjust <= 1e-6) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_count_1e_6 <- get_count(theta_filtered_1e_6, "Theta_1e_6")

#.001
psi3_filtered_001 <- all_uncompiled %>% filter(padjust <= .001) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_count_001 <- get_count(psi3_filtered_001, "Psi3_001")

psi5_filtered_001 <- all_uncompiled %>% filter(padjust <= .001) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_count_001 <- get_count(psi5_filtered_001, "Psi5_001")

theta_filtered_001 <- all_uncompiled %>% filter(padjust <= .001) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_count_001 <- get_count(theta_filtered_001, "Theta_001")

#join
joined <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% 
                  left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% left_join(theta_count_001)

########################################################################
######################-----Significance Plots-----######################
########################################################################

################-----Psi3 Significance Plot-----#################
title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_psi3 <- ggplot(joined, aes(x = reorder(sampleID, Psi3), y = Psi3)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Psi3, x = reorder(sampleID, Psi3))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Psi3_1e_6, x = reorder(sampleID, Psi3_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Psi3_001, x = reorder(sampleID, Psi3_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Psi3)+sd(Psi3)+sd(Psi3), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Psi3)+sd(Psi3)+sd(Psi3)+sd(Psi3), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_psi3

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/significance_psi3.pdf", plot=significance_plot_psi3,  limitsize = FALSE, units = "in")

################-----Psi5 Significance Plot-----#################
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_psi5 <- ggplot(joined, aes(x = reorder(sampleID, Psi5), y = Psi5)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Psi5, x = reorder(sampleID, Psi5))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Psi5_1e_6, x = reorder(sampleID, Psi5_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Psi5_001, x = reorder(sampleID, Psi5_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Psi5)+sd(Psi5)+sd(Psi5), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Psi5)+sd(Psi5)+sd(Psi5)+sd(Psi5), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_psi5

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/significance_psi5.pdf", plot=significance_plot_psi5,  limitsize = FALSE, units = "in")

################-----Theta Significance Plot-----#################
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_theta <- ggplot(joined, aes(x = reorder(sampleID, Theta), y = Theta)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Theta, x = reorder(sampleID, Theta))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Theta_1e_6, x = reorder(sampleID, Theta_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Theta_001, x = reorder(sampleID, Theta_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Theta)+sd(Theta)+sd(Theta), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Theta)+sd(Theta)+sd(Theta)+sd(Theta), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Theta) + sd(joined$Theta) + sd(joined$Theta) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Theta) + sd(joined$Theta) + sd(joined$Theta) + sd(joined$Theta)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_theta

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/significance_theta.pdf", plot=significance_plot_theta,  limitsize = FALSE, units = "in")

########################################################################
#########################-----By Gene Plots-----########################
########################################################################

################-----Psi3 Plot-----#################
joined <- psi3_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/psi5.pdf", plot=function_plot_psi5,  limitsize = FALSE, units = "in")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/theta.pdf", plot=function_plot_theta,  limitsize = FALSE, units = "in")

########################################################################
####################-----By Junction per Person-----####################
########################################################################
################-----Get Number by Gene per Person-----#################
psi3_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_count <- psi3_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi3 <- setdiff(files$sampleID, psi3_count$sampleID)
missing_psi3_df <- data.frame(sampleID = missing_psi3, n = 0)
psi3_count <- bind_rows(psi3_count, missing_psi3_df)
colnames(psi3_count) <- c("sampleID", "Psi3")
psi3_count <- psi3_count %>% filter(! sampleID %in% low_RIN$sampleID)

psi5_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_count <- psi5_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi5 <- setdiff(files$sampleID, psi5_count$sampleID)
missing_psi5_df <- data.frame(sampleID = missing_psi5, n = 0)
psi5_count <- bind_rows(psi5_count, missing_psi5_df)
colnames(psi5_count) <- c("sampleID", "Psi5")
psi5_count <- psi5_count %>% filter(! sampleID %in% low_RIN$sampleID)

theta_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_count <- theta_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_theta <- setdiff(files$sampleID, theta_count$sampleID)
missing_theta_df <- data.frame(sampleID = missing_theta, n = 0)
theta_count <- bind_rows(theta_count, missing_theta_df)
colnames(theta_count) <- c("sampleID", "Theta")
theta_count <- theta_count %>% filter(! sampleID %in% low_RIN$sampleID)

joined <- left_join(psi3_count, psi5_count) %>% left_join(theta_count)

################-----Psi3 Plot-----#################
joined <- psi3_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Junctions Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi3_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
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
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi3_junctions

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/psi3_junctions.pdf", plot=function_plot_psi3_junctions,  limitsize = FALSE, units = "in")

################-----Psi3 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Junctions Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi5_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
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

function_plot_psi5_junctions

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/psi5_junctions.pdf", plot=function_plot_psi5_junctions,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Junctions Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_theta_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
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

function_plot_theta_junctions

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/theta_junctions.pdf", plot=function_plot_theta_junctions,  limitsize = FALSE, units = "in")


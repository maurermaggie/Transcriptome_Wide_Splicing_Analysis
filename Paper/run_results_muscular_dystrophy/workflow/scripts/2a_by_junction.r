library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(tidyverse)

metadata_joined_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv"
metadata_joined <- read_csv(metadata_joined_fp)

output_theta <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/theta_2a.pdf"
output_psi3 <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/psi3_2a.pdf"
output_psi5 <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/psi5_2a.pdf"
output_jaccard <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/jaccard_2a.pdf"

output_stats <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/figure_2a.pdf"
output_venn <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/venn_2a.pdf"
output_upset <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/upset_2a.pdf"

########################################################################
###################-----By Junction per Person-----#####################
########################################################################
joined <- metadata_joined %>% select(sampleID, ends_with("Junctions"))
psi3_count <- joined %>% select(sampleID, Psi3_Junctions)
psi5_count <- joined %>% select(sampleID, Psi5_Junctions)
theta_count <- joined %>% select(sampleID, Theta_Junctions)
jaccard_count <- joined %>% select(sampleID, Jaccard_Junctions)

colnames(psi3_count) <- c("sampleID", "n")
colnames(psi5_count) <- c("sampleID", "n")
colnames(theta_count) <- c("sampleID", "n")
colnames(jaccard_count) <- c("sampleID", "n")

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

ggsave(filename=output_psi3, plot=function_plot_psi3_junctions,  limitsize = FALSE, units = "in")

################-----Psi5 Plot-----#################
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
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=8, hjust = -.00001) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi5_junctions

ggsave(filename=output_psi5, plot=function_plot_psi5_junctions,  limitsize = FALSE, units = "in")

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

ggsave(filename=output_theta, plot=function_plot_theta_junctions,  limitsize = FALSE, units = "in")

################-----Jaccard Plot-----#################
joined <- jaccard_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- "Jaccard Index"
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant Jaccard Index Events"))
y_lab <- expression("Number of Junctions Significant Jaccard Index Events")

#check which distribution this is considering
#give a z-score
function_plot_jaccard_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
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

function_plot_jaccard_junctions

ggsave(filename=output_jaccard, plot=function_plot_jaccard_junctions,  limitsize = FALSE, units = "in")

########################################################################
####################-----By Junction per Person-----####################
########################################################################
joined <- metadata_joined %>% select(sampleID, ends_with("Junctions"))

two_sd_psi3 <- mean(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions)
two_sd_psi5 <- mean(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions)
two_sd_theta <- mean(joined$Theta_Junctions) + sd(joined$Theta_Junctions) + sd(joined$Theta_Junctions)
two_sd_jaccard <- mean(joined$Jaccard_Junctions) + sd(joined$Jaccard_Junctions) + sd(joined$Jaccard_Junctions)

filter(joined, Psi3_Junctions >= two_sd_psi3) %>% nrow
filter(joined, Psi5_Junctions >= two_sd_psi5) %>% nrow
filter(joined, Theta_Junctions >= two_sd_theta) %>% nrow

psi3_outs <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID)
psi5_outs <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID)
theta_outs <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID)
jaccard_outs <- filter(joined, Jaccard_Junctions >= two_sd_jaccard) %>% pull(sampleID)

psi3_outs_l <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID) %>% length
psi5_outs_l <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID) %>% length
theta_outs_l <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID) %>% length
jaccard_outs_l <- filter(joined, Jaccard_Junctions >= two_sd_jaccard) %>% pull(sampleID) %>% length

stats_list <- list(`Psi3 Outliers`= psi3_outs, `Psi5 Outliers`=psi5_outs, `Theta Outliers` =theta_outs )

stats_df <- data.frame(Metrics = c("Psi3", "Psi5", "Theta", "Jaccard"), Number_Outliers=c(psi3_outs_l, psi5_outs_l, theta_outs_l, jaccard_outs_l))

write_csv(stats_df, output_stats)

stats_venn <- ggVennDiagram(stats_list,  order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)
ggsave(filename=output_venn, plot=stats_venn,  limitsize = FALSE, units = "in")

stats <- ggVennDiagram(stats_list, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)
ggsave(filename=output_upset, plot=stats,  limitsize = FALSE, units = "in")

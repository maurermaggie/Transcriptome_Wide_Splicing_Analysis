library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(tidyverse)

args <- commandArgs(TRUE)
metadata_joined_fp <- args[1]
metadata_joined <- read_csv(metadata_joined_fp)

output_theta <- args[2]
output_psi3 <- args[3]
output_psi5 <- args[4]
output_jaccard <- args[5]

output_stats <- args[6]
output_venn <- args[7]
output_upset <- args[8]

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
                     "Significant" ~ psi ~ 3 ~ "Outliers"))
y_lab <- expression("Number of Junctions Significant" ~ Psi ~ 3 ~ "Outliers")

#check which distribution this is considering
#give a z-score
function_plot_psi3_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
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
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi3_junctions

ggsave(filename=output_psi3, plot=function_plot_psi3_junctions,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Psi5 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 5 ~ "Outliers"))
y_lab <- expression("Number of Junctions Significant" ~ Psi ~ 5 ~ "Outliers")

#check which distribution this is considering
#give a z-score
function_plot_psi5_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#DA786C") +
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
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi5_junctions

ggsave(filename=output_psi5, plot=function_plot_psi5_junctions,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ theta ~ "Outliers"))
y_lab <- expression("Number of Junctions Significant" ~ theta ~ "Outliers")

#check which distribution this is considering
#give a z-score
function_plot_theta_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#88ACAA") +
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
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_theta_junctions

ggsave(filename=output_theta, plot=function_plot_theta_junctions,  limitsize = FALSE, units = "in", height=12, width=11)

################-----Jaccard Plot-----#################
joined <- jaccard_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- "Jaccard Index"
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant Jaccard Index Outliers"))
y_lab <- expression("Number of Junctions Significant Jaccard Index Outliers")

#check which distribution this is considering
#give a z-score
function_plot_jaccard_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#F1AE7A") +
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
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_jaccard_junctions

ggsave(filename=output_jaccard, plot=function_plot_jaccard_junctions,  limitsize = FALSE, units = "in", height=12, width=10)

########################################################################
####################-----By Junction per Person-----####################
########################################################################
joined <- metadata_joined %>% select(sampleID, ends_with("Junctions"))

two_sd_psi3 <- mean(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions) + sd(joined$Psi3_Junctions)
two_sd_psi5 <- mean(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions) + sd(joined$Psi5_Junctions)
two_sd_theta <- mean(joined$Theta_Junctions) + sd(joined$Theta_Junctions) + sd(joined$Theta_Junctions)
two_sd_jaccard <- mean(joined$Jaccard_Junctions) + sd(joined$Jaccard_Junctions) + sd(joined$Jaccard_Junctions)
two_sd_all <- mean(joined$All_Junctions) + sd(joined$All_Junctions) + sd(joined$All_Junctions)

filter(joined, Psi3_Junctions >= two_sd_psi3) %>% nrow
filter(joined, Psi5_Junctions >= two_sd_psi5) %>% nrow
filter(joined, Theta_Junctions >= two_sd_theta) %>% nrow

psi3_outs <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID)
psi5_outs <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID)
theta_outs <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID)
jaccard_outs <- filter(joined, Jaccard_Junctions >= two_sd_jaccard) %>% pull(sampleID)
all_outs <- filter(joined, All_Junctions >= two_sd_all) %>% pull(sampleID)

psi3_outs_l <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID) %>% length
psi5_outs_l <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID) %>% length
theta_outs_l <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID) %>% length
jaccard_outs_l <- filter(joined, Jaccard_Junctions >= two_sd_jaccard) %>% pull(sampleID) %>% length
all_outs_l <- filter(joined, All_Junctions >= two_sd_all) %>% pull(sampleID) %>% length

stats_list <- list(`Psi3 Outliers`= psi3_outs, `Psi5 Outliers`=psi5_outs, `Theta Outliers` =theta_outs)

outliers_all <- c(psi3_outs, psi5_outs, theta_outs, jaccard_outs, all_outs) %>% unique %>% length
stats_df <- data.frame(Metrics = c("Psi3", "Psi5", "Theta", "Jaccard", "All", "All Any Type"), Number_Outliers=c(psi3_outs_l, psi5_outs_l, theta_outs_l, jaccard_outs_l, all_outs_l, outliers_all))

write_csv(stats_df, output_stats)

stats_venn <- ggVennDiagram(stats_list,  order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)
ggsave(filename=output_venn, plot=stats_venn,  limitsize = FALSE, units = "in")

stats <- ggVennDiagram(stats_list, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none",   
                top.bar.show.numbers = TRUE,  top.bar.numbers.size = 10, sets.bar.show.numbers = TRUE, set_size = 10, label_size=10)
ggsave(filename=output_upset, plot=stats,  limitsize = FALSE, units = "in")
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

########################################################################
####################-----By Junction per Person-----####################
########################################################################
psi3_count <- joined %>% select(sampleID, Psi3_Junctions)
psi5_count <- joined %>% select(sampleID, Psi5_Junctions)
theta_count <- joined %>% select(sampleID, Theta_Junctions)
all_count <- joined %>% select(sampleID, All_Junctions)
jaccard_count <- joined %>% select(sampleID, Jaccard_Junctions)

################-----Psi3 Plot-----#################
joined <- psi3_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Junctions with Significant" ~ Psi ~ 3 ~ "Events")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_junction/psi3_junctions.pdf", plot=function_plot_psi3_junctions,  limitsize = FALSE, units = "in")

################-----Psi3 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Junctions with Significant" ~ Psi ~ 5 ~ "Events")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_junction/psi5_junctions.pdf", plot=function_plot_psi5_junctions,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Junctions with Significant" ~ theta ~ "Events")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_junction/theta_junctions.pdf", plot=function_plot_theta_junctions,  limitsize = FALSE, units = "in")

################-----All Plot-----#################
joined <- all_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- "All Significant Junctions"
x_lab <- paste("Samples ordered by Number of", "\n",
              "Significant Outlier Junctions")
y_lab <- "Number of Significant Outlier Junctions"

#check which distribution this is considering
#give a z-score
function_plot_all_junctions <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
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

function_plot_all_junctions

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_junction/all_by_junctions.pdf", plot=function_plot_all_junctions,  limitsize = FALSE, units = "in")

two_sd <- mean(all_count$All) + sd(all_count$All) + sd(all_count$All)
two_sd_Psi3 <- mean(psi3_count$Psi3) + sd(psi3_count$Psi3) + sd(psi3_count$Psi3)
two_sd_Psi5 <- mean(psi5_count$Psi5) + sd(psi5_count$Psi5) + sd(psi5_count$Psi5)
two_sd_Theta <- mean(theta_count$Theta) + sd(theta_count$Theta) + sd(theta_count$Theta)

all <- filter(all_count, All >= two_sd) %>% pull(sampleID)
psi3 <- filter(psi3_count, Psi3 >= two_sd_Psi3) %>% pull(sampleID)
psi5 <- filter(psi5_count, Psi5 >= two_sd_Psi5) %>% pull(sampleID)
theta <- filter(theta_count, Theta >= two_sd_Theta) %>% pull(sampleID)

all_outliers <- c(all, psi3, psi5, theta) %>% unique
37/400

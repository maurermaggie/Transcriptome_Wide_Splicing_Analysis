require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

args <- commandArgs(TRUE)
genes_junctions_joined_fp <- args[1]
genes_junctions_joined <- read_csv(genes_junctions_joined_fp)
output_file <- args[2]
output_file2 <- args[3]


########################################################################
#####################-----Get Long Dataframes-----######################
########################################################################
joined_genes <- genes_junctions_joined %>% select(sampleID, ends_with("Genes"))
joined_junctions <- genes_junctions_joined %>% select(sampleID, ends_with("Junctions"))

joined_long_genes <- melt(setDT(joined_genes), id.vars = c("sampleID"), variable.name = "genes_with_aberrant_junction_of_type")

joined_long_junctions <- melt(setDT(joined_junctions), id.vars = c("sampleID"), variable.name = "aberrant_junction_of_type")


########################################################################
##################-----Make Plots Junctions-----########################
########################################################################

################-----Box Plot-----#################
title_lab <- paste("Number of Aberrant Splice Junctions", "\n", "of each FRASER Type", sep="")

psi3 <- expression(psi ~ 3)
psi5 <- expression(psi ~ 5)
theta <- expression(theta)

MyColour <- c("#7D5E62", "#DA786C", "#88ACAA", "#D7DCD8", "#F1AE7A")
names(MyColour) <- c("Psi3_Junctions", "Psi5_Junctions", "Theta_Junctions", "All_Junctions", "Jaccard_Junctions")

FRASER_boxplot_junctions <- ggplot(joined_long_junctions, aes(fill = aberrant_junction_of_type ,x=aberrant_junction_of_type, y=value)) + 
  geom_boxplot(show.legend = FALSE) +
  scale_x_discrete(labels = c('Psi3_Junctions' = psi3,
                                      'Theta_Junctions' = theta,
                                          'Psi5_Junctions' = psi5,
                                              'All_Junctions' = paste(" ", "\n", "\n", "Combined", "\n", "FRASER"),
                                              'Jaccard_Junctions' = paste("Jaccard", "\n", "Index")))+
  xlab("Splicing Outlier Statistic") +
  ylab("Number of Outlier Splice Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30, vjust=-.7)) +
  theme(axis.title.x = element_text(size = 30, vjust=-1.3)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=MyColour)

FRASER_boxplot_junctions
ggsave(filename=output_file, plot=FRASER_boxplot_junctions,  limitsize = FALSE, units = "in", height=10, width=12)

################-----Get Number by Junction per Person-----#################
all_mean <- mean(joined_junctions$All_Junctions)
all_median <- median(joined_junctions$All_Junctions)
all_sd <- sd(joined_junctions$All_Junctions)
all_first <- quantile(joined_junctions$All_Junctions, 0.25)
all_last <- quantile(joined_junctions$All_Junctions, 0.75)

psi3_mean <- mean(joined_junctions$Psi3_Junctions)
psi3_median <- median(joined_junctions$Psi3_Junctions)
psi3_sd <- sd(joined_junctions$Psi3_Junctions)
psi3_first <- quantile(joined_junctions$Psi3_Junctions, 0.25)
psi3_last <- quantile(joined_junctions$Psi3_Junctions, 0.75)

psi5_mean <- mean(joined_junctions$Psi5_Junctions)
psi5_median <- median(joined_junctions$Psi5_Junctions)
psi5_sd <- sd(joined_junctions$Psi5_Junctions)
psi5_first <- quantile(joined_junctions$Psi5_Junctions, 0.25)
psi5_last <- quantile(joined_junctions$Psi5_Junctions, 0.75)

theta_mean <- mean(joined_junctions$Theta_Junctions)
theta_median <- median(joined_junctions$Theta_Junctions)
theta_sd <- sd(joined_junctions$Theta_Junctions)
theta_first <- quantile(joined_junctions$Theta_Junctions, 0.25)
theta_last <- quantile(joined_junctions$Theta_Junctions, 0.75)

jac_mean <- mean(joined_junctions$Jaccard_Junctions)
jac_median <- median(joined_junctions$Jaccard_Junctions)
jac_sd <- sd(joined_junctions$Jaccard_Junctions)
jac_first <- quantile(joined_junctions$Jaccard_Junctions, 0.25)
jac_last <- quantile(joined_junctions$Jaccard_Junctions, 0.75)

b_stats <- data.frame(metric=c("All", "Psi3", "Psi5", "Theta", "Jaccard Index"), mean=c(all_mean, psi3_mean, psi5_mean, theta_mean, jac_mean), 
                      median= c(all_median, psi3_median, psi5_median, theta_median, jac_median), sd= c(all_sd, psi3_sd, psi5_sd, theta_sd, jac_sd),
                      first=c(all_first, psi3_first, psi5_first, theta_first, jac_first), last=c(all_last, psi3_last, psi5_last, theta_last, jac_last))

write_csv(b_stats, output_file2)



require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

genes_junctions_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv")
output_stats2 <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/figure_1b.csv"
########################################################################
#####################-----Get Long Dataframes-----######################
########################################################################
joined_junctions <- genes_junctions_joined %>% select(sampleID, ends_with("Junctions"))

joined_long_junctions <- melt(setDT(joined_junctions), id.vars = c("sampleID"), variable.name = "aberrant_junction_of_type")
joined_long_junctions <- joined_long_junctions %>% filter(aberrant_junction_of_type != "Psi3_Psi5_Junctions")

########################################################################
##################-----Make Plots Junctions-----########################
########################################################################

################-----Violin Plot-----#################
violin_plot <- ggplot(joined_long_junctions, aes(x=aberrant_junction_of_type, y=value)) + 
  geom_violin()

violin_plot

################-----Box Plot-----#################
title_lab <- paste("Number of Aberrant Splice Junctions", "\n", "of each FRASER Type", sep="")

psi3 <- expression(psi ~ 3)
psi5 <- expression(psi ~ 5)
theta <- expression(theta)

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
  theme(axis.text.x = element_text(size = 27, vjust=-.7)) +
  theme(axis.title.x = element_text(size = 30, vjust=-1.3)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(10,5,2,4,1)])

FRASER_boxplot_junctions
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/Figure1b_junctions.pdf", plot=FRASER_boxplot_junctions,  limitsize = FALSE, units = "in")

########################################################################
######################-----Stats Junctions-----#########################
########################################################################

################-----Get Number by Junction per Person-----#################
all_mean <- mean(genes_junctions_joined$All_Junctions)
all_sd <- sd(genes_junctions_joined$All_Junctions)

psi3_mean <- mean(genes_junctions_joined$Psi3_Junctions)
psi3_sd <- sd(genes_junctions_joined$Psi3_Junctions)

psi5_mean <- mean(genes_junctions_joined$Psi5_Junctions)
psi5_sd <- sd(genes_junctions_joined$Psi5_Junctions)

theta_mean <- mean(genes_junctions_joined$Theta_Junctions)
theta_sd <- sd(genes_junctions_joined$Theta_Junctions)

jac_mean <- mean(genes_junctions_joined$Jaccard_Junctions)
jac_sd <- sd(genes_junctions_joined$Jaccard_Junctions)

b_stats <- data.frame(metric=c("All", "Psi3", "Psi5", "Theta", "Jaccard Index"), mean=c(all_mean, psi3_mean, psi5_mean, theta_mean, jac_mean), sd= c(all_sd, psi3_sd, psi5_sd, theta_sd, jac_sd))

write_csv(b_stats, output_stats2)


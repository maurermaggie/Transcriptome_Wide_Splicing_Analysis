require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

genes_junctions_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_filtered_with_sibs.csv")

########################################################################
#####################-----Get Long Dataframes-----######################
########################################################################
joined_genes <- genes_junctions_joined %>% select(sampleID, ends_with("Genes"))
joined_junctions <- genes_junctions_joined %>% select(sampleID, ends_with("Junctions"))

joined_long_genes <- melt(setDT(joined_genes), id.vars = c("sampleID"), variable.name = "genes_with_aberrant_junction_of_type")

joined_long_junctions <- melt(setDT(joined_junctions), id.vars = c("sampleID"), variable.name = "aberrant_junction_of_type")

########################################################################
#####################-----Make Plots Genes-----#########################
########################################################################

################-----Violin Plot-----#################
violin_plot <- ggplot(joined_long_genes, aes(x=genes_with_aberrant_junction_of_type, y=value)) + 
  geom_violin()

violin_plot

################-----Box Plot-----#################
title_lab <- paste("Number of Genes with Aberrant Splice Junctions", "\n", "of each FRASER Type", sep="")

psi3 <- expression(psi ~ 3)
psi5 <- expression(psi ~ 5)
theta <- expression(theta)

FRASER_boxplot_genes <- ggplot(joined_long_genes, aes(fill = genes_with_aberrant_junction_of_type ,x=genes_with_aberrant_junction_of_type, y=value)) + 
  geom_boxplot(show.legend = FALSE) +
  scale_x_discrete(labels = c('Psi3_Genes' = psi3,
                                      'Theta_Genes' = theta,
                                          'Psi5_Genes' = psi5,
                                            'All_Genes' = paste(" ", "\n", "\n", "Combined", "\n", "FRASER"),
                                              'Jaccard_Genes' = paste("Jaccard", "\n", "Index")))+
  xlab("Splicing Outlier Statistic") +
  ylab("Number of Genes with Outlier Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 27, vjust=-.7)) +
  theme(axis.title.x = element_text(size = 30, vjust=-1.3)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(10,5,2,4, 1)])

FRASER_boxplot_genes
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/B_splicing_boxplot/Figure1b_genes.pdf", plot=FRASER_boxplot_genes,  limitsize = FALSE, units = "in")

########################################################################
#####################-----Make Plots Genes-----#########################
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
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/B_splicing_boxplot/Figure1b_junctions.pdf", plot=FRASER_boxplot_junctions,  limitsize = FALSE, units = "in")






require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

########################################################################
################-----Get Number by Gene per Person-----#################
########################################################################
filter_genes <- function(dataframe, filter, column_name){
      dataframe_filtered <- dataframe %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
      if (filter == "psi3") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi3")
        print("filtering by psi3")
      } else if (filter == "psi5") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi5")
        print("filtering by psi5")
      } else if (filter == "theta") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "theta")
        print("filtering by theta")
      } else {
         dataframe_filtered < dataframe_filtered
      }
      dataframe_by_gene <- dataframe_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
      dataframe_count <- dataframe_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

filter_junctions <- function(dataframe, filter, column_name){
      dataframe_filtered <- dataframe %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
      if (filter == "psi3") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi3")
        print("filtering by psi3")
      } else if (filter == "psi5") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi5")
        print("filtering by psi5")
      } else if (filter == "theta") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "theta")
        print("filtering by theta")
      } else {
         dataframe_filtered < dataframe_filtered
      }
      dataframe_count <- dataframe_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

all_count_genes <- filter_genes(all_uncompiled, "all", "All_Genes")
psi3_count_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes")
psi5_count_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes")
theta_count_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes")

all_count_junctions <- filter_junctions(all_uncompiled, "all", "All_Junctions")
psi3_count_junctions <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions")
psi5_count_junctions <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions")
theta_count_junctions <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions")

joined_genes <- left_join(all_count_genes, psi3_count_genes) %>% left_join(psi5_count_genes) %>% left_join(theta_count_genes)
joined_long_genes <- melt(setDT(joined_genes), id.vars = c("sampleID"), variable.name = "genes_with_aberrant_junction_of_type")

joined_junctions <- left_join(all_count_junctions, psi3_count_junctions) %>% left_join(psi5_count_junctions) %>% left_join(theta_count_junctions)
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
                                            'All_Genes' = paste(" ", "\n", "\n", "Combined", "\n", "FRASER")))+
  xlab("Splicing Outlier Type") +
  ylab("Number of Genes with Outlier Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30, vjust=-.7)) +
  theme(axis.title.x = element_text(size = 30, vjust=-.7)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(10,5,2,4)])

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
                                              'All_Junctions' = paste(" ", "\n", "\n", "Combined", "\n", "FRASER")))+
  xlab("Splicing Outlier Type") +
  ylab("Number of Outlier Splice Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30, vjust=-.7)) +
  theme(axis.title.x = element_text(size = 30, vjust=-.7)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(10,5,2,4)])

FRASER_boxplot_junctions
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/B_splicing_boxplot/Figure1b_junctions.pdf", plot=FRASER_boxplot_junctions,  limitsize = FALSE, units = "in")






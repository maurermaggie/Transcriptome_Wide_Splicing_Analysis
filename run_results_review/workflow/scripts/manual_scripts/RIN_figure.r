library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

################-----rm batch 8-----#################
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/DataFrames/metadata_counts_outlier_joined.csv")
joined <- joined %>% select(sampleID, Theta_Genes_Outlier_Status, Jaccard_Genes_Outlier_Status, Psi3_Genes_Outlier_Status, Psi5_Genes_Outlier_Status, All_Genes_Outlier_Status, RIN)

outlier <- joined %>% filter(Theta_Genes_Outlier_Status == 1 | Jaccard_Genes_Outlier_Status == 1 | Psi3_Genes_Outlier_Status == 1 | Psi5_Genes_Outlier_Status == 1 | All_Genes_Outlier_Status == 1)
outliers <- outlier %>% pull(sampleID)

non_outlier <- joined %>% filter(! sampleID %in% outliers)
non_outliers <- non_outlier %>% pull(sampleID)

outliers_RIN <- joined %>% filter(sampleID %in% outliers) %>% select(sampleID, RIN)
non_outliers_RIN <- joined %>% filter(sampleID %in% non_outliers) %>% select(sampleID, RIN)

outliers_RIN$status <- "Outlier"
non_outliers_RIN$status <- "Non Outlier"

joined <- bind_rows(outliers_RIN, non_outliers_RIN)

RIN_boxplot <- ggplot(joined, aes(x = status, y = RIN)) + 
    geom_boxplot(show.legend=FALSE, outlier.shape=21, outlier.alpha=.7, outlier.size=10) +
    theme_classic(base_size = 25)+
    ylab("RIN") +
    xlab("Outlier Status")

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/Plots/RIN/RIN.pdf", plot=RIN_boxplot,  limitsize = FALSE, units = "in", height=12, width=10)
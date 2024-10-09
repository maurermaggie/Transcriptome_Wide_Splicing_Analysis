#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")

#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

################-----Get Significant MIG Intron Retention Events-----#################
#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG")

#get theta dataframe
results_filtered_theta <- all_uncompiled %>% filter(!is.na(hgncSymbol)) %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta")
results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol")
colnames(results_filtered_theta) <- c("sampleID", "gene_symbol")

#combine dataframes
sig_genes_type_all_samples <- left_join(MIG_table_select, results_filtered_theta) %>% filter(!is.na(sampleID))

joined <- sig_genes_type_all_samples %>% select(sampleID, gene_symbol) %>% 
    group_by(sampleID) %>% tally() %>% arrange(desc(n))

missing_all <- setdiff(files$sampleID, joined$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
all_count <- bind_rows(joined, missing_all_df)
all_count <- all_count %>% filter(! sampleID %in% low_RIN$sampleID)

colnames(all_count) <- c("sampleID", "no_MIGs_theta_juncs")
joined <- all_count

################-----Number Minor Introns Retained-----#################
title_lab <- expression(atop("Comparison of the Number of" ~ theta ~ "Events in MIGs",
                     "Per Sample between Samples with and without an Excess",
                     "Number of" ~ theta ~ "Events in MIGs"))
y_lab <- expression("Number of" ~ theta ~ "Events in MIGs")
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_theta_juncs)) + 
  geom_boxplot() +
  ylab(y_lab) +
  xlab("Groups")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text( vjust = 1))


MIG_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_3/B_MIGs_IR_boxplot_RNU4ATAC/plots/MIG_boxplot_RNU4ATAC.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in")

mean(joined_non_MS$no_MIGs_theta_juncs) 
sd(joined_non_MS$no_MIGs_theta_juncs) 
mean(joined_MS$no_MIGs_theta_juncs) 
sd(joined_MS$no_MIGs_theta_juncs) 
312.5/2.813131
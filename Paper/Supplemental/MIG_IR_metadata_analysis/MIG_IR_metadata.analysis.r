#library(EnsDb.Hsapiens.v79)
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
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################

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

################-----Join to Metadata-----#################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch)

names(all_count)[names(all_count) == 'sampleID'] <- 'RDID'
all_filter_RD <- filter(all_count, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)

names(all_count)[names(all_count) == 'RDID'] <- 'GSS_ID'
all_filter_GSS <- filter(all_count, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)

names(all_count)[names(all_count) == 'GSS_ID'] <- 'UDN_ID'
all_filter_UDN <- filter(all_count, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)

names(all_count)[names(all_count) == 'UDN_ID'] <- 'sampleID'

joined_filtered <- bind_rows(all_RDID, all_GSS, all_UDN)

########################################################################
#######################-----Linear Model-----###########################
########################################################################
make_estimate_df <- function(row, dataframe) {
    Esitmate_lm_df <- lm(dataframe[[row]] ~ RIN + age + sex+ batch, data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
    Value <- summary(Esitmate_lm_df)$coefficients[,0]
    df$Value <- c("Intercept", "RIN", "age", "sexM", "batch")
    df$Numbers <- ifelse(df$Value == "Intercept", 5,
                            ifelse(df$Value == "RIN", 4,
                            ifelse(df$Value == "age", 3,
                            ifelse(df$Value == "sexM", 2, 1))))
   df
}

MIGs <- make_estimate_df("no_MIGs_theta_juncs", joined_filtered)

########################################################################
###################-----Variance Explained----##########################
########################################################################
joined$numbers <- rownames(joined)
title <- expression( "MIGs with" ~ theta ~ "Outliers")
x_lab <- expression(atop("Samples ordered by Number of MIGs",
                     "with Significant" ~ theta~ "Outliers"))
y_lab <- expression("Number of MIGs with Significant" ~ theta ~ "Outliers")

linear_variance_explained_gene_MIGs <- ggplot(MIGs, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.9), width = 0.2)+
  geom_text(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
                          ifelse(padjust < 0.001, "***", 
                            ifelse(padjust < 0.05, "*", "")))), 
            position = position_dodge(width = .5), vjust = 2, size = 30 / .pt, angle= 90) +
  labs(title_lab=title)+
  xlab(x_lab) +
  ylab(y_lab) +
  theme(axis.text.x = element_text(angle = 0, 
                                   vjust = 0.6, 
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal[c(1,3,5,7,9)])
  
linear_variance_explained_gene_MIGs
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/MIG_IR_metadata_analysis/plots/MIG_IR_metadata_analysis.pdf", plot=linear_variance_explained_gene_MIGs,  limitsize = FALSE, units = "in")



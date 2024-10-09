#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(GO.db)
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

################-----Get Stats-----#################
RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
RD268_zscore <- RD268/ sd(joined$no_MIGs_theta_juncs)

GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
GSS225379_zscore <- GSS225379/ sd(joined$no_MIGs_theta_juncs)

UDN550488.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_zscore <- UDN550488.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_theta_juncs)

UDN238929.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_zscore <- UDN238929.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_theta_juncs)

RD380 <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
RD380_zscore <- RD380/ sd(joined$no_MIGs_theta_juncs)

mean(RNU4ATAC$no_MIGs_theta_juncs)
sd(RNU4ATAC$no_MIGs_theta_juncs)

mean(non_RNU4ATAC$no_MIGs_theta_juncs)
sd(non_RNU4ATAC$no_MIGs_theta_juncs)

mean(RNU4ATAC$no_MIGs_theta_juncs)/mean(non_RNU4ATAC$no_MIGs_theta_juncs)

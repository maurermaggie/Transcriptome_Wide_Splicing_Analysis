library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(broom)
library(FRASER)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")
genes <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_filtered.rds")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
#####################-----Get Gene Dataframe-----#######################
########################################################################
#str(genes)
#row_ranges_df <- as.data.frame(genes@rowRanges) %>% filter(passed == TRUE)
#row_ranges <- genes@rowRanges
#findOverlaps
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_annotated <- annotateRangesWithTxDb(genes, txdb=txdb, orgDb=orgDb)

########################################################################
##########################-----Get Stats-----###########################
########################################################################

##########################-----Number Genes and Junctions-----###########################
junctions <- fds_annotated@rowRanges %>% as.data.frame
number_junctions <- nrow(junctions)
nubmer_genes <- junctions %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol) %>%  unique %>% length

##########################-----Get Stats-----###########################
psi3_mean <- mean(joined$Psi3_Junctions)
psi5_mean <- mean(joined$Psi3_Junctions)
theta_mean <- mean(joined$Psi3_Junctions)
all_mean <- mean(joined$Psi3_Junctions)

psi3_sd <- sd(joined$Psi3_Junctions)
psi5_sd <- sd(joined$Psi3_Junctions)
theta_sd <- sd(joined$Psi3_Junctions)
all_sd <- sd(joined$Psi3_Junctions)

two_sd_psi3 <- psi3_mean + psi3_sd + psi3_sd
two_sd_psi5 <- psi5_mean + psi5_sd + psi5_sd
two_sd_theta <- theta_mean + theta_sd + theta_sd
two_sd_all <- all_mean + all_sd + all_sd

psi3_outliers <- filter(joined, Psi3_Junctions >= two_sd_psi3) %>% pull(sampleID)
psi5_outliers <- filter(joined, Psi5_Junctions >= two_sd_psi5) %>% pull(sampleID)
theta_outliers <- filter(joined, Theta_Junctions >= two_sd_theta) %>% pull(sampleID)
all_outliers <- filter(joined, All_Junctions >= two_sd_all) %>% pull(sampleID)

psi3_outliers %>% length
psi5_outliers %>% length
theta_outliers %>% length
all_outliers %>% length

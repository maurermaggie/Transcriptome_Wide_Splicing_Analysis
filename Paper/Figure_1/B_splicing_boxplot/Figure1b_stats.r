require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")
genes <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_filtered.rds")


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

filtered_genes <- fds_annotated@rowRanges %>% as.data.frame %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol) %>% unique %>% length
filtered_junctions <- fds_annotated@rowRanges %>% as.data.frame %>% nrow

########################################################################
########################-----Get Stat Info-----#########################
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

filter_junctions <- function(dataframe, filter, column_name, significance){
      dataframe_filtered <- dataframe %>% filter(padjust <= significance) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
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

################-----Get Number by Gene per Person-----#################
mean(all_count_genes$All_Genes)
sd(all_count_genes$All_Genes)

mean(psi3_count_genes$Psi3_Genes)
sd(psi3_count_genes$Psi3_Genes)

mean(psi5_count_genes$Psi5_Genes)
sd(psi5_count_genes$Psi5_Genes)

mean(theta_count_genes$Theta_Genes)
sd(theta_count_genes$Theta_Genes)


################-----Get Number by Junction per Person-----#################
mean(all_count_junctions$All_Junctions)
sd(all_count_junctions$All_Junctions)

mean(psi3_count_junctions$Psi3_Junctions)
sd(psi3_count_junctions$Psi3_Junctions)

mean(psi5_count_junctions$Psi5_Junctions)
sd(psi5_count_junctions$Psi5_Junctions)

mean(theta_count_junctions$Theta_Junctions)
sd(theta_count_junctions$Theta_Junctions)

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

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")

args <- commandArgs(TRUE)
all_uncompiled_fp <- args[1]
low_RIN_fp <- args[2]
genesets_fp <- args[3]
genes_fp <- args[4]
output_file <- args[5]
output_file2 <- args[6]

FRASER <- args[7]

if (FRASER == "FRASER") {
    all_uncompiled <- read_csv(all_uncompiled_fp)
} else {
    all_uncompiled <- read_csv(all_uncompiled_fp, col_names=FALSE)
}

low_RIN <- read_csv(low_RIN_fp) %>% pull(sampleID)

autosomal_dominant_fp <- paste(genesets_fp, "/autosomal_dominant.tsv", sep="")
autosomal_dominant <- read_tsv(autosomal_dominant_fp, col_names = FALSE)

autosomal_recessive_fp <- paste(genesets_fp, "/autosomal_recessive.tsv", sep="")
autosomal_recessive <- read_tsv(autosomal_recessive_fp, col_names = FALSE)

CRISPR_essential_fp <- paste(genesets_fp, "/CRISPR_essential_genes.tsv", sep="")
CRISPR_essential <- read_tsv(CRISPR_essential_fp, col_names = FALSE)

CRISPR_nonessential_fp <- paste(genesets_fp, "/CRISPR_nonessential_genes.tsv", sep="")
CRISPR_nonessential <- read_tsv(CRISPR_nonessential_fp, col_names = FALSE)

OMIM_fp <- paste(genesets_fp, "/genemap2.tsv", sep="")
OMIM <- read_tsv(OMIM_fp)

homozygous_LoF_fp <- paste(genesets_fp, "/homozygous_Lof_genes.tsv", sep="")
homozygous_LoF <- read_tsv(homozygous_LoF_fp, col_names = FALSE)

Developmental_delay_fp <- paste(genesets_fp, "/developmental_delay_genes.csv", sep="")
Developmental_delay <- read_csv(Developmental_delay_fp)

olfactory_receptors_fp <- paste(genesets_fp, "/olfactory_receptors.tsv", sep="")
olfactory_receptors <- read_tsv(olfactory_receptors_fp, col_names = FALSE)

haplo_fp <- paste(genesets_fp, "/haploinsufficient.tsv", sep="")
haplo <- read_tsv(haplo_fp, col_names = FALSE)

display.brewer.all
pal <- brewer.pal(10, "Paired")

genes <- read_rds(genes_fp)

if (FRASER == "FRASER2") {
    colnames(all_uncompiled) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi",
                      "counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")
}

all_uncompiled <- all_uncompiled %>% filter(! sampleID %in% low_RIN)

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

filtered_genes <- fds_annotated@rowRanges %>% as.data.frame %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol)
junctions <- filtered_genes %>% length()
print("number junctions")
print(junctions)

genes <- filtered_genes %>% unique %>% length()
print("number genes")
print(genes)

########################################################################
##################-----Create Outlier Dataframe-----####################
########################################################################

################-----Get Number by Gene per Person-----#################
forest_plot_function <- function(outlier_genes, filtered_gene_list) {
            ################-----In olfactory or Not-----#################
            HAP <- haplo %>% pull(X1)
            in_HAP <- intersect(outlier_genes$hgncSymbol, HAP)
            not_in_HAP <- setdiff(outlier_genes$hgncSymbol, HAP)

            outlier_genes$HAP <- ifelse(outlier_genes$hgncSymbol %in% in_HAP, TRUE, FALSE)
            HAP_df <- outlier_genes

            ################-----In olfactory or Not-----#################
            OL <- olfactory_receptors %>% pull(X1)
            in_OL <- intersect(outlier_genes$hgncSymbol, OL)
            not_in_OL <- setdiff(outlier_genes$hgncSymbol, OL)

            outlier_genes$OL <- ifelse(outlier_genes$hgncSymbol %in% in_OL, TRUE, FALSE)
            OL_df <- outlier_genes %>% select(hgncSymbol, OL)

            ################-----In DD or Not-----#################
            DD <- Developmental_delay %>% pull(`gene symbol`)
            in_DD <- intersect(outlier_genes$hgncSymbol, DD)
            not_in_DD <- setdiff(outlier_genes$hgncSymbol, DD)

            outlier_genes$DD <- ifelse(outlier_genes$hgncSymbol %in% in_DD, TRUE, FALSE)
            DD_df <- outlier_genes %>% select(hgncSymbol, DD)

            ################-----In Disease or Not-----#################
            OMIM_genes <- OMIM %>% pull(`Gene/Locus And Other Related Symbols`)
            OMIM_genes_list <- unlist(strsplit(OMIM_genes, ","))
            OMIM_genes_cleaned <- sapply(OMIM_genes_list, trimws)
            in_disease <- intersect(outlier_genes$hgncSymbol, OMIM_genes_cleaned)
            not_in_disease <- setdiff(outlier_genes$hgncSymbol, OMIM_genes_cleaned)

            outlier_genes$OMIM_disease_gene <- ifelse(outlier_genes$hgncSymbol %in% in_disease, TRUE, FALSE)
            OMIM_df <- outlier_genes %>% select(hgncSymbol, OMIM_disease_gene)

            ################-----CRISPR Essential or Not-----#################
            essential <- intersect(CRISPR_essential$X1, outlier_genes$hgncSymbol)
            non_essential <- intersect(CRISPR_nonessential$X1, outlier_genes$hgncSymbol)

            outlier_genes$CRISPR_Essential <- ifelse(outlier_genes$hgncSymbol %in% essential, TRUE, FALSE)
            outlier_genes$CRISPR_Non_Essential <- ifelse(outlier_genes$hgncSymbol %in% non_essential, TRUE, FALSE)

            Essential <- outlier_genes %>% select(hgncSymbol, CRISPR_Essential)
            Non_Essential <- outlier_genes %>% select(hgncSymbol, CRISPR_Non_Essential)

            ################-----Inheritance Pattern-----#################
            Recessive <- data.frame(Gene = autosomal_recessive, Role = "Recessive")
            Dominant <- data.frame(Gene = autosomal_dominant, Role = "Dominant")

            Pattern <- bind_rows(Recessive, Dominant)

            outlier_in_pattern <- filter(outlier_genes, hgncSymbol %in% Pattern$X1) %>% select(-OMIM_disease_gene)

            recessive <- intersect(Recessive$X1, outlier_in_pattern$hgncSymbol)
            dominant <- intersect(Dominant$X1, outlier_in_pattern$hgncSymbol)

            outlier_in_pattern$Recessive_Pattern <- ifelse(outlier_in_pattern$hgncSymbol %in% recessive, TRUE, FALSE)
            Recessive_Pattern <- outlier_in_pattern %>% select(hgncSymbol, Recessive_Pattern)

            outlier_in_pattern$Dominant_Pattern <- ifelse(outlier_in_pattern$hgncSymbol %in% dominant, TRUE, FALSE)
            Dominant_Pattern <- outlier_in_pattern %>% select(hgncSymbol, Dominant_Pattern)

            ################-----Homozygous LoF-----#################
            #homo_LoF <- intersect(outlier_genes$hgncSymbol, homozygous_LoF$X1)
            #not_LoF <- setdiff(outlier_genes$hgncSymbol, homozygous_LoF$X1)

            #outlier_genes$Homozygous_LoF <- ifelse(outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            #Homozygous_LoF <- outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

            ########################################################################
            ################-----Create Non-Outlier Dataframe-----##################
            ########################################################################

            ################-----Get Number by Gene per Person-----#################
            all_genes <- as.data.frame(filtered_gene_list)
            colnames(all_genes) <- c("gene_symbol")
            non_outlier_genes <- setdiff(all_genes$gene_symbol, outlier_genes$hgncSymbol) %>% unique %>% data.frame
            colnames(non_outlier_genes) <- c("hgncSymbol")

            ################-----In OL or Not-----#################
            in_HAP <- intersect(non_outlier_genes$hgncSymbol, HAP)
            not_in_HAP <- setdiff(non_outlier_genes$hgncSymbol, HAP)

            non_outlier_genes$HAP <- ifelse(non_outlier_genes$hgncSymbol %in% in_HAP, TRUE, FALSE)
            HAP_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, HAP)

            ################-----In OL or Not-----#################
            in_OL <- intersect(non_outlier_genes$hgncSymbol, OL)
            not_in_OL <- setdiff(non_outlier_genes$hgncSymbol, OL)

            non_outlier_genes$OL <- ifelse(non_outlier_genes$hgncSymbol %in% in_OL, TRUE, FALSE)
            OL_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, OL)

            ################-----In DD or Not-----#################
            in_DD <- intersect(non_outlier_genes$hgncSymbol, DD)
            not_in_DD <- setdiff(non_outlier_genes$hgncSymbol, DD)

            non_outlier_genes$DD <- ifelse(non_outlier_genes$hgncSymbol %in% in_DD, TRUE, FALSE)
            DD_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, DD)

            ################-----In Disease or Not-----#################
            in_disease <- intersect(non_outlier_genes$hgncSymbol, OMIM_genes_cleaned)
            not_in_disease <- setdiff(non_outlier_genes$hgncSymbol, OMIM_genes_cleaned)

            non_outlier_genes$OMIM_disease_gene <- ifelse(non_outlier_genes$hgncSymbol %in% in_disease, TRUE, FALSE)
            OMIM_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, OMIM_disease_gene)

            ################-----CRISPR Essential or Not-----#################
            essential <- intersect(non_outlier_genes$hgncSymbol, CRISPR_essential$X1)
            non_essential <- setdiff(non_outlier_genes$hgncSymbol, CRISPR_nonessential$X1)

            non_outlier_genes$CRISPR_Essential <- ifelse(non_outlier_genes$hgncSymbol %in% essential, TRUE, FALSE)
            non_outlier_genes$CRISPR_Non_Essential <- ifelse(non_outlier_genes$hgncSymbol %in% non_essential, TRUE, FALSE)

            Essential_nonoutlier <- non_outlier_genes %>% select(hgncSymbol, CRISPR_Essential)
            Non_Essential_nonoutlier <- non_outlier_genes %>% select(hgncSymbol, CRISPR_Non_Essential)

            ################-----Inheritance Pattern-----#################
            nonoutlier_in_pattern <- filter(non_outlier_genes, hgncSymbol %in% Pattern$X1) %>% select(-OMIM_disease_gene)

            recessive <- intersect(Recessive$X1, nonoutlier_in_pattern$hgncSymbol)
            dominant <- intersect(Dominant$X1, nonoutlier_in_pattern$hgncSymbol)

            nonoutlier_in_pattern$Recessive_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% recessive, TRUE, FALSE)
            Recessive_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Recessive_Pattern)

            nonoutlier_in_pattern$Dominant_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% dominant, TRUE, FALSE)
            Dominant_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Dominant_Pattern)

            ################-----Homozygous LoF-----#################
            #homo_LoF <- intersect(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)
            #not_LoF <- setdiff(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)

            #non_outlier_genes$Homozygous_LoF <- ifelse(non_outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            #Homozygous_LoF_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

            ########################################################################
            #####################-----Combine Dataframes-----#######################
            ########################################################################
            ################-----In Disease or Not-----#################
            HAP_df_nonoutliers$Splicing_Outlier <- FALSE
            HAP_df$Splicing_Outlier <- TRUE

            HAP_joined <- bind_rows(HAP_df, HAP_df_nonoutliers)

            HAP_Splicing <- filter(HAP_joined, Splicing_Outlier == TRUE) %>% filter(HAP == TRUE) %>% nrow
            HAP_Non_Splicing <- filter(HAP_joined, Splicing_Outlier != TRUE) %>% filter(HAP == TRUE) %>% nrow
            Non_HAP_Splicing <- filter(HAP_joined, Splicing_Outlier == TRUE) %>% filter(HAP != TRUE) %>% nrow
            Non_HAP_Non_Splicing <- filter(HAP_joined, Splicing_Outlier != TRUE) %>% filter(HAP != TRUE) %>% nrow

            HAP_ct <- data.frame(
                "Splicing_Outlier" = c(HAP_Splicing, Non_HAP_Splicing),
                "Spilicing_Inlier" = c(HAP_Non_Splicing, Non_HAP_Non_Splicing),
                row.names = c("HAP_Disease_Gene", "Non_HAP_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(HAP_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In Disease or Not-----#################
            OL_df_nonoutliers$Splicing_Outlier <- FALSE
            OL_df$Splicing_Outlier <- TRUE

            OL_joined <- bind_rows(OL_df, OL_df_nonoutliers)

            OL_Splicing <- filter(OL_joined, Splicing_Outlier == TRUE) %>% filter(OL == TRUE) %>% nrow
            OL_Non_Splicing <- filter(OL_joined, Splicing_Outlier != TRUE) %>% filter(OL == TRUE) %>% nrow
            Non_OL_Splicing <- filter(OL_joined, Splicing_Outlier == TRUE) %>% filter(OL != TRUE) %>% nrow
            Non_OL_Non_Splicing <- filter(OL_joined, Splicing_Outlier != TRUE) %>% filter(OL != TRUE) %>% nrow

            OL_ct <- data.frame(
                "Splicing_Outlier" = c(OL_Splicing, Non_OL_Splicing),
                "Spilicing_Inlier" = c(OL_Non_Splicing, Non_OL_Non_Splicing),
                row.names = c("OL_Disease_Gene", "Non_OL_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(OL_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In DD or Not-----#################
            DD_df_nonoutliers$Splicing_Outlier <- FALSE
            DD_df$Splicing_Outlier <- TRUE

            DD_joined <- bind_rows(DD_df, DD_df_nonoutliers)

            DD_Splicing <- filter(DD_joined, Splicing_Outlier == TRUE) %>% filter(DD == TRUE) %>% nrow
            DD_Non_Splicing <- filter(DD_joined, Splicing_Outlier != TRUE) %>% filter(DD == TRUE) %>% nrow
            Non_DD_Splicing <- filter(DD_joined, Splicing_Outlier == TRUE) %>% filter(DD != TRUE) %>% nrow
            Non_DD_Non_Splicing <- filter(DD_joined, Splicing_Outlier != TRUE) %>% filter(DD != TRUE) %>% nrow

            DD_ct <- data.frame(
                "Splicing_Outlier" = c(DD_Splicing, Non_DD_Splicing),
                "Spilicing_Inlier" = c(DD_Non_Splicing, Non_DD_Non_Splicing),
                row.names = c("DD_Disease_Gene", "Non_DD_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(DD_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In Disease or Not-----#################
            OMIM_df_nonoutliers$Splicing_Outlier <- FALSE
            OMIM_df$Splicing_Outlier <- TRUE

            OMIM <- bind_rows(OMIM_df, OMIM_df_nonoutliers)

            OMIM_Splicing <- filter(OMIM, Splicing_Outlier == TRUE) %>% filter(OMIM_disease_gene == TRUE) %>% nrow
            OMIM_Non_Splicing <- filter(OMIM, Splicing_Outlier != TRUE) %>% filter(OMIM_disease_gene == TRUE) %>% nrow
            Non_OMIM_Splicing <- filter(OMIM, Splicing_Outlier == TRUE) %>% filter(OMIM_disease_gene != TRUE) %>% nrow
            Non_OMIM_Non_Splicing <- filter(OMIM, Splicing_Outlier != TRUE) %>% filter(OMIM_disease_gene != TRUE) %>% nrow

            OMIM_ct <- data.frame(
                "Splicing_Outlier" = c(OMIM_Splicing, Non_OMIM_Splicing),
                "Spilicing_Inlier" = c(OMIM_Non_Splicing, Non_OMIM_Non_Splicing),
                row.names = c("OMIM_Disease_Gene", "Non_OMIN_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(OMIM_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----CRISPR Essential or Not-----#################
            Essential_nonoutlier$Splicing_Outlier <- FALSE
            Essential$Splicing_Outlier <- TRUE

            CRISPR_Essentials <- bind_rows(Essential_nonoutlier, Essential)
            colnames(CRISPR_Essentials) <- c("hgncSymbol", "CRISPR_Essential", "Splicing_Outlier")

            Essential_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essential == TRUE) %>% nrow
            Essential_Non_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essential == TRUE) %>% nrow
            Non_Essential_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essential != TRUE) %>% nrow
            Non_Essential_Non_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essential != TRUE) %>% nrow

            Essential_ct <- data.frame(
                "Splicing_Outlier" = c(Essential_Splicing, Non_Essential_Splicing),
                "Spilicing_Inlier" = c(Essential_Non_Splicing, Non_Essential_Non_Splicing),
                row.names = c("CRISPR_Essential_Gene", "CRISPR_Not_Essential_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(Essential_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")


            Non_Essential_nonoutlier$Splicing_Outlier <- FALSE
            Non_Essential$Splicing_Outlier <- TRUE

            CRISPR_Non_Essentials <- bind_rows(Non_Essential_nonoutlier, Non_Essential)
            colnames(CRISPR_Non_Essentials) <- c("hgncSymbol", "CRISPR_Non_Essential", "Splicing_Outlier")

            Non_Essential_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Non_Essential == TRUE) %>% nrow
            Non_Essential_Non_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Non_Essential == TRUE) %>% nrow
            Essential_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Non_Essential != TRUE) %>% nrow
            Essential_Non_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Non_Essential != TRUE) %>% nrow

            Non_Essential_ct <- data.frame(
                "Splicing_Outlier" = c(Non_Essential_Splicing, Essential_Splicing),
                "Spilicing_Inlier" = c(Non_Essential_Non_Splicing, Essential_Non_Splicing),
                row.names = c("CRISPR_Non_Essential_Gene", "CRISPR_Not_Non_Essential_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(Non_Essential_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----Inheritance Pattern-----#################
            Recessive_Pattern_nonoutliers$Splicing_Outliers <- FALSE
            Recessive_Pattern$Splicing_Outliers <- TRUE

            Recessive <- bind_rows(Recessive_Pattern, Recessive_Pattern_nonoutliers)

            Recessive_Splicing <- filter(Recessive, Splicing_Outliers == TRUE) %>% filter(Recessive_Pattern == TRUE) %>% nrow
            Recessive_Non_Splicing <- filter(Recessive, Splicing_Outliers != TRUE) %>% filter(Recessive_Pattern == TRUE) %>% nrow
            Non_Recessive_Splicing <- filter(Recessive, Splicing_Outliers == TRUE) %>% filter(Recessive_Pattern != TRUE) %>% nrow
            Non_Recessive_Non_Splicing <- filter(Recessive, Splicing_Outliers != TRUE) %>% filter(Recessive_Pattern != TRUE) %>% nrow

            Recessive_ct <- data.frame(
                "Splicing_Outlier" = c(Recessive_Splicing, Non_Recessive_Splicing),
                "Spilicing_Inlier" = c(Recessive_Non_Splicing, Non_Recessive_Non_Splicing),
                row.names = c("Recessive_Inheritance_Pattern", "Non_Recessive_Inheritance_Pattern"),
                stringsAsFactors = FALSE
            )

            colnames(Recessive_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            Dominant_Pattern_nonoutliers$Splicing_Outlier <- FALSE
            Dominant_Pattern$Splicing_Outlier <- TRUE 

            Dominant <- bind_rows(Dominant_Pattern, Dominant_Pattern_nonoutliers)

            Dominant_Splicing <- filter(Dominant, Splicing_Outlier == TRUE) %>% filter(Dominant_Pattern == TRUE) %>% nrow
            Dominant_Non_Splicing <- filter(Dominant, Splicing_Outlier != TRUE) %>% filter(Dominant_Pattern == TRUE) %>% nrow
            Non_Dominant_Splicing <- filter(Dominant, Splicing_Outlier == TRUE) %>% filter(Dominant_Pattern != TRUE) %>% nrow
            Non_Dominant_Non_Splicing <- filter(Dominant, Splicing_Outlier != TRUE) %>% filter(Dominant_Pattern != TRUE) %>% nrow

            Dominant_ct <- data.frame(
                "Splicing_Outlier" = c(Dominant_Splicing, Non_Dominant_Splicing),
                "Spilicing_Inlier" = c(Dominant_Non_Splicing, Non_Dominant_Non_Splicing),
                row.names = c("Dominant_Inheritance_Pattern", "Non_Dominant_Inheritance_Pattern"),
                stringsAsFactors = FALSE
            )

            colnames(Dominant_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----Homozygous LoF-----#################
            #Homozygous_LoF_nonoutliers$Splicing_Pattern <- FALSE
            #Homozygous_LoF$Splicing_Pattern <- TRUE

            #Homozygous_Loss_of_Function <- bind_rows(Homozygous_LoF_nonoutliers, Homozygous_LoF)

            #LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            #LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            #Non_LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow
            #Non_LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow

            #LoF_ct <- data.frame(
            #    "Splicing_Outlier" = c(LoF_Splicing, Non_LoF_Splicing),
            #    "Spilicing_Inlier" = c(LoF_Non_Splicing, Non_LoF_Non_Splicing),
            #    row.names = c("Homozygous_Loss_of_Function", "Non_Homozygous_Loss_of_Function"),
            #    stringsAsFactors = FALSE
            #)

            #colnames(LoF_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ########################################################################
            #####################-----Statistical Tests-----########################
            ########################################################################
            ################-----In OL or Not-----#################
            HAP_f <- fisher.test(HAP_ct)
            HAP_tidy <- tidy(HAP_f)

            ################-----In OL or Not-----#################
            OL_f <- fisher.test(OL_ct)
            OL_tidy <- tidy(OL_f)

            ################-----In DD or Not-----#################
            DD_f <- fisher.test(DD_ct)
            DD_tidy <- tidy(DD_f)

            ################-----In Disease or Not-----#################
            OMIM_f <- fisher.test(OMIM_ct)
            OMIM_tidy <- tidy(OMIM_f)

            ################-----CRISPR Essential or Not-----#################
            Essential_f <- fisher.test(Essential_ct)
            Essential_tidy <- tidy(Essential_f)
           
            Non_Essential_f <- fisher.test(Non_Essential_ct)
            Non_Essential_tidy <- tidy(Non_Essential_f)

            ################-----Inheritance Pattern-----#################
            Recessive_f <- fisher.test(Recessive_ct)
            Recessive_tidy <- tidy(Recessive_f)

            Dominant_f <- fisher.test(Dominant_ct)
            Dominant_tidy <- tidy(Dominant_f)

            ################-----Homozygous LoF-----#################
            #LoF_f <- fisher.test(LoF_ct)
            #LoF_tidy <- tidy(LoF_f)

            ########################################################################
            ##########################-----Plots-----###############################
            ########################################################################
            #pvalue <- c(HAP_tidy$p.value, OL_tidy$p.value, DD_tidy$p.value, OMIM_tidy$p.value, Essential_tidy$p.value, Non_Essential_tidy$p.value, Recessive_tidy$p.value, Dominant_tidy$p.value, LoF_tidy$p.value)
            #mean <- c((HAP_tidy$estimate), (OL_tidy$estimate), (DD_tidy$estimate), (OMIM_tidy$estimate), (Essential_tidy$estimate), (Non_Essential_tidy$estimate), (Recessive_tidy$estimate), (Dominant_tidy$estimate), (LoF_tidy$estimate))
            #conf.low <- c((HAP_tidy$conf.low), (OL_tidy$conf.low), (DD_tidy$conf.low), (OMIM_tidy$conf.low), (Essential_tidy$conf.low), (Non_Essential_tidy$estimate), (Recessive_tidy$conf.low), (Dominant_tidy$conf.low), (LoF_tidy$conf.low))
            #conf.high <- c((HAP_tidy$conf.high), (OL_tidy$conf.high), (DD_tidy$conf.high), (OMIM_tidy$conf.high), (Essential_tidy$conf.high), (Non_Essential_tidy$estimate), (Recessive_tidy$conf.high), (Dominant_tidy$conf.high), (LoF_tidy$conf.high))

            #LoF_Tag <- paste("Homozygous Loss of", "\n", "Function Tolerant")
            #CRISPR_Essential_Tag <- paste("CRISPR", "\n", "Essential")
            #CRISPR_Non_Essential_Tag <- paste("CRISPR", "\n", "Non Essential")
            #OL <- paste("Olfactory", "\n", "receptor")
            #DD <- paste("Developmental", "\n", "Delay")
            #values <- c("Haploinsufficient", OL, DD, "OMIM", CRISPR_Essential_Tag, CRISPR_Non_Essential_Tag, "Recessive", "Dominant", LoF_Tag)

            pvalue <- c(OMIM_tidy$p.value, Recessive_tidy$p.value, Dominant_tidy$p.value, HAP_tidy$p.value, DD_tidy$p.value, Non_Essential_tidy$p.value, OL_tidy$p.value)
            mean <- c((OMIM_tidy$estimate), (Recessive_tidy$estimate), (Dominant_tidy$estimate), (HAP_tidy$estimate), (DD_tidy$estimate), (Non_Essential_tidy$estimate), (OL_tidy$estimate))
            conf.low <- c((OMIM_tidy$conf.low), (Recessive_tidy$conf.low), (Dominant_tidy$conf.low), (HAP_tidy$conf.low), (DD_tidy$conf.low), (Non_Essential_tidy$conf.low), (OL_tidy$conf.low))
            conf.high <- c((OMIM_tidy$conf.high), (Recessive_tidy$conf.high), (Dominant_tidy$conf.high), (HAP_tidy$conf.high), (DD_tidy$estimate), (Non_Essential_tidy$conf.high), (OL_tidy$conf.high))

            #LoF_Tag <- paste("Homozygous Loss of", "\n", "Function Tolerant")
            CRISPR_Non_Essential_Tag <- paste("CRISPR", "\n", "Non Essential Gene")
            OL <- paste("Olfactory", "\n", "Receptor Gene")
            DD <- paste("Developmental", "\n", "Delay Gene")
            Rec <- paste("Recessive", "\n", "Disease Gene")
            Dom <- paste("Dominant", "\n", "Disease Gene")
            OMIM <- paste("OMIM", "\n", "Disease Gene")
            Hap <- paste("Haploinsufficient", "\n", "Disease Gene")
            values <- c(OMIM, Rec, Dom, Hap, DD, CRISPR_Non_Essential_Tag, OL)
            number <- c(7,6,5,4,3,2,1)

            forest_plot_data <- data.frame(values, pvalue, mean, conf.low, conf.high)
            forest_plot_data$padjust <- p.adjust(forest_plot_data$pvalue, method="fdr", n=length(forest_plot_data$pvalue))

            forest_plot_odds_ratio <- ggplot(data=forest_plot_data, aes(colour = values, x=reorder(values, number), y=mean, ymin=conf.low, ymax=conf.high)) +
                    geom_pointrange(size=1.5, show.legend=FALSE) + 
                    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
                    coord_flip() +  # flip coordinates (puts labels on y axis)
                    theme_classic(base_size = 30)+
                    scale_colour_manual(values = pal[c(3,2,10,5,1, 8,9)]) +
                    ylab("Odds Ratio") +
                    theme(axis.title.y=element_blank()) +
                    ylim(0, 4)# use a white background
            forest_plot_odds_ratio
}

outlier_genes_05 <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% select(hgncSymbol) %>% filter(!is.na(hgncSymbol)) %>% unique 
outlier_genes_05_fp <- forest_plot_function(outlier_genes_05, filtered_genes)
ggsave(filename=output_file, plot=outlier_genes_05_fp,  limitsize = FALSE, units = "in", height=12, width=10)

fisher_function <- function(outlier_genes, filtered_gene_list) {
            ################-----In olfactory or Not-----#################
            HAP <- haplo %>% pull(X1)
            in_HAP <- intersect(outlier_genes$hgncSymbol, HAP)
            not_in_HAP <- setdiff(outlier_genes$hgncSymbol, HAP)

            outlier_genes$HAP <- ifelse(outlier_genes$hgncSymbol %in% in_HAP, TRUE, FALSE)
            HAP_df <- outlier_genes

            ################-----In olfactory or Not-----#################
            OL <- olfactory_receptors %>% pull(X1)
            in_OL <- intersect(outlier_genes$hgncSymbol, OL)
            not_in_OL <- setdiff(outlier_genes$hgncSymbol, OL)

            outlier_genes$OL <- ifelse(outlier_genes$hgncSymbol %in% in_OL, TRUE, FALSE)
            OL_df <- outlier_genes %>% select(hgncSymbol, OL)

            ################-----In DD or Not-----#################
            DD <- Developmental_delay %>% pull(`gene symbol`)
            in_DD <- intersect(outlier_genes$hgncSymbol, DD)
            not_in_DD <- setdiff(outlier_genes$hgncSymbol, DD)

            outlier_genes$DD <- ifelse(outlier_genes$hgncSymbol %in% in_DD, TRUE, FALSE)
            DD_df <- outlier_genes %>% select(hgncSymbol, DD)

            ################-----In Disease or Not-----#################
            OMIM_genes <- OMIM %>% pull(`Gene/Locus And Other Related Symbols`)
            OMIM_genes_list <- unlist(strsplit(OMIM_genes, ","))
            OMIM_genes_cleaned <- sapply(OMIM_genes_list, trimws)
            in_disease <- intersect(outlier_genes$hgncSymbol, OMIM_genes_cleaned)
            not_in_disease <- setdiff(outlier_genes$hgncSymbol, OMIM_genes_cleaned)

            outlier_genes$OMIM_disease_gene <- ifelse(outlier_genes$hgncSymbol %in% in_disease, TRUE, FALSE)
            OMIM_df <- outlier_genes %>% select(hgncSymbol, OMIM_disease_gene)

            ################-----CRISPR Essential or Not-----#################
            essential <- intersect(CRISPR_essential$X1, outlier_genes$hgncSymbol)
            non_essential <- intersect(CRISPR_nonessential$X1, outlier_genes$hgncSymbol)

            outlier_genes$CRISPR_Essential <- ifelse(outlier_genes$hgncSymbol %in% essential, TRUE, FALSE)
            outlier_genes$CRISPR_Non_Essential <- ifelse(outlier_genes$hgncSymbol %in% non_essential, TRUE, FALSE)

            Essential <- outlier_genes %>% select(hgncSymbol, CRISPR_Essential)
            Non_Essential <- outlier_genes %>% select(hgncSymbol, CRISPR_Non_Essential)

            ################-----Inheritance Pattern-----#################
            Recessive <- data.frame(Gene = autosomal_recessive, Role = "Recessive")
            Dominant <- data.frame(Gene = autosomal_dominant, Role = "Dominant")

            Pattern <- bind_rows(Recessive, Dominant)

            outlier_in_pattern <- filter(outlier_genes, hgncSymbol %in% Pattern$X1) %>% select(-OMIM_disease_gene)

            recessive <- intersect(Recessive$X1, outlier_in_pattern$hgncSymbol)
            dominant <- intersect(Dominant$X1, outlier_in_pattern$hgncSymbol)

            outlier_in_pattern$Recessive_Pattern <- ifelse(outlier_in_pattern$hgncSymbol %in% recessive, TRUE, FALSE)
            Recessive_Pattern <- outlier_in_pattern %>% select(hgncSymbol, Recessive_Pattern)

            outlier_in_pattern$Dominant_Pattern <- ifelse(outlier_in_pattern$hgncSymbol %in% dominant, TRUE, FALSE)
            Dominant_Pattern <- outlier_in_pattern %>% select(hgncSymbol, Dominant_Pattern)

            ################-----Homozygous LoF-----#################
            #homo_LoF <- intersect(outlier_genes$hgncSymbol, homozygous_LoF$X1)
            #not_LoF <- setdiff(outlier_genes$hgncSymbol, homozygous_LoF$X1)

            #outlier_genes$Homozygous_LoF <- ifelse(outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            #Homozygous_LoF <- outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

            ########################################################################
            ################-----Create Non-Outlier Dataframe-----##################
            ########################################################################

            ################-----Get Number by Gene per Person-----#################
            all_genes <- as.data.frame(filtered_gene_list)
            colnames(all_genes) <- c("gene_symbol")
            non_outlier_genes <- setdiff(all_genes$gene_symbol, outlier_genes$hgncSymbol) %>% unique %>% data.frame
            colnames(non_outlier_genes) <- c("hgncSymbol")

            ################-----In OL or Not-----#################
            in_HAP <- intersect(non_outlier_genes$hgncSymbol, HAP)
            not_in_HAP <- setdiff(non_outlier_genes$hgncSymbol, HAP)

            non_outlier_genes$HAP <- ifelse(non_outlier_genes$hgncSymbol %in% in_HAP, TRUE, FALSE)
            HAP_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, HAP)

            ################-----In OL or Not-----#################
            in_OL <- intersect(non_outlier_genes$hgncSymbol, OL)
            not_in_OL <- setdiff(non_outlier_genes$hgncSymbol, OL)

            non_outlier_genes$OL <- ifelse(non_outlier_genes$hgncSymbol %in% in_OL, TRUE, FALSE)
            OL_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, OL)

            ################-----In DD or Not-----#################
            in_DD <- intersect(non_outlier_genes$hgncSymbol, DD)
            not_in_DD <- setdiff(non_outlier_genes$hgncSymbol, DD)

            non_outlier_genes$DD <- ifelse(non_outlier_genes$hgncSymbol %in% in_DD, TRUE, FALSE)
            DD_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, DD)

            ################-----In Disease or Not-----#################
            in_disease <- intersect(non_outlier_genes$hgncSymbol, OMIM_genes_cleaned)
            not_in_disease <- setdiff(non_outlier_genes$hgncSymbol, OMIM_genes_cleaned)

            non_outlier_genes$OMIM_disease_gene <- ifelse(non_outlier_genes$hgncSymbol %in% in_disease, TRUE, FALSE)
            OMIM_df_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, OMIM_disease_gene)

            ################-----CRISPR Essential or Not-----#################
            essential <- intersect(non_outlier_genes$hgncSymbol, CRISPR_essential$X1)
            non_essential <- setdiff(non_outlier_genes$hgncSymbol, CRISPR_nonessential$X1)

            non_outlier_genes$CRISPR_Essential <- ifelse(non_outlier_genes$hgncSymbol %in% essential, TRUE, FALSE)
            non_outlier_genes$CRISPR_Non_Essential <- ifelse(non_outlier_genes$hgncSymbol %in% non_essential, TRUE, FALSE)

            Essential_nonoutlier <- non_outlier_genes %>% select(hgncSymbol, CRISPR_Essential)
            Non_Essential_nonoutlier <- non_outlier_genes %>% select(hgncSymbol, CRISPR_Non_Essential)

            ################-----Inheritance Pattern-----#################
            nonoutlier_in_pattern <- filter(non_outlier_genes, hgncSymbol %in% Pattern$X1) %>% select(-OMIM_disease_gene)

            recessive <- intersect(Recessive$X1, nonoutlier_in_pattern$hgncSymbol)
            dominant <- intersect(Dominant$X1, nonoutlier_in_pattern$hgncSymbol)

            nonoutlier_in_pattern$Recessive_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% recessive, TRUE, FALSE)
            Recessive_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Recessive_Pattern)

            nonoutlier_in_pattern$Dominant_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% dominant, TRUE, FALSE)
            Dominant_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Dominant_Pattern)

            ################-----Homozygous LoF-----#################
            #homo_LoF <- intersect(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)
            #not_LoF <- setdiff(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)

            #non_outlier_genes$Homozygous_LoF <- ifelse(non_outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            #Homozygous_LoF_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

            ########################################################################
            #####################-----Combine Dataframes-----#######################
            ########################################################################
            ################-----In Disease or Not-----#################
            HAP_df_nonoutliers$Splicing_Outlier <- FALSE
            HAP_df$Splicing_Outlier <- TRUE

            HAP_joined <- bind_rows(HAP_df, HAP_df_nonoutliers)

            HAP_Splicing <- filter(HAP_joined, Splicing_Outlier == TRUE) %>% filter(HAP == TRUE) %>% nrow
            HAP_Non_Splicing <- filter(HAP_joined, Splicing_Outlier != TRUE) %>% filter(HAP == TRUE) %>% nrow
            Non_HAP_Splicing <- filter(HAP_joined, Splicing_Outlier == TRUE) %>% filter(HAP != TRUE) %>% nrow
            Non_HAP_Non_Splicing <- filter(HAP_joined, Splicing_Outlier != TRUE) %>% filter(HAP != TRUE) %>% nrow

            HAP_ct <- data.frame(
                "Splicing_Outlier" = c(HAP_Splicing, Non_HAP_Splicing),
                "Spilicing_Inlier" = c(HAP_Non_Splicing, Non_HAP_Non_Splicing),
                row.names = c("HAP_Disease_Gene", "Non_HAP_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(HAP_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In Disease or Not-----#################
            OL_df_nonoutliers$Splicing_Outlier <- FALSE
            OL_df$Splicing_Outlier <- TRUE

            OL_joined <- bind_rows(OL_df, OL_df_nonoutliers)

            OL_Splicing <- filter(OL_joined, Splicing_Outlier == TRUE) %>% filter(OL == TRUE) %>% nrow
            OL_Non_Splicing <- filter(OL_joined, Splicing_Outlier != TRUE) %>% filter(OL == TRUE) %>% nrow
            Non_OL_Splicing <- filter(OL_joined, Splicing_Outlier == TRUE) %>% filter(OL != TRUE) %>% nrow
            Non_OL_Non_Splicing <- filter(OL_joined, Splicing_Outlier != TRUE) %>% filter(OL != TRUE) %>% nrow

            OL_ct <- data.frame(
                "Splicing_Outlier" = c(OL_Splicing, Non_OL_Splicing),
                "Spilicing_Inlier" = c(OL_Non_Splicing, Non_OL_Non_Splicing),
                row.names = c("OL_Disease_Gene", "Non_OL_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(OL_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In DD or Not-----#################
            DD_df_nonoutliers$Splicing_Outlier <- FALSE
            DD_df$Splicing_Outlier <- TRUE

            DD_joined <- bind_rows(DD_df, DD_df_nonoutliers)

            DD_Splicing <- filter(DD_joined, Splicing_Outlier == TRUE) %>% filter(DD == TRUE) %>% nrow
            DD_Non_Splicing <- filter(DD_joined, Splicing_Outlier != TRUE) %>% filter(DD == TRUE) %>% nrow
            Non_DD_Splicing <- filter(DD_joined, Splicing_Outlier == TRUE) %>% filter(DD != TRUE) %>% nrow
            Non_DD_Non_Splicing <- filter(DD_joined, Splicing_Outlier != TRUE) %>% filter(DD != TRUE) %>% nrow

            DD_ct <- data.frame(
                "Splicing_Outlier" = c(DD_Splicing, Non_DD_Splicing),
                "Spilicing_Inlier" = c(DD_Non_Splicing, Non_DD_Non_Splicing),
                row.names = c("DD_Disease_Gene", "Non_DD_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(DD_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----In Disease or Not-----#################
            OMIM_df_nonoutliers$Splicing_Outlier <- FALSE
            OMIM_df$Splicing_Outlier <- TRUE

            OMIM <- bind_rows(OMIM_df, OMIM_df_nonoutliers)

            OMIM_Splicing <- filter(OMIM, Splicing_Outlier == TRUE) %>% filter(OMIM_disease_gene == TRUE) %>% nrow
            OMIM_Non_Splicing <- filter(OMIM, Splicing_Outlier != TRUE) %>% filter(OMIM_disease_gene == TRUE) %>% nrow
            Non_OMIM_Splicing <- filter(OMIM, Splicing_Outlier == TRUE) %>% filter(OMIM_disease_gene != TRUE) %>% nrow
            Non_OMIM_Non_Splicing <- filter(OMIM, Splicing_Outlier != TRUE) %>% filter(OMIM_disease_gene != TRUE) %>% nrow

            OMIM_ct <- data.frame(
                "Splicing_Outlier" = c(OMIM_Splicing, Non_OMIM_Splicing),
                "Spilicing_Inlier" = c(OMIM_Non_Splicing, Non_OMIM_Non_Splicing),
                row.names = c("OMIM_Disease_Gene", "Non_OMIN_Disease_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(OMIM_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----CRISPR Essential or Not-----#################
            Essential_nonoutlier$Splicing_Outlier <- FALSE
            Essential$Splicing_Outlier <- TRUE

            CRISPR_Essentials <- bind_rows(Essential_nonoutlier, Essential)
            colnames(CRISPR_Essentials) <- c("hgncSymbol", "CRISPR_Essential", "Splicing_Outlier")

            Essential_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essential == TRUE) %>% nrow
            Essential_Non_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essential == TRUE) %>% nrow
            Non_Essential_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essential != TRUE) %>% nrow
            Non_Essential_Non_Splicing <- filter(CRISPR_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essential != TRUE) %>% nrow

            Essential_ct <- data.frame(
                "Splicing_Outlier" = c(Essential_Splicing, Non_Essential_Splicing),
                "Spilicing_Inlier" = c(Essential_Non_Splicing, Non_Essential_Non_Splicing),
                row.names = c("CRISPR_Essential_Gene", "CRISPR_Not_Essential_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(Essential_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")


            Non_Essential_nonoutlier$Splicing_Outlier <- FALSE
            Non_Essential$Splicing_Outlier <- TRUE

            CRISPR_Non_Essentials <- bind_rows(Non_Essential_nonoutlier, Non_Essential)
            colnames(CRISPR_Non_Essentials) <- c("hgncSymbol", "CRISPR_Non_Essential", "Splicing_Outlier")

            Non_Essential_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Non_Essential == TRUE) %>% nrow
            Non_Essential_Non_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Non_Essential == TRUE) %>% nrow
            Essential_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier == TRUE) %>% filter(CRISPR_Non_Essential != TRUE) %>% nrow
            Essential_Non_Splicing <- filter(CRISPR_Non_Essentials, Splicing_Outlier != TRUE) %>% filter(CRISPR_Non_Essential != TRUE) %>% nrow

            Non_Essential_ct <- data.frame(
                "Splicing_Outlier" = c(Non_Essential_Splicing, Essential_Splicing),
                "Spilicing_Inlier" = c(Non_Essential_Non_Splicing, Essential_Non_Splicing),
                row.names = c("CRISPR_Non_Essential_Gene", "CRISPR_Not_Non_Essential_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(Non_Essential_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----Inheritance Pattern-----#################
            Recessive_Pattern_nonoutliers$Splicing_Outliers <- FALSE
            Recessive_Pattern$Splicing_Outliers <- TRUE

            Recessive <- bind_rows(Recessive_Pattern, Recessive_Pattern_nonoutliers)

            Recessive_Splicing <- filter(Recessive, Splicing_Outliers == TRUE) %>% filter(Recessive_Pattern == TRUE) %>% nrow
            Recessive_Non_Splicing <- filter(Recessive, Splicing_Outliers != TRUE) %>% filter(Recessive_Pattern == TRUE) %>% nrow
            Non_Recessive_Splicing <- filter(Recessive, Splicing_Outliers == TRUE) %>% filter(Recessive_Pattern != TRUE) %>% nrow
            Non_Recessive_Non_Splicing <- filter(Recessive, Splicing_Outliers != TRUE) %>% filter(Recessive_Pattern != TRUE) %>% nrow

            Recessive_ct <- data.frame(
                "Splicing_Outlier" = c(Recessive_Splicing, Non_Recessive_Splicing),
                "Spilicing_Inlier" = c(Recessive_Non_Splicing, Non_Recessive_Non_Splicing),
                row.names = c("Recessive_Inheritance_Pattern", "Non_Recessive_Inheritance_Pattern"),
                stringsAsFactors = FALSE
            )

            colnames(Recessive_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            Dominant_Pattern_nonoutliers$Splicing_Outlier <- FALSE
            Dominant_Pattern$Splicing_Outlier <- TRUE 

            Dominant <- bind_rows(Dominant_Pattern, Dominant_Pattern_nonoutliers)

            Dominant_Splicing <- filter(Dominant, Splicing_Outlier == TRUE) %>% filter(Dominant_Pattern == TRUE) %>% nrow
            Dominant_Non_Splicing <- filter(Dominant, Splicing_Outlier != TRUE) %>% filter(Dominant_Pattern == TRUE) %>% nrow
            Non_Dominant_Splicing <- filter(Dominant, Splicing_Outlier == TRUE) %>% filter(Dominant_Pattern != TRUE) %>% nrow
            Non_Dominant_Non_Splicing <- filter(Dominant, Splicing_Outlier != TRUE) %>% filter(Dominant_Pattern != TRUE) %>% nrow

            Dominant_ct <- data.frame(
                "Splicing_Outlier" = c(Dominant_Splicing, Non_Dominant_Splicing),
                "Spilicing_Inlier" = c(Dominant_Non_Splicing, Non_Dominant_Non_Splicing),
                row.names = c("Dominant_Inheritance_Pattern", "Non_Dominant_Inheritance_Pattern"),
                stringsAsFactors = FALSE
            )

            colnames(Dominant_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ################-----Homozygous LoF-----#################
            #Homozygous_LoF_nonoutliers$Splicing_Pattern <- FALSE
            #Homozygous_LoF$Splicing_Pattern <- TRUE

            #Homozygous_Loss_of_Function <- bind_rows(Homozygous_LoF_nonoutliers, Homozygous_LoF)

            #LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            #LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            #Non_LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow
            #Non_LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow

            #LoF_ct <- data.frame(
            #    "Splicing_Outlier" = c(LoF_Splicing, Non_LoF_Splicing),
            #    "Spilicing_Inlier" = c(LoF_Non_Splicing, Non_LoF_Non_Splicing),
            #    row.names = c("Homozygous_Loss_of_Function", "Non_Homozygous_Loss_of_Function"),
            #    stringsAsFactors = FALSE
            #)

            #colnames(LoF_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

            ########################################################################
            #####################-----Statistical Tests-----########################
            ########################################################################
            ################-----In OL or Not-----#################
            HAP_f <- fisher.test(HAP_ct)
            HAP_tidy <- tidy(HAP_f)

            ################-----In OL or Not-----#################
            OL_f <- fisher.test(OL_ct)
            OL_tidy <- tidy(OL_f)

            ################-----In DD or Not-----#################
            DD_f <- fisher.test(DD_ct)
            DD_tidy <- tidy(DD_f)

            ################-----In Disease or Not-----#################
            OMIM_f <- fisher.test(OMIM_ct)
            OMIM_tidy <- tidy(OMIM_f)

            ################-----CRISPR Essential or Not-----#################
            Essential_f <- fisher.test(Essential_ct)
            Essential_tidy <- tidy(Essential_f)
           
            Non_Essential_f <- fisher.test(Non_Essential_ct)
            Non_Essential_tidy <- tidy(Non_Essential_f)

            ################-----Inheritance Pattern-----#################
            Recessive_f <- fisher.test(Recessive_ct)
            Recessive_tidy <- tidy(Recessive_f)

            Dominant_f <- fisher.test(Dominant_ct)
            Dominant_tidy <- tidy(Dominant_f)

            ################-----Homozygous LoF-----#################
            #LoF_f <- fisher.test(LoF_ct)
            #LoF_tidy <- tidy(LoF_f)

            ########################################################################
            ##########################-----Plots-----###############################
            ########################################################################
            #pvalue <- c(HAP_tidy$p.value, OL_tidy$p.value, DD_tidy$p.value, OMIM_tidy$p.value, Essential_tidy$p.value, Non_Essential_tidy$p.value, Recessive_tidy$p.value, Dominant_tidy$p.value, LoF_tidy$p.value)
            #mean <- c((HAP_tidy$estimate), (OL_tidy$estimate), (DD_tidy$estimate), (OMIM_tidy$estimate), (Essential_tidy$estimate), (Non_Essential_tidy$estimate), (Recessive_tidy$estimate), (Dominant_tidy$estimate), (LoF_tidy$estimate))
            #conf.low <- c((HAP_tidy$conf.low), (OL_tidy$conf.low), (DD_tidy$conf.low), (OMIM_tidy$conf.low), (Essential_tidy$conf.low), (Non_Essential_tidy$estimate), (Recessive_tidy$conf.low), (Dominant_tidy$conf.low), (LoF_tidy$conf.low))
            #conf.high <- c((HAP_tidy$conf.high), (OL_tidy$conf.high), (DD_tidy$conf.high), (OMIM_tidy$conf.high), (Essential_tidy$conf.high), (Non_Essential_tidy$estimate), (Recessive_tidy$conf.high), (Dominant_tidy$conf.high), (LoF_tidy$conf.high))

            #LoF_Tag <- paste("Homozygous Loss of", "\n", "Function Tolerant")
            #CRISPR_Essential_Tag <- paste("CRISPR", "\n", "Essential")
            #CRISPR_Non_Essential_Tag <- paste("CRISPR", "\n", "Non Essential")
            #OL <- paste("Olfactory", "\n", "receptor")
            #DD <- paste("Developmental", "\n", "Delay")
            #values <- c("Haploinsufficient", OL, DD, "OMIM", CRISPR_Essential_Tag, CRISPR_Non_Essential_Tag, "Recessive", "Dominant", LoF_Tag)

            pvalue <- c(OMIM_tidy$p.value, Recessive_tidy$p.value, Dominant_tidy$p.value, HAP_tidy$p.value, DD_tidy$p.value, Non_Essential_tidy$p.value, OL_tidy$p.value)
            mean <- c((OMIM_tidy$estimate), (Recessive_tidy$estimate), (Dominant_tidy$estimate), (HAP_tidy$estimate), (DD_tidy$estimate), (Non_Essential_tidy$estimate), (OL_tidy$estimate))
            conf.low <- c((OMIM_tidy$conf.low), (Recessive_tidy$conf.low), (Dominant_tidy$conf.low), (HAP_tidy$conf.low), (DD_tidy$estimate), (Non_Essential_tidy$conf.low), (OL_tidy$conf.low))
            conf.high <- c((OMIM_tidy$conf.high), (Recessive_tidy$conf.high), (Dominant_tidy$conf.high), (HAP_tidy$conf.high), (DD_tidy$estimate), (Non_Essential_tidy$conf.high), (OL_tidy$conf.high))

            #LoF_Tag <- paste("Homozygous Loss of", "\n", "Function Tolerant")
            CRISPR_Non_Essential_Tag <- paste("CRISPR", "\n", "Non Essential Gene")
            OL <- paste("Olfactory", "\n", "Receptor Gene")
            DD <- paste("Developmental", "\n", "Delay Gene")
            Rec <- paste("Recessive", "\n", "Disease Gene")
            Dom <- paste("Dominant", "\n", "Disease Gene")
            OMIM <- paste("OMIM", "\n", "Disease Gene")
            Hap <- paste("Haploinsufficient", "\n", "Disease Gene")
            values <- c(OMIM, Rec, Dom, Hap, DD, CRISPR_Non_Essential_Tag, OL)

            forest_plot_data <- data.frame(values, pvalue, mean, conf.low, conf.high)
            forest_plot_data$padjust <- p.adjust(forest_plot_data$pvalue, method="fdr", n=length(forest_plot_data$pvalue))

            forest_plot_data
}

outlier_genes_05 <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% select(hgncSymbol) %>% filter(!is.na(hgncSymbol)) %>% unique
outlier_genes_05_fp <- fisher_function(outlier_genes_05, filtered_genes)

outlier_genes_05_fp$p_adjust <- p.adjust(outlier_genes_05_fp$pvalue, method="fdr", n=length(outlier_genes_05_fp$pvalue))

write_csv(outlier_genes_05_fp, output_file2)
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

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")

autosomal_dominant <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/autosomal_dominant.tsv", col_names = FALSE)
autosomal_recessive <- read_tsv("Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/autosomal_recessive.tsv", col_names = FALSE)
CRISPR_essential <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/CRISPR_essential_genes.tsv", col_names = FALSE)
CRISPR_nonessential <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/CRISPR_nonessential_genes.tsv", col_names = FALSE)
OMIM <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/genemap2.tsv")
homozygous_LoF <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/homozygous_Lof_genes.tsv", col_names = FALSE)
Developmental_delay <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/developmental_delay_genes.csv")
olfactory_receptors <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/olfactory_receptors.tsv", col_names=FALSE)
haplo <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_1/Gene_Information/haploinsufficient.tsv", col_names=FALSE)

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

display.brewer.all
pal <- brewer.pal(10, "Paired")

genes <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_filtered.rds")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

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

########################################################################
##################-----Create Outlier Dataframe-----####################
########################################################################

################-----Get Number by Gene per Person-----#################
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
            CRISPR_essential_df <- data.frame(Gene = CRISPR_essential, Role = "Essential")
            CRISPR_nonessential_df <- data.frame(Gene = CRISPR_nonessential, Role = "Nonessential")

            CRISPR <- bind_rows(CRISPR_essential_df, CRISPR_nonessential_df)

            outlier_in_CRISPR <- filter(outlier_genes, hgncSymbol %in% CRISPR$X1) %>% select(-OMIM_disease_gene)

            essential <- intersect(CRISPR_essential_df$X1, outlier_in_CRISPR$hgncSymbol)
            non_essential <- intersect(CRISPR_nonessential_df$X1, outlier_in_CRISPR$hgncSymbol)

            outlier_in_CRISPR$CRISPR_Role <- ifelse(outlier_in_CRISPR$hgncSymbol %in% essential, TRUE, FALSE)
            CRISPR_Role <- outlier_in_CRISPR %>% select(hgncSymbol, CRISPR_Role)

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
            homo_LoF <- intersect(outlier_genes$hgncSymbol, homozygous_LoF$X1)
            not_LoF <- setdiff(outlier_genes$hgncSymbol, homozygous_LoF$X1)

            outlier_genes$Homozygous_LoF <- ifelse(outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            Homozygous_LoF <- outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

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
            nonoutlier_in_CRISPR <- filter(non_outlier_genes, hgncSymbol %in% CRISPR$X1) %>% select(-OMIM_disease_gene)

            essential <- intersect(CRISPR_essential_df$X1, nonoutlier_in_CRISPR$hgncSymbol)
            non_essential <- intersect(CRISPR_nonessential_df$X1, nonoutlier_in_CRISPR$hgncSymbol)

            nonoutlier_in_CRISPR$CRISPR_Role <- ifelse(nonoutlier_in_CRISPR$hgncSymbol %in% essential, TRUE, FALSE)
            CRISPR_Role_nonoutliers <- nonoutlier_in_CRISPR %>% select(hgncSymbol, CRISPR_Role)

            ################-----Inheritance Pattern-----#################
            nonoutlier_in_pattern <- filter(non_outlier_genes, hgncSymbol %in% Pattern$X1) %>% select(-OMIM_disease_gene)

            recessive <- intersect(Recessive$X1, nonoutlier_in_pattern$hgncSymbol)
            dominant <- intersect(Dominant$X1, nonoutlier_in_pattern$hgncSymbol)

            nonoutlier_in_pattern$Recessive_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% recessive, TRUE, FALSE)
            Recessive_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Recessive_Pattern)

            nonoutlier_in_pattern$Dominant_Pattern <- ifelse(nonoutlier_in_pattern$hgncSymbol %in% dominant, TRUE, FALSE)
            Dominant_Pattern_nonoutliers <- nonoutlier_in_pattern %>% select(hgncSymbol, Dominant_Pattern)

            ################-----Homozygous LoF-----#################
            homo_LoF <- intersect(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)
            not_LoF <- setdiff(non_outlier_genes$hgncSymbol, homozygous_LoF$X1)

            non_outlier_genes$Homozygous_LoF <- ifelse(non_outlier_genes$hgncSymbol %in% homo_LoF, TRUE, FALSE)
            Homozygous_LoF_nonoutliers <- non_outlier_genes %>% select(hgncSymbol, Homozygous_LoF)

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
            CRISPR_Role_nonoutliers$Splicing_Outlier <- FALSE
            CRISPR_Role$Splicing_Outlier <- TRUE

            CRISPR_Essential <- bind_rows(CRISPR_Role, CRISPR_Role_nonoutliers)
            colnames(CRISPR_Essential) <- c("hgncSymbol", "CRISPR_Essentials", "Splicing_Outlier")

            CRISPR_Splicing <- filter(CRISPR_Essential, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essentials == TRUE) %>% nrow
            CRISPR_Non_Splicing <- filter(CRISPR_Essential, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essentials == TRUE) %>% nrow
            Non_CRISPR_Splicing <- filter(CRISPR_Essential, Splicing_Outlier == TRUE) %>% filter(CRISPR_Essentials != TRUE) %>% nrow
            Non_CRISPR_Non_Splicing <- filter(CRISPR_Essential, Splicing_Outlier != TRUE) %>% filter(CRISPR_Essentials != TRUE) %>% nrow

            CRISPR_ct <- data.frame(
                "Splicing_Outlier" = c(CRISPR_Splicing, Non_CRISPR_Splicing),
                "Spilicing_Inlier" = c(CRISPR_Non_Splicing, Non_CRISPR_Non_Splicing),
                row.names = c("CRISPR_Essential_Gene", "CRISPR_Non_Essential_Gene"),
                stringsAsFactors = FALSE
            )

            colnames(CRISPR_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

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
            Homozygous_LoF_nonoutliers$Splicing_Pattern <- FALSE
            Homozygous_LoF$Splicing_Pattern <- TRUE

            Homozygous_Loss_of_Function <- bind_rows(Homozygous_LoF_nonoutliers, Homozygous_LoF)

            LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF == TRUE) %>% nrow
            Non_LoF_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern == TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow
            Non_LoF_Non_Splicing <- filter(Homozygous_Loss_of_Function, Splicing_Pattern != TRUE) %>% filter(Homozygous_LoF != TRUE) %>% nrow

            LoF_ct <- data.frame(
                "Splicing_Outlier" = c(LoF_Splicing, Non_LoF_Splicing),
                "Spilicing_Inlier" = c(LoF_Non_Splicing, Non_LoF_Non_Splicing),
                row.names = c("Homozygous_Loss_of_Function", "Non_Homozygous_Loss_of_Function"),
                stringsAsFactors = FALSE
            )

            colnames(LoF_ct) <- c("Splicing_Outlier", "Spilicing_Inlier")

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
            CRISPR_f <- fisher.test(CRISPR_ct)
            CRISPR_tidy <- tidy(CRISPR_f)

            ################-----Inheritance Pattern-----#################
            Recessive_f <- fisher.test(Recessive_ct)
            Recessive_tidy <- tidy(Recessive_f)

            Dominant_f <- fisher.test(Dominant_ct)
            Dominant_tidy <- tidy(Dominant_f)

            ################-----Homozygous LoF-----#################
            LoF_f <- fisher.test(LoF_ct)
            LoF_tidy <- tidy(LoF_f)

            ########################################################################
            ##########################-----Plots-----###############################
            ########################################################################
            pvalue <- c(HAP_tidy$p.value, OL_tidy$p.value, DD_tidy$p.value, OMIM_tidy$p.value, CRISPR_tidy$p.value, Recessive_tidy$p.value, Dominant_tidy$p.value, LoF_tidy$p.value)
            mean <- c((HAP_tidy$estimate), (OL_tidy$estimate), (DD_tidy$estimate), (OMIM_tidy$estimate), (CRISPR_tidy$estimate), (Recessive_tidy$estimate), (Dominant_tidy$estimate), (LoF_tidy$estimate))
            conf.low <- c((HAP_tidy$conf.low), (OL_tidy$conf.low), (DD_tidy$conf.low), (OMIM_tidy$conf.low), (CRISPR_tidy$conf.low), (Recessive_tidy$conf.low), (Dominant_tidy$conf.low), (LoF_tidy$conf.low))
            conf.high <- c((HAP_tidy$conf.high), (OL_tidy$conf.high), (DD_tidy$conf.high), (OMIM_tidy$conf.high), (CRISPR_tidy$conf.high), (Recessive_tidy$conf.high), (Dominant_tidy$conf.high), (LoF_tidy$conf.high))

            LoF_Tag <- paste("Homozygous", "\n", "Loss of Function")
            CRISPR_Tag <- paste("CRISPR", "\n", "Essential")
            values <- c("HAP", "OL", "DD", "OMIM", CRISPR_Tag, "Recessive", "Dominant", LoF_Tag)

            forest_plot_data <- data.frame(values, pvalue, mean, conf.low, conf.high)
            forest_plot_data
}


outlier_genes_1e_6 <- all_uncompiled %>% filter(padjust <= 1e-6) %>% filter(abs(zScore) >= 2) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(! sampleID %in% low_RIN$sampleID) %>% select(hgncSymbol) %>% filter(!is.na(hgncSymbol)) %>% unique
outlier_genes_1e_6_fp <- fisher_function(outlier_genes_1e_6, filtered_genes)

outlier_genes_001 <- all_uncompiled %>% filter(padjust <= 0.001) %>% filter(abs(zScore) >= 2) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(! sampleID %in% low_RIN$sampleID) %>% select(hgncSymbol) %>% filter(!is.na(hgncSymbol)) %>% unique
outlier_genes_001_fp <- fisher_function(outlier_genes_001, filtered_genes)

outlier_genes_05 <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(zScore) >= 2) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(! sampleID %in% low_RIN$sampleID) %>% select(hgncSymbol) %>% filter(!is.na(hgncSymbol)) %>% unique
outlier_genes_05_fp <- fisher_function(outlier_genes_05, filtered_genes)

outlier_genes_05_fp$p_adjust <- p.adjust(outlier_genes_05_fp$pvalue, method="fdr", n=length(outlier_genes_05_fp$pvalue))


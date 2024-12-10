library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

twenty5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_25_filtered_stanford/FRASER_output.csv")
twenty5_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_25_filtered_stanford/FRASER_input.csv")
fifty <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_50_filtered_stanford/FRASER_output.csv")
fifty_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_50_filtered_stanford/FRASER_input.csv")
one00  <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_100_filtered_stanford/FRASER_output.csv")
one00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_100_filtered_stanford/FRASER_input.csv")
one50 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_150_filtered_stanford/FRASER_output.csv")
one50_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_150_filtered_stanford/FRASER_input.csv")
two00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_200_filtered_stanford/FRASER_output.csv")
two00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_200_filtered_stanford/FRASER_input.csv")
three00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_300_filtered_stanford/FRASER_output.csv")
three00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_300_filtered_stanford/FRASER_input.csv")
four00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_400_filtered_stanford/FRASER_output.csv")
four00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_400_filtered_stanford/FRASER_input.csv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

########################################################################
######################-----Make Dataframe-----##########################
########################################################################

######################-----25-----##########################
get_number_MIGs_theta <- function(dataframe, files, col_name){
        MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
        colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
        MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG")

        #get theta dataframe
        results_filtered_theta <- dataframe %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta")
        results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol") %>% unique
        colnames(results_filtered_theta) <- c("sampleID", "gene_symbol")

        #combine dataframes
        sig_genes_type_all_samples <- left_join(MIG_table_select, results_filtered_theta) %>% filter(!is.na(sampleID))

        joined <- sig_genes_type_all_samples %>% select(sampleID, gene_symbol) %>% 
            group_by(sampleID) %>% tally() %>% arrange(desc(n))

        missing_all <- setdiff(files$sampleID, joined$sampleID)
        missing_all_df <- data.frame(sampleID = missing_all, n = 0)
        MIG_count <- bind_rows(joined, missing_all_df)

        MIG_colname <- paste("no_MIGs_theta_juncs", col_name, sep="_")
        colnames(MIG_count) <- c("sampleID", MIG_colname)
        MIG_count

}

twenty5_MIGs <- get_number_MIGs_theta(twenty5, twenty5_files, "25") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_25)
fifty_MIGs <- get_number_MIGs_theta(fifty, fifty_files, "50") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_50)
one00_MIGs <- get_number_MIGs_theta(one00, one00_files, "100") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_100)
one50_MIGs <- get_number_MIGs_theta(one50, one50_files, "150") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_150)
two00_MIGs <- get_number_MIGs_theta(two00, two00_files, "200") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_200)
three00_MIGs <- get_number_MIGs_theta(three00, three00_files, "300") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_300)
four00_MIGs <- get_number_MIGs_theta(four00, four00_files, "400") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_400)

number_MIGs_with_IR <- data.frame(sample_size = c(25, 50, 100, 150, 200, 300, 400), 
                            number_MIGs_theta=c(twenty5_MIGs, fifty_MIGs, one00_MIGs, one50_MIGs, two00_MIGs, three00_MIGs, four00_MIGs))


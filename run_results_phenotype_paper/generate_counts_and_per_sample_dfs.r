library(readxl)
library(magrittr)
library(biomaRt)
library(tidyverse)
require(data.table)

############################################################################################################
###################################----Read in dataframes----###############################################
############################################################################################################
#Fibroblast data - input and output
    #you need the input because FRASER output doesn't include samples with 0 outliers
fibro_theta_FRASER_output <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibroblasts_theta.csv")
fibro_FRASER_input <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibroblast_input.xlsx")

#Blood data - input and output
    #you need the input because FRASER output doesn't include samples with 0 outliers
blood_FRASER_output <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/output/output_drafts/FRASER_no_missing_rm_seqbatch_8/FRASER_output.csv")
blood_FRASER_input <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/output/output_drafts/FRASER_no_missing_rm_seqbatch_8/FRASER_input.csv")

#MIDB file used for filtering on Dec 2nd, 2025
MIDB_all <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/Homo_sapiens_intron.csv") 

############################################################################################################
########################################----Functions----###################################################
############################################################################################################
make_outlier_df_gene <- function(filtered_tables, filepath) {   
    if (!dir.exists(filepath)) {
        dir.create(filepath, recursive = TRUE)
    }

    filenumber <- length(list.files(filepath))

    ids_w_outliers <- filtered_tables$sampleID %>% unique
    ids_w_outliers_count <- filtered_tables %>% pull(sampleID) %>% unique %>% length()

    if (filenumber < ids_w_outliers_count) {
        for (i in ids_w_outliers) {
            ind_gene_df_sample <- filtered_tables %>% filter(sampleID == i)
            
            ind_gene_df <- ind_gene_df_sample %>% select(gene_symbol) %>% unique()
            
            pro_ensg <- mean(startsWith(ind_gene_df$gene_symbol, "ENSG"), na.rm = TRUE)

            if (pro_ensg > 0.5) {
                ind_gene_df <- convert_to_gene_names(ind_gene_df)
            }
            
            filepath_sample <- paste0(filepath, "/", i, ".csv")
            print(filepath_sample)
            write_csv(ind_gene_df, filepath_sample)
        }
    }
} 

convert_to_gene_names <- function(ensg_df) {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    ensg_ids <- unique(gsub("\\.[0-9]+$", "", ensg_df$gene_symbol))

    # Query BioMart and map ENSG ids to hgnc/ gene symbols
    mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ensg_ids,
    mart       = ensembl
    )

    # Remove version suffix from the returned ENSG ids so they match my input
    mapping$ensembl_gene_id <- gsub("\\.[0-9]+$", "", mapping$ensembl_gene_id)

    colnames(ensg_df) <- c("ensembl_gene_id")

    # Merge back to original data frame
    filtered_mapped <- left_join(ensg_df, mapping)
    filtered_mapped
}

make_outlier_df_junction <- function(filtered_tables, filepath) {  
    if (!dir.exists(filepath)) {
        dir.create(filepath, recursive = TRUE)
    }

    filenumber <- length(list.files(filepath))

    ids_w_outliers <- filtered_tables$sampleID %>% unique
    ids_w_outliers_count <- filtered_tables %>% pull(sampleID) %>% unique() %>% length()

    if (filenumber < ids_w_outliers_count) {
        for (i in ids_w_outliers) {
            ind_junction_df_sample <- filtered_tables %>% filter(sampleID == i) %>% unique()
            
            filepath_sample <- paste0(filepath, "/", i, ".csv")
            write_csv(ind_junction_df_sample, filepath_sample)
        }
    }
}

get_genes_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(dir_path, level_name, "/unique_gene_counts.csv")
    if (!file.exists(count_output_fp)) {
        genes <- filtered_table %>% select(sampleID, gene_symbol) %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
        zero_count_samples <- input_sid %>% filter(! sampleID %in% filtered_table$sampleID)
        zero_count_samples$n <- 0
        count_joined <- bind_rows(genes, zero_count_samples)
        count_joined_sdv <- count_joined %>% mutate(Z_score = (n - mean(n)) / sd(n)) %>%
            arrange(sampleID)
        write_csv(count_joined_sdv, count_output_fp)
        count_joined_sdv
    }
}

get_junctions_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(dir_path, level_name, "/unique_junction_counts.csv")
    if (!file.exists(count_output_fp)) {
        junctions <- filtered_table %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
        zero_count_samples <- input_sid %>% filter(! sampleID %in% filtered_table$sampleID)
        zero_count_samples$n <- 0
        count_joined <- bind_rows(junctions, zero_count_samples)
        count_joined_sdv <- count_joined %>% mutate(Z_score = (n - mean(n)) / sd(n)) %>%
            arrange(sampleID)
        write_csv(count_joined_sdv, count_output_fp)
        count_joined_sdv
    } 
}

main_function <- function(filtered_FRASER_results, FRASER_input_files, minor_intron_file, dir_path, FRASER_output_fp, FRASER_input_fp, MIG_fp) {
    ###############----Define all levels----##############
    gene_minor_juncfiltered <- "MIGs_w_theta_outliers_in_minor_introns"
    gene_minor_NONfiltered <- "MIGs_w_theta_outliers"
    gene_nonmajor_juncfiltered <- "non_major_genes_w_theta_outliers_in_non_major_introns"
    gene_nonmajor_NONfiltered <- "non_major_genes_w_theta_outliers"

    junc_minor_juncfiltered <- "theta_outliers_in_minor_introns"
    junc_minor_NONfiltered <- "theta_outliers_in_MIGs"
    junc_nonmajor_juncfiltered <- "theta_outliers_in_non_major_introns"
    junc_nonmajor_NONfiltered <- "theta_outliers_in_non_major_genes"

    level <- list(gene_minor_juncfiltered, gene_minor_NONfiltered, 
                        gene_nonmajor_juncfiltered, gene_nonmajor_NONfiltered,
                        junc_minor_juncfiltered, junc_minor_NONfiltered,
                        junc_nonmajor_juncfiltered, junc_nonmajor_NONfiltered)
    
    ###############----Make output directories----##############
    for (i in level) {
            dir_path_level <- paste0(dir_path, i, "/")
            if (!dir.exists(dir_path_level)) {
                dir.create(dir_path_level, recursive = TRUE)
            }
    }

    ###############----Join and filter dataframes----##############
    pro_ensg <- mean(startsWith(filtered_FRASER_results$hgncSymbol, "ENSG"), na.rm = TRUE)

    if (pro_ensg > 0.5) {
        names(minor_intron_file)[names(minor_intron_file) == "ensembl_gene_id"] <- "hgncSymbol"
        
        ####Filter FRASER results to MIGs
        MIG_list <- minor_intron_file %>% filter(intron_class == "minor") %>% pull(hgncSymbol) %>% unique
        filtered_FRASER_results_MIGs <- filtered_FRASER_results %>% filter(hgncSymbol %in% MIG_list) #this is an input for the counts function, ect.

        ####Filter FRASER results to non-major
        non_major_list <- minor_intron_file %>% pull(hgncSymbol) %>% unique
        filtered_FRASER_results_non_major <- filtered_FRASER_results %>% filter(hgncSymbol %in% non_major_list) #this is an input for the counts function, ect.
     
        ####Join minor
        names(minor_intron_file)[names(minor_intron_file) == "gene_symbol"] <- "gene_name"
        MIGs_joined <- left_join(minor_intron_file, filtered_FRASER_results_MIGs, by="hgncSymbol")
        names(MIGs_joined)[names(MIGs_joined) == "hgncSymbol"] <- "gene_symbol"
        names(filtered_FRASER_results_MIGs)[names(filtered_FRASER_results_MIGs) == "hgncSymbol"] <- "gene_symbol"

        ####Join non-major
        names(minor_intron_file)[names(minor_intron_file) == "gene_symbol"] <- "gene_name"
        non_major_joined <- left_join(minor_intron_file, filtered_FRASER_results_non_major, by="hgncSymbol")
        names(non_major_joined)[names(non_major_joined) == "hgncSymbol"] <- "gene_symbol"
        names(filtered_FRASER_results_non_major)[names(filtered_FRASER_results_non_major) == "hgncSymbol"] <- "gene_symbol"
    } else {
        ####Filter FRASER results to MIGs
        MIG_list <- minor_intron_file %>% filter(intron_class == "minor") %>% pull(gene_symbol) %>% unique
        filtered_FRASER_results_MIGs <- filtered_FRASER_results %>% filter(hgncSymbol %in% MIG_list) #this is an input for the counts function, ect.

        ####Filter FRASER results to non-major
        non_major_list <- minor_intron_file %>% pull(gene_symbol) %>% unique
        filtered_FRASER_results_non_major <- filtered_FRASER_results %>% filter(hgncSymbol %in% non_major_list) #this is an input for the counts function, ect.

        ####Join minor
        names(filtered_FRASER_results_MIGs)[names(filtered_FRASER_results_MIGs) == "hgncSymbol"] <- "gene_symbol"
        MIGs_joined <- left_join(minor_intron_file, filtered_FRASER_results_MIGs, by="gene_symbol")
        names(filtered_FRASER_results_MIGs)[names(filtered_FRASER_results_MIGs) == "hgncSymbol"] <- "gene_symbol"

        ####Join non-major
        names(filtered_FRASER_results_non_major)[names(filtered_FRASER_results_non_major) == "hgncSymbol"] <- "gene_symbol"
        non_major_joined <- left_join(minor_intron_file, filtered_FRASER_results_non_major, by="gene_symbol")
        names(filtered_FRASER_results_non_major)[names(filtered_FRASER_results_non_major) == "hgncSymbol"] <- "gene_symbol"        
    }

    ####Filter FRASER results to minor
    MIGs_joined_filter_on_start_end <- MIGs_joined %>% filter(intron_class == "minor") %>% 
                                       filter(start == intron_start | start == intron_start-1 | start == intron_start+1 |
                                       end == intron_end|  end == intron_end-1 | end == intron_end+1) %>% 
                                       select(sampleID, gene_symbol, intron_name, intron_start, intron_end, start, end, intron_class, pValue, padjust, 
                                       zScore, psiValue, deltaPsi, meanCounts, meanTotalCounts, counts, totalCounts) %>% unique #this is an input for the counts function, ect.

    ####Filter FRASER results to non-major
    non_major_joined_filter_on_start_end <- non_major_joined %>% 
                                       filter(start == intron_start | start == intron_start-1 | start == intron_start+1 |
                                       end == intron_end|  end == intron_end-1 | end == intron_end+1) %>% 
                                       select(sampleID, gene_symbol, intron_name, intron_start, intron_end, start, end, intron_class, pValue, padjust, 
                                       zScore, psiValue, deltaPsi, meanCounts, meanTotalCounts, counts, totalCounts) %>% unique #this is an input for the counts function, ect.

    #########################----get gene/ junction/ transcript counts----#######################
    get_genes_count_df(MIGs_joined_filter_on_start_end, FRASER_input_files, level[1], dir_path)
    get_genes_count_df(filtered_FRASER_results_MIGs, FRASER_input_files, level[2], dir_path)
    get_genes_count_df(non_major_joined_filter_on_start_end, FRASER_input_files, level[3], dir_path)
    get_genes_count_df(filtered_FRASER_results_non_major, FRASER_input_files, level[4], dir_path)

    get_junctions_count_df(MIGs_joined_filter_on_start_end, FRASER_input_files, level[5], dir_path)
    get_junctions_count_df(filtered_FRASER_results_MIGs, FRASER_input_files, level[6], dir_path)
    get_junctions_count_df(non_major_joined_filter_on_start_end, FRASER_input_files, level[7], dir_path)
    get_junctions_count_df(filtered_FRASER_results_non_major, FRASER_input_files, level[8], dir_path)

    #########################----get gene counts----#######################
    make_outlier_df_gene(MIGs_joined_filter_on_start_end, paste0(dir_path, level[1], "/unique_genes_per_sample"))
    make_outlier_df_gene(filtered_FRASER_results_MIGs, paste0(dir_path, level[2], "/unique_genes_per_sample"))
    make_outlier_df_gene(non_major_joined_filter_on_start_end, paste0(dir_path, level[3], "/unique_genes_per_sample"))
    make_outlier_df_gene(filtered_FRASER_results_non_major, paste0(dir_path, level[4], "/unique_genes_per_sample"))

    make_outlier_df_junction(MIGs_joined_filter_on_start_end, paste0(dir_path, level[5], "/unique_junctions_per_sample"))
    make_outlier_df_junction(filtered_FRASER_results_MIGs, paste0(dir_path, level[6], "/unique_junctions_per_sample"))
    make_outlier_df_junction(non_major_joined_filter_on_start_end, paste0(dir_path, level[7], "/unique_junctions_per_sample"))
    make_outlier_df_junction(filtered_FRASER_results_non_major, paste0(dir_path, level[8], "/unique_junctions_per_sample"))

    readme_simple <- c(
        "# Run completed on:",
            paste("**", format(Sys.Date()), "**"),
        "## Filepath MIG file used:",
            paste("`", MIG_fp, "`"),
        "## Filepath of FRASER output used",
            paste("`", FRASER_output_fp, "`"),
        "## Filepath of FRASER input used",
            paste("`", FRASER_input_fp, "`")
    )
    
    writeLines(readme_simple, paste0(dir_path, "README.md"))
}

############################################################################################################
#######################################----Get results----##################################################
############################################################################################################
fibro_filtered_05_abs3 <- fibro_theta_FRASER_output %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)

blood_filtered_05_abs3 <- blood_FRASER_output %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")

#######################################----Fibroblasts----##################################################
dirpath_fibro <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibro_filtered_05_abs3_Jan_28/"
main_function(fibro_filtered_05_abs3, fibro_FRASER_input, MIDB_all, dirpath_fibro, fibro_theta_FRASER_output_fp, fibro_FRASER_input_fp, MIDB_all_fp)
test %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))

#######################################----Fibroblasts----##################################################
dirpath_blood <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/blood_filtered_05_abs3_try2/"
blood_FRASER_input <- blood_FRASER_input %>% select(sampleID) %>% unique
main_function(blood_filtered_05_abs3, blood_FRASER_input, MIDB_all, dirpath_blood, blood_FRASER_output_fp, blood_FRASER_input_fp, MIDB_all_fp)

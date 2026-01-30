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
fibro_theta_FRASER_output_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibroblasts_theta.csv"
fibro_theta_FRASER_output <- read_csv(fibro_theta_FRASER_output_fp)

fibro_FRASER_input_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibroblast_input.xlsx"
fibro_FRASER_input <- read_excel(fibro_FRASER_input_fp)

#Blood data - input and output
    #you need the input because FRASER output doesn't include samples with 0 outliers
blood_FRASER_output_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/output/output_drafts/FRASER_no_missing_rm_seqbatch_8/FRASER_output.csv"
blood_FRASER_output <- read_csv(blood_FRASER_output_fp)

blood_FRASER_input_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/output/output_drafts/FRASER_no_missing_rm_seqbatch_8/FRASER_input.csv"
blood_FRASER_input <- read_csv(blood_FRASER_input_fp)

#MIDB file used for filtering on Dec 2nd, 2025
MIDB_all_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/Homo_sapiens_intron.csv" 
MIDB_all <- read_csv(MIDB_all_fp) 

############################################################################################################
##########################################----Functions----#################################################
############################################################################################################

#####################################----convert_to_gene_names----##########################################
RATIONAL:
    #the fibroblast FRASER output from the Broad included ENSGs as its output in the hgncSymbol column INSTEAD of hgnc symbols/ gene names
        #so, this function converts ENSGs to hgnc symbols/ gene names
        #it is called in make_outlier_df_gene after checking if the rows in the hgncSymbol column contain the string "ENSG" > 50% of the time
#INPUT: 
    #ensg_df: the filtered FRASER table given to make_outlier_df_gene from change_colnames_bind_rows from the main_function
        #this table is sometimes altered in change_colnames_bind_rows before being passed to make_outlier_df_gene to only include introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
#OUTPUT: the same dataframe given in the INPUT with the values in the "hgncSymbol" column changed from ENSG IDs to hgnc symbols/ gene names

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

####################################----make_outlier_df_gene----############################################
#INPUT: 
    #filtered_table: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #input_sid: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #dir_path: the directory name where the output will be saved under the "unique_genes_per_sample" folder
#OUTPUT:
    #a folder (unique_genes_per_sample) containing one file per sample named [sampleID].csv with two columns
        #ensembl_gene_id: ENSG of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)
        #hgnc_symbol: gene names of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)

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

######################################----get_genes_count_df----############################################
#INPUT: 
    #filtered_tables: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #filepath: the filepath where you would like to save the csv output
#OUTPUT:
    #a csv with three columns:
        #sampleID: the sampleID as found in the FRASER output
        #n: the number of genes containing outliers of the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)
        #Z_score: the Z-score for that sample taken from the mean and sd of the n column

get_genes_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(level_name, "unique_gene_counts.csv")
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

####################################----get_junctions_count_df----###########################################
#INPUT: 
    #filtered_tables: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #filepath: the filepath where you would like to save the csv output
#OUTPUT:
    #a csv with three columns:
        #sampleID: the sampleID as found in the FRASER output
        #n: the number of junctions with outliers of the specified intron type (as found in the folder name) which are EITHER
            #outliers of any intron type (ex: number of intron retention events of any type within minor intron containing genes, or MIGs)
            #outliers of only the specified intron type (ex: number of junctions with minor intron retention)
        #Z_score: the Z-score for that sample taken from the mean and sd of the n column

get_junctions_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(level_name, "unique_junction_counts.csv")
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

#################################----change_colnames_bind_rows----###########################################
#INPUT: 
    #intron: the intron type you want to filter the MIDB file by 
        #this will do through every intron type in the MIDB database in for loops in main_function
    #MIDB_file: the file from the MIDB database read in earlier via read_csv
    #FRASER_results_filtered: the results of FRASER given to the main function
        #this should ALREADY be filtered for the padjust, deltaPsi, and type of your choice
    #FRASER_input: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #dir_path: 
        #the path you want the outputs of get_genes_count_df and, if focusing on a per-gene level, make_outlier_df_gene
    #index: the number in the list given to the function by main_funcion
        #here, the odd numbers indicate that we are looking for either only:
            #genes with outliers of ONLY the specified intron type (as found in the folder name) OR
            #junctions with ONLY outliers of the specified intron type (as found in the folder name)
                #as opposed to, for example, junctions with outliers of ANY intron type in a gene containing at least one intron of 
                    #the specified intron type (as found in the folder name)
#OUTPUT:
    #if looking on a per-gene basis:
        #the output of get_genes_count_df AND make_outlier_df_gene
    #if looking on a per-junction basis:
        #the output of get_junctions_count_df

change_colnames_bind_rows <- function(intron, MIDB_file, FRASER_results_filtered, FRASER_input, level, index, dir_path) {
    pro_ensg <- mean(startsWith(FRASER_results_filtered$hgncSymbol, "ENSG"), na.rm = TRUE)

    if (pro_ensg > 0.5) {
        intron_list <- MIDB_file %>% filter(intron_class == intron) %>% pull(ensembl_gene_id) %>% unique

        names(MIDB_file)[names(MIDB_file) == "ensembl_gene_id"] <- "hgncSymbol"
        
        ####Filter FRASER results to MIGs
        FRASER_results_filtered <- FRASER_results_filtered %>% filter(hgncSymbol %in% intron_list) #this is an input for the counts function, ect.

        ####Join to intron file
        names(MIDB_file)[names(MIDB_file) == "gene_symbol"] <- "gene_name"
        MIDB_joined <- left_join(MIDB_file, FRASER_results_filtered, by="hgncSymbol")
        names(MIDB_joined)[names(MIDB_joined) == "hgncSymbol"] <- "gene_symbol"
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"

    } else {
        intron_list <- MIDB_file %>% filter(intron_class == intron) %>% pull(gene_symbol) %>% unique

        ####Filter FRASER results to MIGs
        FRASER_results_filtered <- FRASER_results_filtered %>% filter(hgncSymbol %in% intron_list) #this is an input for the counts function, ect.

        ####Join to intron file
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"
        MIDB_joined <- left_join(MIDB_file, FRASER_results_filtered, by="gene_symbol")
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"
    }

    if (index %% 2 == 1) {
    FRASER_results_filtered <- MIDB_joined %>% filter(intron_class == "minor") %>% 
                                       filter(start == intron_start | start == intron_start-1 | start == intron_start+1 |
                                       end == intron_end|  end == intron_end-1 | end == intron_end+1) %>% 
                                       select(sampleID, gene_symbol, intron_name, intron_start, intron_end, start, end, intron_class, pValue, padjust, 
                                       zScore, psiValue, deltaPsi, meanCounts, meanTotalCounts, counts, totalCounts) %>% unique #this is an input for the counts function, ect.

    }

    if (grepl("genes_w", level)) {
        get_genes_count_df(FRASER_results_filtered, FRASER_input, level, dir_path)
        make_outlier_df_gene(FRASER_results_filtered, paste0(level, "unique_genes_per_sample"))
    } else {
        get_junctions_count_df(FRASER_results_filtered, FRASER_input, level, dir_path)
    }
    
}

#################################----change_colnames_bind_rows----###########################################
#INPUT: 
    #filtered_FRASER_results: the results of FRASER given to the main function
        #this should ALREADY be filtered for the padjust, deltaPsi, and type of your choice
    #FRASER_input_files: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #minor_intron_file: the file from the MIDB database read in earlier via read_csv
    #dir_path: the directory where you want all results saved
    #FRASER_output_fp: the filepath where filtered_FRASER_results originates from
        #this will be recorded in a readme
    #FRASER_input_fp: the filepath where FRASER_input_files originates from
        #this will be recorded in a readme
    #MIG_fp: the filepath where minor_intron_file originates from
        #this will be recorded in a readme
#OUTPUT:
    #output of get_junctions_count_df for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical): 
        #csv with three columns:
            #sampleID: the sampleID as found in the FRASER output
            #n: the number of junctions with outliers of the specified intron type (as found in the folder name) which are EITHER
                #outliers of any intron type (ex: number of intron retention events of any type within minor intron containing genes, or MIGs)
                #outliers of only the specified intron type (ex: number of junctions with minor intron retention)
            #Z_score: the Z-score for that sample taken from the mean and sd of the n column
    #output of get_genes_count_df for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical):
        #a csv with three columns:
            #sampleID: the sampleID as found in the FRASER output
            #n: the number of genes containing outliers of the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
            #Z_score: the Z-score for that sample taken from the mean and sd of the n column
    #the output of make_outlier_df_gene for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical):
         #a folder (unique_genes_per_sample) containing one file per sample named [sampleID].csv with two columns
            #ensembl_gene_id: ENSG of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
            #hgnc_symbol: gene names of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
    #a readme containing:
            #the date this was run
            #filepath of MIDB file used
            #filepath of the FRASER output used
            #filepath of the FRASER input used

main_function <- function(filtered_FRASER_results, FRASER_input_files, minor_intron_file, dir_path, FRASER_output_fp, FRASER_input_fp, MIG_fp) {
    ###############----Define all levels----##############
    gene_minor_juncfiltered <- "genes_w_minor_intron_retention"
    gene_minor_NONfiltered <- "genes_w_minor_introns_w_intron_retention"
    junc_minor_juncfiltered <- "juncs_w_minor_intron_retention"
    junc_minor_NONfiltered <- "juncs_w_intron_retention_in_genes_w_minor_intron"

    gene_majorlike_juncfiltered <- "genes_w_major_like_intron_retention"
    gene_majorlike_NONfiltered <- "genes_w_major_like_introns_w_intron_retention"
    junc_majorlike_juncfiltered <- "juncs_w_major_like_intron_retention"
    junc_majorlike_NONfiltered <- "juncs_w_intron_retention_in_genes_w_major_like_intron"

    gene_minorlike_juncfiltered <- "genes_w_minor_like_intron_retention"
    gene_minorlike_NONfiltered <- "genes_w_minor_like_introns_w_intron_retention"
    junc_minorlike_juncfiltered <- "juncs_w_minor_like_intron_retention"
    junc_minorlike_NONfiltered <- "juncs_w_intron_retention_in_genes_w_minor_like_intron"

    gene_noncanonical_juncfiltered <- "genes_w_non_canonical_intron_retention"
    gene_noncanonical_NONfiltered <- "genes_w_non_canonical_introns_w_intron_retention"
    junc_noncanonical_juncfiltered <- "juncs_w_non_canonical_intron_retention"
    junc_noncanonical_NONfiltered <- "juncs_w_intron_retention_in_genes_w_non_canonical_genes"

    gene_majorhybrid_juncfiltered <- "genes_w_major_hybrid_intron_retention"
    gene_majorhybrid_NONfiltered <- "genes_w_major_hybrid_introns_w_intron_retention"
    junc_majorhybrid_juncfiltered <- "juncs_w_major_hybrid_intron_retention"
    junc_majorhybrid_NONfiltered <- "juncs_w_major_hybrid_in_genes_w_major_hybrid_genes"

    gene_minorhybrid_juncfiltered <- "genes_w_minor_hybrid_intron_retention"
    gene_minorhybrid_NONfiltered <- "genes_w_minor_hybrid_introns_w_intron_retention"
    junc_minorhybrid_juncfiltered <- "juncs_w_minor_hybrid_intron_retention"
    junc_minorhybrid_NONfiltered <- "juncs_w_minor_hybrid_in_genes_w_minor_hybrid_genes"

    level <- list(gene_minor_juncfiltered, gene_minor_NONfiltered, 
                        junc_minor_juncfiltered, junc_minor_NONfiltered, 
                        gene_majorlike_juncfiltered, gene_majorlike_NONfiltered,
                        junc_majorlike_juncfiltered, junc_majorlike_NONfiltered,
                        gene_minorlike_juncfiltered, gene_minorlike_NONfiltered,
                        junc_minorlike_juncfiltered, junc_minorlike_NONfiltered,
                        gene_noncanonical_juncfiltered, gene_noncanonical_NONfiltered,
                        junc_noncanonical_juncfiltered, junc_noncanonical_NONfiltered,
                        gene_majorhybrid_juncfiltered, gene_majorhybrid_NONfiltered,
                        junc_majorhybrid_juncfiltered, junc_majorhybrid_NONfiltered,
                        gene_minorhybrid_juncfiltered, gene_minorhybrid_NONfiltered,
                        junc_minorhybrid_juncfiltered,junc_minorhybrid_NONfiltered)
    
    ###############----Make output directories----##############
    level_value <- list()
    for (i in level) {
         if (grepl("genes_w", i)) {
            type <- "gene_level/"         
        }

        if (grepl("juncs_w", i)) {
            type <- "junction_level/"
        }
                
        if (grepl("minor_intron|minor_gene", i)) {
            intron <- "minor/"
        }

        if (grepl("minor_hybrid", i)) {
            intron <- "minor_hybrid/"
        }

        if (grepl("major_hybrid", i)) {
            intron <- "major_hybrid/"
        }

        if (grepl("minor_like", i)) {
            intron <- "minor_like/"
        }

        if (grepl("major_like", i)) {
            intron <- "major_like/"
        }

        if (grepl("non_canonical", i)) {
            intron <- "non_canonical/"
        }

        dir_path_level <- paste0(dir_path, intron, type, i, "/")

        if (!dir.exists(dir_path_level)) {
                print(dir_path_level)
                dir.create(dir_path_level, recursive = TRUE)
        }

        level_value[[length(level_value) + 1]] <- dir_path_level
        level_value <- unlist(level_value)
    }

    ###############----Join and filter dataframes----##############
    for (i in 1:4) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 5:8) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("major-like", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 9:12) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor-like", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 13:16) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("non-canonical", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 17:20) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("major_hybrid", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 21:24) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor_hybrid", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)

    }

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
fibro_filtered_05_abs3 <- fibro_theta_FRASER_output %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")

blood_filtered_05_abs3 <- blood_FRASER_output %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")

#######################################----Fibroblasts----##################################################
dirpath_fibro <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibro_filtered_05_abs3_Jan_30/"
main_function(fibro_filtered_05_abs3, fibro_FRASER_input, MIDB_all, dirpath_fibro, fibro_theta_FRASER_output_fp, fibro_FRASER_input_fp, MIDB_all_fp)
test %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))

#######################################----Fibroblasts----##################################################
dirpath_blood <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/blood_filtered_05_abs3_try2/"
blood_FRASER_input <- blood_FRASER_input %>% select(sampleID) %>% unique
main_function(blood_filtered_05_abs3, blood_FRASER_input, MIDB_all, dirpath_blood, blood_FRASER_output_fp, blood_FRASER_input_fp, MIDB_all_fp)

dir_path <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/fibro_filtered_05_abs3_Jan_28/"
filtered_FRASER_results <- fibro_filtered_05_abs3
FRASER_input_files <- fibro_FRASER_input
minor_intron_file <- MIDB_all
MIDB_file <- MIDB_all

FRASER_results_filtered <- fibro_filtered_05_abs3

library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

args <- commandArgs(TRUE)
FRASER_sample_table_file <- args[1]
FRASER_outputs_file <- args[2]
metadata_file <- args[3]
mig_file <- args[4]
output_directory <- args[5]

#make parameters flexible

get_junction_count <- function(data_type, filter_on){
  results_all_samples <- FRASER %>% filter(type == data_type) 
  
  if (filter_on == "z_score"){
        results_filtered <- results_all_samples %>% filter(abs(zScore) > 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "z_score_over") {
        results_filtered <- results_all_samples %>% filter(zScore > 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "z_score_under") {
        results_filtered <- results_all_samples %>% filter(zScore < 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "p_value") {
        results_filtered <- results_all_samples %>% filter(padjust < 1E-6) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "p_value_delta_psi") {
        results_filtered <- results_all_samples %>% filter(padjust < 1E-6) %>% filter(abs(deltaPsi) > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi") {
        results_filtered <- results_all_samples %>% filter(abs(deltaPsi) > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi_over") {
        results_filtered <- results_all_samples %>% filter(deltaPsi > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi_under") {
        results_filtered <- results_all_samples %>% filter(deltaPsi < 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "all") {
        results_filtered <- results_all_samples %>% select(sampleID, hgncSymbol)
        print("results_all_samples")
        print(results_all_samples)
  }  else {
        print("did not select valid 'filter_on' type")
        #quit()
  }

  results_grouped <- results_filtered %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
  results_grouped
}

get_MIGs_count <- function(data_type, filter_on){
  results_all_samples <- FRASER %>% filter(type == data_type) 
  
  if (filter_on == "z_score"){
        results_filtered <- results_all_samples %>% filter(abs(zScore) > 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "z_score_over") {
        results_filtered <- results_all_samples %>% filter(zScore > 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "z_score_under") {
        results_filtered <- results_all_samples %>% filter(zScore < 2) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "p_value") {
        results_filtered <- results_all_samples %>% filter(padjust < 1E-6) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "p_value_delta_psi") {
        results_filtered <- results_all_samples %>% filter(padjust < 1E-6) %>% filter(abs(deltaPsi) > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi") {
        results_filtered <- results_all_samples %>% filter(abs(deltaPsi) > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi_over") {
        results_filtered <- results_all_samples %>% filter(deltaPsi > 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "delta_psi_under") {
        results_filtered <- results_all_samples %>% filter(deltaPsi < 0.3) %>% select(sampleID, hgncSymbol)
  } else if (filter_on == "all") {
        results_filtered <- results_all_samples %>% select(sampleID, hgncSymbol)
  } else {
        print("did not select valid 'filter_on' type")
        #quit()
  }
  
  MIG_table_select <- MIG_table %>% select(gene_symbol, gene_class) %>% unique
  colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
  colnames(results_filtered) <- c("sampleID", "gene_symbol")
  results_filtered <- results_filtered %>% filter(!is.na(gene_symbol))

  sig_genes_type_all_samples <- left_join(results_filtered, MIG_table_select)
  gene_type_all_samples <- sig_genes_type_all_samples %>% select(sampleID, gene_class) %>% 
    filter(gene_class == "MIG") %>%
    group_by(sampleID) %>% tally() %>% arrange(desc(n))

  MIG_sample_list <- gene_type_all_samples %>% pull(sampleID) %>% unique
  no_events_MIGs <- setdiff(fraser_samples, MIG_sample_list)
  no_MIGS <- data.frame(sampleID = no_events_MIGs, n = 0)

  MIGs_all_samples <- bind_rows(no_MIGS, gene_type_all_samples)
  MIGs_all_samples
}

###########################################################################
##########################-----Read in Files-----##########################
###########################################################################
FRASER <- read_csv(FRASER_outputs_file)
FRASER_sample_table <- read_csv(FRASER_sample_table_file)
metadata <- read_tsv(metadata_file)
MIG_table <- read_csv(mig_file)

###########################################################################
#########################-----FRASER Results-----##########################
###########################################################################

#############################-----FRASER Count All-----##############################
theta_all <- get_junction_count("theta", "all")
psi5_all <- get_junction_count("psi5", "all")
psi3_all <- get_junction_count("psi3", "all")

colnames(theta_all) <- c("RDID", "no_theta_juncs")
colnames(psi5_all) <- c("RDID", "no_psi5_juncs")
colnames(psi3_all) <- c("RDID", "no_psi3_juncs")

#############################-----FRASER Z Score Theta-----##############################
theta_sig_z_score <- get_junction_count("theta", "z_score")
theta_under_z_score <- get_junction_count("theta", "z_score_under")
theta_over_z_score <- get_junction_count("theta", "z_score_over")

colnames(theta_sig_z_score) <- c("RDID", "no_theta_abs_z_score")
colnames(theta_under_z_score) <- c("RDID", "no_theta_z_score_under")
colnames(theta_over_z_score) <- c("RDID", "no_theta_z_score_over")

#############################-----FRASER Z Score Psi5-----##############################
psi5_sig_z_score <- get_junction_count("psi5", "z_score")
psi5_under_z_score <- get_junction_count("psi5", "z_score_under")
psi5_over_z_score <- get_junction_count("psi5", "z_score_over")

colnames(psi5_sig_z_score) <- c("RDID", "no_psi5_abs_z_score")
colnames(psi5_under_z_score) <- c("RDID", "no_psi5_z_score_under")
colnames(psi5_over_z_score) <- c("RDID", "no_psi5_z_score_over")

#############################-----FRASER Z Score Psi3-----##############################
psi3_sig_z_score <- get_junction_count("psi3", "z_score")
psi3_under_z_score <- get_junction_count("psi3", "z_score_under")
psi3_over_z_score <- get_junction_count("psi3", "z_score_over")

colnames(psi3_sig_z_score) <- c("RDID", "no_psi3_abs_z_score")
colnames(psi3_under_z_score) <- c("RDID", "no_psi3_z_score_under")
colnames(psi3_over_z_score) <- c("RDID", "no_psi3_z_score_over")

#############################-----FRASER p value Broad Theta-----##############################
theta_1E6_significant <- get_junction_count("theta", "p_value")
colnames(theta_1E6_significant) <- c("RDID", "no_theta_pval_1E_neg6")

#############################-----FRASER p value, delta psi Broad Theta-----##############################
theta_1E6_deltaPsi_significant <- get_junction_count("theta", "p_value_delta_psi")
colnames(theta_1E6_deltaPsi_significant) <- c("RDID", "no_theta_pval_1E_neg6_abs_deltaPsi")

#############################-----FRASER deltaPsi Theta-----##############################
theta_deltaPsi <- get_junction_count("theta", "delta_psi")
theta_deltaPsi_over <- get_junction_count("theta", "delta_psi_over")
theta_deltaPsi_under <- get_junction_count("theta", "delta_psi_under")

colnames(theta_deltaPsi) <- c("RDID", "no_theta_abs_deltaPsi")
colnames(theta_deltaPsi_over) <- c("RDID", "no_theta_deltaPsi_over")
colnames(theta_deltaPsi_under) <- c("RDID", "no_theta_deltaPsi_under")

#############################-----FRASER deltaPsi Psi3-----##############################
psi3_deltaPsi <- get_junction_count("psi3", "delta_psi")
psi3_deltaPsi_over <- get_junction_count("psi3", "delta_psi_over")
psi3_deltaPsi_under <- get_junction_count("psi3", "delta_psi_under")

colnames(psi3_deltaPsi) <- c("RDID", "no_psi3_abs_deltaPsi")
colnames(psi3_deltaPsi_over) <- c("RDID", "no_psi3_deltaPsi_over")
colnames(psi3_deltaPsi_under) <- c("RDID", "no_psi3_deltaPsi_under")

#############################-----FRASER deltaPsi Psi5-----##############################
psi5_deltaPsi <- get_junction_count("psi5", "delta_psi")
psi5_deltaPsi_over <- get_junction_count("psi5", "delta_psi_over")
psi5_deltaPsi_under <- get_junction_count("psi5", "delta_psi_under")

colnames(psi5_deltaPsi) <- c("RDID", "no_psi5_abs_deltaPsi")
colnames(psi5_deltaPsi_over) <- c("RDID", "no_psi5_deltaPsi_over")
colnames(psi5_deltaPsi_under) <- c("RDID", "no_psi5_deltaPsi_under")

#############################-----Combine FRASER-----##############################
fraser <- left_join(theta_all, psi5_all) %>% left_join(psi3_all) %>% left_join(theta_sig_z_score) %>%
    left_join(theta_under_z_score) %>% left_join(theta_over_z_score) %>% left_join(psi5_sig_z_score) %>%
    left_join(psi5_under_z_score) %>% left_join(psi5_over_z_score) %>% left_join(psi3_sig_z_score) %>%
    left_join(psi3_under_z_score) %>% left_join(psi3_over_z_score) %>% left_join(theta_1E6_significant) %>%
    left_join(theta_1E6_deltaPsi_significant) %>% left_join(theta_deltaPsi) %>% left_join(theta_deltaPsi_over) %>%
    left_join(theta_deltaPsi_under) %>% left_join(psi3_deltaPsi) %>% left_join(psi3_deltaPsi_over) %>%
    left_join(psi3_deltaPsi_under) %>% left_join(psi5_deltaPsi) %>% left_join(psi5_deltaPsi_over) %>%
    left_join(psi5_deltaPsi_under)

     
fraser[is.na(fraser)] <- 0

no_events <- setdiff(FRASER_sample_table$sampleID, fraser$RDID)
no_events_df <- data_frame(RDID=no_events, no_theta_juncs = 0, no_psi3_juncs = 0, no_psi5_juncs = 0,
                           no_theta_abs_z_score = 0, no_theta_z_score_under = 0, no_theta_z_score_over = 0,
                            no_psi5_abs_z_score = 0, no_psi5_z_score_under = 0, no_psi5_z_score_over = 0,
                            no_psi3_abs_z_score = 0, no_psi3_z_score_under = 0, no_psi3_z_score_over = 0,
                            no_theta_pval_1E_neg6 = 0, no_theta_pval_1E_neg6_abs_deltaPsi = 0, 
                            no_theta_abs_deltaPsi = 0, no_theta_deltaPsi_over = 0, no_theta_deltaPsi_under = 0,
                            no_psi5_abs_deltaPsi = 0, no_psi5_deltaPsi_over = 0, no_psi5_deltaPsi_under = 0,
                            no_psi3_abs_deltaPsi = 0, no_psi3_deltaPsi_over = 0, no_psi3_deltaPsi_under = 0)

fraser_all <- bind_rows(no_events_df, fraser)
fraser_samples <- fraser_all %>% pull(RDID)

###########################################################################
########################-----Metadata Results-----#########################
###########################################################################

###########################-----metadata-----##############################
metadata_select <- metadata %>% select(sample_id, indv_id) %>%
                    set_colnames(c("RDID", "Individual_ID"))
metadata_in_fraser <- metadata_select %>% filter(RDID %in% fraser_samples)

###########################################################################
##########################-----MIGs Results-----###########################
###########################################################################

#############################-----FRASER Count All-----##############################
MIGs_theta_all <- get_MIGs_count("theta", "all")
MIGs_psi5_all <- get_MIGs_count("psi5", "all")
MIGs_psi3_all <- get_MIGs_count("psi3", "all")

colnames(MIGs_theta_all) <- c("RDID", "no_MIGs_theta_juncs")
colnames(MIGs_psi5_all) <- c("RDID", "no_MIGs_psi5_juncs")
colnames(MIGs_psi3_all) <- c("RDID", "no_MIGs_psi3_juncs")

###############################-----MIGs theta zScore-----##############################
MIGs_theta_z_score <- get_MIGs_count("theta", "z_score")
MIGs_theta_z_score_over <- get_MIGs_count("theta", "z_score_over")
MIGs_theta_z_score_under <- get_MIGs_count("theta", "z_score_under")

colnames(MIGs_theta_z_score) <- c("RDID", "no_MIGs_theta_abs_z_score")
colnames(MIGs_theta_z_score_over) <- c("RDID", "no_MIGs_theta_z_score_over")
colnames(MIGs_theta_z_score_under) <- c("RDID", "no_MIGs_theta_z_score_under")

###############################-----MIGs psi3 zScore-----##############################
MIGs_psi3_z_score <- get_MIGs_count("psi3", "z_score")
MIGs_psi3_z_score_over <- get_MIGs_count("psi3", "z_score_over")
MIGs_psi3_z_score_under <- get_MIGs_count("psi3", "z_score_under")

colnames(MIGs_psi3_z_score) <- c("RDID", "no_MIGs_psi3_abs_z_score")
colnames(MIGs_psi3_z_score_over) <- c("RDID", "no_MIGs_psi3_z_score_over")
colnames(MIGs_psi3_z_score_under) <- c("RDID", "no_MIGs_psi3_z_score_under")

###############################-----MIGs psi5 zScore-----##############################
MIGs_psi5_z_score <- get_MIGs_count("psi5", "z_score")
MIGs_psi5_z_score_over <- get_MIGs_count("psi5", "z_score_over")
MIGs_psi5_z_score_under <- get_MIGs_count("psi5", "z_score_under")

colnames(MIGs_psi5_z_score) <- c("RDID", "no_MIGs_psi5_abs_z_score")
colnames(MIGs_psi5_z_score_over) <- c("RDID", "no_MIGs_psi5_z_score_over")
colnames(MIGs_psi5_z_score_under) <- c("RDID", "no_MIGs_psi5_z_score_under")

#############################-----MIGs p value Broad Theta-----##############################
MIGs_theta_1E6_significant <- get_MIGs_count("theta", "p_value")
colnames(MIGs_theta_1E6_significant) <- c("RDID", "no_MIGs_theta_pval_1E_neg6")

#############################-----MIGs p value, delta psi Broad Theta-----##############################
MIGs_theta_1E6_deltaPsi_significant <- get_MIGs_count("theta", "p_value_delta_psi")
colnames(MIGs_theta_1E6_deltaPsi_significant) <- c("RDID", "no_MIGs_theta_pval_1E_neg6_abs_deltaPsi")

###############################-----MIGs psi5 deltaPsi-----##############################
MIGs_theta_deltaPsi <- get_MIGs_count("theta", "delta_psi")
MIGs_theta_deltaPsi_over <- get_MIGs_count("theta", "delta_psi_over")
MIGs_theta_deltaPsi_under <- get_MIGs_count("theta", "delta_psi_under")

colnames(MIGs_theta_deltaPsi) <- c("RDID", "no_MIGs_theta_abs_deltaPsi")
colnames(MIGs_theta_deltaPsi_over) <- c("RDID", "no_MIGs_theta_deltaPsi_over")
colnames(MIGs_theta_deltaPsi_under) <- c("RDID", "no_MIGs_theta_deltaPsi_under")

###############################-----MIGs psi3 zScore-----##############################
MIGs_psi3_deltaPsi <- get_MIGs_count("psi3", "delta_psi")
MIGs_psi3_deltaPsi_over <- get_MIGs_count("psi3", "delta_psi_over")
MIGs_psi3_deltaPsi_under <- get_MIGs_count("psi3", "delta_psi_under")

colnames(MIGs_psi3_deltaPsi) <- c("RDID", "no_MIGs_psi3_abs_deltaPsi")
colnames(MIGs_psi3_deltaPsi_over) <- c("RDID", "no_MIGs_psi3_deltaPsi_over")
colnames(MIGs_psi3_deltaPsi_under) <- c("RDID", "no_MIGs_psi3_deltaPsi_under")

###############################-----MIGs psi5 zScore-----##############################
MIGs_psi5_deltaPsi <- get_MIGs_count("psi5", "delta_psi")
MIGs_psi5_deltaPsi_over <- get_MIGs_count("psi5", "delta_psi_over")
MIGs_psi5_deltaPsi_under <- get_MIGs_count("psi5", "delta_psi_under")

colnames(MIGs_psi5_deltaPsi) <- c("RDID", "no_MIGs_psi5_abs_deltaPsi")
colnames(MIGs_psi5_deltaPsi_over) <- c("RDID", "no_MIGs_psi5_deltaPsi_over")
colnames(MIGs_psi5_deltaPsi_under) <- c("RDID", "no_MIGs_psi5_deltaPsi_under")

#############################-----Combine MIGs-----##############################
MIGs <- left_join(MIGs_theta_all, MIGs_psi5_all) %>% left_join(MIGs_psi3_all) %>% left_join(MIGs_theta_z_score) %>%
    left_join(MIGs_theta_z_score_under) %>% left_join(MIGs_theta_z_score_over) %>% left_join(MIGs_psi5_z_score) %>%
    left_join(MIGs_psi5_z_score_under) %>% left_join(MIGs_psi5_z_score_over) %>% left_join(MIGs_psi3_z_score) %>%
    left_join(MIGs_psi3_z_score_under) %>% left_join(MIGs_psi3_z_score_over) %>% left_join(MIGs_theta_1E6_significant) %>%
    left_join(MIGs_theta_1E6_deltaPsi_significant) %>% left_join(MIGs_theta_deltaPsi) %>% left_join(MIGs_theta_deltaPsi_over) %>%
    left_join(MIGs_theta_deltaPsi_under) %>% left_join(MIGs_psi3_deltaPsi) %>% left_join(MIGs_psi3_deltaPsi_over) %>%
    left_join(MIGs_psi3_deltaPsi_under) %>% left_join(MIGs_psi5_deltaPsi) %>% left_join(MIGs_psi5_deltaPsi_over) %>%
    left_join(MIGs_psi5_deltaPsi_under)

     
MIGs[is.na(MIGs)] <- 0

no_events <- setdiff(FRASER_sample_table$sampleID, MIGs$RDID)
no_events_df <- data_frame(RDID=no_events, no_MIGs_theta_juncs = 0, no_MIGs_psi3_juncs = 0, no_MIGs_psi5_juncs = 0,
                           no_MIGs_theta_abs_z_score = 0, no_MIGs_theta_z_score_under = 0, no_MIGs_theta_z_score_over = 0,
                            no_MIGs_psi5_abs_z_score = 0, no_MIGs_psi5_z_score_under = 0, no_MIGs_psi5_z_score_over = 0,
                            no_MIGs_psi3_abs_z_score = 0, no_MIGs_psi3_z_score_under = 0, no_MIGs_psi3_z_score_over = 0,
                            no_MIGs_theta_pval_1E_neg6 = 0, no_MIGs_theta_pval_1E_neg6_abs_deltaPsi = 0, 
                            no_MIGs_theta_abs_deltaPsi = 0, no_MIGs_theta_deltaPsi_over = 0, no_MIGs_theta_deltaPsi_under = 0,
                            no_MIGs_psi5_abs_deltaPsi = 0, no_MIGs_psi5_deltaPsi_over = 0, no_MIGs_psi5_deltaPsi_under = 0,
                            no_MIGs_psi3_abs_deltaPsi = 0, no_MIGs_psi3_deltaPsi_over = 0, no_MIGs_psi3_deltaPsi_under = 0)

MIGs_all <- bind_rows(no_events_df, MIGs)

###########################################################################
##########################-----Other Metrics-----##########################
###########################################################################

###########################-----psi3/ psi5-----############################
fraser_all$psi3_psi5 <- fraser_all$no_psi3_juncs + fraser_all$no_psi5_juncs

###########################################################################
########################-----Join Dataframes-----##########################
###########################################################################

#########################-----joined_raw-----##############################
joined_raw <- left_join(fraser_all, MIGs_all) %>% left_join(metadata_in_fraser)
names(joined_raw)[names(joined_raw) == 'RDID'] <- 'Sample_ID'

filename_raw <- paste0(output_directory, "/FRASER_results_raw.csv", sep="")
write.csv(joined_raw, filename_raw, row.names=FALSE)

#######################-----joined_zscores-----############################
metadata_z <- joined_raw %>% select(c(Sample_ID, Individual_ID))
joined_z <- joined_raw %>% select(-c(Sample_ID, Individual_ID))
z_scores <- joined_z %>% mutate_all(~scale(.))

joined_z_select <- bind_cols(metadata_z, z_scores) %>% as.data.frame()
joined_z_select_colnames <- colnames(joined_z_select)

column_names_zscore <- paste("zscore_", joined_z_select_colnames, sep="")
colnames(joined_z_select) <- column_names_zscore
names(joined_z_select)[names(joined_z_select) == 'zscore_Sample_ID'] <- 'Sample_ID'
names(joined_z_select)[names(joined_z_select) == 'zscore_Individual_ID'] <- 'Individual_ID'
                           
filename_zscores <- paste0(output_directory, "/FRASER_results_zscores.csv", sep="")
write.csv(as.matrix(joined_z_select), filename_zscores, row.names=FALSE)

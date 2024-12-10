######################################################
##############----Split Dataframes----################
######################################################
theta <- FRASER_metadata_important %>% select(no_theta_juncs, no_theta_abs_z_score, no_theta_z_score_under, no_theta_z_score_over,
            no_theta_pval_1E_neg6, no_theta_pval_1E_neg6_abs_deltaPsi, no_theta_abs_deltaPsi, no_theta_deltaPsi_over, no_theta_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

theta_zscore <- FRASER_metadata_important %>% select(zscore_no_theta_juncs, zscore_no_theta_abs_z_score, zscore_no_theta_z_score_under, zscore_no_theta_z_score_over,
            zscore_no_theta_pval_1E_neg6, zscore_no_theta_pval_1E_neg6_abs_deltaPsi, zscore_no_theta_abs_deltaPsi, zscore_no_theta_deltaPsi_over, zscore_no_theta_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

psi3 <- FRASER_metadata_important %>% select(no_psi3_juncs, no_psi3_abs_z_score, no_psi3_z_score_under, no_psi3_z_score_over,
            no_psi3_abs_deltaPsi, no_psi3_deltaPsi_over, no_psi3_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

psi3_zscore<- FRASER_metadata_important %>% select(zscore_no_psi3_juncs, zscore_no_psi3_abs_z_score, zscore_no_psi3_z_score_under, zscore_no_psi3_z_score_over,
            zscore_no_psi3_abs_deltaPsi, zscore_no_psi3_deltaPsi_over, zscore_no_psi3_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

psi5 <- FRASER_metadata_important %>% select(no_psi5_juncs, no_psi5_abs_z_score, no_psi5_z_score_under, no_psi5_z_score_over,
            no_psi5_abs_deltaPsi, no_psi5_deltaPsi_over, no_psi5_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

psi5_zscore <- FRASER_metadata_important %>% select(zscore_no_psi5_juncs, zscore_no_psi5_abs_z_score, zscore_no_psi5_z_score_under, zscore_no_psi5_z_score_over,
            zscore_no_psi5_abs_deltaPsi, zscore_no_psi5_deltaPsi_over, zscore_no_psi5_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

FRASER_metadata_important$all <- FRASER_metadata_important$no_theta_juncs + FRASER_metadata_important$no_psi3_juncs + FRASER_metadata_important$no_psi5_juncs
combined <- FRASER_metadata_important %>% select(psi3_psi5, all,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

MIGs_zscore <- FRASER_metadata_important %>% select(zscore_no_MIGs_theta_juncs, zscore_no_MIGs_psi3_juncs, zscore_no_MIGs_psi5_juncs,                   
              zscore_no_MIGs_theta_abs_z_score, zscore_no_MIGs_theta_z_score_under, zscore_no_MIGs_theta_z_score_over,             
              zscore_no_MIGs_psi5_abs_z_score, zscore_no_MIGs_psi5_z_score_under, zscore_no_MIGs_psi5_z_score_over,              
              zscore_no_MIGs_psi3_abs_z_score, zscore_no_MIGs_psi3_z_score_under, zscore_no_MIGs_psi3_z_score_over,              
              zscore_no_MIGs_theta_pval_1E_neg6, zscore_no_MIGs_theta_pval_1E_neg6_abs_deltaPsi, zscore_no_MIGs_theta_abs_deltaPsi,             
              zscore_no_MIGs_theta_deltaPsi_over, zscore_no_MIGs_theta_deltaPsi_under, zscore_no_MIGs_psi5_abs_deltaPsi,              
              zscore_no_MIGs_psi5_deltaPsi_over, zscore_no_MIGs_psi5_deltaPsi_under, zscore_no_MIGs_psi3_abs_deltaPsi,              
              zscore_no_MIGs_psi3_deltaPsi_over, zscore_no_MIGs_psi3_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

MIGs <- FRASER_metadata_important %>% select(no_MIGs_theta_juncs, no_MIGs_psi3_juncs, no_MIGs_psi5_juncs,                   
              no_MIGs_theta_abs_z_score, no_MIGs_theta_z_score_under, no_MIGs_theta_z_score_over,             
              no_MIGs_psi5_abs_z_score, no_MIGs_psi5_z_score_under, no_MIGs_psi5_z_score_over,              
              no_MIGs_psi3_abs_z_score, no_MIGs_psi3_z_score_under, no_MIGs_psi3_z_score_over,              
              no_MIGs_theta_pval_1E_neg6, no_MIGs_theta_pval_1E_neg6_abs_deltaPsi, no_MIGs_theta_abs_deltaPsi,             
              no_MIGs_theta_deltaPsi_over, no_MIGs_theta_deltaPsi_under, no_MIGs_psi5_abs_deltaPsi,              
              no_MIGs_psi5_deltaPsi_over, no_MIGs_psi5_deltaPsi_under, no_MIGs_psi3_abs_deltaPsi,              
              no_MIGs_psi3_deltaPsi_over, no_MIGs_psi3_deltaPsi_under,
            wetlab_id, RDID, UDN_ID, GSS_ID,
            affected_status, sex, age, solved_status, patient_description, HPO_terms, HPO_IDs, 
            candidate_gene, candidate_variants, disease_category, condition_id, condition_inheritance, UDN_causal_gene, gene_known_for_phenotype, 
            RIN, batch, seq_batch, notes,
            family_identifier, GREGoR_maternal_id, GREGoR_paternal_id, proband_relationship, GREGoR_twin_id, additional_family_members_with_variant, UDN_family_type, relationship_details, UDN_related_proband, family_identifier, family_sequenced)

write_csv(theta, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/theta.csv")
write_csv(theta_zscore, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/zscore_theta.csv")
write_csv(psi3, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/psi3.csv")
write_csv(psi3_zscore, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/zscore_psi3.csv")
write_csv(psi5, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/psi5.csv")
write_csv(psi5_zscore, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/zscore_psi5.csv")
write_csv(combined, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/combined.csv")
write_csv(MIGs, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/MIGs.csv")
write_csv(MIGs_zscore, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/zscore_MIGs.csv")

######################################################
############----Get Z-score combined----##############
######################################################
combined_mean_psi3_psi5 <- mean(combined$psi3_psi5)
combined_sd_psi3_psi5 <- sd(combined$psi3_psi5)
combined$zscore_psi3_psi_5 <- ((combined$psi3_psi5 - combined_mean_psi3_psi5)/ combined_sd_psi3_psi5)

combined_mean_all <- mean(combined$all)
combined_sd_all <- sd(combined$all)
combined$zscore_all <- ((combined$all - combined_mean_all)/ combined_sd_all)

combined_zscore <- combined %>% select(-psi3_psi5, -all)
write_csv(combined_zscore, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/zscore_combined.csv")

######################################################
##############----Get Top Z-scores----################
######################################################

##################----Theta----#######################
theta_zscore$ID <- paste(theta_zscore$GSS_ID, theta_zscore$UDN_ID, theta_zscore$RDID, theta_zscore$wetlab_id, sep="_")
no_theta_juncs <- theta_zscore %>% filter(zscore_no_theta_juncs > 2) %>% select(ID) %>% as.list
no_theta_abs_z_score <- theta_zscore %>% filter(zscore_no_theta_abs_z_score > 2) %>% select(ID) %>% as.list
no_theta_z_score_under <- theta_zscore %>% filter(zscore_no_theta_z_score_under > 2) %>% select(ID) %>% as.list
no_theta_z_score_over <- theta_zscore %>% filter(zscore_no_theta_z_score_over > 2) %>% select(ID) %>% as.list
no_theta_pval_1E_neg6 <- theta_zscore %>% filter(zscore_no_theta_pval_1E_neg6 > 2) %>% select(ID) %>% as.list
no_theta_pval_1E_neg6_abs_deltaPsi <- theta_zscore %>% filter(zscore_no_theta_pval_1E_neg6_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_theta_abs_deltaPsi <- theta_zscore %>% filter(zscore_no_theta_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_theta_deltaPsi_over <- theta_zscore %>% filter(zscore_no_theta_deltaPsi_over > 2) %>% select(ID) %>% as.list
no_theta_deltaPsi_under <- theta_zscore %>% filter(zscore_no_theta_deltaPsi_under > 2) %>% select(ID) %>% as.list

####################----psi3----######################
psi3_zscore$ID <- paste(psi3_zscore$GSS_ID, psi3_zscore$UDN_ID, psi3_zscore$RDID, psi3_zscore$wetlab_id,  sep="_")
no_psi3_juncs <- psi3_zscore %>% filter(zscore_no_psi3_juncs > 2) %>% select(ID) %>% as.list
no_psi3_z_score_under <- psi3_zscore %>% filter(zscore_no_psi3_z_score_under > 2) %>% select(ID) %>% as.list
no_psi3_z_score_over <- psi3_zscore %>% filter(zscore_no_psi3_z_score_over > 2)%>% select(ID) %>% as.list
no_psi3_abs_deltaPsi <- psi3_zscore %>% filter(zscore_no_psi3_abs_deltaPsi > 2)%>% select(ID) %>% as.list
no_psi3_deltaPsi_over <- psi3_zscore %>% filter(zscore_no_psi3_deltaPsi_over > 2) %>% select(ID) %>% as.list
no_psi3_deltaPsi_under <- psi3_zscore %>% filter(zscore_no_psi3_deltaPsi_under > 2) %>% select(ID) %>% as.list

####################----psi5----######################
psi5_zscore$ID <- paste(psi5_zscore$GSS_ID, psi5_zscore$UDN_ID, psi5_zscore$RDID, psi5_zscore$wetlab_id,  sep="_")
no_psi5_juncs <- psi5_zscore %>% filter(zscore_no_psi5_juncs > 2) %>% select(ID) %>% as.list
no_psi5_z_score_under <- psi5_zscore %>% filter(zscore_no_psi5_z_score_under > 2) %>% select(ID) %>% as.list
no_psi5_z_score_over <- psi5_zscore %>% filter(zscore_no_psi5_z_score_over > 2) %>% select(ID) %>% as.list
no_psi5_abs_deltaPsi <- psi5_zscore %>% filter(zscore_no_psi5_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_psi5_deltaPsi_over <- psi5_zscore %>% filter(zscore_no_psi5_deltaPsi_over > 2) %>% select(ID) %>% as.list
no_psi5_deltaPsi_under <- psi5_zscore %>% filter(zscore_no_psi5_deltaPsi_under > 2) %>% select(ID) %>% as.list

#################----combined----#####################
combined_zscore$ID <- paste(combined_zscore$GSS_ID, combined_zscore$UDN_ID, combined_zscore$RDID, combined_zscore$wetlab_id,  sep="_")
no_psi3_psi5 <- combined_zscore %>% filter(zscore_psi3_psi_5 > 3) %>% select(ID) %>% as.list
no_all <- combined_zscore %>% filter(zscore_all > 3) %>% select(ID) %>% as.list

###################----MIGs----#######################
MIGs_zscore$ID <- paste(MIGs_zscore$GSS_ID, MIGs_zscore$UDN_ID, MIGs_zscore$RDID, MIGs_zscore$wetlab_id,  sep="_")
no_MIGs_theta_juncs <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_juncs > 2) %>% select(ID) %>% as.list
no_MIGs_psi3_juncs <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_juncs > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_juncs <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_juncs > 2) %>% select(ID) %>% as.list
no_MIGs_theta_abs_z_score <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_abs_z_score > 2) %>% select(ID) %>% as.list
no_MIGs_theta_z_score_under <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_z_score_under > 2) %>% select(ID) %>% as.list
no_MIGs_theta_z_score_over <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_z_score_over > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_abs_z_score <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_abs_z_score > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_z_score_under <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_z_score_under > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_z_score_over <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_z_score_over > 2) %>% select(ID) %>% as.list
no_MIGs_psi3_abs_z_score <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_abs_z_score > 2) %>% select(ID) %>% as.list
no_MIGs_psi3_z_score_under <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_z_score_under > 2) %>% select(ID) %>% as.list
no_MIGs_psi3_z_score_over <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_z_score_over > 2) %>% select(ID) %>% as.list
no_MIGs_theta_pval_1E_neg6 <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_pval_1E_neg6 > 2) %>% select(ID) %>% as.list
no_MIGs_theta_pval_1E_neg6_abs_deltaPsi <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_pval_1E_neg6_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_MIGs_theta_abs_deltaPsi <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_MIGs_theta_deltaPsi_over <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_deltaPsi_over > 2) %>% select(ID) %>% as.list
no_MIGs_theta_deltaPsi_under <- MIGs_zscore %>% filter(zscore_no_MIGs_theta_deltaPsi_under > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_abs_deltaPsi <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_abs_deltaPsi > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_deltaPsi_over <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_deltaPsi_over > 2) %>% select(ID) %>% as.list
no_MIGs_psi5_deltaPsi_under <- MIGs_zscore %>% filter(zscore_no_MIGs_psi5_deltaPsi_under > 2) %>% select(ID) %>% as.list
no_MIGs_psi3_abs_deltaPsi <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_abs_deltaPsi > 2)%>% select(ID) %>% as.list
no_MIGs_psi3_deltaPsi_under <- MIGs_zscore %>% filter(zscore_no_MIGs_psi3_deltaPsi_under > 2)%>% select(ID) %>% as.list

######################################################
#########----Create Dataframe of Conern----###########
######################################################
FRASER_metadata_important$ID <- paste(FRASER_metadata_important$GSS_ID, FRASER_metadata_important$UDN_ID, FRASER_metadata_important$RDID, FRASER_metadata_important$wetlab_id,  sep="_")

##################----Theta----#######################
theta_outlier <- c(no_theta_juncs, no_theta_abs_z_score, no_theta_z_score_under, no_theta_z_score_over, no_theta_pval_1E_neg6,
                    no_theta_pval_1E_neg6_abs_deltaPsi, no_theta_abs_deltaPsi, no_theta_deltaPsi_over, no_theta_deltaPsi_under) %>% unlist %>% unique()

FRASER_metadata_theta <- FRASER_metadata_important %>% filter(is_in(ID, theta_outlier))
FRASER_metadata_theta$outlier_type <- "theta"

###################----Psi3----#######################
psi3_outlier <- c(no_psi3_juncs, no_psi3_z_score_under, no_psi3_z_score_over, no_psi3_abs_deltaPsi, no_psi3_deltaPsi_over,
                    no_psi3_deltaPsi_under) %>% unlist %>% unique()

FRASER_metadata_psi3 <- FRASER_metadata_important %>% filter(is_in(ID, psi3_outlier))
FRASER_metadata_psi3$outlier_type<- "psi3"

###################----Psi5----#######################
psi5_outlier <- c(no_psi5_juncs, no_psi5_z_score_under, no_psi5_z_score_over, no_psi5_abs_deltaPsi, no_psi5_deltaPsi_over,
                    no_psi5_deltaPsi_under) %>% unlist %>% unique()

FRASER_metadata_psi5 <- FRASER_metadata_important %>% filter(is_in(ID, psi5_outlier))
FRASER_metadata_psi5$outlier_type <- "psi5"

#################----combined----#####################
combined_outlier <- c(no_psi3_psi5, no_all) %>% unlist %>% unique()

FRASER_metadata_combined <- FRASER_metadata_important %>% filter(is_in(ID, combined_outlier))
FRASER_metadata_combined$outlier_type <- "combined"

###################----MIGs----#######################
MIGs_outlier <- c(no_MIGs_theta_juncs, no_MIGs_psi3_juncs, no_MIGs_psi5_juncs,no_MIGs_theta_abs_z_score,
                        no_MIGs_theta_z_score_under, no_MIGs_theta_z_score_over, no_MIGs_psi5_abs_z_score,
                        no_MIGs_psi5_z_score_under, no_MIGs_psi5_z_score_over, no_MIGs_psi3_abs_z_score,
                        no_MIGs_psi3_z_score_under, no_MIGs_psi3_z_score_over, no_MIGs_theta_pval_1E_neg6,
                        no_MIGs_theta_pval_1E_neg6_abs_deltaPsi,no_MIGs_theta_abs_deltaPsi, no_MIGs_theta_deltaPsi_over,
                        no_MIGs_theta_deltaPsi_under, no_MIGs_psi5_abs_deltaPsi, no_MIGs_psi5_deltaPsi_over, 
                        no_MIGs_psi5_deltaPsi_under, no_MIGs_psi3_abs_deltaPsi, no_MIGs_psi3_deltaPsi_under) %>% unlist %>% unique()

FRASER_metadata_MIG <- FRASER_metadata_important %>% filter(is_in(ID, MIGs_outlier))
FRASER_metadata_MIG$outlier_type <- "MIG"

#############----combine dataframes----###############
all <- list(FRASER_metadata_theta$ID, FRASER_metadata_psi3$ID, FRASER_metadata_psi5$ID,
                 FRASER_metadata_combined$ID, FRASER_metadata_MIG$ID) %>% unlist %>% unique
none <- setdiff(FRASER_metadata_important$ID, all)
FRASER_metadata_none <- filter(FRASER_metadata_important, is_in(ID, none))

FRASER_metadata_all <- bind_rows(FRASER_metadata_MIG, FRASER_metadata_combined, FRASER_metadata_psi5, FRASER_metadata_psi3, FRASER_metadata_theta)

FRASER_metadata_all <- FRASER_metadata_all %>% group_by(ID) %>% 
       mutate(outlier_type = paste(outlier_type, collapse=", ")) %>% unique

FRASER_metadata_none$outlier_type <- "NA"

FRASER_metadata_important <- bind_rows(FRASER_metadata_all, FRASER_metadata_none)

#############----combine dataframes----###############
#this part isn't working
#FRASER_metadata_important <- FRASER_metadata_important %>% mutate(outlier_type = ifelse(outlier_type == "NA", NA, FRASER_metadata_important$outlier_type)) 
concerning_splicing_patients <- FRASER_metadata_important %>% filter(outlier_type != "NA") %>%
                                   select(ID, GSS_ID, UDN_ID, wetlab_id, RDID, affected_status, outlier_type, patient_description,
                                   solved_status, condition_id, UDN_causal_gene, candidate_variants, candidate_gene,
                                   disease_category, HPO_terms, HPO_IDs, RIN, family_identifier, age, sex)

#############----find missing RIN----###############


write_csv(concerning_splicing_patients, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/concerning_splicing_patterns.csv")

#need to add number of genes outliers
#look at one a gene level

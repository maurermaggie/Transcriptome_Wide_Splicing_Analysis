library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

metadata_UDN <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/2017to2022_metadata.tsv")
metadata_GREGoR<- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/Rare_Disease_Metadata.xlsx") 
shruti <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/stanford_participant_metadata2.csv") %>%
        select(-ethnicity, -race, -otherRace, -primaryLanguage, state) %>% filter()
GSS_participant <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/participant.tsv")
GSS_phenotype <- read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/phenotype.tsv")
GSS_genetics <-read_tsv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/genetic_findings.tsv")
missing_1 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/Missing_Metadata_1.csv")
#27
missing_2 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/Missing_Metadata_2.csv")
#68
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_no_missing/FRASER_output.csv") %>% select(sampleID) %>% unique
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_no_missing/FRASER_input.csv")
missing_Miami <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_metadata_miami_NIH.csv") %>% select(-`...6`)
GSS_UDN_MAP <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/GSS_UDN_MAP.csv")
RNA_seq_info <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/RIN_sequencing_info.csv")
missing_age_df <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_age_new.csv") %>% select(GSS_ID, RDID, UDN_ID, `AGE AT COLLECTION`)
batchz <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchz.txt")
batchy <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchy.txt")
batchx <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/batchx.txt")
new_batch_info <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_batch_new.csv")
batch_info <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/Kevins_Batch_plot.csv") %>% select(PPID, `RD - ID`, `SEQ RUN #`)
rachels_table <- fread("/home/maurertm/smontgom/shared/UDN/Data/SampleTables/sample_table_UDN_sep2023.txt")

######################################################
###########----Join and Clean Metadata----############
######################################################

#################----UDN Metadata----#################
metadata_UDN_no_NA_columns <- metadata_UDN[,colSums(is.na(metadata_UDN))<nrow(metadata_UDN)]
metadata_UDN_no_NA_columns <- metadata_UDN_no_NA_columns[vapply(metadata_UDN_no_NA_columns, function(x) length(unique(x)) > 1, logical(1L))]

#they are not always equal- will need to check matches to both
ifelse(metadata_UDN_no_NA_columns$wetlab_id == metadata_UDN_no_NA_columns$sample_id, "Equal", "Not Equal")

#################----Batch vs Seq Batch----#################
identical(metadata_UDN_no_NA_columns[['seq_batch...8']],metadata_UDN_no_NA_columns[['seq_batch...56']])
identical(metadata_UDN_no_NA_columns[['batch']],metadata_UDN_no_NA_columns[['seq_batch...8']])

metadata_UDN_no_NA_columns %>% select(batch, seq_batch...8, indv_id) %>% print(n=600)
metadata_UDN_no_NA_columns$seq_batch...8 <- metadata_UDN_no_NA_columns$batch

metadata_UDN_no_NA_columns <- metadata_UDN_no_NA_columns %>% select(-seq_batch...56)

#################----GREGoR Metadata----##############
metadata_GREGoR_no_NA_columns <- metadata_GREGoR[,colSums(is.na(metadata_GREGoR))<nrow(metadata_GREGoR)]
metadata_GREGoR_no_NA_columns <- metadata_GREGoR_no_NA_columns[vapply(metadata_GREGoR_no_NA_columns, function(x) length(unique(x)) > 1, logical(1L))]

#they are all equal (except cases with NA in both, so I kept only one column)
ifelse(metadata_GREGoR_no_NA_columns$indv_id == metadata_GREGoR_no_NA_columns$institution_id, "Equal", "Not Equal")
ifelse(metadata_GREGoR_no_NA_columns$wetlab_id == metadata_GREGoR_no_NA_columns$sample_id, "Equal", "Not Equal")
#keeping wetlab_id and indv_id to differentiate which one have RNA-seq? And so it matches with metadata_UDN

################----Combine Metadata----##############
names(metadata_UDN_no_NA_columns)[names(metadata_UDN_no_NA_columns) == 'seq_batch...8'] <- 'seq_batch'
setdiff(colnames(metadata_GREGoR_no_NA_columns), colnames(metadata_UDN_no_NA_columns))
setdiff(colnames(metadata_UDN_no_NA_columns), colnames(metadata_GREGoR_no_NA_columns))


metadata_GREGoR_pruned <- metadata_GREGoR_no_NA_columns %>% select(-is_RD)
metadata_UDN_pruned <- metadata_UDN_no_NA_columns %>% select(-hg19_vcf_wasp, -hg38_vcf_wasp, -fastq_dir, -bcl_dir,
                                -index1, -index2, -genome_vcf, -exome_vcf, -library_index_I5, -library_index_I7, -library_index_identifier,
                                -library_conc_nM, -library_conc_source, -library_conc_ng_per_ul, -library_avg_bp, -institution,
                                -initials, -prioritize, -demultiplex_status, -analysis_notes, -origin, -dbGap_sharing_restriction,
                                -date_collected, -date_prepped, -date_sequenced, -"260_280", -"260_230", rna_conc_qubit_ng_per_ul,
                                -paxgene_tube_num, -bloodcenter_barcode, -globin_depletion, -library_notes, -library_avg_bp, -library_conc_ng_per_ul,
                                -library_conc_source, -library_conc_nM, -library_index_identifier, -seq_machine, -sample_notes,
                                -rna_conc_qubit_ng_per_ul)

metadata_GREGoR_pruned$age <- as.character(metadata_GREGoR_pruned$age)
joined_metadata <- bind_rows(metadata_GREGoR_pruned, metadata_UDN_pruned)

##############----Combine with Shruti----############
setdiff(colnames(shruti), colnames(joined_metadata))
setdiff(colnames(joined_metadata), colnames(shruti))
names(shruti)[names(shruti) == 'symptom_category'] <- 'disease_category'

shruti <- shruti %>% group_by(UDN_ID) %>% 
       mutate(diagnosis_name = paste(diagnosis_name, collapse=",")  %>% unique) %>%
       mutate(diagnosis_certainty = paste(diagnosis_certainty, collapse=",")) %>%
       mutate(causative_gene = paste(causative_gene, collapse=",")) %>%
       mutate(diagnosis_method = paste(diagnosis_method, collapse=",")) %>%
       mutate(diagnosis_method_type = paste(diagnosis_method_type, collapse=",")) %>%
       mutate(phenotype_degree_explained = paste(phenotype_degree_explained, collapse=",")) %>%
       mutate(diagnosis_consequences = paste(diagnosis_consequences, collapse=",")) %>%
       mutate(diagnosis_method_explain = paste(diagnosis_method_explain, collapse=",")) %>%
       mutate(variant_type = paste(variant_type, collapse=",")) %>% unique

overlap_institution_id <- intersect(shruti$UDN_ID, joined_metadata$institution_id)
in_instit_id <- filter(joined_metadata, is_in(institution_id, overlap_institution_id))
#452
not_in_instit_id <- filter(joined_metadata, !is_in(institution_id, overlap_institution_id))
#205

overlap_indv_id <- intersect(shruti$UDN_ID, not_in_instit_id$indv_id)
in_indv_id <- filter(not_in_instit_id, is_in(indv_id, overlap_indv_id))
#49
not_in_indv_id_or_instit_id <- filter(not_in_instit_id, !is_in(indv_id, overlap_indv_id))
#156

#total 657
#156+49+452=657

names(shruti)[names(shruti) == 'UDN_ID'] <- 'institution_id'
in_instit_id_joined <- left_join(in_instit_id, shruti, by ="institution_id")

names(shruti)[names(shruti) == 'institution_id'] <- 'indv_id'
in_indv_id_joined <- left_join(in_indv_id, shruti, by ="indv_id")

joined_shruti <- bind_rows(in_instit_id_joined, in_indv_id_joined, not_in_indv_id_or_instit_id)
#657 long

##############----Combine with GSS Info----############
GSS_phenotype_joined <- GSS_phenotype %>% group_by(participant_id) %>% 
       mutate(additional_details = paste(additional_details, collapse=";")  %>% unique) %>%
       mutate(additional_modifiers = paste(additional_modifiers, collapse=";")) %>%
       mutate(onset_age_range = paste(onset_age_range, collapse=";")  %>% unique) %>%
       mutate(ontology = paste(ontology, collapse=";")) %>%
       mutate(phenotype_id = paste(phenotype_id, collapse=";")  %>% unique) %>%
       mutate(presence = paste(presence, collapse=";")) %>%
       mutate(term_id = paste(term_id, collapse=";")) %>% unique

setdiff(GSS_participant$participant_id, GSS_phenotype$participant_id)
setdiff(GSS_phenotype$participant_id, GSS_participant$participant_id)

GSS_genetics_joined <- GSS_genetics %>% group_by(participant_id) %>% 
       mutate(additional_family_members_with_variant = paste(additional_family_members_with_variant, collapse=";")  %>% unique) %>%
       mutate(condition_id = paste(condition_id, collapse=";")) %>%
       mutate(condition_inheritance = paste(condition_inheritance, collapse=";")  %>% unique) %>%
       mutate(gene = paste(gene, collapse=";")) %>%
       mutate(known_condition_name = paste(known_condition_name, collapse=";")  %>% unique) %>%
       mutate(gene_known_for_phenotype = paste(gene_known_for_phenotype, collapse=";")) %>%
       mutate(partial_contribution_explained = paste(partial_contribution_explained, collapse=";")) %>%
       select(additional_family_members_with_variant, condition_id, condition_inheritance, gene, known_condition_name,
       gene_known_for_phenotype, partial_contribution_explained) %>% unique

setdiff(GSS_participant$participant_id, GSS_genetics_joined$participant_id)
setdiff(GSS_genetics_joined$participant_id, GSS_participant$participant_id)

setdiff(GSS_phenotype$participant_id, GSS_genetics_joined$participant_id)
setdiff(GSS_genetics_joined$participant_id, GSS_phenotype$participant_id)

GSS_results <- left_join(GSS_participant, GSS_phenotype_joined) %>% left_join(GSS_genetics_joined)

setdiff(colnames(GSS_results), colnames(joined_metadata))
setdiff(colnames(joined_metadata), colnames(GSS_results))
#need to get RIN for GREGoR cases

intersect(GSS_results$participant_id, joined_metadata$institution_id)
intersect(GSS_results$participant_id, joined_metadata$indv_id)
intersect(GSS_results$participant_id, joined_metadata$wetlab_id)

joined_metadata$family_id <- as.character(joined_metadata$family_id)
names(GSS_results)[names(GSS_results) == 'sex'] <- 'sex_GREGoR'
names(GSS_results)[names(GSS_results) == 'affected_status'] <- 'affected_status_GREGoR'
names(GSS_results)[names(GSS_results) == 'family_id'] <- 'family_id_GREGoR'

joined_GSS <- bind_rows(GSS_results, joined_metadata)
#327+657=984

#########----Combine with Missing Metadata----########
setdiff(missing_1$institution_id, missing_2$institution_id)
setdiff(missing_2$institution_id, missing_1$institution_id)
setdiff(colnames(missing_1), colnames(missing_2))
setdiff(colnames(missing_2), colnames(missing_1))

missing_1 <- missing_1 %>% select(-Relationship)
names(missing_1)[names(missing_1) == 'notes (list candidate gene name if entered in Gateway, otherwise no genetic candidates)'] <- 'notes'

missing_2 <- missing_2 %>% select(-Test, -`Both Parents?`)

intersect(missing_1$institution_id, missing_2$institution_id)

missing_combined <- bind_rows(missing_1, missing_2)
#27+68=95

names(missing_combined)[names(missing_combined) == 'institution_id'] <- 'ID'
intersect(joined_GSS$institution_id, missing_combined$ID)
intersect(joined_GSS$wetlab_id, missing_combined$ID)
intersect(joined_GSS$sample_id, missing_combined$ID)
intersect(joined_GSS$indv_id, missing_combined$ID)

in_both <- intersect(joined_GSS$institution_id, missing_combined$ID)
setdiff(colnames(joined_GSS), colnames(missing_combined))
setdiff(colnames(missing_combined), colnames(joined_GSS))

missing_combined <- missing_combined %>% select(-mutation_chr, -mutation_bp, -genome_vcf_on_durga, -number_of_reads_mapped_to_transcriptome,
-exome_vcf, -number_of_uniquely_mapped_reads, -percent_duplicates, -number_after_filtering_MAPQ_30, -number_after_deduplication, -number_of_mapped_reads,
-number_of_trimmed_reads, -number_of_raw_reads, -sequencing_status, -"date sequenced", -passage, -bioanalyzer_cDNA_conc_pM,
-"260_230", -"260_280", -avg_library_size_, -globin_depleted, -dbGap_sharing_restriction, -INDEX_ID, -"Collection Date")

joined_GSS <- joined_GSS %>% select(-consent_code, -internal_project_id, -pmid_id, -recontactable)
#974
missing_filtered <- missing_combined %>% filter(!is_in(ID, in_both))
#91

missing_filtered$age <- as.character(missing_filtered$age)
joined_missing <- bind_rows(missing_filtered, joined_GSS)
#1075 = 984+91

######################################################
################----Join FRASER data----##############
######################################################
missing_dataframe <- setdiff(files$sampleID, all_uncompiled$sampleID)
missing_dataframe_df <- data.frame(sampleID = missing_dataframe)
FRASER_result_combined <- bind_rows(all_uncompiled, missing_dataframe_df)
colnames(FRASER_result_combined) <- c("Sample_ID")

######################################################
#############----Simplify joined_missing----##########
######################################################
#~~~~~~Split joined_missing into three parts~~~~~~#
joined_missing$number <- rownames(joined_missing)
joined_missing_wetlab_id_doubles <- joined_missing %>% filter(!is.na(wetlab_id)) %>% group_by(wetlab_id) %>% tally() %>%
                                   arrange(desc(n)) %>% filter(n >= 2) %>% select(wetlab_id) %>% as.list %>% unlist
joined_missing_wetlab_id_singles <- joined_missing %>% filter(!is.na(wetlab_id)) %>% group_by(wetlab_id) %>% tally() %>%
                                   arrange(desc(n)) %>% filter(n == 1) %>% select(wetlab_id) %>% as.list %>% unlist
joined_missing_wetlab_id_singles <- joined_missing %>% filter(wetlab_id %in% joined_missing_wetlab_id_singles)                                   
joined_missing_wetlab_id_na <- joined_missing %>% filter(is.na(wetlab_id))

#~~~~~~Filter doubles~~~~~~#
doubles <- filter(joined_missing, wetlab_id %in% joined_missing_wetlab_id_doubles)

doubles_filtered <- doubles %>% filter(sample_type == "PAXGene RNA tube")
doubles_filtered_inv <- doubles %>% filter(is.na(sample_type))
setdiff(doubles_filtered_inv$wetlab_id, doubles_filtered$wetlab_id)
#"RD320" "RD338"
doubles_filtered_inv <- doubles_filtered_inv %>% filter(!is_in(wetlab_id, c("RD320", "RD338")))
doubles_filtered_unique <- doubles_filtered_inv %>% filter(is_in(wetlab_id, c("RD320", "RD338")))

#~~~~~~Test for similarities between the doubles groups~~~~~~#
dataframe_test <- data.frame(inv = (doubles_filtered_inv %>% arrange(desc(wetlab_id)) %>% select(HPO_terms)), pax = (doubles_filtered %>% arrange(desc(wetlab_id)) %>% select(HPO_terms)))
dataframe_test <- dataframe_test %>% mutate(HPO_terms.1 = ifelse(HPO_terms.1 == "N/A", NA, dataframe_test$HPO_terms.1))
dataframe_test %>% filter(is.na(HPO_terms)) %>% filter(!is.na(HPO_terms.1))
dataframe_test %>% filter(is.na(HPO_terms.1)) %>% filter(!is.na(HPO_terms))

#~~~~~~Filter doubles_filtered_unique~~~~~~#
RD320 <- doubles_filtered_unique %>% filter(wetlab_id == "RD320")
RD320_1 <- RD320[1,] %>% as.list
RD320_2 <- RD320[2,] %>% as.list
setdiff(RD320_1, RD320_2)

RD320 %>% select(batch, resolved_case, causal_gene, seq_batch, number)
RD320 <- RD320[1,]

RD338 <- doubles_filtered_unique %>% filter(wetlab_id == "RD338")
RD338_1 <- RD338[1,] %>% as.list
RD338_2 <- RD338[2,] %>% as.list
setdiff(RD338_1, RD338_2)

RD338 %>% select(read_length, variant_data, seq_batch, family_id, number)
RD338 <- RD338[1,]

#~~~~~~Join doubles together~~~~~~#
doubles <- bind_rows(RD338, RD320, doubles_filtered)

#~~~~~~Join all together~~~~~~#
joined_missing <- bind_rows(doubles, joined_missing_wetlab_id_singles, joined_missing_wetlab_id_na)

######################################################
###########----Join Metadata and FRASER----###########
#################----on Sample_ID----#################
######################################################

##############----wetlab_id x Sample_ID----###########
wetlab_id_na <- joined_missing %>% filter(is.na(wetlab_id)) %>% filter(wetlab_id == "N/A")
wetlab_id_not_na <- joined_missing %>% filter(!is.na(wetlab_id)) %>% filter(wetlab_id != "N/A")

names(FRASER_result_combined)[names(FRASER_result_combined) == 'Sample_ID'] <- 'wetlab_id'
wet_lab_matched_list <- intersect(FRASER_result_combined$wetlab_id, wetlab_id_not_na$wetlab_id)
wet_lab_matched <- filter(FRASER_result_combined, is_in(wetlab_id, wet_lab_matched_list))
#length 263

wet_lab_matched_joined <- left_join(wet_lab_matched, joined_missing, by = "wetlab_id")
#263

wet_lab_matched_joined %>% group_by(wetlab_id) %>% tally() %>% arrange(desc(n))

wet_lab_left_over_list <- setdiff(FRASER_result_combined$wetlab_id, wetlab_id_not_na$wetlab_id)
wet_lab_left_over <- filter(FRASER_result_combined, is_in(wetlab_id, wet_lab_left_over_list))
#length 159

#length FRASER_result_combined = 422
#159+263=422

#############----sample_id x Sample_ID----############
sample_id_na <- joined_missing %>% filter(is.na(sample_id)) %>% filter(sample_id == "N/A")
sample_id_not_na <- joined_missing %>% filter(!is.na(sample_id)) %>% filter(sample_id != "N/A")

names(wet_lab_left_over)[names(wet_lab_left_over) == 'wetlab_id'] <- 'sample_id'
sample_id_matched_list <- intersect(wet_lab_left_over$sample_id, sample_id_not_na$sample_id)

sample_id_matched <- filter(wet_lab_left_over, is_in(sample_id, sample_id_matched_list))
#length 24

sample_id_matched_joined <- left_join(sample_id_matched, joined_missing, by = "sample_id")

sample_id_left_over_list <- setdiff(wet_lab_left_over$sample_id, wetlab_id_not_na$sample_id)
sample_id_left_over <- filter(wet_lab_left_over, is_in(sample_id, sample_id_left_over_list))
#length 135

#263+24+135=422

###########----institution_id x Sample_ID----#########
institution_id_na <- joined_missing %>% filter(is.na(institution_id)) %>% filter(institution_id == "N/A")
institution_id_not_na <- joined_missing %>% filter(!is.na(institution_id)) %>% filter(institution_id != "N/A")

names(sample_id_left_over)[names(sample_id_left_over) == 'sample_id'] <- 'institution_id'
institution_id_matched_list <- intersect(sample_id_left_over$institution_id, institution_id_not_na$institution_id)
institution_id_matched <- filter(sample_id_left_over, is_in(institution_id, institution_id_matched_list))
#length 0

institution_id_left_over_list <- setdiff(sample_id_left_over$institution_id, institution_id_not_na$institution_id)
institution_id_left_over <- filter(sample_id_left_over, is_in(institution_id, institution_id_left_over_list))
#length 135

#263+24+135+0=422

###############----indv_id x Sample_ID----############
indv_id_na <- joined_missing %>% filter(is.na(indv_id)) %>% filter(indv_id == "N/A")
indv_id_not_na <- joined_missing %>% filter(!is.na(indv_id)) %>% filter(indv_id != "N/A")

names(institution_id_left_over)[names(institution_id_left_over) == 'institution_id'] <- 'indv_id'
indv_id_matched_list <- intersect(institution_id_left_over$indv_id, indv_id_not_na$indv_id)
indv_id_Ind_matched <- filter(institution_id_left_over, is_in(indv_id, indv_id_matched_list))
#length 0

indv_id_left_over_list <- setdiff(institution_id_left_over$indv_id, indv_id_not_na$indv_id)
indv_id_left_over <- filter(institution_id_left_over, is_in(indv_id, indv_id_left_over_list))
#length 135

#263+24+135+0+0=422

##########----participant_id x Sample_ID----##########
participant_id_na <- joined_missing %>% filter(is.na(participant_id)) %>% filter(participant_id == "N/A")
participant_id_not_na <- joined_missing %>% filter(!is.na(participant_id)) %>% filter(participant_id != "N/A")

names(indv_id_left_over)[names(indv_id_left_over) == 'indv_id'] <- 'participant_id'
participant_id_matched_list <- intersect(indv_id_left_over$participant_id, participant_id_not_na$participant_id)
participant_id_matched <- filter(indv_id_left_over, is_in(participant_id, participant_id_matched_list))
#length 58

participant_id_matched_joined <- left_join(participant_id_matched, joined_missing, by = "participant_id")

participant_id_left_over_list <- setdiff(indv_id_left_over$participant_id, participant_id_not_na$participant_id)
participant_id_left_over <- filter(indv_id_left_over, is_in(participant_id, participant_id_left_over_list))
#length 77

#263+24+77+58=422

################----ID x Sample_ID----################
ID_na <- joined_missing %>% filter(is.na(ID)) %>% filter(ID == "N/A")
ID_not_na <- joined_missing %>% filter(!is.na(ID)) %>% filter(ID != "N/A")

names(participant_id_left_over)[names(participant_id_left_over) == 'participant_id'] <- 'ID'
ID_matched_list <- intersect(participant_id_left_over$ID, ID_not_na$ID)
ID_matched <- filter(participant_id_left_over, is_in(ID, ID_matched_list))
#length 70

ID_matched_joined <- left_join(ID_matched, joined_missing, by = "ID")

ID_left_over_list <- setdiff(participant_id_left_over$ID, ID_not_na$ID)
ID_left_over <- filter(participant_id_left_over, is_in(ID, ID_left_over_list))
#length 7

#263+24+58+7+70=422

################----fixed_ID x Sample_ID----################
#ID_left_over$ID <- gsub(".Aligned.sortedByCoord.out.bam", "", ID_left_over$ID)
names(missing_Miami)[names(missing_Miami) == 'sampleID'] <- 'ID'
missing_Miami$ID <- paste(missing_Miami$ID, ".Aligned.sortedByCoord.out.bam", sep="") 

ID_missing_matched_joined <- left_join(ID_left_over, missing_Miami)
ID_missing_matched_joined$batch <- 67
ID_missing_matched_joined$seq_batch <- 67

ID_colnames <- ID_matched_joined %>% select(-RIN, -sex, -age, -affected_status, -batch, -seq_batch) %>% colnames
ID_missing_matched_joined <- merge(ID_missing_matched_joined, ID_matched_joined[,ID_colnames], by="ID",all.x=T)
ID_missing_matched_joined$age <- as.character(ID_missing_matched_joined$age)
#length 7

#255+24+70+58+7=414
##################----joined----###############
ID_matched_joined %>% nrow #70
participant_id_matched_joined %>% nrow #58
participant_id_matched_joined$ID <- participant_id_matched_joined$participant_id
sample_id_matched_joined %>% nrow #24
sample_id_matched_joined$ID <- sample_id_matched_joined$sample_id
wet_lab_matched_joined %>% nrow #255
wet_lab_matched_joined$ID <- wet_lab_matched_joined$wetlab_id
ID_missing_matched_joined %>% nrow #7

FRASER_metadata_joined <- bind_rows(ID_matched_joined, participant_id_matched_joined) %>%
              bind_rows(sample_id_matched_joined) %>% bind_rows(wet_lab_matched_joined) %>% bind_rows(ID_missing_matched_joined)
#n=422
FRASER_metadata_joined_missing <- data.frame(participant_id= c("GSS189457", "GSS153634"))
FRASER_colnames <- FRASER_metadata_joined %>% colnames
FRASER_missing <- merge(FRASER_metadata_joined_missing, FRASER_metadata_joined[,FRASER_colnames], by="participant_id",all.x=T)

FRASER_metadata_joined <- bind_rows(FRASER_metadata_joined, FRASER_missing)
#n=424

FRASER_metadata_joined$sample_ID <- FRASER_metadata_joined$ID

#######################################################
############----Combine GSS and UDN IDs----###########
######################################################
#indv_id (GSS and UDN), participant_id (GSS), wetlab_id (RDID), sample_id (RDID), institution_id (UDN_ID)
colnames(GSS_UDN_MAP) <- c("GSS", "UDN")
intersect(GSS_UDN_MAP$UDN, FRASER_metadata_joined$indv_id)
intersect(GSS_UDN_MAP$UDN, FRASER_metadata_joined$institution_id)
intersect(GSS_UDN_MAP$UDN, FRASER_metadata_joined$ID)

intersect(GSS_UDN_MAP$GSS, FRASER_metadata_joined$indv_id)
#character(0)
intersect(GSS_UDN_MAP$GSS, FRASER_metadata_joined$participant_id)
#character(0)

names(GSS_UDN_MAP)[names(GSS_UDN_MAP) == 'UDN'] <- 'indv_id'
FRASER_metadata1 <- left_join(FRASER_metadata_joined, GSS_UDN_MAP)

names(GSS_UDN_MAP)[names(GSS_UDN_MAP) == 'indv_id'] <- 'ID'
FRASER_metadata <- left_join(FRASER_metadata1, GSS_UDN_MAP)

######################################################
##########----Clean Metadata FRASER file----##########
######################################################
FRASER_metadata_remove_allNA_cols <- FRASER_metadata %>% select(-institution, -in_freeze, -is_RD, -case_id, -status,
                                        -genome_vcf, -exome_vcf_on_durga, -institution, -partial_contribution_explained)

FRASER_metadata_important <- FRASER_metadata_remove_allNA_cols %>% select(-country_of_origin,
                                        -is_duplicated, -sample_type, -source_of_RNA, -read_length, -variant_data,
                                        -age_at_last_observation, -ancestry_detail, -gregor_center, -missing_variant_case,
                                        -prior_testing, -reported_ethnicity, -reported_race, -phenotype_id, -ontology)

################----IDENTIFIERS----###################
#indv_id (UDN), participant_id (GSS), wetlab_id (RDID), sample_id (RDID), institution_id (UDN), GSS, ID (UDN)
FRASER_metadata_identifiers <- FRASER_metadata_important %>% select(indv_id, participant_id, wetlab_id, sample_ID, institution_id, GSS, ID)

FRASER_metadata_identifiers$numbers <- rownames(FRASER_metadata_identifiers)

#~~~~~~ID x instit_id = UDN_ID~~~~~~#
intersect(FRASER_metadata_identifiers$sample_ID, FRASER_metadata_identifiers$institution_id)

FRASER_metadata_identifiers_ID <- FRASER_metadata_identifiers %>% filter(grepl("UDN", ID))
FRASER_metadata_identifiers_institution_id <- FRASER_metadata_identifiers %>% filter(!is.na(institution_id))
FRASER_metadata_identifiers_NA <- FRASER_metadata_identifiers %>% filter(is.na(institution_id)) %>% filter(! grepl("UDN", ID))

FRASER_metadata_identifiers_ID$UDN_ID <- FRASER_metadata_identifiers_ID$ID
FRASER_metadata_identifiers_institution_id$UDN_ID <- FRASER_metadata_identifiers_institution_id$institution_id
FRASER_metadata_identifiers_NA$UDN_ID <- "NA"

FRASER_metadata_identifiers <- bind_rows(FRASER_metadata_identifiers_ID, FRASER_metadata_identifiers_institution_id, FRASER_metadata_identifiers_NA) %>%
                                select(-ID, -institution_id)
FRASER_metadata_identifiers$UDN_ID<-gsub("-","",as.character(FRASER_metadata_identifiers$UDN_ID))
FRASER_metadata_identifiers$UDN_ID<-gsub(" ","",as.character(FRASER_metadata_identifiers$UDN_ID))
FRASER_metadata_identifiers$indv_id<-gsub(" ","",as.character(FRASER_metadata_identifiers$indv_id))
#424

#~~~~~~UDN_ID x indv_id = UDN_ID~~~~~~#
inUDN_not_indv <- setdiff(FRASER_metadata_identifiers$UDN_ID, FRASER_metadata_identifiers$indv_id)
setdiff(FRASER_metadata_identifiers$indv_id, FRASER_metadata_identifiers$UDN_ID)
#character(0)

filter(FRASER_metadata_identifiers, is_in(UDN_ID, inUDN_not_indv)) %>% print(n=100)
FRASER_metadata_identifiers <- FRASER_metadata_identifiers %>% select(-indv_id)

#~~~~~~participant_id x GSS = UDN_ID~~~~~~#
FRASER_metadata_identifiers %>% select(participant_id, GSS) %>% print(n=450)

FRASER_metadata_identifiers_GSS <- FRASER_metadata_identifiers %>% filter(!is.na(GSS))
FRASER_metadata_identifiers_participant_id <- FRASER_metadata_identifiers %>% filter(!is.na(participant_id))
FRASER_metadata_identifiers_NA_GSS <- FRASER_metadata_identifiers %>% filter(is.na(participant_id)) %>% filter(is.na(GSS))

FRASER_metadata_identifiers_GSS$GSS_ID <- FRASER_metadata_identifiers_GSS$GSS
FRASER_metadata_identifiers_participant_id$GSS_ID <- FRASER_metadata_identifiers_participant_id$participant_id
FRASER_metadata_identifiers_NA_GSS$GSS_ID <- NA

FRASER_metadata_identifiers <- bind_rows(FRASER_metadata_identifiers_GSS, FRASER_metadata_identifiers_participant_id, FRASER_metadata_identifiers_NA_GSS) %>%
                                select(-GSS, -participant_id)

#~~~~~~sample_id = RDID~~~~~~#
#names(FRASER_metadata_identifiers)[names(FRASER_metadata_identifiers) == 'sample_id'] <- 'RDID'
FRASER_metadata_identifiers_wetlab_RDID <- FRASER_metadata_identifiers %>% filter(grepl("RD", FRASER_metadata_identifiers$wetlab_id))
FRASER_metadata_identifiers_wetlab_non_RDID <- FRASER_metadata_identifiers %>% filter(! grepl("RD", FRASER_metadata_identifiers$wetlab_id))

FRASER_metadata_identifiers_wetlab_RDID$RDID <- FRASER_metadata_identifiers_wetlab_RDID$wetlab_id
FRASER_metadata_identifiers_wetlab_non_RDID$RDID <- NA 

FRASER_metadata_identifiers <- bind_rows(FRASER_metadata_identifiers_wetlab_RDID, FRASER_metadata_identifiers_wetlab_non_RDID)

FRASER_metadata_identifiers_no_RDID <- FRASER_metadata_identifiers %>% filter(is.na(RDID))
FRASER_metadata_identifiers_RDID <- FRASER_metadata_identifiers %>% filter(!is.na(RDID))

FRASER_metadata_identifiers_no_RDID_sample_ID_RDID <- FRASER_metadata_identifiers_no_RDID %>% filter(grepl("RD", FRASER_metadata_identifiers_no_RDID$sample_ID))
FRASER_metadata_identifiers_no_RDID_sample_ID_non_RDID <- FRASER_metadata_identifiers_no_RDID %>% filter(! grepl("RD", FRASER_metadata_identifiers_no_RDID$sample_ID))

FRASER_metadata_identifiers_no_RDID_sample_ID_RDID$RDID <- FRASER_metadata_identifiers_no_RDID_sample_ID_RDID$sample_ID
FRASER_metadata_identifiers_no_RDID_sample_ID_non_RDID$RDID <- NA 

FRASER_metadata_identifiers_no_RDID <- bind_rows(FRASER_metadata_identifiers_no_RDID_sample_ID_RDID, FRASER_metadata_identifiers_no_RDID_sample_ID_non_RDID)

FRASER_metadata_identifiers <- bind_rows(FRASER_metadata_identifiers_RDID, FRASER_metadata_identifiers_no_RDID)

#~~~~~~rejoin dataframes~~~~~~#
FRASER_metadata_important$numbers <- rownames(FRASER_metadata_important)
FRASER_metadata_important <- FRASER_metadata_important %>% left_join(FRASER_metadata_identifiers) %>% 
                                select(-indv_id, -participant_id, -sample_id, -institution_id, -GSS, -ID)

###############----CASE v CONTROL----#################
#proband, affected_status, affected_status_GREGoR
FRASER_metadata_important <- FRASER_metadata_important %>% mutate(affected_status_GREGoR = ifelse(affected_status_GREGoR == "Unaffected", "Control",
                                                                                                ifelse(affected_status_GREGoR == "Affected", "Case", NA)))

FRASER_metadata_important <- FRASER_metadata_important %>% mutate(proband = ifelse(proband == 1, "Case",
                                                                                                ifelse(proband == 0, "Control", NA)))
FRASER_metadata_important$case_control <- paste(FRASER_metadata_important$proband, FRASER_metadata_important$affected_status, FRASER_metadata_important$affected_status_GREGoR, sep="")

FRASER_metadata_important <- FRASER_metadata_important %>% mutate(case_control = ifelse(case_control == "CaseCaseN", "Case",
                                                                                                ifelse(case_control == "NAControlNA", "Control", 
                                                                                                ifelse(case_control == "ControlControlNA", "Control",
                                                                                                ifelse(case_control == "NACaseNA", "Case",
                                                                                                ifelse(case_control == "NANANA", NA,
                                                                                                ifelse(case_control == "NANACase", "Case",
                                                                                                ifelse(case_control == "NANAControl", "Control", NA))))))))
FRASER_metadata_important$affected_status <- FRASER_metadata_important$case_control

FRASER_metadata_important <- FRASER_metadata_important %>% select(-proband, -case_control, -affected_status_GREGoR)

missing_affected_status <- FRASER_metadata_important %>% filter(is.na(affected_status))
not_missing_affected_status <- FRASER_metadata_important %>% filter(!is.na(affected_status))

not_na <- missing_affected_status %>% filter(!is.na(clinical_description))
na <- missing_affected_status %>% filter(is.na(clinical_description))

na$phenotype_description
na$term_id

GSS_participant %>% filter(participant_id %in% na$GSS_ID)

na[1, 3] <- "Unknown"
na[2, 3] <- "Case"
na[3, 3] <- "Unknown"
na[4, 3] <- "Case"
na[5, 3] <- "Case"

not_na$affected_status <- "Case"

FRASER_metadata_important <- bind_rows(not_missing_affected_status, not_na, na)
FRASER_metadata_important %>% filter(solve_status == "Unaffected") %>% filter(affected_status == "Case") %>% select(number)

####################----SEX----#######################
#sex, sex_GREGoR, inferred_sex, sex_detail
#~~~~~~sex x sex_GREGoR = RDID~~~~~~#
FRASER_metadata_important_GREGoR <- FRASER_metadata_important %>% filter(!is.na(sex_GREGoR))
FRASER_metadata_important_GREGoR <- FRASER_metadata_important_GREGoR %>% mutate(sex_GREGoR = ifelse(sex_GREGoR == "Female", "F", "M"))
FRASER_metadata_important_GREGoR$sex <- FRASER_metadata_important_GREGoR$sex_GREGoR
FRASER_metadata_important_GREGoR <- FRASER_metadata_important_GREGoR %>% select(-sex_GREGoR)

FRASER_metadata_important_UDN <- FRASER_metadata_important %>% filter(is.na(sex_GREGoR)) %>% select(-sex_GREGoR)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_UDN, FRASER_metadata_important_GREGoR)

#~~~~~~inferred_sex & sex_detail = RDID~~~~~~#
FRASER_metadata_important <- FRASER_metadata_important %>% select(-inferred_sex, -sex_detail)

####################----AGE----#######################
#onset_age_range, age, age_at_enrollment
#~~~~~~onset_age_range~~~~~~#
FRASER_metadata_important <- FRASER_metadata_important %>% select(-onset_age_range)
#no meaningful info (as far as I can tell)

#~~~~~~age & age_at_enrollment = RDID~~~~~~#
FRASER_metadata_important_age <- FRASER_metadata_important %>% filter(!is.na(age))
FRASER_metadata_important_enrollment <- FRASER_metadata_important %>% filter(!is.na(age_at_enrollment))
FRASER_metadata_important_NA_age <- FRASER_metadata_important %>% filter(is.na(age)) %>% filter(is.na(age_at_enrollment))

FRASER_metadata_important_age$age_of <- FRASER_metadata_important_age$age
FRASER_metadata_important_enrollment$age_of <- FRASER_metadata_important_enrollment$age_at_enrollment
FRASER_metadata_important_NA_age$age_of <- NA

FRASER_metadata_important_enrollment$age_of <- as.character(FRASER_metadata_important_enrollment$age_of)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_age, FRASER_metadata_important_enrollment, FRASER_metadata_important_NA_age) %>%
                                select(-age, -age_at_enrollment)

names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'age_of'] <- 'age'

no_age <- FRASER_metadata_important %>% filter(is.na(age))
#get age for these patients

################----CLINICAL DS----###################
#clinical_description, phenotype_description, additional_details, 
#known_condition_name, gene_known_for_phenotype
#~~~~~~clinical_description & phenotype_description = patient_description~~~~~~#
FRASER_metadata_important_clinical_description <- FRASER_metadata_important %>% filter(!is.na(clinical_description))
FRASER_metadata_important_phenotype_description <- FRASER_metadata_important %>% filter(!is.na(phenotype_description))
FRASER_metadata_important_NA_clin <- FRASER_metadata_important %>% filter(is.na(phenotype_description)) %>% filter(is.na(clinical_description))

FRASER_metadata_important_clinical_description$patient_description <- FRASER_metadata_important_clinical_description$clinical_description
FRASER_metadata_important_phenotype_description$patient_description <- FRASER_metadata_important_phenotype_description$phenotype_description
FRASER_metadata_important_NA_clin$patient_description <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_clinical_description, FRASER_metadata_important_phenotype_description, FRASER_metadata_important_NA_clin) %>%
                                select(-clinical_description, -phenotype_description)

#~~~~~~add additional_details to patient_description~~~~~~#
FRASER_metadata_important$additional_details <-gsub("NA","",as.character(FRASER_metadata_important$additional_details))
FRASER_metadata_important$additional_details <-gsub(";"," ",as.character(FRASER_metadata_important$additional_details))
FRASER_metadata_important$additional_details <-gsub("    ","",as.character(FRASER_metadata_important$additional_details))
FRASER_metadata_important$additional_details <-gsub("   ","",as.character(FRASER_metadata_important$additional_details))
FRASER_metadata_important$additional_details <-gsub("  ","",as.character(FRASER_metadata_important$additional_details))
FRASER_metadata_important$additional_details <-trimws(FRASER_metadata_important$additional_details)

nav <- c('', ' ', "history of")
FRASER_metadata_important <- transform(FRASER_metadata_important, additional_details=replace(additional_details, additional_details %in% nav, NA))

FRASER_metadata_important_add <- FRASER_metadata_important %>% filter(!is.na(additional_details)) %>% filter(is.na(patient_description))
FRASER_metadata_important_patient_description <- FRASER_metadata_important %>% filter(!is.na(patient_description)) %>% filter(is.na(additional_details))
FRASER_metadata_important_NA_add <- FRASER_metadata_important %>% filter(is.na(patient_description)) %>% filter(is.na(additional_details))
FRASER_metadata_important_NA_no_add <- FRASER_metadata_important %>% filter(!is.na(patient_description)) %>% filter(!is.na(additional_details))

FRASER_metadata_important_add$patient_description <- FRASER_metadata_important_add$additional_details
FRASER_metadata_important_patient_description$patient_description <- FRASER_metadata_important_patient_description$patient_description
FRASER_metadata_important_NA_add$patient_description <- NA
FRASER_metadata_important_NA_no_add$patient_description <- paste(FRASER_metadata_important_NA_no_add$patient_description, "ADDITIONAL DETAILS: ", FRASER_metadata_important_NA_no_add$additional_details)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_add, FRASER_metadata_important_patient_description, FRASER_metadata_important_NA_add, FRASER_metadata_important_NA_no_add) %>%
                                select(-additional_details)

#~~~~~~add known_condition_name to patient_description~~~~~~#
FRASER_metadata_important_patient_description_cond <- FRASER_metadata_important %>% filter(!is.na(patient_description)) %>% filter(is.na(known_condition_name))
FRASER_metadata_important_NA_cond <- FRASER_metadata_important %>% filter(is.na(patient_description)) %>% filter(is.na(known_condition_name))
FRASER_metadata_important_NA_no_cond <- FRASER_metadata_important %>% filter(!is.na(patient_description)) %>% filter(!is.na(known_condition_name))

FRASER_metadata_important_patient_description_cond$patient_description <- FRASER_metadata_important_patient_description_cond$patient_description
FRASER_metadata_important_NA_cond$patient_description <- NA
FRASER_metadata_important_NA_no_cond$patient_description <- paste(FRASER_metadata_important_NA_no_cond$patient_description, " KNOWN CONDITION NAME: ", FRASER_metadata_important_NA_no_cond$additional_details)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_patient_description_cond, FRASER_metadata_important_NA_cond, FRASER_metadata_important_NA_no_cond) %>%
                                select(-known_condition_name)

#################----HPO TERMS----####################
#HPO_terms, HPO_terms_ID, presence, term_id, additional_modifiers
#~~~~~~delete presence~~~~~~#
unknown <- FRASER_metadata_important %>% filter(presence == "Present;Present;Present;Present;Present;Present;Present;Present;Unknown")
unknown$notes <- "unknown if has HP:0032046"
unknown$term_id <- "HP:0002069;HP:0001344;HP:0001250;HP:0000729;HP:0001249;HP:0000722;HP:0001051;HP:0000545"

absent <- FRASER_metadata_important %>% filter(presence == "Present;Present;Present;Present;Absent;Present")
absent$notes <- "absent HP:0000159"
absent$term_id <- "Unicoronalsynostosis;HP:0000587;Unicoronal synostosis;HP:0001347;HP:0004442"

not_present_absent <- FRASER_metadata_important %>% filter(!is_in(numbers, c(unknown$numbers, absent$numbers)))

FRASER_metadata_important <- bind_rows(unknown, absent, not_present_absent) %>% select(-presence)

#~~~~~~add HPO_terms_ID x term_id = HPO_IDs~~~~~~#
FRASER_metadata_important_term <- FRASER_metadata_important %>% filter(!is.na(term_id)) %>% filter(is.na(HPO_terms_ID))
FRASER_metadata_important_hpo <- FRASER_metadata_important %>% filter(!is.na(HPO_terms_ID)) %>% filter(is.na(term_id))
FRASER_metadata_important_NA_term <- FRASER_metadata_important %>% filter(is.na(HPO_terms_ID)) %>% filter(is.na(term_id))

FRASER_metadata_important_term$HPO_IDs <- FRASER_metadata_important_term$term_id
FRASER_metadata_important_hpo$HPO_IDs <- FRASER_metadata_important_hpo$HPO_terms_ID
FRASER_metadata_important_NA_term$HPO_IDs <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_term, FRASER_metadata_important_hpo, FRASER_metadata_important_NA_term) %>%
                                select(-HPO_terms_ID, -term_id)

#~~~~~~add additional_modifiers~~~~~~#
FRASER_metadata_important$additional_modifiers <-gsub("NA;","",as.character(FRASER_metadata_important$additional_modifiers))
FRASER_metadata_important$additional_modifiers <-gsub(";NA"," ",as.character(FRASER_metadata_important$additional_modifiers))
FRASER_metadata_important$additional_modifiers <-trimws(FRASER_metadata_important$additional_modifiers)

FRASER_metadata_important <- FRASER_metadata_important %>% mutate(additional_modifiers = ifelse(additional_modifiers == "NA", NA, FRASER_metadata_important$additional_modifiers))

FRASER_metadata_important_HPO <- FRASER_metadata_important %>% filter(!is.na(HPO_IDs)) %>% filter(is.na(additional_modifiers))
FRASER_metadata_important_NA_mod <- FRASER_metadata_important %>% filter(is.na(HPO_IDs)) %>% filter(is.na(additional_modifiers))
FRASER_metadata_important_NA_no_mod <- FRASER_metadata_important %>% filter(!is.na(HPO_IDs)) %>% filter(!is.na(additional_modifiers))

FRASER_metadata_important_HPO$HPO_IDs <- FRASER_metadata_important_HPO$HPO_IDs
FRASER_metadata_important_NA_mod$HPO_IDs <- NA
FRASER_metadata_important_NA_no_mod$HPO_IDs <- paste(FRASER_metadata_important_NA_no_mod$HPO_IDs, "ADDITIONAL MODIFIERS: ", FRASER_metadata_important_NA_no_mod$additional_modifiers)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_HPO, FRASER_metadata_important_NA_mod, FRASER_metadata_important_NA_no_mod) %>%
                                select(-additional_modifiers)

#################----HPO TERMS----####################
#candidate_gene_list, candidate_variant_list\nCoordinates (hg38): | Gene: |Isoform: |c. |p. | Zyg: | Inh:
#causal_gene, gene, candidate_variant_list

#~~~~~~candidate_gene_list x gene = candidate_genes~~~~~~#
FRASER_metadata_important$candidate_gene_list %>% unique
FRASER_metadata_important %>% filter(!is.na(gene)) %>% select(RDID, wetlab_id, UDN_ID, GSS_ID, resolved_case, solve_status, gene)

FRASER_metadata_important_gene <- FRASER_metadata_important %>% filter(!is.na(gene)) %>% filter(is.na(candidate_gene_list))
FRASER_metadata_important_cand <- FRASER_metadata_important %>% filter(!is.na(candidate_gene_list)) %>% filter(is.na(gene))
FRASER_metadata_important_NA_gene <- FRASER_metadata_important %>% filter(is.na(gene)) %>% filter(is.na(candidate_gene_list))

FRASER_metadata_important_gene$candidate_gene <- FRASER_metadata_important_gene$gene
FRASER_metadata_important_cand$candidate_gene <- FRASER_metadata_important_cand$candidate_gene_list
FRASER_metadata_important_NA_gene$candidate_gene <- NA

FRASER_metadata_important %>% filter(!is.na(gene)) %>% select(RDID, wetlab_id, UDN_ID, GSS_ID, resolved_case, solve_status, gene)
FRASER_metadata_important <- bind_rows(FRASER_metadata_important_gene, FRASER_metadata_important_cand, FRASER_metadata_important_NA_gene) %>%
                                select(-gene, -candidate_gene_list)

#~~~~~~candidate_variant_list.Coordinates (hg38): | Gene: |Isoform: |c. |p. | Zyg: | Inh: x candidate_variant_list = candidate_variants~~~~~~#
FRASER_metadata_important$candidate_gene_list %>% unique

FRASER_metadata_important_long <- FRASER_metadata_important %>% filter(!is.na(`candidate_variant_list.Coordinates..hg38.....Gene...Isoform...c...p....Zyg....Inh.`)) %>% filter(is.na(candidate_variant_list))
FRASER_metadata_important_var <- FRASER_metadata_important %>% filter(!is.na(candidate_variant_list)) %>% filter(is.na(`candidate_variant_list.Coordinates..hg38.....Gene...Isoform...c...p....Zyg....Inh.`))
FRASER_metadata_important_NA_long <- FRASER_metadata_important %>% filter(is.na(candidate_variant_list)) %>% filter(is.na(`candidate_variant_list.Coordinates..hg38.....Gene...Isoform...c...p....Zyg....Inh.`))

FRASER_metadata_important_long$candidate_variants <- FRASER_metadata_important_long$`candidate_variant_list.Coordinates..hg38.....Gene...Isoform...c...p....Zyg....Inh.`
FRASER_metadata_important_var$candidate_variants <- FRASER_metadata_important_var$candidate_variant_list
FRASER_metadata_important_NA_long$candidate_variants <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_long, FRASER_metadata_important_var, FRASER_metadata_important_NA_long) %>%
                                select(-"candidate_variant_list.Coordinates..hg38.....Gene...Isoform...c...p....Zyg....Inh.", -candidate_variant_list)

#~~~~~~fix causal_gene = UDN_causal_gene~~~~~~#
FRASER_metadata_important <- FRASER_metadata_important %>% mutate(causal_gene = ifelse(causal_gene == "#N/A", NA, FRASER_metadata_important$causal_gene))
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'causal_gene'] <- 'UDN_causal_gene'

##############----SOLVED STATUS----###################
#resolved_case, solve_status

#~~~~~~resolved_case x solve_status = solved_status~~~~~~#
FRASER_metadata_important <- FRASER_metadata_important %>% mutate(resolved_case = ifelse(resolved_case == "#N/A", NA, 
                                                                    ifelse(resolved_case == "No", "Unsolved",
                                                                    ifelse(resolved_case == "no", "Unsolved",
                                                                    ifelse(resolved_case == "Yes", "Solved",
                                                                    ifelse(resolved_case == "yes", "Solved",
                                                                    ifelse(resolved_case == "N/A", NA, "Other")))))))

FRASER_metadata_important_res <- FRASER_metadata_important %>% filter(!is.na(resolved_case)) %>% filter(is.na(solve_status))
FRASER_metadata_important_solve <- FRASER_metadata_important %>% filter(!is.na(solve_status)) %>% filter(is.na(resolved_case))
FRASER_metadata_important_NA_res <- FRASER_metadata_important %>% filter(is.na(solve_status)) %>% filter(is.na(resolved_case))

FRASER_metadata_important_res$solved_status <- FRASER_metadata_important_res$resolved_case
FRASER_metadata_important_solve$solved_status <- FRASER_metadata_important_solve$solve_status
FRASER_metadata_important_NA_res$solved_status <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_res, FRASER_metadata_important_solve, FRASER_metadata_important_NA_res) %>%
                                select(-resolved_case, -solve_status)

unclear_status <- filter(FRASER_metadata_important, is.na(solved_status)) %>% select(GSS_ID, UDN_ID, RDID, wetlab_id)

############----FAMILY IDENTIFIER----#################
#family_id_GREGoR, family_id, relationship, proband_relationship_detail
#maternal_id, paternal_id, twin_id, family, family_type, relationship, related_proband
#~~~~~~family_id_GREGoR x family_id = family_identifier~~~~~~#
FRASER_metadata_important$family_id %>% unique
FRASER_metadata_important$family_id_GREGoR %>% unique
intersect(FRASER_metadata_important$family_id_GREGoR, FRASER_metadata_important$GSS_ID)

FRASER_metadata_important_fam <- FRASER_metadata_important %>% filter(!is.na(family_id)) %>% filter(is.na(family_id_GREGoR))
FRASER_metadata_important_GREG <- FRASER_metadata_important %>% filter(!is.na(family_id_GREGoR)) %>% filter(is.na(family_id))
FRASER_metadata_important_no_NA_fam <- FRASER_metadata_important %>% filter(is.na(family_id)) %>% filter(is.na(family_id_GREGoR))

FRASER_metadata_important_fam$family_identifier <- as.character(FRASER_metadata_important_fam$family_id)
FRASER_metadata_important_GREG$family_identifier <- FRASER_metadata_important_GREG$family_id_GREGoR
FRASER_metadata_important_no_NA_fam$family_identifier <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_fam, FRASER_metadata_important_GREG, FRASER_metadata_important_no_NA_fam) %>%
                                select(-family_id, -family_id_GREGoR)

no_family_id <- FRASER_metadata_important %>% filter(is.na(family_identifier)) %>% select(UDN_ID, GSS_ID, RDID, wetlab_id, affected_status)

#~~~~~~relationship x proband_relationship_detail = relationship_details~~~~~~#
FRASER_metadata_important <- FRASER_metadata_important %>% mutate(relationship = ifelse(relationship == "sister", "Sister", 
                                                                    ifelse(relationship == "brother", "Brother",
                                                                    ifelse(relationship == "mother", "Mother",
                                                                    ifelse(relationship == "father", "Father", FRASER_metadata_important$relationship)))))

#FRASER_metadata_important_half <- FRASER_metadata_important %>% filter(!is.na(proband_relationship_detail))
#FRASER_metadata_important_half$relationship <- "Sister (unclear if half or full)"

#FRASER_metadata_important_not_half <- FRASER_metadata_important %>% filter(is.na(proband_relationship_detail))

#FRASER_metadata_important <- bind_rows(FRASER_metadata_important_half, FRASER_metadata_important_not_half) %>% select(-proband_relationship_detail)
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'relationship'] <- 'relationship_details'

#~~~~~~rename maternal_id, paternal_id, twin_id~~~~~~#
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'maternal_id'] <- 'GREGoR_maternal_id'
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'paternal_id'] <- 'GREGoR_paternal_id'
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'twin_id'] <- 'GREGoR_twin_id'

#~~~~~~rename maternal_id, paternal_id, twin_id~~~~~~#
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'family_type'] <- 'UDN_family_type'
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'related_proband'] <- 'UDN_related_proband'

#~~~~~~family to family_sequenced~~~~~~#
FRASER_metadata_important$family %>% unique
proband_only <- FRASER_metadata_important %>% filter(UDN_family_type == "Proband")
proband_only$family_sequenced <- "No"

NA_family <- FRASER_metadata_important %>% filter(is.na(UDN_family_type))
Other <- FRASER_metadata_important %>% filter(is_in(UDN_family_type, c("Trio", "Quad")))

not_proband_only <- bind_rows(NA_family, Other)

no_family <- not_proband_only %>% filter(is.na(UDN_related_proband)) %>% filter(is.na(GREGoR_twin_id)) %>% filter(is.na(GREGoR_paternal_id)) %>%
                filter(is.na(GREGoR_maternal_id)) %>% filter(is.na(relationship_details)) %>% filter(is.na(family_identifier))
no_family$family_sequenced <- "No"
no_family_numbers <- select(no_family, numbers) %>% as.list() %>% unlist %>% as.double()

family <- not_proband_only %>% filter(!is_in(numbers, no_family_numbers))
family$family_sequenced <- "Yes"

FRASER_metadata_important <- bind_rows(proband_only, no_family, family) %>% select(-family)

#################----BATCH----########################
#seq_batch, batch, batch_sample_number, seq_batch_2

#~~~~~~remove batch_sample_number~~~~~~#
FRASER_metadata_important$batch_sample_number %>% unique
FRASER_metadata_important <- FRASER_metadata_important %>% select(-batch_sample_number)

#~~~~~~remove seq_batch_2~~~~~~#
missing_batch_information <- FRASER_metadata_important %>% filter(is.na(seq_batch)) %>% filter(is.na(batch))

FRASER_metadata_important %>% group_by(seq_batch) %>% tally
FRASER_metadata_important %>% group_by(batch) %>% tally

#difference between batch and seq_batch?

##################----Misc----########################
#notes
#disease_category
#RIN
#additional_family_members_with_variant
#condition_id -> use to say if solved or not
#condition_inheritance
FRASER_metadata_important <- FRASER_metadata_important %>% select(-origin, -contact)

######################################################
###############----Add in RIN info----################
######################################################
#~~~~~~combine RIN and FRASER dataframe~~~~~~#
colnames(RNA_seq_info) <- c("GSS_ID", "Date_drawn", "Data_received", "Tissue_type", "My_prep_date", "Aliquots", "RNA_RIN", "260_280", "260_230", "QUBIT_RNA_ng_ul",	"uL_VOLUME_USED_FOR_LIB",	
"uL_WATER_ADDED",	"LIBRARY_DATE", "LIBRARY_TECH",	"AVG_BP",	"QUBIT_DNA_ng_ul", "nM_CONC", "PLATE_POSITION", "TECAN_BARCODE", "NOTES")

GSS_overlap <- intersect(FRASER_metadata_important$GSS_ID, RNA_seq_info$GSS_ID)
UDN_overlap <- intersect(FRASER_metadata_important$UDN_ID, RNA_seq_info$GSS_ID)

GSS_overlap_df <- FRASER_metadata_important %>% filter(is_in(GSS_ID, GSS_overlap))
UDN_overlap_df <- FRASER_metadata_important %>% filter(is_in(UDN_ID, UDN_overlap))

GSS_UDN_overlap <- intersect(GSS_overlap_df$numbers, UDN_overlap_df$numbers)

Non_overlap_df <- FRASER_metadata_important %>% filter(!is_in(GSS_ID, GSS_overlap)) %>%
                    filter(!is_in(UDN_ID, UDN_overlap))

RNA_seq_info <- RNA_seq_info %>% select(GSS_ID, RNA_RIN)

GSS_overlap_df <- GSS_overlap_df %>% left_join(RNA_seq_info)

names(RNA_seq_info)[names(RNA_seq_info) == 'GSS_ID'] <- 'UDN_ID'
UDN_overlap_df <- UDN_overlap_df %>% left_join(RNA_seq_info)

Non_overlap_df$RNA_RIN <- NA

FRASER_metadata_important <- bind_rows(Non_overlap_df, UDN_overlap_df, GSS_overlap_df)

#~~~~~~combine RNA_RIN with RIN~~~~~~#
FRASER_metadata_important_RIN <- FRASER_metadata_important %>% filter(!is.na(RIN)) %>% filter(is.na(RNA_RIN))
FRASER_metadata_important_RNA <- FRASER_metadata_important %>% filter(!is.na(RNA_RIN)) %>% filter(is.na(RIN))
FRASER_metadata_important_NA_RIN <- FRASER_metadata_important %>% filter(is.na(RIN)) %>% filter(is.na(RNA_RIN))

FRASER_metadata_important_RIN$RIN <- as.character(FRASER_metadata_important_RIN$RIN)
FRASER_metadata_important_RNA$RIN <- FRASER_metadata_important_RNA$RNA_RIN
FRASER_metadata_important_NA_RIN$RIN <- NA

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_RIN, FRASER_metadata_important_RNA, FRASER_metadata_important_NA_RIN) %>%
                                select(-RNA_RIN)

no_RIN <- FRASER_metadata_important %>% filter(is.na(RIN)) %>% select(GSS_ID, UDN_ID, RDID, wetlab_id)

RD199 <- FRASER_metadata_important %>% filter(RDID == "RD199")
not_RD199 <- FRASER_metadata_important %>% filter(!is_in(RDID, c("RD199")))

RD199$RIN <- "8.9"
FRASER_metadata_important <- bind_rows(RD199, not_RD199)

no_RIN <- FRASER_metadata_important %>% filter(is.na(RIN)) %>% select(GSS_ID, UDN_ID, RDID, wetlab_id)

######################################################
#############----Find Mising Values----###############
######################################################

#############----Age----###############
FRASER_metadata_important_with_age <- FRASER_metadata_important %>% filter(!is.na(age))
FRASER_metadata_important_missing_age <- FRASER_metadata_important %>% filter(is.na(age))
FRASER_metadata_important_missing_age <- FRASER_metadata_important_missing_age %>% select(-age)

FRASER_missing_age_df_GSS <- filter(FRASER_metadata_important_missing_age, grepl("GSS", GSS_ID))
missing_age_df_GSS <- missing_age_df %>% select(GSS_ID, `AGE AT COLLECTION`)
colnames(missing_age_df_GSS) <- c("GSS_ID", "age")
missing_age_df_joined_GSS <- left_join(FRASER_missing_age_df_GSS, missing_age_df_GSS)

FRASER_missing_age_df_RDID <- filter(FRASER_metadata_important_missing_age, grepl("RD", RDID))
missing_age_df_RDID <- missing_age_df %>% select(RDID, `AGE AT COLLECTION`)
colnames(missing_age_df_RDID) <- c("RDID", "age")
missing_age_df_joined_RDID <- left_join(FRASER_missing_age_df_RDID, missing_age_df_RDID)

FRASER_missing_age_df_UDN <- filter(FRASER_metadata_important_missing_age, grepl("UDN", UDN_ID))
missing_age_df_UDN <- missing_age_df %>% select(UDN_ID, `AGE AT COLLECTION`)
colnames(missing_age_df_UDN) <- c("UDN_ID", "age")
missing_age_df_joined_UDN <- left_join(FRASER_missing_age_df_UDN, missing_age_df_UDN)

FRASER_missing_age <- rbind(missing_age_df_joined_GSS, missing_age_df_joined_RDID, missing_age_df_joined_UDN)

FRASER_metadata_important <- rbind(FRASER_metadata_important_with_age, FRASER_missing_age)

FRASER_metadata_important %>% filter(is.na(age))

#############----Batch----###############
FRASER_metadata_important_with_batch <- FRASER_metadata_important %>% filter(!is.na(batch))
FRASER_metadata_important_missing_batch <- FRASER_metadata_important %>% filter(is.na(batch))

GSS_z <- intersect(batchz$Batch, FRASER_metadata_important_missing_batch$GSS_ID)
GSS_y <- intersect(batchy$Batch, FRASER_metadata_important_missing_batch$GSS_ID)
GSS_x <- intersect(batchx$Batch, FRASER_metadata_important_missing_batch$GSS_ID)

UDN_z <- intersect(batchz$Batch, FRASER_metadata_important_missing_batch$UDN_ID)
UDN_y <- intersect(batchy$Batch, FRASER_metadata_important_missing_batch$UDN_ID)
UDN_x <- intersect(batchx$Batch, FRASER_metadata_important_missing_batch$UDN_ID)

batchz_combined <- c(GSS_z, UDN_z)
batchy_combined <- c(GSS_y, UDN_y)

FRASER_metadata_important_missing_batch$batch <- ifelse(FRASER_metadata_important_missing_batch$GSS_ID %in% batchz_combined, 21, FRASER_metadata_important_missing_batch$batch)
FRASER_metadata_important_missing_batch$batch <- ifelse(FRASER_metadata_important_missing_batch$GSS_ID %in% batchy_combined, 22, FRASER_metadata_important_missing_batch$batch)
FRASER_metadata_important_missing_batch$batch <- ifelse(FRASER_metadata_important_missing_batch$UDN_ID %in% batchz_combined, 21, FRASER_metadata_important_missing_batch$batch)
FRASER_metadata_important_missing_batch$batch <- ifelse(FRASER_metadata_important_missing_batch$UDN_ID %in% batchy_combined, 22, FRASER_metadata_important_missing_batch$batch)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_missing_batch, FRASER_metadata_important_with_batch)

FRASER_metadata_important$batch %>% unique

#############----Batch AGAIN----###############
FRASER_metadata_important_with_batch <- FRASER_metadata_important %>% filter(!is.na(batch))
FRASER_metadata_important_missing_batch <- FRASER_metadata_important %>% filter(is.na(batch)) %>% select(-batch)

new_batch_info <- new_batch_info %>% select(RDID, `SEQ RUN #`)
colnames(new_batch_info) <- c("RDID", "batch")
missing_age_df_joined_RDID <- left_join(FRASER_metadata_important_missing_batch, new_batch_info)

FRASER_metadata_important <- rbind(FRASER_metadata_important_with_batch, missing_age_df_joined_RDID)

FRASER_metadata_important_missing_batch <- FRASER_metadata_important %>% filter(UDN_ID %in% c("UDN479881", "UDN969133"))
FRASER_metadata_important_ok_batch <- FRASER_metadata_important %>% filter(!UDN_ID %in% c("UDN479881", "UDN969133"))

FRASER_metadata_important_missing_batch$batch <- NA

FRASER_metadata_important <- rbind(FRASER_metadata_important_missing_batch, FRASER_metadata_important_ok_batch)

FRASER_metadata_important %>% filter(is.na(batch))
FRASER_metadata_important %>% group_by(batch) %>% tally()

#############----Batch AGAIN----###############
batch_info_5 <- batch_info %>% filter(`SEQ RUN #` == 5)
colnames(batch_info_5) <- c("UDN_ID", "RDID", "batch")

batch_5_UDNID <- intersect(batch_info_5$UDN_ID, FRASER_metadata_important$UDN_ID)
batch_5_RDID <- intersect(batch_info_5$RDID, FRASER_metadata_important$RDID)

filter(FRASER_metadata_important, UDN_ID %in% batch_5_UDNID) %>% select(batch)
#batch 5 = 15

batch_info_6 <- batch_info %>% filter(`SEQ RUN #` == 6)
colnames(batch_info_6) <- c("UDN_ID", "RDID", "batch")
batch_6_UDNID <- intersect(batch_info_6$UDN_ID, FRASER_metadata_important$UDN_ID)
batch_6_RDID <- intersect(batch_info_6$RDID, FRASER_metadata_important$RDID)

filter(FRASER_metadata_important, UDN_ID %in% batch_6_UDNID) %>% select(batch)
#batch 6 = 16

FRASER_metadata_important_wrong_batch <- FRASER_metadata_important %>% filter(UDN_ID %in% c("UDN479881", "UDN969133"))
FRASER_metadata_important_ok_batch <- FRASER_metadata_important %>% filter(!UDN_ID %in% c("UDN479881", "UDN969133"))

FRASER_metadata_important_wrong_batch$batch <- c(15, 16)

FRASER_metadata_important <- rbind(FRASER_metadata_important_wrong_batch, FRASER_metadata_important_ok_batch)
FRASER_metadata_important %>% filter(is.na(batch))

#############----Add ID----###############
FRASER_metadata_important$ID <- paste(FRASER_metadata_important$RDID, FRASER_metadata_important$GSS_ID, FRASER_metadata_important$UDN_ID, sep="_")
bad_IDs <- FRASER_metadata_important %>% group_by(sample_ID) %>% tally() %>% arrange(desc(n)) %>% filter(n ==2) %>% pull(sample_ID)

FRASER_metadata_important_bad <- FRASER_metadata_important %>% filter(sample_ID %in% bad_IDs) %>% arrange(sample_ID) 
FRASER_metadata_important_bad_cases <- FRASER_metadata_important_bad %>% filter(affected_status == "Case") %>% filter(!is.na(notes))
FRASER_metadata_important_bad_controls <- FRASER_metadata_important_bad %>% filter(affected_status == "Control") %>% select(-number, -numbers) %>% unique

number<- FRASER_metadata_important_bad %>% filter(affected_status == "Control") %>% select(number) %>% slice(1) %>% pull(number)
numbers<- FRASER_metadata_important_bad %>% filter(affected_status == "Control") %>% select(numbers) %>% slice(1) %>% pull(numbers)

FRASER_metadata_important_bad_controls$numbers <- c(numbers)
FRASER_metadata_important_bad_controls$number <- c(number)

FRASER_metadata_important_bad <- bind_rows(FRASER_metadata_important_bad_cases, FRASER_metadata_important_bad_controls)
FRASER_metadata_important_good <- FRASER_metadata_important %>% filter(! sample_ID %in% bad_IDs)

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_bad, FRASER_metadata_important_good)

names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'sample_ID'] <- 'sampleID'
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'sampleID'] <- 'sample_ID'

write_csv(FRASER_metadata_important, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")

#############----Fix Batch----###############
batch_info_1_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 1) %>% pull(PPID)
one <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_1_UDN_ID)
wrong_UDN_one <- one %>% filter(batch != 1) %>% pull(UDN_ID)
wrong_UDN_one_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_one)
wrong_UDN_one_df$batch <- 1
right_UDN_one_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_one)
FRASER_metadata_important <- rbind(right_UDN_one_df, wrong_UDN_one_df)

batch_info_2_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 2) %>% pull(PPID)
two <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_2_UDN_ID)
wrong_UDN_two <- two %>% filter(batch != 2) %>% pull(UDN_ID)
wrong_UDN_two_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_two)
wrong_UDN_two_df$batch <- 2
right_UDN_two_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_two)
FRASER_metadata_important <- rbind(right_UDN_two_df, wrong_UDN_two_df)

batch_info_3_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 3) %>% pull(PPID)
three <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_3_UDN_ID)
wrong_UDN_three <- three %>% filter(batch != 3) %>% pull(UDN_ID)
wrong_UDN_three_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_three)
wrong_UDN_three_df$batch <- 3
right_UDN_three_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_three)
FRASER_metadata_important <- rbind(right_UDN_three_df, wrong_UDN_three_df)

batch_info_4_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 4) %>% pull(PPID)
four <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_4_UDN_ID)
wrong_UDN_four <- four %>% filter(batch != 4) %>% pull(UDN_ID)
wrong_UDN_four_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_four)
wrong_UDN_four_df$batch <- 4
right_UDN_four_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_four)
FRASER_metadata_important <- rbind(right_UDN_four_df, wrong_UDN_four_df)

batch_info_5_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 5) %>% pull(PPID)
five <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_5_UDN_ID)
wrong_UDN_five <- five %>% filter(batch != 5) %>% pull(UDN_ID)
wrong_UDN_five_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_five)
wrong_UDN_five_df$batch <- 5
right_UDN_five_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_five)
FRASER_metadata_important <- rbind(right_UDN_five_df, wrong_UDN_five_df)

batch_info_6_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 6) %>% pull(PPID)
six <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_6_UDN_ID)
wrong_UDN_six <- six %>% filter(batch != 6) %>% pull(UDN_ID)
wrong_UDN_six_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_six)
wrong_UDN_six_df$batch <- 6
right_UDN_six_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_six)
FRASER_metadata_important <- rbind(right_UDN_six_df, wrong_UDN_six_df)

batch_info_7_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 7) %>% pull(PPID)
seven <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_7_UDN_ID)
wrong_UDN_seven <- seven %>% filter(batch != 7) %>% pull(UDN_ID)
wrong_UDN_seven_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_seven)
wrong_UDN_seven_df$batch <- 7
right_UDN_seven_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_seven)
FRASER_metadata_important <- rbind(right_UDN_seven_df, wrong_UDN_seven_df)

batch_info_8_UDN_ID <- batch_info %>% filter(`SEQ RUN #` == 8) %>% pull(PPID)
eight <- FRASER_metadata_important %>% filter(UDN_ID %in% batch_info_8_UDN_ID)
wrong_UDN_eight <- eight %>% filter(batch != 8) %>% pull(UDN_ID)
wrong_UDN_eight_df <- FRASER_metadata_important %>% filter(UDN_ID %in% wrong_UDN_eight)
wrong_UDN_eight_df$batch <- 8
right_UDN_eight_df <- FRASER_metadata_important %>% filter(! UDN_ID %in% wrong_UDN_eight)
FRASER_metadata_important <- rbind(right_UDN_eight_df, wrong_UDN_eight_df)

#############----Add Enrollment Type----###############
metadata_GSS_only <- FRASER_metadata_important %>% filter(!is.na(GSS_ID)) %>% filter(is.na(RDID)) %>% filter(is.na(UDN_ID)) %>% pull(ID)
metadata_GSS_RDID <- FRASER_metadata_important %>% filter(!is.na(GSS_ID)) %>% filter(!is.na(RDID))  %>% pull(ID)
metadata_GSS_UDN <- FRASER_metadata_important %>% filter(!is.na(GSS_ID)) %>% filter(!is.na(UDN_ID)) %>% pull(ID)

GSS_only <- metadata_GSS_only
GSS_UDN <- c(metadata_GSS_RDID, metadata_GSS_UDN) %>% unique

UDN_ID_only <- FRASER_metadata_important %>% filter(!is.na(UDN_ID)) %>% filter(is.na(GSS_ID))  %>% pull(ID)
RDID_only <- FRASER_metadata_important %>% filter(!is.na(RDID)) %>% filter(is.na(GSS_ID)) %>% pull(ID)

UDN_only <- c(UDN_ID_only, RDID_only) %>% unique

FRASER_metadata_important$Enrollment_Type <- ifelse(FRASER_metadata_important$ID %in% GSS_only, "GSS",
                                ifelse(FRASER_metadata_important$ID %in% GSS_UDN, "GSS_UDN",
                                ifelse(FRASER_metadata_important$ID %in% UDN_only, "UDN", NA)))
FRASER_metadata_important %>% filter(is.na(Enrollment_Type))

######################################################
##############----Add Rachel's Data----###############
######################################################
colnames(rachels_table) <- c("UDN_ID", "sampleID", "RDID", "Institution", "Tissue", "affected_status_rachel", "batch_rachel", "family_id", "relationship", "related_proband",
                                   "age", "sex", "origin", "disease_category", "resolved_case", "causal_gene", "RIN_rachel", "INDEX1", "INDEX2", "READ_LENGTH", "PIXEL_DIST",
                                   "BCL_DIR", "FASTQ_DIR")
rachels_table_select <- rachels_table %>% select(UDN_ID, RDID, BCL_DIR, FASTQ_DIR, Institution)
FRASER_metadata_important <- left_join(FRASER_metadata_important, rachels_table_select)
FRASER_metadata_important <- FRASER_metadata_important %>% select(-number, -numbers)

######################################################
################----Add Sample ID----#################
######################################################
setdiff(FRASER_result_combined$wetlab_id, FRASER_metadata_important$sample_ID)
#colnames(FRASER_result_combined) <- c("sample_ID")

#FRASER_result_combined_RD <- filter(FRASER_result_combined, grepl("RD", sampleID))
#colnames(FRASER_result_combined_RD) <- c("RDID")
#FRASER_RD <- left_join(FRASER_result_combined_RD, FRASER_metadata_important)
#FRASER_RD$sampleID <- FRASER_RD$RDID

#FRASER_result_combined_GSS <- filter(FRASER_result_combined, grepl("GSS", sampleID))
#colnames(FRASER_result_combined_GSS) <- c("GSS_ID")
#FRASER_GSS <- left_join(FRASER_result_combined_GSS, FRASER_metadata_important)
#FRASER_GSS$sampleID <- FRASER_GSS$GSS_ID

#FRASER_result_combined_UDN <- filter(FRASER_result_combined, grepl("UDN", sampleID))
#colnames(FRASER_result_combined_UDN) <- c("UDN_ID")
#FRASER_UDN <- left_join(FRASER_result_combined_UDN, FRASER_metadata_important)
#FRASER_UDN$sampleID <- FRASER_UDN$UDN_ID

#FRASER_metadata_important_with_missing <- bind_rows(FRASER_RD, FRASER_GSS, FRASER_UDN) %>% filter(!is.na(sampleID)) %>% unique
#missing <- FRASER_metadata_important_with_missing %>% filter(is.na(affected_status)) %>% pull(GSS_ID)
#found <- joined_missing %>% filter(participant_id %in% missing) 
#found <- found[, !colSums(is.na(found)), drop = FALSE] %>% select(-ancestry_detail, -gregor_center, -maternal_id, -missing_variant_case, -paternal_id, 
 #                    -reported_ethnicity, -reported_race, -twin_id, -phenotype_id, -onset_age_range, -number)
#colnames(found) <- c("GSS_ID", "affected_status", "family_identifier", "notes", "proband_relationship", "sex", "solved_status", "additional_details", 
  #                   "additional_modifiers", "ontology", "presence", "term_id")
#found$HPO_terms <- paste(found$term_id, ";", found$additional_modifiers)
#found$HPO_terms <- gsub("NA;", "", found$HPO_terms)
#found$HPO_terms <- gsub(";NA", "", found$HPO_terms)
#found$HPO_terms <- gsub(" ; NA", "", found$HPO_terms)

#found$additional_details <- c("chronic congestion", "Accretions on teeth, Maxillary lip tie;single central incisor behind two normal upper incisors")
#found$notes <- paste(found$notes, "|", found$additional_details)

#found$affected_status <- c("Case", "Case")
#found$sex <- c("M", "F")

#found <- found %>% select(-additional_details, -additional_modifiers, -ontology, -presence, -term_id)

#missing_cols <- setdiff(colnames(FRASER_metadata_important_with_missing), colnames(found))
#found$sampleID <- found$GSS_ID
#found$age <- c(NA, "2")

#FRASER_metadata_important_without_missing <-  FRASER_metadata_important_with_missing %>% filter(! GSS_ID %in% missing) %>% unique

#FRASER_metadata_important <- bind_rows(FRASER_metadata_important_without_missing, found)

#n=422

################----Low RIN----#################
FRASER_metadata_important$RIN <- as.numeric(FRASER_metadata_important$RIN)
low_RIN <- FRASER_metadata_important %>% filter(RIN < 7) %>% select(sample_ID, RIN)
colnames(low_RIN) <- c("sampleID", "RIN")

write_csv(low_RIN, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

##############----Missing RIN found by Devon----##############
found_RIN <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/no_RIN.xlsx")
colnames(found_RIN) <- c("GSS_ID", "UDN_ID", "RDID", "RIN")

FRASER_metadata_important_with_RIN <- FRASER_metadata_important %>% filter(!is.na(RIN))

FRASER_metadata_important_missing_RIN <- FRASER_metadata_important %>% filter(is.na(RIN)) %>% select(-RIN)
FRASER_metadata_important_missing_RIN_sample_ID <- FRASER_metadata_important_missing_RIN %>% pull(sample_ID)

UDN_RIN <- filter(found_RIN, UDN_ID %in% FRASER_metadata_important_missing_RIN_sample_ID) %>% select(UDN_ID, RIN)
UDN_RIN$RIN <- as.numeric(UDN_RIN$RIN)
GSS_RIN <- filter(found_RIN, GSS_ID %in% FRASER_metadata_important_missing_RIN_sample_ID) %>% select(GSS_ID, RIN)
GSS_RIN$RIN <- as.numeric(GSS_RIN$RIN)

FRASER_metadata_important_missing_RIN_UDN <- left_join(FRASER_metadata_important_missing_RIN, UDN_RIN) %>% filter(!is.na(RIN))
FRASER_metadata_important_missing_RIN_GSS <- left_join(FRASER_metadata_important_missing_RIN, GSS_RIN) %>% filter(!is.na(RIN))

FRASER_metadata_important_missing_RIN <- bind_rows(FRASER_metadata_important_missing_RIN_UDN, FRASER_metadata_important_missing_RIN_GSS)
FRASER_metadata_important_with_RIN$RIN <- as.numeric(FRASER_metadata_important_with_RIN$RIN )

FRASER_metadata_important <- bind_rows(FRASER_metadata_important_with_RIN, FRASER_metadata_important_missing_RIN)
names(FRASER_metadata_important)[names(FRASER_metadata_important) == 'sample_ID'] <- 'sampleID'
write_csv(FRASER_metadata_important, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")

################----Missing Samples----#################
not_missing <- FRASER_metadata_important %>% filter(!is.na(RIN)) %>% filter(!is.na(age)) %>% filter(!is.na(sex)) %>% filter(!is.na(batch)) %>% filter(!is.na(affected_status)) %>% filter(affected_status != "Unknown") %>% pull(sampleID)
missing <- FRASER_metadata_important %>% filter(! sampleID %in% not_missing) %>% select(sampleID)

write_csv(missing, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_metadata.csv")

################----Analysis Only----#################
low_RIN <- low_RIN %>% pull(sampleID)
missing <- missing %>% pull(sampleID)

removed <- c(low_RIN, missing) %>% unique

FRASER_metadata_important_analysis <- FRASER_metadata_important %>% filter(! sampleID %in% missing)
write_csv(FRASER_metadata_important_analysis, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/analyzed_FRASER_metadata.csv")

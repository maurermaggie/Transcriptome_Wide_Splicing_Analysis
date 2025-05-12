library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Filter and Add Group to Joined-----#################
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/DataFrames/metadata_counts_outlier_joined.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")
Devon <- read_excel("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches_dbedit .xlsx") %>% select(sampleID, disease_category)
colnames(Devon) <- c("sampleID", "category")

sampleIDs <- joined %>% pull(sampleID)
run_FRASER <- metadata %>% filter(sampleID %in% sampleIDs)

run_FRASER %>% group_by(Enrollment_Type) %>% tally

run_FRASER <- run_FRASER %>% filter(affected_status != "Unknown")

metadata_GSS_only <- run_FRASER %>% filter(!is.na(GSS_ID)) %>% filter(is.na(RDID)) %>% filter(is.na(UDN_ID)) %>% pull(ID) %>% length
metadata_GSS_RDID <- run_FRASER %>% filter(!is.na(GSS_ID)) %>% filter(!is.na(RDID))  %>% pull(ID) %>% length
metadata_GSS_UDN <- run_FRASER %>% filter(!is.na(GSS_ID)) %>% filter(!is.na(UDN_ID)) %>% pull(ID) %>% length
UDN_ID_only <- run_FRASER %>% filter(!is.na(UDN_ID)) %>% filter(is.na(GSS_ID))  %>% pull(ID) %>% length

metadata_cleaned_case <- run_FRASER %>% filter(affected_status == "Case")
metadata_cleaned_case <- left_join(metadata_cleaned_case, Devon)

metadata_cleaned_case %>% group_by(affected_status) %>% tally

metadata_cleaned_case <- metadata_cleaned_case %>% mutate(disease_category = ifelse(metadata_cleaned_case$disease_category == "Nephrology (kidney diseases and disorders)" , "Nephrology",
               ifelse(metadata_cleaned_case$disease_category == "N/A", NA, 
               ifelse(metadata_cleaned_case$disease_category == "Musculoskeletal and Ortho", "Musculoskeletal/orthopedics",
               ifelse(metadata_cleaned_case$disease_category == "other", "Other", 
               ifelse(metadata_cleaned_case$disease_category ==  "Cardiology and vascular", "Cardiology/vascular conditions",
               ifelse(metadata_cleaned_case$disease_category == "Neurologic", "Neurology",
               ifelse(metadata_cleaned_case$disease_category == "Neuromuscular", "Neurology", 
               ifelse(metadata_cleaned_case$disease_category == "Neurodevelopmental", "Neurology", 
               ifelse(metadata_cleaned_case$disease_category == "Gastroenterology (disorder of the stomach and intestines)", "Gastroenterology", 
               ifelse(metadata_cleaned_case$disease_category == "Allergies and disorders of the immune system", "Immunology/Allergies", 
               ifelse(metadata_cleaned_case$disease_category == "#N/A", NA, 
               ifelse(metadata_cleaned_case$disease_category == "Musculoskeletal", "Musculoskeletal/orthopedics", 
               ifelse(metadata_cleaned_case$disease_category == "Allergy/Immune", "Immunology/Allergies",
               ifelse(metadata_cleaned_case$disease_category == "Musculoskeletal and orthopedics (structural and functional disorders of muscles, bones, and joints)", "Musculoskeletal/orthopedics",
               ifelse(metadata_cleaned_case$disease_category == "neurology", "Neurology",
               ifelse(metadata_cleaned_case$disease_category == "Cardiology and vascular conditions (heart, artery, vein, and lymph disorders)","Cardiology/vascular conditions",
               ifelse(metadata_cleaned_case$disease_category == "Neurology (disorders of the nervous system, including brain and spinal cord)", "Neurology",
               ifelse(metadata_cleaned_case$disease_category == "Musculoskeletal and orthopedics", "Musculoskeletal/orthopedics", 
               ifelse(metadata_cleaned_case$disease_category =="Pulmonary", "Pulmonology", 
               ifelse(metadata_cleaned_case$disease_category =="Nephrology (kidney diseases and disorders)", "Nephrology",
               ifelse(metadata_cleaned_case$disease_category == "Cardiolgy and other vascular", "Cardiology/vascular conditions", 
               ifelse(metadata_cleaned_case$disease_category == "Gastrointestinal","Gastroenterology", 
               ifelse(metadata_cleaned_case$disease_category == "Infectious diseases","Immunology/Allergies", 
               ifelse(is.na(metadata_cleaned_case$disease_category), "Other", metadata_cleaned_case$disease_category)
               ))))))))))))))))))))))))

metadata_cleaned_case <- metadata_cleaned_case %>% mutate(disease_category = ifelse(is.na(metadata_cleaned_case$disease_category), "Other", metadata_cleaned_case$disease_category))
metadata_cleaned_case %>% group_by(disease_category) %>% tally() %>% arrange(desc(n))
run_FRASER %>% group_by(affected_status) %>% tally() %>% arrange(desc(n))

metadata_cleaned_case <- metadata_cleaned_case %>% mutate(category = ifelse(is.na(metadata_cleaned_case$category), "Other", 
               ifelse(metadata_cleaned_case$category == "Allergy/Immune", "Immunology/Allergies", 
               ifelse(metadata_cleaned_case$category == "NA", "Other", 
               ifelse(metadata_cleaned_case$category == "Infectious diseases", "Immunology/Allergies", 
               ifelse(metadata_cleaned_case$category == "Neurodevelopmental", "Neurology", metadata_cleaned_case$category)
               )))))

metadata_cleaned_case %>% group_by(category) %>% tally() %>% arrange(desc(n)) %>% print(n=40)

run_FRASER <- metadata %>% filter(sampleID %in% sampleIDs)

rachels_table <- fread("/home/maurertm/smontgom/shared/UDN/Data/SampleTables/sample_table_UDN_sep2023.txt")
colnames(rachels_table) <- c("UDN_ID", "sampleID", "RDID", "Institution", "Tissue", "affected_status_rachel", "batch_rachel", "family_id", "relationship", "related_proband",
                                   "age", "sex", "origin", "disease_category", "resolved_case", "causal_gene", "RIN_rachel", "INDEX1", "INDEX2", "READ_LENGTH", "PIXEL_DIST",
                                   "BCL_DIR", "FASTQ_DIR")
rachels_table_select <- rachels_table %>% select(UDN_ID, RDID, BCL_DIR, FASTQ_DIR, Institution)

run_FRASER <- run_FRASER %>% filter(affected_status != "Unknown")
run_FRASER_RDID_in_Rachels_table <- run_FRASER %>% filter(RDID %in% rachels_table_select$RDID)
run_FRASER_RDID_NOT_in_Rachels_table <- run_FRASER %>% filter(! RDID %in% rachels_table_select$RDID)

run_FRASER_RDID_NOT_UDN_in_Rachels_table <- run_FRASER_RDID_NOT_in_Rachels_table %>% filter(UDN_ID %in% rachels_table_select$RDID)
run_FRASER_RDID_NOT_OR_UDN_in_Rachels_table <- run_FRASER_RDID_NOT_in_Rachels_table %>% filter(! UDN_ID %in% rachels_table_select$RDID)

run_FRASER_RDID_in_Rachels_table %>% nrow
run_FRASER_RDID_NOT_OR_UDN_in_Rachels_table %>% nrow

leafcutter <- fread("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/LeafCutterMD/SplitByInd/all_annotated_adjusted_zscores_collapsedbygenes.hg38.sorted_dedupOptical_minMQ255.bed.gz")
metadata_old <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_minor_spliceosome/output/output_Nov_18_try5/DataFrames/metadata_counts_outlier_joined.csv")
metadata_old <- metadata_old %>% filter(!is.na(sex)) %>% filter(!is.na(age)) %>% filter(!is.na(affected_status)) %>% filter(!is.na(RIN)) %>% filter(RIN>=7) %>% filter(affected_status!="Unknown")

intersect(metadata_old$RDID, leafcutter$V5) %>% unique %>% length #270 in length
intersect(run_FRASER$RDID, leafcutter$V5) %>% unique %>% length
#252

novel_new <- run_FRASER %>% filter(is.na(RDID)) %>% select(GSS_ID, UDN_ID) #128
old_new <- run_FRASER %>% filter(!is.na(RDID)) %>% select(GSS_ID, UDN_ID) #257

novel_old <- metadata_old %>% filter(is.na(RDID)) %>% select(GSS_ID, UDN_ID) #120
old_old <- metadata_old %>% filter(!is.na(RDID)) %>% select(GSS_ID, UDN_ID) #270

setdiff(metadata_old$sampleID, run_FRASER$sampleID)
setdiff(run_FRASER$sampleID, metadata_old$sampleID)

FRASER1 <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_no_missing_rm_seqbatch_8/FRASER_filtered.rds")
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_annotated <- annotateRangesWithTxDb(FRASER1, txdb=txdb, orgDb=orgDb)

filtered_genes <- fds_annotated@rowRanges %>% as.data.frame %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol) %>% unique %>% length

FRASER2 <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake_old_filters/output/FRASER_no_missing_rm_seqbatch_8/FRASER2/FRASER_filtered.rds")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_annotated <- annotateRangesWithTxDb(FRASER2, txdb=txdb, orgDb=orgDb)

filtered_genes <- fds_annotated@rowRanges %>% as.data.frame %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol) %>% unique %>% length

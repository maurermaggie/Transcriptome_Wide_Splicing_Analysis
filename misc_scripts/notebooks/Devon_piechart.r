library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")
missing <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_metadata.csv")

remove <- c(low_RIN$sampleID, missing$sampleID)

metadata_cleaned <- metadata %>% filter(! sampleID %in% remove)

#affected_status
metadata_cleaned %>% group_by(affected_status) %>% tally()
217+173

390-287
103-7

#disease_category
metadata_cleaned$disease_category %>% unique
metadata_cleaned_case <- metadata_cleaned %>% filter(affected_status == "Case")

metadata_cleaned_case <- metadata_cleaned_case %>% mutate(disease_category = ifelse(metadata_cleaned_case$disease_category == "Nephrology (kidney diseases and disorders)" , "Nephrology",
               ifelse(metadata_cleaned_case$disease_category == "N/A", NA, 
               ifelse(metadata_cleaned_case$disease_category == "Musculoskeletal and Ortho", "Musculoskeletal/orthopedics",
               ifelse(metadata_cleaned_case$disease_category == "other", "Other", 
               ifelse(metadata_cleaned_case$disease_category ==  "Cardiology and vascular", "Cardiology/vascular conditions",
               ifelse(metadata_cleaned_case$disease_category == "Neurologic", "Neurology",
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
               ifelse(is.na(metadata_cleaned_case$disease_category), "Other", metadata_cleaned_case$disease_category)
               )))))))))))))))))))

metadata_pie <- metadata_cleaned_case %>% group_by(disease_category) %>% tally() 

pie(metadata_pie$n, labels = metadata_pie$disease_category, main="Pie Chart of Countries", radius = 1, cex = 0.5)

#disease_category
metadata_cleaned %>% group_by(Enrollment_Type) %>% tally()

#get number genes and junctions
rds <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/output/FRASER_Stanford_Filters/FRASER_filtered.rds")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_annotated <- annotateRangesWithTxDb(rds, txdb=txdb, orgDb=orgDb)

filtered_genes <- fds_annotated@rowRanges %>% as.data.frame %>% filter(!is.na(hgnc_symbol)) %>% pull(hgnc_symbol)

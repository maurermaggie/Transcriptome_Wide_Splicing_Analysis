library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

########################################################################
####################-----Remove missing metadata-----###################
########################################################################
#metadata %>% filter(is.na(batch))
#metadata %>% filter(is.na(sex))
#no_age <- metadata %>% filter(is.na(age)) %>% pull(sampleID)
#no_RIN <- metadata %>% filter(is.na(RIN)) %>% pull(sampleID)
#no_affected_status <- metadata %>% filter(affected_status == "Unknown") %>% pull(sampleID)

#missing_fields <- c(no_age, no_RIN, no_affected_status)

#metadata_filtered <- metadata %>% filter(! sampleID %in% missing_fields)
#metadata_filtered %>% group_by(affected_status) %>% tally()

#missing_metadata <- setdiff(metadata$sampleID, metadata_filtered$sampleID) %>% data.frame
#colnames(missing_metadata) <- c("sampleID")

#write_csv(missing_metadata, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_metadata.csv")

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
RNU6ATAC <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs)
RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs)
GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs)
UDN550488 <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)
UDN238929 <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)

non_RNU4ATAC <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")) %>% pull(no_MIGs_with_theta_juncs)

non_RNU4ATAC_6ATAC <- joined %>% filter(! sampleID %in% c("RD380", "RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")) %>% pull(no_MIGs_with_theta_juncs)

########################################################################
########################-----Get Stats----##############################
########################################################################
RNU4ATAC <- c(RD268, GSS225379, UDN550488, UDN238929)
mean_RNU4ATAC <- mean(RNU4ATAC)
sd_RNU4ATAC <- sd(RNU4ATAC)

mean_non_RNU4ATAC <- mean(non_RNU4ATAC)
sd_non_RNU4ATAC <- sd(non_RNU4ATAC)

mean_non_RNU4ATAC_6ATAC <- mean(non_RNU4ATAC_6ATAC)
sd_non_RNU4ATAC_6ATAC <- sd(non_RNU4ATAC_6ATAC)

mean_RNU4ATAC/mean_non_RNU4ATAC

RNU6ATAC/mean_non_RNU4ATAC_6ATAC 

all <- c(RNU6ATAC, RNU4ATAC, non_RNU4ATAC_6ATAC)
mean(all)

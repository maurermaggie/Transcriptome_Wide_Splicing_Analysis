library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Get Stats-----#################
RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
RD268_zscore <- RD268/ sd(joined$no_MIGs_theta_juncs)

GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
GSS225379_zscore <- GSS225379/ sd(joined$no_MIGs_theta_juncs)

UDN550488.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_zscore <- UDN550488.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_theta_juncs)

UDN238929.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_zscore <- UDN238929.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_theta_juncs)

RD380 <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_theta_juncs) - mean(joined$no_MIGs_theta_juncs)
RD380_zscore <- RD380/ sd(joined$no_MIGs_theta_juncs)

mean(RNU4ATAC$no_MIGs_theta_juncs)
sd(RNU4ATAC$no_MIGs_theta_juncs)

mean(non_RNU4ATAC$no_MIGs_theta_juncs)
sd(non_RNU4ATAC$no_MIGs_theta_juncs)

mean(joined$no_MIGs_theta_juncs)
sd(joined$no_MIGs_theta_juncs)

mean(RNU4ATAC$no_MIGs_theta_juncs)/mean(non_RNU4ATAC$no_MIGs_theta_juncs)


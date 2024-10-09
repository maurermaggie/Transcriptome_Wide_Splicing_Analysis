library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Introns Retained-----#################
joined_RNU4ATAC <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_RNU4ATAC$MS <- "RNU4ATAC-opathies"

joined_RNU6ATAC <- joined %>% filter(sampleID == "RD380")
joined_RNU6ATAC$MS <- "RNU6ATAC-opathy"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_non_MS$MS <- "All Other Samples"

joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))

joined <- bind_rows(joined_RNU4ATAC, joined_non_MS, joined_RNU6ATAC)

joined_RNU6ATAC$no_MIGs_theta_juncs/mean(joined_non_MS$no_MIGs_theta_juncs)

mean(joined_non_MS$no_MIGs_theta_juncs)
sd(joined_non_MS$no_MIGs_theta_juncs)

mean(joined_MS$no_MIGs_theta_juncs)
sd(joined_MS$no_MIGs_theta_juncs)

mean(joined_RNU6ATAC$no_MIGs_theta_juncs)/mean(joined_non_MS$no_MIGs_theta_juncs)
mean(joined_MS$no_MIGs_theta_juncs)/mean(joined_non_MS$no_MIGs_theta_juncs)

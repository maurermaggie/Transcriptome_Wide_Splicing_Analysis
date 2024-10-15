library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
joined_fp <- args[1]
joined <- read_csv(joined_fp)

number_MIGs_affected_boxplot <- args[2]
number_MIGs_with_theta <- args[3]
means_number_MIGs_with_theta <- args[4]

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Introns Retained-----#################
y_lab <- expression("Number of MIGs with Outlier" ~ theta ~ "Events")
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_with_theta_juncs)) + 
  geom_boxplot() +
  ylab(y_lab) +
  xlab("Groups")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text( vjust = 1))


MIG_boxplot
ggsave(filename=number_MIGs_affected_boxplot, plot=MIG_boxplot,  limitsize = FALSE, units = "in")

################-----Get Stats-----#################
RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
RD268_zscore <- RD268/ sd(joined$no_MIGs_with_theta_juncs)
RD268_no <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs) %>% length

GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
GSS225379_zscore <- GSS225379/ sd(joined$no_MIGs_with_theta_juncs)
GSS225379_no <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs) %>% length

UDN550488.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_zscore <- UDN550488.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_with_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_no <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) %>% length

UDN238929.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_zscore <- UDN238929.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_with_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_no <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) %>% length

RD380 <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
RD380_zscore <- RD380/ sd(joined$no_MIGs_with_theta_juncs)
RD380_no <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs) %>% length

zscores <- data.frame(values=c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"),
                      zscores_theta_in_MIGs=c(RD268_zscore, GSS225379_zscore, UDN550488.Aligned.sortedByCoord.out.bam_zscore, UDN238929.Aligned.sortedByCoord.out.bam_zscore, RD380_zscore),
                      number_theta_juncs_in_MIGs = c(RD268_no, GSS225379_no, UDN550488.Aligned.sortedByCoord.out.bam_no, UDN550488.Aligned.sortedByCoord.out.bam_no, RD380_no))

write_csv(zscores, number_MIGs_with_theta)

RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
non_RNU4ATAC$type <- NA

RNU4ATAC_mean <- mean(RNU4ATAC$no_MIGs_with_theta_juncs)
RNU4ATA_sd <- sd(RNU4ATAC$no_MIGs_with_theta_juncs)

non_RNU4ATAC_mean <- mean(non_RNU4ATAC$no_MIGs_with_theta_juncs)
non_RNU4ATAC_sd <- sd(non_RNU4ATAC$no_MIGs_with_theta_juncs)

overall_mean <- mean(joined$no_MIGs_with_theta_juncs)
overall_sd <- sd(joined$no_MIGs_with_theta_juncs)

means_df <- data.frame(values= c("RNU4ATAC", "non-RNU4ATAC", "overall"), mean=c(RNU4ATAC_mean, non_RNU4ATAC_mean, overall_mean), sd=c(RNU4ATA_sd, non_RNU4ATAC_sd, overall_sd))

write_csv(means_df, means_number_MIGs_with_theta)

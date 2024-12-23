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
y_lab <- expression("Number of MIGs with" ~ theta ~ "Outlier Junctions")
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)
title_lab <- expression(atop("Number of Minor Intron containing Genes",
                     "with Significant" ~theta~ "Outliers"))

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c("All Other Samples", "RNU4ATAC-opathies")
MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_with_theta_juncs, fill=MS)) + 
  geom_boxplot(show.legend=FALSE, outlier.shape=21, outlier.alpha=.7, outlier.size=10) +
  ylab(y_lab) +
  xlab("Groups")+
  labs(title=title_lab)  +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text( vjust = 1))+
  scale_fill_manual(values=colours)


MIG_boxplot
ggsave(filename=number_MIGs_affected_boxplot, plot=MIG_boxplot,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Get Stats-----#################
RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
RD268_zscore <- RD268/ sd(joined$no_MIGs_with_theta_juncs)
RD268_no <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs)

GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
GSS225379_zscore <- GSS225379/ sd(joined$no_MIGs_with_theta_juncs)
GSS225379_no <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs)

UDN550488.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_zscore <- UDN550488.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_with_theta_juncs)
UDN550488.Aligned.sortedByCoord.out.bam_no <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)

UDN238929.Aligned.sortedByCoord.out.bam <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_zscore <- UDN238929.Aligned.sortedByCoord.out.bam/ sd(joined$no_MIGs_with_theta_juncs)
UDN238929.Aligned.sortedByCoord.out.bam_no <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)

RD380 <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs) - mean(joined$no_MIGs_with_theta_juncs)
RD380_zscore <- RD380/ sd(joined$no_MIGs_with_theta_juncs)
RD380_no <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs)

zscores <- data.frame(values=c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"),
                      zscores_theta_in_MIGs=c(RD268_zscore, GSS225379_zscore, UDN550488.Aligned.sortedByCoord.out.bam_zscore, UDN238929.Aligned.sortedByCoord.out.bam_zscore, RD380_zscore),
                      number_theta_juncs_in_MIGs = c(RD268_no, GSS225379_no, UDN550488.Aligned.sortedByCoord.out.bam_no, UDN550488.Aligned.sortedByCoord.out.bam_no, RD380_no))

write_csv(zscores, number_MIGs_with_theta)

RNU4ATAC <- joined %>% arrange(desc(no_MIGs_with_theta_juncs)) %>% select(sampleID, no_MIGs_with_theta_juncs) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_MIGs_with_theta_juncs)) %>% select(sampleID, no_MIGs_with_theta_juncs) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
non_RNU4ATAC$type <- NA
RNU6ATAC <- joined %>% arrange(desc(no_MIGs_with_theta_juncs)) %>% select(sampleID, no_MIGs_with_theta_juncs) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "RNU6ATAC"
non_MS <- joined %>% arrange(desc(no_MIGs_with_theta_juncs)) %>% select(sampleID, no_MIGs_with_theta_juncs) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
non_MS$type <- "Non Minor Spliceosome"

RNU4ATAC_mean <- mean(RNU4ATAC$no_MIGs_with_theta_juncs)
RNU4ATA_sd <- sd(RNU4ATAC$no_MIGs_with_theta_juncs)
RNU4ATA_first <- quantile(RNU4ATAC$no_MIGs_with_theta_juncs, 0.25)
RNU4ATA_last <- quantile(RNU4ATAC$no_MIGs_with_theta_juncs, 0.75)

non_RNU4ATAC_mean <- mean(non_RNU4ATAC$no_MIGs_with_theta_juncs)
non_RNU4ATAC_sd <- sd(non_RNU4ATAC$no_MIGs_with_theta_juncs)
non_RNU4ATA_first <- quantile(non_RNU4ATAC$no_MIGs_with_theta_juncs, 0.25)
non_RNU4ATA_last <- quantile(non_RNU4ATAC$no_MIGs_with_theta_juncs, 0.75)

overall_mean <- mean(joined$no_MIGs_with_theta_juncs)
overall_sd <- sd(joined$no_MIGs_with_theta_juncs)
overall_first <- quantile(joined$no_MIGs_with_theta_juncs, 0.25)
overall_last <- quantile(joined$no_MIGs_with_theta_juncs, 0.75)

RNU6ATAC_mean <- mean(RNU6ATAC$no_MIGs_with_theta_juncs)
RNU6ATAC_sd <- sd(RNU6ATAC$no_MIGs_with_theta_juncs)
RNU6ATAC_first <- quantile(RNU6ATAC$no_MIGs_with_theta_juncs, 0.25)
RNU6ATAC_last <- quantile(RNU6ATAC$no_MIGs_with_theta_juncs, 0.75)

non_MS_mean <- mean(non_MS$no_MIGs_with_theta_juncs)
non_MS_sd <- sd(non_MS$no_MIGs_with_theta_juncs)
non_MS_first <- quantile(non_MS$no_MIGs_with_theta_juncs, 0.25)
non_MS_last <- quantile(non_MS$no_MIGs_with_theta_juncs, 0.75)

means_df <- data.frame(values= c("RNU4ATAC", "non-RNU4ATAC", "overall", "RNU6ATAC", "non-Minor_Spliceosome-opathy"), 
mean=c(RNU4ATAC_mean, non_RNU4ATAC_mean, overall_mean, RNU6ATAC_mean, non_MS_mean), 
sd=c(RNU4ATA_sd, non_RNU4ATAC_sd, overall_sd, RNU6ATAC_sd, non_MS_sd), 
first_quantile=c(RNU4ATA_first, non_RNU4ATA_first, overall_first, RNU6ATAC_first, non_MS_first), 
last_quantile=c(RNU4ATA_last, non_RNU4ATA_last, overall_last, RNU6ATAC_last, non_MS_last))

write_csv(means_df, means_number_MIGs_with_theta)

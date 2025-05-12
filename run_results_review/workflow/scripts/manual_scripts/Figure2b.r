library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Filter and Add Group to Joined-----#################
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/DataFrames/metadata_counts_outlier_joined.csv")

joined <- joined %>% select(sampleID, no_MIGs_with_theta_juncs, no_MIGs_with_psi3_juncs, no_MIGs_with_psi5_juncs, no_MIGs_with_jaccard_juncs)
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)
joined_long <- melt(setDT(joined), id.vars = c("sampleID", "MS"), variable.name = "outlier_type")
joined_long$code <- paste(joined_long$MS, joined_long$outlier_type, sep="")

joined_MS <- joined_long %>% filter(MS=="RNU4ATAC-opathies")

colours <- c("#D7DCD8", "#DA786C", "#D7DCD8", "#7D5E62", "#D7DCD8", "#88ACAA", "#D7DCD8", "#F1AE7A")
names(colours) <- c("All Other Samplesno_MIGs_with_theta_juncs", "RNU4ATAC-opathiesno_MIGs_with_theta_juncs", "All Other Samplesno_MIGs_with_psi3_juncs",
"RNU4ATAC-opathiesno_MIGs_with_psi3_juncs", "All Other Samplesno_MIGs_with_psi5_juncs", "RNU4ATAC-opathiesno_MIGs_with_psi5_juncs",
"All Other Samplesno_MIGs_with_jaccard_juncs", "RNU4ATAC-opathiesno_MIGs_with_jaccard_juncs")

MIG_boxplot <- ggplot(joined_long, aes(x = outlier_type, y = value, fill=code)) + 
      geom_boxplot(show.legend=FALSE, outlier.shape=21, outlier.alpha=.7, outlier.size=10) +
      scale_fill_manual(values=colours)+
      ylab("Number of MIGs with Outlier Junctions") +
      xlab("Groups")+
  #labs(title=title_lab)  +
      scale_fill_manual(values = colours)+
      theme_update(plot.title = element_text(hjust = 0.5)) +
      theme_bw(base_size = 25)+
      theme(axis.text.x = element_text(vjust = 1))+
      geom_point(aes(1.2, 145), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(1.1, 194), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(1.2, 167), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(1.15, 179), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.25, 10), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.15, 15), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(2.1, 12), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.15, 6), colour="black", show.legend=FALSE, size=3) +      
      geom_point(aes(3.2, 9), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(3.1, 13), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(3.25, 10), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(3.15, 12), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.2, 93), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.1, 144), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(4.2, 152), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.15, 198), colour="black", show.legend=FALSE, size=3)

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/rm_batch_8/Plots/Figure2/boxplox_RNU4ATAC.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Get Stats RNU6ATAC-----#################
RNU6ATAC_filtered <- RNU6ATAC %>% select(sampleID, no_MIGs_with_theta_juncs, no_MIGs_with_psi3_juncs, no_MIGs_with_psi5_juncs, no_MIGs_with_jaccard_juncs)
write_csv(RNU6ATAC_filtered, stats_RD380)

################-----Filter and Add Group to Joined-----#################
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/NOT_rm_batch_8/DataFrames/metadata_counts_outlier_joined.csv")

joined <- joined %>% select(sampleID, no_MIGs_with_theta_juncs, no_MIGs_with_psi3_juncs, no_MIGs_with_psi5_juncs, no_MIGs_with_jaccard_juncs)
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)
joined_long <- melt(setDT(joined), id.vars = c("sampleID", "MS"), variable.name = "outlier_type")
joined_long$code <- paste(joined_long$MS, joined_long$outlier_type, sep="")

joined_MS <- joined_long %>% filter(MS=="RNU4ATAC-opathies")

colours <- c("#D7DCD8", "#DA786C", "#D7DCD8", "#7D5E62", "#D7DCD8", "#88ACAA", "#D7DCD8", "#F1AE7A")
names(colours) <- c("All Other Samplesno_MIGs_with_theta_juncs", "RNU4ATAC-opathiesno_MIGs_with_theta_juncs", "All Other Samplesno_MIGs_with_psi3_juncs",
"RNU4ATAC-opathiesno_MIGs_with_psi3_juncs", "All Other Samplesno_MIGs_with_psi5_juncs", "RNU4ATAC-opathiesno_MIGs_with_psi5_juncs",
"All Other Samplesno_MIGs_with_jaccard_juncs", "RNU4ATAC-opathiesno_MIGs_with_jaccard_juncs")

joined_MS$y_value <- c(1.2,1.1,1.2,1.15, 2.2,2.1,1.2,1.15, 3.2,3.1,3.2,3.15, 4.2,4.1,4.2,4.15)

MIG_boxplot <- ggplot(joined_long, aes(x = outlier_type, y = value, fill=code)) + 
      geom_boxplot(show.legend=FALSE, outlier.shape=21, outlier.alpha=.7, outlier.size=10) +
      scale_fill_manual(values=colours)+
      ylab("Number of MIGs with Outlier Junctions") +
      xlab("Groups")+
  #labs(title=title_lab)  +
      scale_fill_manual(values = colours)+
      theme_update(plot.title = element_text(hjust = 0.5)) +
      theme_bw(base_size = 25)+
      theme(axis.text.x = element_text(vjust = 1))+
      geom_point(aes(1.2, 194), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(1.1, 168), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(1.2, 141), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(1.15, 179), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.25, 15), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.2, 12), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(2.1, 10), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(2.25, 7), colour="black", show.legend=FALSE, size=3) +      
      geom_point(aes(3.3, 12), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(3.15, 9), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(3.25, 11), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(3.2, 11), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.2, 137), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.1, 149), colour="black", show.legend=FALSE, size=3) + 
      geom_point(aes(4.2, 88), colour="black", show.legend=FALSE, size=3) +
      geom_point(aes(4.15, 197), colour="black", show.legend=FALSE, size=3)

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/Plots/Figure2/boxplox_RNU4ATAC.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Get Stats RNU6ATAC-----#################
RNU6ATAC_filtered <- joined %>% select(sampleID, no_MIGs_with_theta_juncs, no_MIGs_with_psi3_juncs, no_MIGs_with_psi5_juncs, no_MIGs_with_jaccard_juncs) %>% filter(sampleID == "RD380")
write_csv(RNU6ATAC_filtered, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/Stats/RNU6ATAC.csv")

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
                      number_theta_juncs_in_MIGs = c(RD268_no, GSS225379_no, UDN550488.Aligned.sortedByCoord.out.bam_no, UDN238929.Aligned.sortedByCoord.out.bam_no, RD380_no))

write_csv(zscores, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/Stats/RNU4ATAC_zscores.csv")

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

write_csv(means_df, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/run_results_review/output/output_rm_seqbatch_8/Stats/means.csv")

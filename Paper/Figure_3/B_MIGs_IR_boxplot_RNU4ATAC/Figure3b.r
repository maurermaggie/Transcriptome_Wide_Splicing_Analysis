library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Introns Retained-----#################
title_lab <- expression(atop("Comparison of the Number of" ~ theta ~ "Events in MIGs",
                     "Per Sample between Samples with and without an Excess",
                     "Number of" ~ theta ~ "Events in MIGs"))
y_lab <- expression("Number of" ~ theta ~ "Events in MIGs")
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_theta_juncs)) + 
  geom_boxplot() +
  ylab(y_lab) +
  xlab("Groups")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text( vjust = 1))


MIG_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_3/B_MIGs_IR_boxplot_RNU4ATAC/plots/MIG_boxplot_RNU4ATAC.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in")

mean(joined_non_MS$no_MIGs_theta_juncs) 
sd(joined_non_MS$no_MIGs_theta_juncs) 
mean(joined_MS$no_MIGs_theta_juncs) 
sd(joined_MS$no_MIGs_theta_juncs) 
312.5/2.813131

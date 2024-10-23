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
joined_RNU4ATAC <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_RNU4ATAC$MS <- "RNU4ATAC-opathies"

joined_RNU6ATAC <- joined %>% filter(sampleID == "RD380")
joined_RNU6ATAC$MS <- "RNU6ATAC-opathy"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_RNU4ATAC, joined_non_MS, joined_RNU6ATAC)

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_theta_juncs)) + 
  geom_boxplot() +
  ylab(y_lab) +
  xlab("Groups")+
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) 


MIG_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_4/B_MIGs_IR_boxplot/plots/MIG_boxplot_RNU6ATAC.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in")

mean(joined_non_MS$no_MIGs_theta_juncs) 
sd(joined_non_MS$no_MIGs_theta_juncs) 
mean(joined_RNU4ATAC$no_MIGs_theta_juncs) 
sd(joined_RNU4ATAC$no_MIGs_theta_juncs) 
312.5/2.813131

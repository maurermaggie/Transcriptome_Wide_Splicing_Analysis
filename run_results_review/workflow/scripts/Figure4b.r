library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
joined_fp <- args[1]
joined <- read_csv(joined_fp)

number_MIGs_affected_boxplot_RNU6ATAC_RNU4ATAC <- args[2]

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Introns Retained-----#################
title_lab <- expression(atop("Comparison of the Number of" ~ theta ~ "Outlier Junctions in MIGs",
                     "Per Sample between Samples with and without an Excess",
                     "Number of" ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta ~ "Outlier Junctions in MIGs")
joined_RNU4ATAC <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
joined_RNU4ATAC$MS <- "RNU4ATAC-opathies"

joined_RNU6ATAC <- joined %>% filter(sampleID == "RD380")
joined_RNU6ATAC$MS <- "RNU6ATAC-opathy"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_RNU4ATAC, joined_non_MS, joined_RNU6ATAC)

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_with_theta_juncs)) + 
  geom_boxplot() +
  ylab(y_lab) +
  xlab("Groups")+
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) 


MIG_boxplot
ggsave(filename=number_MIGs_affected_boxplot_RNU6ATAC_RNU4ATAC, plot=MIG_boxplot,  limitsize = FALSE, units = "in", height=10, width=10)

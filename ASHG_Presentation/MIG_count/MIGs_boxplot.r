#library(EnsDb.Hsapiens.v79)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

################-----Number Minor Introns Retained-----#################
title_lab <- expression(atop("Comparison of the Number of MIGs with" ~theta~ "Events",
                     "Per Sample between Samples with and without an Excess",
                     "Number of MIGs with" ~theta~ "Events"))
y_lab <- expression("Number of MIGs with" ~ theta ~ "Events")
joined_MS <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_MS$MS <- "Samples 1-5"

joined_non_MS <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
joined_non_MS$MS <- "All Other Samples"

joined <- bind_rows(joined_MS, joined_non_MS)

dark_colors = c("white", "blue") %>% 
  col2rgb %>% #convert HEX colors to RGB Matrix
  "*"(0.7) %>% # make each component "darker"
  apply(2, lift_dv(rgb, maxColorValue = 255))
color <- c("white", "blue")
alpha_color <- 10
alpha_fill <- 0.5

MIG_boxplot <- ggplot(joined, aes(x=MS, y=no_MIGs_with_theta_juncs, fill=MS)) + 
  geom_boxplot(show.legend=FALSE) +
  ylab(y_lab) +
  xlab("Groups")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 40)+
  scale_fill_manual(values = alpha(color, alpha_fill)) + 
  scale_color_manual(values = alpha(color, alpha_color))

MIG_boxplot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/MIG_count/MIG_boxplot.pdf", plot=MIG_boxplot,  limitsize = FALSE, units = "in")

mean(joined_non_MS$no_MIGs_with_theta_juncs) 
sd(joined_non_MS$no_MIGs_with_theta_juncs) 
mean(joined_MS$no_MIGs_with_theta_juncs) 
sd(joined_MS$no_MIGs_with_theta_juncs) 
mean(joined_MS$no_MIGs_with_theta_juncs) /mean(joined_non_MS$no_MIGs_with_theta_juncs) 


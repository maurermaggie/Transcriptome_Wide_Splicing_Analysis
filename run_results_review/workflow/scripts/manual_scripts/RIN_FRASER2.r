library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

################-----rm batch 8-----#################
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Arriaga_2025/run_results_review_cleaned_github/output/response_to_review2/DataFrames/metadata_counts_outlier_joined.csv")
joined <- joined %>% select(sampleID, Jaccard_Junctions, RIN, age, sex, batch)
joined <- joined %>% filter(RIN > 7)

Esitmate_lm <- lm(Jaccard_Junctions ~ RIN + age + sex + batch, data = joined)
p_values <- summary(Esitmate_lm)$coefficients[,4] 

p.adjust(p_values, method="fdr", n=length(p_values))

RIN_plot <- ggplot(joined,aes(RIN, Jaccard_Junctions)) + 
  #labs(title="RIN vs # of Significant Outlier Junctions per Person")+
  geom_point(method='lm', formula= Jaccard_Junctions~RIN) +
  xlab("RIN") +
  ylab("Number of Genes with Significant Jaccard Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 25)

RIN_plot

ggsave(filename=RIN_fp, plot=RIN_plot,  limitsize = FALSE, units = "in", height=10, width=10)

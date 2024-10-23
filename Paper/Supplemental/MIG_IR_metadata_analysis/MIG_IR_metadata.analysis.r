library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(FRASER)
library(tidyverse)

joined_filtered <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
#######################-----Linear Model-----###########################
########################################################################
make_estimate_df <- function(row, dataframe) {
    Esitmate_lm_df <- lm(dataframe[[row]] ~ RIN + age + sex+ batch, data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="bonferroni", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
    Value <- summary(Esitmate_lm_df)$coefficients[,0]
    df$Value <- c("Intercept", "RIN", "age", "sexM", "batch")
    df$Numbers <- ifelse(df$Value == "Intercept", 5,
                            ifelse(df$Value == "RIN", 4,
                            ifelse(df$Value == "age", 3,
                            ifelse(df$Value == "sexM", 2, 1))))
   df
}

MIGs <- make_estimate_df("no_MIGs_theta_juncs", joined_filtered)

########################################################################
###################-----Variance Explained----##########################
########################################################################
title <- expression( "MIGs with" ~ theta ~ "Outliers")
x_lab <- expression(atop("Samples ordered by Number of MIGs",
                     "with Significant" ~ theta~ "Outliers"))
y_lab <- expression("Number of MIGs with Significant" ~ theta ~ "Outliers")

linear_variance_explained_gene_MIGs <- ggplot(MIGs, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.9), width = 0.2)+
  geom_text(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
                          ifelse(padjust < 0.001, "***", 
                            ifelse(padjust < 0.05, "*", "")))), 
            position = position_dodge(width = .5), vjust = 2, size = 30 / .pt, angle= 90) +
  labs(title_lab=title)+
  xlab(x_lab) +
  ylab(y_lab) +
  theme(axis.text.x = element_text(angle = 0, 
                                   vjust = 0.6, 
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal[c(1,3,5,7,9)])
  
linear_variance_explained_gene_MIGs
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/MIG_IR_metadata_analysis/plots/MIG_IR_metadata_analysis.pdf", plot=linear_variance_explained_gene_MIGs,  limitsize = FALSE, units = "in")



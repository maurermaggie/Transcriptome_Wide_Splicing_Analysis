library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
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

################-----Linear Model Genes-----#################
#.005
All_05_Genes <- make_estimate_df("All_Genes", joined_filtered)
Psi3_05_Genes <- make_estimate_df("Psi3_Genes", joined_filtered)
Psi5_05_Genes <- make_estimate_df("Psi5_Genes", joined_filtered)
Theta_05_Genes <- make_estimate_df("Theta_Genes", joined_filtered)
Jaccard_05_Genes <- make_estimate_df("Jaccard_Genes", joined_filtered)


################-----Linear Model Junctions-----#################
#.005
All_05_Junctions <- make_estimate_df("All_Junctions", joined_filtered)
Psi3_05_Junctions <- make_estimate_df("Psi3_Junctions", joined_filtered)
Psi5_05_Junctions <- make_estimate_df("Psi5_Junctions", joined_filtered)
Theta_05_Junctions <- make_estimate_df("Theta_Junctions", joined_filtered)
Jaccard_05_Genes <- make_estimate_df("Jaccard_Genes", joined_filtered)



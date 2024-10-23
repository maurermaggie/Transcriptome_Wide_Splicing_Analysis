library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined_filtered <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))
########################################################################
######################-----Logistic Model-----##########################
########################################################################
make_estimate_df <- function(dataframe, column) {
    Esitmate_lm_df <- glm(dataframe[[column]] ~ RIN + age + sex+ batch, family=binomial(link='logit'), data = dataframe)
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
All_05_Genes <- make_estimate_df(joined_filtered, "All_Genes_Outlier_Status")
All_01_Genes <- make_estimate_df(joined_filtered, "All_Genes_01_Outlier_Status")
All_001_Genes <- make_estimate_df(joined_filtered, "All_Genes_001_Outlier_Status")
All_1e_6_Genes <- make_estimate_df(joined_filtered, "All_Genes_1e_6_Outlier_Status")

Psi3_05_Genes <- make_estimate_df(joined_filtered, "Psi3_Genes_Outlier_Status")
Psi3_01_Genes <- make_estimate_df(joined_filtered, "Psi3_Genes_01_Outlier_Status")
Psi3_001_Genes <- make_estimate_df(joined_filtered, "Psi3_Genes_001_Outlier_Status")
Psi3_1e_6_Genes <- make_estimate_df(joined_filtered, "Psi3_Genes_1e_6_Outlier_Status")

Psi5_05_Genes <- make_estimate_df(joined_filtered, "Psi5_Genes_Outlier_Status")
Psi5_01_Genes <- make_estimate_df(joined_filtered, "Psi5_Genes_01_Outlier_Status")
Psi5_001_Genes <- make_estimate_df(joined_filtered, "Psi5_Genes_001_Outlier_Status")
Psi5_1e_6_Genes <- make_estimate_df(joined_filtered, "Psi5_Genes_1e_6_Outlier_Status")

Theta_05_Genes <- make_estimate_df(joined_filtered, "Theta_Genes_Outlier_Status")
Theta_01_Genes <- make_estimate_df(joined_filtered, "Theta_Genes_01_Outlier_Status")
Theta_001_Genes <- make_estimate_df(joined_filtered, "Theta_Genes_001_Outlier_Status")
Theta_1e_6_Genes <- make_estimate_df(joined_filtered, "Theta_Genes_1e_6_Outlier_Status")

Jaccard_05_Genes <- make_estimate_df(joined_filtered, "Jaccard_Genes_Outlier_Status")
Jaccard_01_Genes <- make_estimate_df(joined_filtered, "Jaccard_Genes_01_Outlier_Status")
Jaccard_001_Genes <- make_estimate_df(joined_filtered, "Jaccard_Genes_001_Outlier_Status")
Jaccard_1e_6_Genes <- make_estimate_df(joined_filtered, "Jaccard_Genes_1e_6_Outlier_Status")

################-----Linear Model Junctions-----#################
#.005
All_05_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_Outlier_Status")
All_01_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_01_Outlier_Status")
All_001_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_001_Outlier_Status")
All_1e_6_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_1e_6_Outlier_Status")

Psi3_05_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_Outlier_Status")
Psi3_01_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_01_Outlier_Status")
Psi3_001_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_001_Outlier_Status")
Psi3_1e_6_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_1e_6_Outlier_Status")

Psi5_05_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_Outlier_Status")
Psi5_01_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_01_Outlier_Status")
Psi5_001_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_001_Outlier_Status")
Psi5_1e_6_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_1e_6_Outlier_Status")

Theta_05_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_Outlier_Status")
Theta_01_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_01_Outlier_Status")
Theta_001_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_001_Outlier_Status")
Theta_1e_6_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_1e_6_Outlier_Status")

Jaccard_05_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_Outlier_Status")
Jaccard_01_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_01_Outlier_Status")
Jaccard_001_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_001_Outlier_Status")
Jaccard_1e_6_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_1e_6_Outlier_Status")

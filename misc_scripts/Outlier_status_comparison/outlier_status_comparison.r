library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

all_outliers <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_minor_spliceosome/output/output_Nov_18_try5/DataFrames/outliers.csv")

Theta_lm_df <- glm(Theta_Junctions_Outlier_Status ~ Jaccard_Junctions_Outlier_Status + Psi5_Junctions_Outlier_Status + Psi3_Junctions_Outlier_Status+ All_Junctions_Outlier_Status, family=binomial(link='logit'), data = all_outliers)
Theta <- data.frame(matrix(nrow=5, ncol=0))
Theta$Estimate <- summary(Theta_lm_df)$coefficients[,1] 
Theta$Absolute_Effect_Estimate <- abs(summary(Theta_lm_df)$coefficients[,1]) 
Theta$pvalue <- summary(Theta_lm_df)$coefficients[,4] 
Theta$padjust <- p.adjust(Theta$pvalue, method="fdr", n=length(Theta$pvalue))
total <- sum(Theta_lm_df$Absolute_Effect_Estimate)
Theta$variance_explained <- Theta$Absolute_Effect_Estimate/total
Theta$std_error <- summary(Theta_lm_df)$coefficients[,2]/total
Value <- summary(Theta_lm_df)$coefficients[,0]
Theta$Value <- c("Intercept", "Jaccard_Junctions_Outlier_Status", "Psi5_Junctions_Outlier_Status", "Psi3_Junctions_Outlier_Status", "All_Junctions_Outlier_Status")

Jaccard_lm_df <- glm(Jaccard_Junctions_Outlier_Status ~ Theta_Junctions_Outlier_Status + Psi5_Junctions_Outlier_Status + Psi3_Junctions_Outlier_Status+ All_Junctions_Outlier_Status, family=binomial(link='logit'), data = all_outliers)
Jaccard <- data.frame(matrix(nrow=5, ncol=0))
Jaccard$Estimate <- summary(Jaccard_lm_df)$coefficients[,1] 
Jaccard$Absolute_Effect_Estimate <- abs(summary(Jaccard_lm_df)$coefficients[,1]) 
Jaccard$pvalue <- summary(Jaccard_lm_df)$coefficients[,4] 
Jaccard$padjust <- p.adjust(Jaccard$pvalue, method="fdr", n=length(Jaccard$pvalue))
total <- sum(Jaccard_lm_df$Absolute_Effect_Estimate)
Jaccard$variance_explained <- Jaccard$Absolute_Effect_Estimate/total
Jaccard$std_error <- summary(Jaccard_lm_df)$coefficients[,2]/total
Value <- summary(Jaccard_lm_df)$coefficients[,0]
Jaccard$Value <- c("Intercept", "Theta_Junctions_Outlier_Status", "Psi5_Junctions_Outlier_Status", "Psi3_Junctions_Outlier_Status", "All_Junctions_Outlier_Status")

Psi5_lm_df <- glm(Psi5_Junctions_Outlier_Status ~ Theta_Junctions_Outlier_Status + Jaccard_Junctions_Outlier_Status + Psi3_Junctions_Outlier_Status+ All_Junctions_Outlier_Status, family=binomial(link='logit'), data = all_outliers)
Psi5 <- data.frame(matrix(nrow=5, ncol=0))
Psi5$Estimate <- summary(Psi5_lm_df)$coefficients[,1] 
Psi5$Absolute_Effect_Estimate <- abs(summary(Psi5_lm_df)$coefficients[,1]) 
Psi5$pvalue <- summary(Psi5_lm_df)$coefficients[,4] 
Psi5$padjust <- p.adjust(Psi5$pvalue, method="fdr", n=length(Psi5$pvalue))
total <- sum(Psi5_lm_df$Absolute_Effect_Estimate)
Psi5$variance_explained <- Psi5$Absolute_Effect_Estimate/total
Psi5$std_error <- summary(Psi5_lm_df)$coefficients[,2]/total
Value <- summary(Psi5_lm_df)$coefficients[,0]
Psi5$Value <- c("Intercept", "Theta_Junctions_Outlier_Status", "Psi5_Junctions_Outlier_Status", "Psi3_Junctions_Outlier_Status", "All_Junctions_Outlier_Status")

Psi3_lm_df <- glm(Psi3_Junctions_Outlier_Status ~ Theta_Junctions_Outlier_Status + Jaccard_Junctions_Outlier_Status + Psi5_Junctions_Outlier_Status+ All_Junctions_Outlier_Status, family=binomial(link='logit'), data = all_outliers)
Psi3 <- data.frame(matrix(nrow=5, ncol=0))
Psi3$Estimate <- summary(Psi3_lm_df)$coefficients[,1] 
Psi3$Absolute_Effect_Estimate <- abs(summary(Psi3_lm_df)$coefficients[,1]) 
Psi3$pvalue <- summary(Psi3_lm_df)$coefficients[,4] 
Psi3$padjust <- p.adjust(Psi3$pvalue, method="fdr", n=length(Psi3$pvalue))
total <- sum(Psi3_lm_df$Absolute_Effect_Estimate)
Psi3$variance_explained <- Psi3$Absolute_Effect_Estimate/total
Psi3$std_error <- summary(Psi3_lm_df)$coefficients[,2]/total
Value <- summary(Psi3_lm_df)$coefficients[,0]
Psi3$Value <- c("Intercept", "Theta_Junctions_Outlier_Status", "Psi5_Junctions_Outlier_Status", "Psi3_Junctions_Outlier_Status", "All_Junctions_Outlier_Status")

All_lm_df <- glm(All_Junctions_Outlier_Status ~ Theta_Junctions_Outlier_Status + Jaccard_Junctions_Outlier_Status + Psi5_Junctions_Outlier_Status+ Psi3_Junctions_Outlier_Status, family=binomial(link='logit'), data = all_outliers)
All <- data.frame(matrix(nrow=5, ncol=0))
All$Estimate <- summary(All_lm_df)$coefficients[,1] 
All$Absolute_Effect_Estimate <- abs(summary(All_lm_df)$coefficients[,1]) 
All$pvalue <- summary(All_lm_df)$coefficients[,4] 
All$padjust <- p.adjust(All$pvalue, method="fdr", n=length(All$pvalue))
total <- sum(All_lm_df$Absolute_Effect_Estimate)
All$variance_explained <- All$Absolute_Effect_Estimate/total
All$std_error <- summary(All_lm_df)$coefficients[,2]/total
Value <- summary(All_lm_df)$coefficients[,0]
All$Value <- c("Intercept", "Theta_Junctions_Outlier_Status", "Psi5_Junctions_Outlier_Status", "Psi3_Junctions_Outlier_Status", "All_Junctions_Outlier_Status")

gene_list <- all_outliers %>% select(sampleID, All_Junctions_Outlier_Status, Theta_Junctions_Outlier_Status, Psi3_Junctions_Outlier_Status, Psi5_Junctions_Outlier_Status, Jaccard_Junctions_Outlier_Status) 
genes <- data.frame(gene_list)
genes$`All Junctions` <- ifelse(genes$All_Junctions_Outlier_Status == 1, TRUE, FALSE)
genes$`Theta Junctions` <- ifelse(genes$Theta_Junctions_Outlier_Status == 1, TRUE, FALSE)
genes$`Psi3 Junctions` <- ifelse(genes$Psi3_Junctions_Outlier_Status == 1, TRUE, FALSE)
genes$`Psi5 Junctions` <- ifelse(genes$Psi5_Junctions_Outlier_Status == 1, TRUE, FALSE)
genes$`Jaccard Junctions` <- ifelse(genes$Jaccard_Junctions_Outlier_Status == 1, TRUE, FALSE)

rownames(genes) <- genes$sampleID
genes <- genes %>% select(-All_Junctions_Outlier_Status, -Psi3_Junctions_Outlier_Status, -Psi5_Junctions_Outlier_Status, -Theta_Junctions_Outlier_Status, -Jaccard_Junctions_Outlier_Status)
#genes <- genes %>% select(-gene_list)
y_lab <- expression("Shared Excess Outlier Samples of Type:")


upset <- upset(genes, colnames(genes)[2:6], name='sampleID', width_ratio=0.1,
    matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=3))
        + scale_color_manual(
            values=c('All Junctions'='#DA786C', 'Theta Junctions'='#88ACAA', 'Psi3 Junctions'='#F2D177', 'Psi5 Junctions'= '#D6BDBF', 'Jaccard Junctions' = '#9F5E64'),
            guide=guide_legend(override.aes=list(shape='circle'))
        )
    ),
    queries=list(
        upset_query(set='All Junctions', fill='#DA786C'),
        upset_query(set='Theta Junctions', fill='#88ACAA'),
        upset_query(set='Psi3 Junctions', fill='#F2D177'),
        upset_query(set='Psi5 Junctions', fill='#D6BDBF'),
        upset_query(set='Jaccard Junctions', fill='#9F5E64')
    ),
    themes=upset_modify_themes(list('intersections_matrix'=theme(axis.text.y=element_text(size=20))))
) +
theme_bw(base_size=25)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  xlab("FRASER and FRASER2 Metrics")+
  ylab(y_lab)
  

upset
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/Outlier_status_comparison/outlier_status_comparison.pdf", plot=upset,  limitsize = FALSE, units = "in", height=15, width=12)

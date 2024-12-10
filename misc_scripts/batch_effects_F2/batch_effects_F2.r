library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

joined_filtered <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_minor_spliceosome/output/output_Nov_18_try5/DataFrames/metadata_counts_outlier_joined.csv")

########################################################################
#######################-----Linear Model-----###########################
########################################################################
make_estimate_df <- function(row, dataframe) {
    Esitmate_lm_df <- lm(dataframe[[row]] ~ RIN + age + sex+ batch, data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="fdr", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
    Value <- summary(Esitmate_lm_df)$coefficients[,0]
    df$Value <- c("Intercept", "RIN", "Age", "Sex", "Batch")
    df$Numbers <- ifelse(df$Value == "Intercept", 5,
                            ifelse(df$Value == "RIN", 4,
                            ifelse(df$Value == "Age", 3,
                            ifelse(df$Value == "Sex", 2, 1))))
   df
}

################-----Linear Model Junctions-----#################
#.005
Jaccard_05_Junctions <- make_estimate_df("Jaccard_Junctions", joined_filtered)
Esitmate <- Jaccard_05_Junctions
Esitmate$Numbers <- ifelse(Esitmate$Value == "Intercept", 5,
                        ifelse(Esitmate$Value == "RIN", 4,
                        ifelse(Esitmate$Value == "Age", 3,
                        ifelse(Esitmate$Value == "Sex", 2, 1))))

total <- sum(Esitmate$Absolute_Effect_Estimate)                 
Esitmate$std_error_percent <- (Esitmate$Absolute_Effect_Estimate/total)*100
Esitmate$variance_explained <- (Esitmate$variance_explained)*100

title=paste("Linear Regression")
Esitmate <- Esitmate %>% filter(Value != "Intercept") %>% arrange(Numbers)

MyColors <- c("#D6BDBF", "#F2D177", "#FFBEAB", "#B4B689")
names(MyColors) <- c("RIN", "Age", "Sex", "Batch")

linear_variance_explained_junction <- ggplot(Esitmate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error_percent,xmax = variance_explained+std_error_percent), position = position_dodge(0.9), width = 0.2)+
  geom_text(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
                          ifelse(padjust < 0.001, "***", 
                            ifelse(padjust < 0.05, "*", "")))), 
            position = position_dodge(width = .5), vjust = 2, size = 30 / .pt, angle= 90) +
  labs(title_lab=title)+
  xlab('Percent Variance Explained') +
  ylab('Metadata Values') +
  theme(axis.text.x = element_text(angle = 0, 
                                   vjust = 0.6, 
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = MyColors)
  
linear_variance_explained_junction
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/batch_effects_F2/linear_regression.pdf", plot=linear_variance_explained_junction,  limitsize = FALSE, units = "in", height=10, width=10)

########################################################################
#####################-----Logisitic Model-----##########################
########################################################################
make_estimate_df <- function(dataframe, column) {
    Esitmate_lm_df <- glm(dataframe[[column]] ~ RIN + age + sex+ batch, family=binomial(link='logit'), data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="fdr", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total
    Value <- summary(Esitmate_lm_df)$coefficients[,0]
    df$Value <- c("Intercept", "RIN", "Age", "Sex", "Batch")
    df$Numbers <- ifelse(df$Value == "Intercept", 5,
                            ifelse(df$Value == "RIN", 4,
                            ifelse(df$Value == "Age", 3,
                            ifelse(df$Value == "Sex", 2, 1))))
   df
}


################-----Linear Model Junctions-----#################
#.005
Jaccard_05_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_Outlier_Status")
Esitmate <- Jaccard_05_Junctions
Esitmate$Numbers <- ifelse(Esitmate$Value == "Intercept", 5,
                        ifelse(Esitmate$Value == "RIN", 4,
                        ifelse(Esitmate$Value == "Age", 3,
                        ifelse(Esitmate$Value == "Sex", 2, 1))))

total <- sum(Esitmate$Absolute_Effect_Estimate)                 
Esitmate$std_error_percent <- (Esitmate$Absolute_Effect_Estimate/total)*100
Esitmate$variance_explained <- (Esitmate$variance_explained)*100

title=paste("Linear Regression")
Esitmate <- Esitmate %>% filter(Value != "Intercept") %>% arrange(Numbers)

MyColors <- c("#D6BDBF", "#F2D177", "#FFBEAB", "#B4B689")
names(MyColors) <- c("RIN", "Age", "Sex", "Batch")

logisitic_variance_explained_junction <- ggplot(Esitmate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error_percent,xmax = variance_explained+std_error_percent), position = position_dodge(0.9), width = 0.2)+
  geom_text(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
                          ifelse(padjust < 0.001, "***", 
                            ifelse(padjust < 0.05, "*", "")))), 
            position = position_dodge(width = .5), vjust = 2, size = 30 / .pt, angle= 90) +
  labs(title_lab=title)+
  xlab('Percent Variance Explained') +
  ylab('Metadata Values') +
  theme(axis.text.x = element_text(angle = 0, 
                                   vjust = 0.6, 
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = MyColors)
  
logisitic_variance_explained_junction
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/batch_effects_F2/logisitic_regression.pdf", plot=logisitic_variance_explained_junction,  limitsize = FALSE, units = "in", height=10, width=10)


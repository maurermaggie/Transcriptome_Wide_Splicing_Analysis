library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
joined_filtered_fp <- args[1]
joined_filtered <- read_csv(joined_filtered_fp)

All_05_stats <- args[2]
Psi3_05_stats <- args[3]
Psi5_05_stats <- args[4]
Theta_05_stats <- args[5]
Jaccard_05_stats <- args[6]
output_var_explained <- args[7]

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

################-----Logistic Model Junctions-----#################
#.005
#All_05_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_Outlier_Status")

Psi3_05_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_Outlier_Status")

Psi5_05_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_Outlier_Status")

Theta_05_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_Outlier_Status")

Jaccard_05_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_Outlier_Status")

write_csv(Psi3_05_Junctions, Psi3_05_stats)
write_csv(Psi5_05_Junctions, Psi5_05_stats)
write_csv(Theta_05_Junctions, Theta_05_stats)
write_csv(Jaccard_05_Junctions, Jaccard_05_stats)

joined_filtered$outlier <- ifelse(joined_filtered$All_Junctions_Outlier_Status == 1, 1,
                            ifelse(joined_filtered$Psi3_Junctions_Outlier_Status == 1, 1,
                             ifelse(joined_filtered$Psi5_Junctions_Outlier_Status == 1, 1,
                             ifelse(joined_filtered$Theta_Junctions_Outlier_Status == 1, 1, 0))))

All_05_Junctions <- make_estimate_df(joined_filtered, "outlier")
write_csv(All_05_Junctions, All_05_stats)

########################################################################
####################-----Variance Explained-----########################
########################################################################

####################-----Variance Explained Junctions-----########################
Estimate <- All_05_Junctions
Estimate$Numbers <- ifelse(Estimate$Value == "Intercept", 5,
                        ifelse(Estimate$Value == "RIN", 4,
                        ifelse(Estimate$Value == "Age", 3,
                        ifelse(Estimate$Value == "Sex", 2, 1))))

Estimate$std_error <- (Estimate$std_error)*100
Estimate$variance_explained <- (Estimate$variance_explained)*100
title=paste("Logistic Regression")
Estimate <- Estimate %>% filter(Value != "Intercept") %>% arrange(Numbers)

####################-----Variance Explained Junctions-----########################
Estimate_j <- Jaccard_05_Junctions
Estimate_j$Numbers <- ifelse(Estimate_j$Value == "Intercept", 5,
                        ifelse(Estimate_j$Value == "RIN", 4,
                        ifelse(Estimate_j$Value == "Age", 3,
                        ifelse(Estimate_j$Value == "Sex", 2, 1))))

Estimate_j$std_error <- (Estimate_j$std_error)*100
Estimate_j$variance_explained <- (Estimate_j$variance_explained)*100
title=paste("Logistic Regression")
Estimate_j <- Estimate_j %>% filter(Value != "Intercept") %>% arrange(Numbers)

Estimate$type <- "FRASER"
Estimate_j$type <- "FRASER2"

Est <- bind_rows(Estimate, Estimate_j)

MyColors <- c("#D6BDBF", "#F2D177")
names(MyColors) <- c("FRASER", "FRASER2")

logistical_variance_explained_junctions <- ggplot(Est, aes(fill = type, x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", position=position_dodge(width = 0.5), width = 0.5, size=0.7) + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.5), width = 0.2)+
  labs(title_lab=title)+
  xlab('Percent Variance Explained') +
  ylab('Metadata Values') +
  theme(axis.text.x = element_text(angle = 0, 
                                   vjust = 0.6, 
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  geom_text(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
                          ifelse(padjust < 0.001, "***", 
                            ifelse(padjust < 0.05, "*", "")))), 
            position = position_dodge(width = .5), vjust = 2, size = 30 / .pt, angle= 90) +
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = MyColors)
  
logistical_variance_explained_junctions
ggsave(filename=output_var_explained, plot=logistical_variance_explained_junctions,  limitsize = FALSE, units = "in", height=10, width=10)
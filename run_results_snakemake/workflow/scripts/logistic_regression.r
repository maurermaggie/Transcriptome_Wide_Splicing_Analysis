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
Specific_Metrics_Var_Explained_fp <- args[8]

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

################-----Logistic Model Junctions-----#################
#.005
All_05_Junctions <- make_estimate_df(joined_filtered, "All_Junctions_Outlier_Status")

Psi3_05_Junctions <- make_estimate_df(joined_filtered, "Psi3_Junctions_Outlier_Status")

Psi5_05_Junctions <- make_estimate_df(joined_filtered, "Psi5_Junctions_Outlier_Status")

Theta_05_Junctions <- make_estimate_df(joined_filtered, "Theta_Junctions_Outlier_Status")

Jaccard_05_Junctions <- make_estimate_df(joined_filtered, "Jaccard_Junctions_Outlier_Status")

write_csv(All_05_Junctions, All_05_stats)
write_csv(Psi3_05_Junctions, Psi3_05_stats)
write_csv(Psi5_05_Junctions, Psi5_05_stats)
write_csv(Theta_05_Junctions, Theta_05_stats)
write_csv(Jaccard_05_Junctions, Jaccard_05_stats)

########################################################################
####################-----Variance Explained-----########################
########################################################################

####################-----Variance Explained Junctions-----########################
Estimate <- All_05_Junctions
Estimate$Numbers <- ifelse(Estimate$Value == "Intercept", 5,
                        ifelse(Estimate$Value == "RIN", 4,
                        ifelse(Estimate$Value == "Age", 3,
                        ifelse(Estimate$Value == "Sex", 2, 1))))

title=paste("Outlier Status Variance Explained", "\n", "(Logistic Regression)")

logistical_variance_explained_junctions <- ggplot(Estimate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.9), width = 0.2)+
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
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal[c(1,3,5,7,9)])
  
logistical_variance_explained_junctions
ggsave(filename=output_var_explained, plot=logistical_variance_explained_junctions,  limitsize = FALSE, units = "in")

################---By junctions type----###################
psi3 <- make_estimate_df(joined_filtered, "Psi3_Junctions_Outlier_Status")
psi3$type <- "Psi3"
psi5 <- make_estimate_df(joined_filtered, "Psi5_Junctions_Outlier_Status")
psi5$type <- "Psi5"
theta <- make_estimate_df(joined_filtered, "Theta_Junctions_Outlier_Status")
theta$type <- "Theta"
jaccard <- make_estimate_df(joined_filtered, "Jaccard_Junctions_Outlier_Status")
jaccard$type <- "Jaccard Index"

all <- bind_rows(psi3, psi5, theta, jaccard)
names(all)[names(all) == 'type'] <- 'Outlier Type'

title=paste("Outlier Status Variance Explained", "\n", "(Logistic Regression)")

logistic_variance_explained_junctions <- ggplot(all, aes(x= variance_explained, y=reorder(Value, Numbers), fill=`Outlier Type`)) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, aes(fill=`Outlier Type`), position="dodge") + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.5), width = 0.2)+
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
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal)+
  theme(legend.title = element_text(face = "bold"))+
  guides(fill=guide_legend(title="Outlier Type"))
  
  
logistic_variance_explained_junctions
ggsave(filename=Specific_Metrics_Var_Explained_fp, plot=logistic_variance_explained_junctions,  limitsize = FALSE, units = "in")
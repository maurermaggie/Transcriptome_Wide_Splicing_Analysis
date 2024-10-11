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
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Logistic_Variance_Explained_Junctions.pdf", plot=logistical_variance_explained_junctions,  limitsize = FALSE, units = "in")

####################-----Variance Explained Genes-----########################
Estimate <- All_05_Genes
Estimate$Numbers <- ifelse(Estimate$Value == "Intercept", 5,
                        ifelse(Estimate$Value == "RIN", 4,
                        ifelse(Estimate$Value == "Age", 3,
                        ifelse(Estimate$Value == "Sex", 2, 1))))

title=paste("Outlier Status Variance Explained", "\n", "(Logistic Regression)")

logistical_variance_explained_genes <- ggplot(Estimate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
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
  
logistical_variance_explained_genes
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Logistic_Variance_Explained_Gene.pdf", plot=logistical_variance_explained_genes,  limitsize = FALSE, units = "in")

################---By gene type----###################
psi3 <- make_estimate_df(joined_filtered, "Psi3_Genes_Outlier_Status")
psi3$type <- "Psi3"
psi5 <- make_estimate_df(joined_filtered, "Psi5_Genes_Outlier_Status")
psi5$type <- "Psi5"
theta <- make_estimate_df(joined_filtered, "Theta_Genes_Outlier_Status")
theta$type <- "Theta"
jaccard <- make_estimate_df(joined_filtered, "Jaccard_Genes_Outlier_Status")
jaccard$type <- "Jaccard Index"

all <- bind_rows(psi3, psi5, theta, jaccard)
names(all)[names(all) == 'type'] <- 'Outlier Type'

title=paste("Outlier Status Variance Explained", "\n", "(Logistic Regression)")

logistic_variance_explained_gene <- ggplot(all, aes(x= variance_explained, y=reorder(Value, Numbers), fill=`Outlier Type`)) +
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
  
  
logistic_variance_explained_gene
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Logistic_Variance_Explained_Gene_Types.pdf", plot=logistic_variance_explained_gene,  limitsize = FALSE, units = "in")

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
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Logistic_Variance_Explained_Junction_Types.pdf", plot=logistic_variance_explained_junctions,  limitsize = FALSE, units = "in")

########################################################################
##########################-----Plots-----###############################
########################################################################

################-----RIN Plot-----#################
RIN_plot_glm <- ggplot(joined_filtered,aes(RIN, All_Junctions_Outlier_Status)) + 
  #labs(title="RIN vs # of Significant Outlier Junctions per Person")+
  geom_point(method='glm', formula= All_Junctions_Outlier_Status~RIN) +
  xlab("RIN") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "glm", se = FALSE)+
  theme_classic(base_size = 25)

RIN_plot_glm

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/RIN/RIN_glm.pdf", plot=RIN_plot_glm,  limitsize = FALSE, units = "in")

################-----Age Plot-----#################
Age_plot_glm <- ggplot(joined_filtered,aes(age, All_Junctions_Outlier_Status)) + 
  geom_point(method='lm', formula= All_Junctions_Outlier_Status~age) +
  xlab("Age") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)

Age_plot_glm

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Age/Age_glm.pdf", plot=Age_plot_glm,  limitsize = FALSE, units = "in")

################-----Sex Plot-----#################
sex_boxplot_glm <- ggplot(joined_filtered, aes(fill = sex ,x=sex, y=All_Junctions_Outlier_Status)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Sex") +
  ylab("Number of Genes with Aberrant Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal[c(5,2)])

sex_boxplot_glm

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Sex/Sex_glm.pdf", plot=sex_boxplot_glm,  limitsize = FALSE, units = "in")

################-----Batch Plot-----#################
batch_boxplot_glm <- ggplot(joined_filtered, aes(fill = batch, x=batch, y=All_Junctions_Outlier_Status)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Batch") +
  ylab("Number of Genes with Aberrant Junctions") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30)) +
  scale_fill_manual(values=pal)

batch_boxplot_glm 

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Batch/Batch_glm.pdf", plot=batch_boxplot_glm,  limitsize = FALSE, units = "in")



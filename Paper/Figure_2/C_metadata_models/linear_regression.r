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
get_p_adjust <- function(dataframe, variable){
            column<-eval(substitute(variable),dataframe, parent.frame())
            get_lm <- lm(column ~ RIN + age + sex+ batch, data = joined_filtered)
            print(summary(get_lm))
            p_values <- summary(get_lm)$coefficients[,4] 
            p_adj <- p.adjust(p_values, method="fdr", n=length(p_values))
            p_adj
}

################-----Linear Model Genes-----#################
#.005
All_05_Genes <- get_p_adjust(joined_filtered, All_Genes)
All_01_Genes <- get_p_adjust(joined_filtered, All_Genes_01)
All_001_Genes <- get_p_adjust(joined_filtered, All_Genes_001)
All_1e_6_Genes <- get_p_adjust(joined_filtered, All_Genes_1e_6)

Psi3_05_Genes <- get_p_adjust(joined_filtered, Psi3_Genes)
Psi3_01_Genes <- get_p_adjust(joined_filtered, Psi3_Genes_01)
Psi3_001_Genes <- get_p_adjust(joined_filtered, Psi3_Genes_001)
Psi3_1e_6_Genes <- get_p_adjust(joined_filtered, Psi3_Genes_1e_6)

Psi5_05_Genes <- get_p_adjust(joined_filtered, Psi5_Genes)
Psi5_01_Genes <- get_p_adjust(joined_filtered, Psi5_Genes_01)
Psi5_001_Genes <- get_p_adjust(joined_filtered, Psi5_Genes_001)
Psi5_1e_6_Genes <- get_p_adjust(joined_filtered, Psi5_Genes_1e_6)

Theta_05_Genes <- get_p_adjust(joined_filtered, Theta_Genes)
Theta_01_Genes <- get_p_adjust(joined_filtered, Theta_Genes_01)
Theta_001_Genes <- get_p_adjust(joined_filtered, Theta_Genes_001)
Theta_1e_6_Genes <- get_p_adjust(joined_filtered, Theta_Genes_1e_6)

Jaccard_05_Genes <- get_p_adjust(joined_filtered, Jaccard_Genes)
Jaccard_01_Genes <- get_p_adjust(joined_filtered, Jaccard_Genes_01)
Jaccard_001_Genes <- get_p_adjust(joined_filtered, Jaccard_Genes_001)
Jaccard_1e_6_Genes <- get_p_adjust(joined_filtered, Jaccard_Genes_1e_6)

################-----Linear Model Junctions-----#################
#.005
All_05_Junctions <- get_p_adjust(joined_filtered, All_Junctions)
All_01_Junctions <- get_p_adjust(joined_filtered, All_Junctions_01)
All_001_Junctions <- get_p_adjust(joined_filtered, All_Junctions_001)
All_1e_6_Junctions <- get_p_adjust(joined_filtered, All_Junctions_1e_6)

Psi3_05_Junctions <- get_p_adjust(joined_filtered, Psi3_Junctions)
Psi3_01_Junctions <- get_p_adjust(joined_filtered, Psi3_Junctions_01)
Psi3_001_Junctions <- get_p_adjust(joined_filtered, Psi3_Junctions_001)
Psi3_1e_6_Junctions <- get_p_adjust(joined_filtered, Psi3_Junctions_1e_6)

Psi5_05_Junctions <- get_p_adjust(joined_filtered, Psi5_Junctions)
Psi5_01_Junctions <- get_p_adjust(joined_filtered, Psi5_Junctions_01)
Psi5_001_Junctions <- get_p_adjust(joined_filtered, Psi5_Junctions_001)
Psi5_1e_6_Junctions <- get_p_adjust(joined_filtered, Psi5_Junctions_1e_6)

Theta_05_Junctions <- get_p_adjust(joined_filtered, Theta_Junctions)
Theta_01_Junctions <- get_p_adjust(joined_filtered, Theta_Junctions_01)
Theta_001_Junctions <- get_p_adjust(joined_filtered, Theta_Junctions_001)
Theta_1e_6_Junctions <- get_p_adjust(joined_filtered, Theta_Junctions_1e_6)

Jaccard_05_Junctions <- get_p_adjust(joined_filtered, Jaccard_Junctions)
Jaccard_01_Junctions <- get_p_adjust(joined_filtered, Jaccard_Junctions_01)
Jaccard_001_Junctions <- get_p_adjust(joined_filtered, Jaccard_Junctions_001)
Jaccard_1e_6_Junctions <- get_p_adjust(joined_filtered, Jaccard_Junctions_1e_6)

########################################################################
##########################-----Plots-----###############################
########################################################################
################-----RIN Plot-----#################
RIN_plot <- ggplot(joined_filtered,aes(RIN, All_Junctions)) + 
  #labs(title="RIN vs # of Significant Outlier Junctions per Person")+
  geom_point(method='lm', formula= All_Junctions~RIN) +
  xlab("RIN") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 25)

RIN_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/RIN/RIN_lm_all.pdf", plot=RIN_plot,  limitsize = FALSE, units = "in")

################-----Age Plot-----#################
Age_plot <- ggplot(joined_filtered,aes(age, All_Junctions)) + 
  geom_point(method='lm', formula= All_Junctions~age) +
  xlab("Age") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Age_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Age_lm_all.pdf", plot=Age_plot,  limitsize = FALSE, units = "in")

################-----Sex Plot-----#################
sex_boxplot <- ggplot(joined_filtered, aes(fill = sex ,x=sex, y=All_Junctions)) + 
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

sex_boxplot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Sex/Sex_lm_all.pdf", plot=sex_boxplot,  limitsize = FALSE, units = "in")

################-----Batch Plot-----#################
joined_filtered$batch_letters <- ifelse(joined_filtered$batch == 20, 'A', 
                            ifelse(joined_filtered$batch == 18, "B",
                            ifelse(joined_filtered$batch == 8, "C",
                            ifelse(joined_filtered$batch == 15, "D",
                            ifelse(joined_filtered$batch == 14, "E",
                            ifelse(joined_filtered$batch == 4, "F",
                            ifelse(joined_filtered$batch == 13, "G",
                            ifelse(joined_filtered$batch == 16, "H",
                            ifelse(joined_filtered$batch == 9, "I",
                            ifelse(joined_filtered$batch ==3, "J",
                            ifelse(joined_filtered$batch == 10, "K", 
                            ifelse(joined_filtered$batch == 21, "L",
                            ifelse(joined_filtered$batch == 22, "M",
                            ifelse(joined_filtered$batch == 23, "N", NA))))))))))))))

batch_boxplot <- ggplot(joined_filtered, aes(fill = batch_letters, x=batch_letters, y=All_Junctions)) + 
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

batch_boxplot 

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Batch_lm_all.pdf", plot=batch_boxplot,  limitsize = FALSE, units = "in")

joined_filtered %>% group_by(batch_letters) %>% tally()
################-----Batch vs RIN Plot-----#################
batch_RIN_boxplot <- ggplot(joined_filtered, aes(fill = batch_letters,x=batch_letters, y=RIN)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Batch") +
  ylab("RIN") +
  #labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30))+
  scale_fill_manual(values=pal)

batch_RIN_boxplot 

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Batch/Batch_RIN.pdf", plot=batch_RIN_boxplot,  limitsize = FALSE, units = "in")

################-----Batch Violin Plot-----#################
batch_violin_plot <- ggplot(joined_filtered, aes(fill = batch_letters, x=batch_letters, y=All_Junctions)) + 
  geom_violin(show.legend = FALSE) +
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

batch_violin_plot 

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Batch/Batch_Violin_Plot.pdf", plot=batch_violin_plot,  limitsize = FALSE, units = "in")

########################################################################
###################-----Variance Explained----##########################
########################################################################

###################-----By Junction----##########################
Esitmate_lm <- lm(All_Junctions ~ RIN + age + sex+ batch, data = joined_filtered)
Esitmate <- summary(Esitmate_lm)$coefficients[,1] %>% as.data.frame
Esitmate <- transpose(Esitmate) %>% unlist
Esitmate <- data.frame(Estiamate=abs(Esitmate), Values = c("Intercept", "RIN", "Age", "Sex", "Batch"))
colnames(Esitmate) <- c("Absolute_Effect_Estimate", "Value")

total <- sum(Esitmate$Absolute_Effect_Estimate)
Esitmate$variance_explained <- Esitmate$Absolute_Effect/total
Esitmate$std_error <- summary(Esitmate_lm)$coefficients[,2]/total
Esitmate$Numbers <- ifelse(Esitmate$Value == "Intercept", 5,
                        ifelse(Esitmate$Value == "RIN", 4,
                        ifelse(Esitmate$Value == "Age", 3,
                        ifelse(Esitmate$Value == "Sex", 2, 1))))
                        
p_values <- summary(Esitmate_lm)$coefficients[,4] 
Esitmate$padjust <- p.adjust(p_values, method="bonferroni", n=length(p_values))
Esitmate$pvalue <- summary(Esitmate_lm)$coefficients[,4]
Esitmate$std_error_percent <- Esitmate$std_error/total

title=paste("Number of Outlier Junctions", "\n", "Variance Explained (Linear Regression)")

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
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal[c(1,3,5,7,9)])
  
linear_variance_explained_junction
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Linear_Variance_Explained_Junction.pdf", plot=linear_variance_explained_junction,  limitsize = FALSE, units = "in")

##############-----By Gene-----##################
Esitmate_lm <- lm(All_Genes ~ RIN + age + sex+ batch, data = joined_filtered)
Esitmate <- summary(Esitmate_lm)$coefficients[,1] %>% as.data.frame
Esitmate <- transpose(Esitmate) %>% unlist
Esitmate <- data.frame(Estiamate=abs(Esitmate), Values = c("Intercept", "RIN", "Age", "Sex", "Batch"))
colnames(Esitmate) <- c("Absolute_Effect_Estimate", "Value")

total <- sum(Esitmate$Absolute_Effect_Estimate)
Esitmate$variance_explained <- Esitmate$Absolute_Effect/total
Esitmate$std_error <- summary(Esitmate_lm)$coefficients[,2]/total
Esitmate$Numbers <- ifelse(Esitmate$Value == "Intercept", 5,
                        ifelse(Esitmate$Value == "RIN", 4,
                        ifelse(Esitmate$Value == "Age", 3,
                        ifelse(Esitmate$Value == "Sex", 2, 1))))
                        
p_values <- summary(Esitmate_lm)$coefficients[,4] 
Esitmate$padjust <- p.adjust(p_values, method="bonferroni", n=length(p_values))
Esitmate$pvalue <- summary(Esitmate_lm)$coefficients[,4]
Esitmate$std_error_percent <- Esitmate$std_error/total

title=paste("Number of Genes with Outlier Junctions", "\n", "Variance Explained (Linear Regression)")

linear_variance_explained_gene <- ggplot(Esitmate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
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
  theme_gray(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal[c(1,3,5,7,9)])
  
linear_variance_explained_gene
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Linear_Variance_Explained_Gene.pdf", plot=linear_variance_explained_gene,  limitsize = FALSE, units = "in")

########################################################################
################---Variance Explained of Specfics----###################
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

################---By gene----###################
psi3 <- make_estimate_df("Psi3_Genes", joined_filtered)
psi3$type <- "Psi3"
psi5 <- make_estimate_df("Psi5_Genes", joined_filtered)
psi5$type <- "Psi5"
theta <- make_estimate_df("Theta_Genes", joined_filtered)
theta$type <- "Theta"
jaccard <- make_estimate_df("Jaccard_Genes", joined_filtered)
jaccard$type <- "Jaccard Index"

all <- bind_rows(psi3, psi5, theta, jaccard)
names(all)[names(all) == 'type'] <- 'Outlier Type'

title=paste("Number of Genes with Outlier Junctions", "\n", "Variance Explained (Linear Regression)")

linear_variance_explained_gene <- ggplot(all, aes(x= variance_explained, y=reorder(Value, Numbers), fill=`Outlier Type`)) +
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
  
  
linear_variance_explained_gene
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Linear_Variance_Explained_Gene_Types.pdf", plot=linear_variance_explained_gene,  limitsize = FALSE, units = "in")

################---By JUNCTION----###################
psi3 <- make_estimate_df("Psi3_Junctions", joined_filtered)
psi3$type <- "Psi3"
psi5 <- make_estimate_df("Psi5_Junctions", joined_filtered)
psi5$type <- "Psi5"
theta <- make_estimate_df("Theta_Junctions", joined_filtered)
theta$type <- "Theta"
jaccard <- make_estimate_df("Jaccard_Junctions", joined_filtered)
jaccard$type <- "Jaccard Index"

all <- bind_rows(psi3, psi5, theta, jaccard)
names(all)[names(all) == 'type'] <- 'Outlier Type'

title=paste("Number of Outlier Junctions", "\n", "Variance Explained (Linear Regression)")

linear_variance_explained_junctions <- ggplot(all, aes(x= variance_explained, y=reorder(Value, Numbers), fill=`Outlier Type`)) +
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
  
  
linear_variance_explained_junctions
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Variance_Explained/Linear_Variance_Explained_Junction_Types.pdf", plot=linear_variance_explained_junctions,  limitsize = FALSE, units = "in")



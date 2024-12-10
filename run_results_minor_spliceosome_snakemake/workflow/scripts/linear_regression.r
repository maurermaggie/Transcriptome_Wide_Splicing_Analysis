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
RIN_fp <- args[7]
Age_fp <- args[8]
Sex_fp <- args[9]
Batch_fp <- args[10]
Batch_v_RIN_fp <- args[11]
Batch_Violin_fp <- args[12]
output_var_explained <- args[13]
Specific_Metrics_Var_Explained_fp <- args[14]

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
All_05_Junctions <- make_estimate_df("All_Junctions", joined_filtered)
Psi3_05_Junctions <- make_estimate_df("Psi3_Junctions", joined_filtered)
Psi5_05_Junctions <- make_estimate_df("Psi5_Junctions", joined_filtered)
Theta_05_Junctions <- make_estimate_df("Theta_Junctions", joined_filtered)
Jaccard_05_Junctions <- make_estimate_df("Jaccard_Junctions", joined_filtered)

write_csv(All_05_Junctions, All_05_stats)
write_csv(Psi3_05_Junctions, Psi3_05_stats)
write_csv(Psi5_05_Junctions, Psi5_05_stats)
write_csv(Theta_05_Junctions, Theta_05_stats)
write_csv(Jaccard_05_Junctions, Jaccard_05_stats)

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

ggsave(filename=RIN_fp, plot=RIN_plot,  limitsize = FALSE, units = "in")

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

ggsave(filename=Age_fp, plot=Age_plot,  limitsize = FALSE, units = "in")

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

ggsave(filename=Sex_fp, plot=sex_boxplot,  limitsize = FALSE, units = "in")

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

ggsave(filename=Batch_fp, plot=batch_boxplot,  limitsize = FALSE, units = "in")

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

ggsave(filename=Batch_v_RIN_fp, plot=batch_RIN_boxplot,  limitsize = FALSE, units = "in")

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

ggsave(filename=Batch_Violin_fp, plot=batch_violin_plot,  limitsize = FALSE, units = "in")

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
Esitmate$padjust <- p.adjust(p_values, method="fdr", n=length(p_values))
Esitmate$pvalue <- summary(Esitmate_lm)$coefficients[,4]
Esitmate$std_error_percent <- Esitmate$std_error/total

Esitmate$std_error_percent <- (Esitmate$std_error_percent)*100
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
ggsave(filename=output_var_explained, plot=linear_variance_explained_junction,  limitsize = FALSE, units = "in", height=10, width=10)

########################################################################
################---Variance Explained of Specfics----###################
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
  geom_text_repel(aes(label = ifelse(between(padjust, 0.001, 0.01), "**", 
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
  theme_gray(base_size = 25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual("legend", values = pal)+
  theme(legend.title = element_text(face = "bold"))+
  guides(fill=guide_legend(title="Outlier Type"))
  
  
linear_variance_explained_junctions
ggsave(filename=Specific_Metrics_Var_Explained_fp, plot=linear_variance_explained_junctions,  limitsize = FALSE, units = "in", height=10, width=12)



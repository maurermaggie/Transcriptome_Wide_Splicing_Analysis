library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(broom)
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

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

#joined_filtered <- joined_filtered %>% filter(affected_status != "Unknown")

########################################################################
#######################-----Linear Model-----###########################
########################################################################
make_estimate_df <- function(row, dataframe) {
    RIN_df <- lm(dataframe[[row]] ~ RIN, data = dataframe) %>% summary
    age_df <- lm(dataframe[[row]] ~ age, data = dataframe) %>% summary
    sex_df <- lm(dataframe[[row]] ~ sex, data = dataframe) %>% summary
    batch_df <- lm(dataframe[[row]] ~ batch, data = dataframe) %>% summary
    #AS_df <- lm(dataframe[[row]] ~ affected_status, data = dataframe) %>% summary

    RIN_r <- RIN_df$r.square 
    age_r <- age_df$r.square 
    sex_r <- sex_df$r.square 
    batch_r <- batch_df$r.square 
    #AS_r <- AS_df$r.squared

    RIN_p <- glance(RIN_df)$p.value
    age_p <- glance(age_df)$p.value
    sex_p <- glance(sex_df)$p.value
    batch_p <- glance(batch_df)$p.value
    #AS_p <- glance(AS_df)$p.value

    df <- data.frame(value = c("RIN", "Age", "Sex", "Batch"), rsq = c(RIN_r, age_r, sex_r, batch_r), p_value=c(RIN_p, age_p, sex_p, batch_p))
    df$p_adjust <- p.adjust(df$p_value, method="fdr")
    df
}

make_estimate_df <- function(row, dataframe) {
    Esitmate_lm_df <- lm(dataframe[[row]] ~ RIN + age + sex+ batch, data = dataframe)
    df <- data.frame(matrix(nrow=5, ncol=0))
    df$Estimate <- summary(Esitmate_lm_df)$coefficients[,1] 
    df$Absolute_Effect_Estimate <- abs(summary(Esitmate_lm_df)$coefficients[,1]) 
    df$pvalue <- summary(Esitmate_lm_df)$coefficients[,4] 
    df$padjust <- p.adjust(df$pvalue, method="fdr", n=length(df$pvalue))
    total <- sum(df$Absolute_Effect_Estimate)
    df$variance_explained <- df$Absolute_Effect_Estimate/total * 100
    df$std_error <- summary(Esitmate_lm_df)$coefficients[,2]/total * 100
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

ggsave(filename=RIN_fp, plot=RIN_plot,  limitsize = FALSE, units = "in", height=10, width=10)

################-----Age Plot-----#################
Age_plot <- ggplot(joined_filtered,aes(age, All_Junctions)) + 
  geom_point(method='lm', formula= All_Junctions~age) +
  xlab("Age") +
  ylab("Number of Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Age_plot

ggsave(filename=Age_fp, plot=Age_plot,  limitsize = FALSE, units = "in", height=10, width=10)

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

ggsave(filename=Sex_fp, plot=sex_boxplot,  limitsize = FALSE, units = "in", height=10, width=10)

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

ggsave(filename=Batch_fp, plot=batch_boxplot,  limitsize = FALSE, units = "in", height=10, width=10)

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

ggsave(filename=Batch_v_RIN_fp, plot=batch_RIN_boxplot,  limitsize = FALSE, units = "in", height=10, width=10)

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

ggsave(filename=Batch_Violin_fp, plot=batch_violin_plot,  limitsize = FALSE, units = "in", height=10, width=10)

########################################################################
###################-----Variance Explained----##########################
########################################################################

###################-----By Junction FRASER----##########################
Esitmate_lm <- lm(All_Junctions ~ RIN + age + sex+ batch, data = joined_filtered)
Esitmate <- summary(Esitmate_lm)$coefficients[,1] %>% as.data.frame
Esitmate <- transpose(Esitmate) %>% unlist
Esitmate <- data.frame(Estiamate=abs(Esitmate), Values = c("Intercept", "RIN", "Age", "Sex", "Batch"))
colnames(Esitmate) <- c("Absolute_Effect_Estimate", "Value")

print("XXXXXXXXXXXXXXX")

Values <- Esitmate %>% select(Value)
Esitmate <- Esitmate %>% select(Absolute_Effect_Estimate)

print("YYYYYYYYYYYYYY")
total <- sum(Esitmate$Absolute_Effect_Estimate)
Esitmate$variance_explained <- Esitmate$Absolute_Effect/total
Esitmate$std_error <- summary(Esitmate_lm)$coefficients[,2]/total

print("ZZZZZZZZZZZZZZZZ")
Esitmate <- bind_cols(Esitmate, Values)
Esitmate$Numbers <- ifelse(Esitmate$Value == "Intercept", 5,
                        ifelse(Esitmate$Value == "RIN", 4,
                        ifelse(Esitmate$Value == "Age", 3,
                        ifelse(Esitmate$Value == "Sex", 2, 1))))

print("ZZZZZZZZZZZZZZZZ111111")                       
p_values <- summary(Esitmate_lm)$coefficients[,4] 
Esitmate$padjust <- p.adjust(p_values, method="fdr", n=length(p_values))
Esitmate$pvalue <- summary(Esitmate_lm)$coefficients[,4]
Esitmate$std_error_percent <- Esitmate$std_error

print("ZZZZZZZZZZZZZZZZ222222")
Esitmate$std_error_percent <- (Esitmate$std_error_percent)*100
Esitmate$variance_explained <- (Esitmate$variance_explained)*100

print("ZZZZZZZZZZZZZZZZ3333333")
title=paste("Linear Regression")
Esitmate <- Esitmate %>% filter(Value != "Intercept") %>% arrange(Numbers)

print("ZZZZZZZZZZZZZZZZ44444AAA")

###################-----By Junction FRASER2----##########################
Esitmate_lm_j <- lm(Jaccard_Junctions ~ RIN + age + sex+ batch, data = joined_filtered)
Esitmate_j <- summary(Esitmate_lm_j)$coefficients[,1] %>% as.data.frame
Esitmate_j <- transpose(Esitmate_j) %>% unlist
Esitmate_j <- data.frame(Estiamate_j=abs(Esitmate_j), Values = c("Intercept", "RIN", "Age", "Sex", "Batch"))
colnames(Esitmate_j) <- c("Absolute_Effect_Estimate", "Value")

print("ZZZZZZZZZZZZZZZZ44444")

total <- sum(Esitmate_j$Absolute_Effect_Estimate)
Esitmate_j$variance_explained <- Esitmate_j$Absolute_Effect/total
Esitmate_j$std_error <- summary(Esitmate_lm_j)$coefficients[,2]/total
Esitmate_j$Numbers <- ifelse(Esitmate_j$Value == "Intercept", 5,
                        ifelse(Esitmate_j$Value == "RIN", 4,
                        ifelse(Esitmate_j$Value == "Age", 3,
                        ifelse(Esitmate_j$Value == "Sex", 2, 1))))

print("ZZZZZZZZZZZZZZZZ5555555")          
p_values <- summary(Esitmate_lm_j)$coefficients[,4] 
Esitmate_j$padjust <- p.adjust(p_values, method="fdr", n=length(p_values))
Esitmate_j$pvalue <- summary(Esitmate_lm_j)$coefficients[,4]
Esitmate_j$std_error_percent <- Esitmate_j$std_error

print("ZZZZZZZZZZZZZZZZ6666666")
Esitmate_j$std_error_percent <- (Esitmate_j$std_error_percent)*100
Esitmate_j$variance_explained <- (Esitmate_j$variance_explained)*100

Esitmate_j <- Esitmate_j %>% filter(Value != "Intercept") %>% arrange(Numbers)

print("ZZZZZZZZZZZZZZZZ7777777")

Esitmate$type <- "FRASER"
Esitmate_j$type <- "FRASER2"

Est <- bind_rows(Esitmate, Esitmate_j)

print("ZZZZZZZZZZZZZZZZ8888888")

#(data = df, aes(x = decades, y = value, fill = variable)
MyColors <- c("#D6BDBF", "#F2D177")
names(MyColors) <- c("FRASER", "FRASER2")

linear_variance_explained_junction <- ggplot(Est, aes(x= variance_explained, y=reorder(Value, Numbers), fill=type)) +
  geom_bar(data=Est, position=position_dodge(width = 0.5), stat = "identity", width = 0.5, size=0.7) + 
  geom_errorbar(aes(xmin = variance_explained-std_error_percent,xmax = variance_explained+std_error_percent), position = position_dodge(0.5), width = 0.2)+
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
  scale_fill_manual(values = MyColors) + 
  xlim(-10, 40)
  
linear_variance_explained_junction
ggsave(filename=output_var_explained, plot=linear_variance_explained_junction,  limitsize = FALSE, units = "in", height=10, width=10)
library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(reshape)
library(GO.db)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")

display.brewer.all
pal <- brewer.pal(12, "Paired")

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################

################-----Get Number by Gene per Person-----#################
all_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2)
all_by_gene <- all_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
all_count <- all_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_all <- setdiff(files$sampleID, all_count$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
all_count <- bind_rows(all_count, missing_all_df)
colnames(all_count) <- c("sampleID", "All Junctions")

################-----Join to Metadata-----#################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch)

colnames(all_count) <- c("RDID", "All Junctions")
all_filter_RD <- filter(all_count, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)
colnames(all_RDID) <- c("sampleID", "All_Junctions", "RIN", "sex", "age", "notes", "affected_status", "batch")

colnames(all_count) <- c("GSS_ID", "All Junctions")
all_filter_GSS <- filter(all_count, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)
colnames(all_GSS) <- c("sampleID", "All_Junctions", "RIN", "sex", "age", "notes", "affected_status", "batch")

colnames(all_count) <- c("UDN_ID", "All Junctions")
all_filter_UDN <- filter(all_count, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)
colnames(all_UDN) <- c("sampleID", "All_Junctions", "RIN", "sex", "age", "notes", "affected_status", "batch")

joined <- bind_rows(all_RDID, all_GSS, all_UDN)

################-----Get Missing RIN-----#################
RIN_missing <- joined %>% filter(is.na(RIN)) %>% pull(sampleID)
metadata_RDID_RIN <- metadata %>% filter(RDID %in% RIN_missing)
metadata_GSS_RIN <- metadata %>% filter(GSS_ID %in% RIN_missing)
metadata_UDN_RIN <- metadata %>% filter(UDN_ID %in% RIN_missing)
missing_RIN <- bind_rows(metadata_GSS_RIN, metadata_RDID_RIN, metadata_UDN_RIN) %>% select(RDID, UDN_ID, GSS_ID)

################-----Get Missing Age-----#################
age_missing <- joined %>% filter(is.na(age)) %>% pull(sampleID)
metadata_RDID_age <- metadata %>% filter(RDID %in% age_missing)
metadata_GSS_age <- metadata %>% filter(GSS_ID %in% age_missing)
metadata_UDN_age <- metadata %>% filter(UDN_ID %in% age_missing)
missing_age<- bind_rows(metadata_RDID_age, metadata_GSS_age, metadata_UDN_age) %>% select(RDID, UDN_ID, GSS_ID)

################-----Get Missing Batch-----#################
batch_missing <- joined %>% filter(is.na(batch)) %>% pull(sampleID)
metadata_RDID_batch <- metadata %>% filter(RDID %in% batch_missing)
metadata_GSS_batch <- metadata %>% filter(GSS_ID %in% batch_missing)
metadata_UDN_batch <- metadata %>% filter(UDN_ID %in% batch_missing)
missing_batch<- bind_rows(metadata_RDID_batch, metadata_GSS_batch, metadata_UDN_batch) %>% select(RDID, UDN_ID, GSS_ID)

################-----Remove Low RIN-----#################
low_RIN <- joined %>% filter(RIN <= 7) %>% select(sampleID, RIN)
write_csv(low_RIN, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

joined_filtered <- joined %>% filter(! sampleID %in% low_RIN$sampleID)

################-----Make New Outlier Category-----#################
three_sd <- mean(joined_filtered$All_Junctions) + sd(joined_filtered$All_Junctions) + sd(joined_filtered$All_Junctions) + sd(joined_filtered$All_Junctions)
two_sd <- mean(joined_filtered$All_Junctions) + sd(joined_filtered$All_Junctions) + sd(joined_filtered$All_Junctions)

joined_filtered$outlier_3sd <- ifelse(joined_filtered$All_Junctions >=three_sd, 1, 0)
joined_filtered$outlier_2sd <- ifelse(joined_filtered$All_Junctions >=two_sd, 1, 0)

########################################################################
##################-----Linear Model and Plots-----######################
########################################################################

################-----Linear Model-----#################
RIN_lm <- lm(All_Junctions ~ RIN + age + sex+ batch, data = joined_filtered)
summary(RIN_lm)
p_values <- summary(RIN_lm)$coefficients[2:5,4] 

p.adjust(p_values, method="bonferroni", n=length(p_values))

################-----RIN Plot-----#################
RIN_plot <- ggplot(joined_filtered,aes(RIN, All_Junctions)) + 
  #labs(title="RIN vs # of Significant Outlier Junctions per Person")+
  geom_point(method='lm', formula= All_Junctions~RIN) +
  xlab("RIN") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 25)

RIN_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/RIN.pdf", plot=RIN_plot,  limitsize = FALSE, units = "in")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Age.pdf", plot=Age_plot,  limitsize = FALSE, units = "in")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Sex.pdf", plot=sex_boxplot,  limitsize = FALSE, units = "in")

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
                            ifelse(joined_filtered$batch == 10, "K", NA)))))))))))

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Batch.pdf", plot=batch_boxplot,  limitsize = FALSE, units = "in")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Batch_RIN.pdf", plot=batch_RIN_boxplot,  limitsize = FALSE, units = "in")

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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Batch_Violin_Plot.pdf", plot=batch_violin_plot,  limitsize = FALSE, units = "in")

################-----Batch Plot 3sd-----#################
number_people <- joined_filtered %>% group_by(batch_letters) %>% tally()
number_outliers_2 <- joined_filtered %>% filter(All_Junctions >= two_sd) %>% group_by(batch_letters) %>% tally
number_outliers <- joined_filtered %>% filter(All_Junctions >= three_sd) %>% group_by(batch_letters) %>% tally
colnames(number_outliers) <- c("batch_letters", "n_outliers_3")
colnames(number_outliers_2) <- c("batch_letters", "n_outliers_2")
batch_stats <- left_join(number_people, number_outliers) %>% left_join(number_outliers_2)
batch_stats$batch_letters <- ifelse(is.na(batch_stats$batch_letters), "Undefined", batch_stats$batch_letters)
batch_stats[is.na(batch_stats)] <- 0
batch_stats$pro_outliers_3sd <- batch_stats$n_outliers_3/batch_stats$n
batch_stats$pro_outliers_2sd <- batch_stats$n_outliers_2/batch_stats$n

batch_box_plot_3sd <- ggplot(batch_stats, aes(fill = batch_letters, x=batch_letters, y=pro_outliers_3sd)) + 
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

batch_box_plot_3sd 

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Batch_3sd.pdf", plot=batch_box_plot_3sd,  limitsize = FALSE, units = "in")

########################################################################
##################-----Linear Model and Plots-----######################
########################################################################

################-----Logistic Model-----#################
Batch_lm_three <- glm(outlier_3sd ~ RIN + age + sex+ batch_letters, data = joined_filtered)
summary(Batch_lm_three)
p_values <- summary(Batch_lm_three)$coefficients[2:13,4] 

p.adjust(p_values, method="bonferroni", n=length(p_values))

Batch_lm_two <- glm(outlier_2sd ~ RIN + age + sex+ batch, data = joined_filtered)
summary(Batch_lm_two)
p_values <- summary(Batch_lm_two)$coefficients[2:5,4] 

p.adjust(p_values, method="bonferroni", n=length(p_values))

################-----RIN Plot-----#################
RIN_plot_glm <- ggplot(joined_filtered,aes(RIN, outlier_2sd)) + 
  #labs(title="RIN vs # of Significant Outlier Junctions per Person")+
  geom_point(method='glm', formula= outlier_2sd~RIN) +
  xlab("RIN") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "glm", se = FALSE)+
  theme_classic(base_size = 25)

RIN_plot_glm

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/RIN_glm.pdf", plot=RIN_plot_glm,  limitsize = FALSE, units = "in")

################-----Age Plot-----#################
Age_plot_glm <- ggplot(joined_filtered,aes(age, outlier_2sd)) + 
  geom_point(method='lm', formula= outlier_2sd~age) +
  xlab("Age") +
  ylab("Number of Genes with Significant Splicing Outliers") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)

Age_plot_glm

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Age_glm.pdf", plot=Age_plot_glm,  limitsize = FALSE, units = "in")

################-----Sex Plot-----#################
sex_boxplot_glm <- ggplot(joined_filtered, aes(fill = sex ,x=sex, y=outlier_2sd)) + 
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Sex_glm.pdf", plot=sex_boxplot_glm,  limitsize = FALSE, units = "in")

################-----Batch Plot-----#################
batch_boxplot_glm <- ggplot(joined_filtered, aes(fill = batch_letters, x=batch_letters, y=outlier_2sd)) + 
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure2/Batch_glm.pdf", plot=batch_boxplot_glm,  limitsize = FALSE, units = "in")
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
#####################-----Arrange Dataframe-----########################
########################################################################
filter_genes <- function(dataframe, filter, column_name, significance){
      dataframe_filtered <- dataframe %>% filter(padjust <= significance) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
      if (filter == "psi3") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi3")
        print("filtering by psi3")
      } else if (filter == "psi5") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi5")
        print("filtering by psi5")
      } else if (filter == "theta") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "theta")
        print("filtering by theta")
      } else {
         dataframe_filtered < dataframe_filtered
      }
      dataframe_by_gene <- dataframe_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
      dataframe_count <- dataframe_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

filter_junctions <- function(dataframe, filter, column_name, significance){
      dataframe_filtered <- dataframe %>% filter(padjust <= significance) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
      if (filter == "psi3") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi3")
        print("filtering by psi3")
      } else if (filter == "psi5") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "psi5")
        print("filtering by psi5")
      } else if (filter == "theta") {
        dataframe_filtered <- dataframe_filtered %>% filter(type == "theta")
        print("filtering by theta")
      } else {
         dataframe_filtered < dataframe_filtered
      }
      dataframe_count <- dataframe_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
      missing_dataframe <- setdiff(files$sampleID, dataframe_count$sampleID)
      missing_dataframe_df <- data.frame(sampleID = missing_dataframe, n = 0)
      dataframe_count <- bind_rows(dataframe_count, missing_dataframe_df)
      colnames(dataframe_count) <- c("sampleID", column_name)
      dataframe_count <- dataframe_count %>% filter(! sampleID %in% low_RIN$sampleID)
}

################-----Get Number by Junction per Person-----#################
psi3_count <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions", 0.05)
psi5_count <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions", 0.05)
theta_count <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions", 0.05)
all_count <- filter_junctions(all_uncompiled, "all", "All_Junctions", 0.05)

psi3_count_001 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_001", .001)
psi5_count_001 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_001", .001)
theta_count_001 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_001", .001)
all_count_001 <- filter_junctions(all_uncompiled, "all", "All_Junctions_001", .001)

psi3_count_01 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_01", .01)
psi5_count_01 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_01", .01)
theta_count_01 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_01", .01)
all_count_01 <- filter_junctions(all_uncompiled, "all", "All_Junctions_01", .01)

psi3_count_1e_6 <- filter_junctions(all_uncompiled, "psi3", "Psi3_Junctions_1e_6", 1e-6)
psi5_count_1e_6 <- filter_junctions(all_uncompiled, "psi5", "Psi5_Junctions_1e_6", 1e-6)
theta_count_1e_6 <- filter_junctions(all_uncompiled, "theta", "Theta_Junctions_1e_6", 1e-6)
all_count_1e_6 <- filter_junctions(all_uncompiled, "all", "All_Junctions_1e_6", 1e-6)

all_count_junctions <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(all_count) %>% left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% 
              left_join(theta_count_001) %>% left_join(all_count_001) %>% left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% 
              left_join(all_count_1e_6) %>% left_join(psi3_count_01) %>% left_join(psi5_count_01) %>% left_join(theta_count_01) %>% left_join(all_count_01)

################-----Get Number by Genes per Person-----#################
psi3_count_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes", 0.05)
psi5_count_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes", 0.05)
theta_count_genes<- filter_genes(all_uncompiled, "theta", "Theta_Genes", 0.05)
all_count_genes <- filter_genes(all_uncompiled, "all", "All_Genes", 0.05)

psi3_count_001_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_001", .001)
psi5_count_001_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_001", .001)
theta_count_001_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_001", .001)
all_count_001_genes <- filter_genes(all_uncompiled, "all", "All_Genes_001", .001)

psi3_count_01_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_01", .01)
psi5_count_01_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_01", .01)
theta_count_01_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_01", .01)
all_count_01_genes <- filter_genes(all_uncompiled, "all", "All_Genes_01", .01)

psi3_count_1e_6_genes <- filter_genes(all_uncompiled, "psi3", "Psi3_Genes_1e_6", 1e-6)
psi5_count_1e_6_genes <- filter_genes(all_uncompiled, "psi5", "Psi5_Genes_1e_6", 1e-6)
theta_count_1e_6_genes <- filter_genes(all_uncompiled, "theta", "Theta_Genes_1e_6", 1e-6)
all_count_1e_6_genes <- filter_genes(all_uncompiled, "all", "All_Genes_1e_6", 1e-6)

all_count_genes <- left_join(psi3_count_genes, psi5_count_genes) %>% left_join(theta_count_genes) %>% left_join(all_count_genes) %>% left_join(psi3_count_001_genes) %>% left_join(psi5_count_001_genes) %>% 
              left_join(theta_count_001_genes) %>% left_join(all_count_001_genes) %>% left_join(psi3_count_1e_6_genes) %>% left_join(psi5_count_1e_6_genes) %>% left_join(theta_count_1e_6_genes) %>% 
              left_join(all_count_1e_6_genes) %>% left_join(psi3_count_01_genes) %>% left_join(psi5_count_01_genes) %>% left_join(theta_count_01_genes) %>% left_join(all_count_01_genes)

################-----Combined Junctions to Genes-----#################
all_count <- left_join(all_count_junctions, all_count_genes)

################-----Join to Metadata-----#################
metadata <- metadata %>% select(GSS_ID, UDN_ID, RDID, RIN, sex, age, notes, affected_status, batch)

names(all_count)[names(all_count) == 'sampleID'] <- 'RDID'
all_filter_RD <- filter(all_count, grepl("RD", RDID))
all_RDID <- left_join(all_filter_RD, metadata) %>% filter(!is.na(RDID)) %>% select(-GSS_ID, -UDN_ID)

names(all_count)[names(all_count) == 'RDID'] <- 'GSS_ID'
all_filter_GSS <- filter(all_count, grepl("GSS", GSS_ID))
all_GSS <- left_join(all_filter_GSS, metadata) %>% filter(!is.na(GSS_ID)) %>% select(-RDID, -UDN_ID)

names(all_count)[names(all_count) == 'GSS_ID'] <- 'UDN_ID'
all_filter_UDN <- filter(all_count, grepl("UDN", UDN_ID))
all_UDN <- left_join(all_filter_UDN, metadata) %>% filter(!is.na(UDN_ID)) %>% select(-RDID, -GSS_ID)

names(all_count)[names(all_count) == 'UDN_ID'] <- 'sampleID'

joined_filtered <- bind_rows(all_RDID, all_GSS, all_UDN)

########################################################################
#######################-----Linear Model-----###########################
########################################################################
get_p_adjust <- function(dataframe, variable){
            column<-eval(substitute(variable),dataframe, parent.frame())
            get_lm <- lm(column ~ RIN + age + sex+ batch, data = joined_filtered)
            print(summary(get_lm))
            p_values <- summary(get_lm)$coefficients[,4] 
            p_adj <- p.adjust(p_values, method="bonferroni", n=length(p_values))
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
####################-----Variance Explained-----########################
########################################################################
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
p.adjust(p_values, method="bonferroni", n=length(p_values))

title=paste("Number of Genes with Outlier Junctions", "\n", "Variance Explained (Linear Regression)")

linear_variance_explained <- ggplot(Esitmate, aes(fill = reorder(Value, Numbers), x= variance_explained, y=reorder(Value, Numbers))) +
  geom_bar(stat = "identity", width = 0.5, size=0.7, show.legend=FALSE) + 
  geom_errorbar(aes(xmin = variance_explained-std_error,xmax = variance_explained+std_error), position = position_dodge(0.9), width = 0.2)+
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
  
linear_variance_explained
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/C_metadata_models/Linear_Variance_Explained.pdf", plot=logistical_variance_explained,  limitsize = FALSE, units = "in")


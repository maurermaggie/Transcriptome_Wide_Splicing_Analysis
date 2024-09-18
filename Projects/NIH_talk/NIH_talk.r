#library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
filepath_blood <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/theta_juncs_blood.pdf"
filepath_fibro <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/theta_juncs_fibro.pdf"
filepath_subset <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/theta_juncs_subset.pdf"
filepath_counts <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/total_counts.pdf"
filepath_counts_per <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/per_person_counts.pdf"
filepath_RD268 <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/RD268_counts.pdf"

fibro <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/fibroblasts/outputs_6_26/compiled_counts/FRASER_results_raw.csv")
subset <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/FRASER_results_raw.csv")

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
subset_uncompiled <- read_csv("/home/maurertm/smontgom/maurertm/FRASER_analysis/outputs_old/FRASER_Output_Apr_1/FRASER_output_2024-04-01_14:51.csv")
################-----Number Theta Blood-----#################
y_lab <- "Number of MIGs with intron retention events"
title_lab <- paste("Number of MIGs with intron retention events in the", "\n", "UDN and GREGoR Stanford Sites", sep=" ")

joined$numbers <- rownames(joined)

top_10 <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(numbers, no_MIGs_theta_juncs) %>% head(3) %>% pull(numbers)
g1 <- subset(joined, numbers %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood <- ggplot(joined, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample_ID") +
  geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)-sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_MIGs_theta_juncs)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) - sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_plot_blood

ggsave(filename=filepath_blood, plot=function_plot_blood,  limitsize = FALSE, units = "in")

################-----Number Theta Fibro-----#################
y_lab <- "Number of MIGs with intron retention events"
title_lab <- paste("Number of MIGs with intron retention events in the", "\n", "UDN and GREGoR Stanford Sites", sep=" ")

fibro$numbers <- rownames(fibro)

g1 <- subset(fibro, Sample_ID == "UDN916111")

#check which distribution this is considering
#give a z-score

function_plot_fibro <- ggplot(fibro, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample_ID") +
  geom_hline(data= fibro, aes(yintercept=mean(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=fibro, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=fibro, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=fibro, aes(yintercept=mean(no_MIGs_theta_juncs)-sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=fibro, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = fibro, aes(label = "Mean", x = 25, y=mean(no_MIGs_theta_juncs)+0.03),  size=2) + 
  geom_text(data = fibro, aes(label = "Mean + SD", x = 25, y=mean(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs)+0.03),  size=2) +
  geom_text(data = fibro, aes(label = "Mean + 2*SD", x = 25, y=mean(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs)+0.03),  size=2) +
  geom_text(data = fibro, aes(label = "Mean + 3*SD", x = 25, y=mean(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs) + sd(no_MIGs_theta_juncs)+0.03),  size=2) +
  geom_text(data = fibro, aes(label = "Mean - SD", x = 25, y=mean(no_MIGs_theta_juncs) - sd(no_MIGs_theta_juncs)+0.03),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_plot_fibro

ggsave(filename=filepath_fibro, plot=function_plot_fibro,  limitsize = FALSE, units = "in")

################-----Number Subset Blood-----#################
y_lab <- "Number of MIGs with intron retention events"
title_lab <- paste("Number of MIGs with intron retention events in a subset of the", "\n", "UDN and GREGoR Stanford Sites", sep=" ")

subset$numbers <- rownames(subset)

g1 <- subset(subset, Sample_ID %in% c("RD268"))

#check which distribution this is considering
#give a z-score

function_plot_subset <- ggplot(subset, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Sample_ID") +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample_ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_plot_subset

ggsave(filename=filepath_subset, plot=function_plot_subset,  limitsize = FALSE, units = "in")

################-----Compare Subset vs Full RD268-----#################
total_subset <- nrow(filter(subset_uncompiled, sampleID == "RD268"))
total_full <- nrow(filter(all_uncompiled, sampleID == "RD268"))

number_subset <- nrow(subset)
number_all <- nrow(joined)

data_RD268 <- data.frame(Number_of_Aberrant_Junctions = c(as.numeric(total_subset), as.numeric(total_full)), Data_Type = c("Proband 44 in Sample Size of 32 Particpants", "Proband 44 in Sample Size of 407 participants"),
                Size_of_Data = c(number_subset, number_all))

title_lab <- paste("Comparing the Total Number of Aberrat Junctions", "\n", "for Proband 44 when analysed in full and subsetted datasets", sep= " ")
function_RD268 <- ggplot(data_RD268, aes(Data_Type, Number_of_Aberrant_Junctions, fill = Data_Type)) +     
                    geom_col(position = 'dodge') +
                    geom_line()+
                    geom_point()+
                    xlab("Sample Sizes") +
                    ylab("Number of Total Aberrant Junctions Detected in Proband 44") +
                    labs(title= title_lab) +
                    theme_update(plot.title = element_text(hjust = 0.5))

ggsave(filename=filepath_RD268, plot=function_RD268,  limitsize = FALSE, units = "in")

################-----Compare Subset vs Full-----#################
total_subset <- nrow(subset_uncompiled)
total_full <- nrow(all_uncompiled)

number_subset <- nrow(subset)
number_all <- nrow(joined)

data_summary <- data.frame(Number_of_Aberrant_Junctions = c(as.numeric(total_subset), as.numeric(total_full)), Data_Type = c("32 Particpants", "407 participants"),
                Size_of_Data = c(number_subset, number_all))

title_lab <- paste("The Total Number of Aberrant Junctions for a Sample Size of 407 Participants", "\n", "is over 2920 Times Larger than the", "\n", "Total Number of Aberrant Junction for a Sample Size of 32 Participants", sep= " ")
function_counts <- ggplot(data_summary, aes(Data_Type, Number_of_Aberrant_Junctions, fill = Data_Type)) +     
                    geom_col(position = 'dodge') +
                    geom_line()+
                    geom_point()+
                    xlab("Sample Populations") +
                    ylab("Number of Total Aberrant Junctions Across All Samples") +
                    labs(title= title_lab) +
                    theme_update(plot.title = element_text(hjust = 0.5))

ggsave(filename=filepath_counts, plot=function_counts,  limitsize = FALSE, units = "in")

################-----Compare Number of Events per Individual-----#################
title_lab <- paste("Comparison of the Number of Aberrant Splicing Events", "\n", "Per Sample in a Full vs Subsetted Dataset", sep=" ")

total_subset <- subset_uncompiled %>% group_by(sampleID) %>% tally()
total_subset$Type <- "Subset_32_Probands"
total_full <- all_uncompiled %>% group_by(sampleID) %>% tally()
total_full$Type <- "Full_407_Probands"

combined <- bind_rows(total_subset, total_full)
colnames(combined) <- c("sampleID", "Number_of_Aberrant_Junctions", "Type")

function_counts_per <- ggplot(combined, aes(x=Type, y=Number_of_Aberrant_Junctions)) + 
  geom_boxplot() +
  xlab("Sample Populations") +
  ylab("Number of Total Aberrant Junctions Per Sample") +
  labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5))

ggsave(filename=filepath_counts_per, plot=function_counts_per,  limitsize = FALSE, units = "in")

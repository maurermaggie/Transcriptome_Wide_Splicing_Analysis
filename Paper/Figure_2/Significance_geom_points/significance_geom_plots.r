library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

########################################################################
######################-----By Gene per Person-----######################
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


################-----Get Number by Gene per Person-----#################
#0.05
psi3_count <- filter_genes(all_uncompiled, "psi3", "Psi3", 0.05)
psi5_count <- filter_genes(all_uncompiled, "psi5", "Psi5", 0.05)
theta_count <- filter_genes(all_uncompiled, "theta", "Theta", 0.05)

#1e-6
psi3_count_1e_6 <- filter_genes(all_uncompiled, "psi3", "Psi3_1e_6", 1e-6)
psi5_count_1e_6 <- filter_genes(all_uncompiled, "psi5", "Psi5_1e_6", 1e-6)
theta_count_1e_6 <- filter_genes(all_uncompiled, "theta", "Theta_1e_6", 1e-6)

#.001
psi3_count_001 <- filter_genes(all_uncompiled, "psi3", "Psi3_001", .001)
psi5_count_001 <- filter_genes(all_uncompiled, "psi5", "Psi5_001", .001)
theta_count_001 <- filter_genes(all_uncompiled, "theta", "Theta_001", .001)

#join
joined <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(psi3_count_1e_6) %>% left_join(psi5_count_1e_6) %>% left_join(theta_count_1e_6) %>% 
                  left_join(psi3_count_001) %>% left_join(psi5_count_001) %>% left_join(theta_count_001)

########################################################################
######################-----Significance Plots-----######################
########################################################################

################-----Psi3 Significance Plot-----#################
title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_psi3 <- ggplot(joined, aes(x = reorder(sampleID, Psi3), y = Psi3)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Psi3, x = reorder(sampleID, Psi3))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Psi3_1e_6, x = reorder(sampleID, Psi3_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Psi3_001, x = reorder(sampleID, Psi3_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Psi3)+sd(Psi3)+sd(Psi3), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Psi3)+sd(Psi3)+sd(Psi3)+sd(Psi3), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3) + sd(joined$Psi3)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_psi3

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/Significance_geom_points/significance_psi3.pdf", plot=significance_plot_psi3,  limitsize = FALSE, units = "in")

################-----Psi5 Significance Plot-----#################
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_psi5 <- ggplot(joined, aes(x = reorder(sampleID, Psi5), y = Psi5)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Psi5, x = reorder(sampleID, Psi5))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Psi5_1e_6, x = reorder(sampleID, Psi5_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Psi5_001, x = reorder(sampleID, Psi5_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Psi5)+sd(Psi5)+sd(Psi5), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Psi5)+sd(Psi5)+sd(Psi5)+sd(Psi5), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5) + sd(joined$Psi5)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_psi5

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/Significance_geom_points/significance_psi5.pdf", plot=significance_plot_psi5,  limitsize = FALSE, units = "in")

################-----Theta Significance Plot-----#################
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
significance_plot_theta <- ggplot(joined, aes(x = reorder(sampleID, Theta), y = Theta)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red", aes(y=Theta, x = reorder(sampleID, Theta))) +
  geom_point(shape=21, alpha=.6, size=5, fill="blue", aes(y=Theta_1e_6, x = reorder(sampleID, Theta_1e_6))) +
  geom_point(shape=21, alpha=.6, size=5, fill="green", aes(y=Theta_001, x = reorder(sampleID, Theta_001))) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(Theta)+sd(Theta)+sd(Theta), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(Theta)+sd(Theta)+sd(Theta)+sd(Theta), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$Theta) + sd(joined$Theta) + sd(joined$Theta) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$Theta) + sd(joined$Theta) + sd(joined$Theta) + sd(joined$Theta)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

significance_plot_theta

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_2/Significance_geom_points/significance_theta.pdf", plot=significance_plot_theta,  limitsize = FALSE, units = "in")

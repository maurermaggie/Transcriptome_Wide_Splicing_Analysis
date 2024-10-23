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
      dataframe_filtered <- dataframe %>% filter(padjust < significance) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
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
all_count <- filter_genes(all_uncompiled, "psi3", "All", 0.05)
psi3_count <- filter_genes(all_uncompiled, "psi3", "Psi3", 0.05)
psi5_count <- filter_genes(all_uncompiled, "psi5", "Psi5", 0.05)
theta_count <- filter_genes(all_uncompiled, "theta", "Theta", 0.05)

#join
joined <- left_join(psi3_count, psi5_count) %>% left_join(theta_count) %>% left_join(all_count)

########################################################################
#########################-----By Gene Plots-----########################
########################################################################

################-----Psi3 Plot-----#################
joined <- psi3_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi3 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 3),  size=8, hjust = -.00001, vjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+3),  size=8, hjust = -.00001, vjust=0.000025) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi3

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_gene/psi3_by_gene.pdf", plot=function_plot_psi3,  limitsize = FALSE, units = "in")

################-----Psi3 Plot-----#################
joined <- psi5_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_psi5 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_psi5

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_gene/psi5_by_gene.pdf", plot=function_plot_psi5,  limitsize = FALSE, units = "in")

################-----Theta Plot-----#################
joined <- theta_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- expression(theta)
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
function_plot_theta <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001, vjust = -.01) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001, vjust =-.01) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_theta

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_gene/theta_by_gene.pdf", plot=function_plot_theta,  limitsize = FALSE, units = "in")

################-----All Plot-----#################
joined <- all_count
colnames(joined) <- c("sampleID", "n")

joined$numbers <- rownames(joined)
title_lab <- "All Genes with Significant Junctions"
x_lab <- paste("Samples ordered by Number of Genes with", "\n",
                     "Significant Outlier Events")
y_lab <- "Number of Genes with Significant Outlier Events"

#check which distribution this is considering
#give a z-score
function_plot_all <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + 6),  size=8, hjust = -.00001, vjust = -.01) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+6),  size=8, hjust = -.00001, vjust =-.01) +
  theme_classic(base_size = 27)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_all

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/excess_junctions/by_gene/all_by_gene.pdf", plot=function_plot_all,  limitsize = FALSE, units = "in")

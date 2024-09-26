library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(GO.db)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

################-----Get Number by Gene per Person-----#################
all_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2)
all_by_gene <- all_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
all_count <- all_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_all <- setdiff(files$sampleID, all_count$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
all_count <- bind_rows(all_count, missing_all_df)
all_count <- all_count %>% filter(! sampleID %in% low_RIN$sampleID)

psi3_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi3")
psi3_by_gene <- psi3_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
psi3_count <- psi3_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi3 <- setdiff(files$sampleID, psi3_count$sampleID)
missing_psi3_df <- data.frame(sampleID = missing_psi3, n = 0)
psi3_count <- bind_rows(psi3_count, missing_psi3_df)
psi3_count <- psi3_count %>% filter(! sampleID %in% low_RIN$sampleID)

psi5_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "psi5")
psi5_by_gene <- psi5_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
psi5_count <- psi5_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi5 <- setdiff(files$sampleID, psi5_count$sampleID)
missing_psi5_df <- data.frame(sampleID = missing_psi5, n = 0)
psi5_count <- bind_rows(psi5_count, missing_psi5_df)
psi5_count <- psi5_count %>% filter(! sampleID %in% low_RIN$sampleID)

theta_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) > 0.3) %>% filter(abs(zScore) > 2) %>% filter(type == "theta")
theta_by_gene <- theta_filtered %>% select(sampleID, hgncSymbol) %>% unique %>% filter(!is.na(hgncSymbol))
theta_count <- theta_by_gene %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_theta <- setdiff(files$sampleID, theta_count$sampleID)
missing_theta_df <- data.frame(sampleID = missing_theta, n = 0)
theta_count <- bind_rows(theta_count, missing_theta_df)
theta_count <- theta_count %>% filter(! sampleID %in% low_RIN$sampleID)

################-----Get Number by Gene per Person-----#################
joined <- all_count

title_lab <- "Number of Genes with Aberrant Splicing"
x_lab <- paste("Samples Ordered by Number of", "\n", "Genes with Significant Outlier Events")

all_by_gene_plot <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
  xlab(x_lab) +
  ylab("Number of Genes with Significant Outlier Events") +
  labs(title=title_lab, size=12)  +
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

all_by_gene_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/excess_junctions/by_gene_measures/all_by_gene.pdf", plot=all_by_gene_plot,  limitsize = FALSE, units = "in")

################-----Get Number by Gene per Person Psi3-----#################
joined <- psi3_count

title_lab <- expression("Number of Genes with" ~ psi ~ 3 ~ "Aberrant Splicing")
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 3 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 3 ~ "Events")

#check which distribution this is considering
#give a z-score
psi3_by_gene_plot <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
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

psi3_by_gene_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/excess_junctions/by_gene_measures/psi3_by_gene.pdf", plot=psi3_by_gene_plot,  limitsize = FALSE, units = "in")

################-----Get Number by Gene per Person Psi5-----#################
joined <- psi5_count

title_lab <- expression("Number of Genes with" ~ psi ~ 5 ~ "Aberrant Splicing")
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ psi ~ 5 ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ Psi ~ 5 ~ "Events")

#check which distribution this is considering
#give a z-score
psi5_by_gene_plot <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
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

psi5_by_gene_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/excess_junctions/by_gene_measures/psi5_by_gene.pdf", plot=psi5_by_gene_plot,  limitsize = FALSE, units = "in")

################-----Get Number by Gene per Person Theta-----#################
joined <- theta_count

title_lab <- expression("Number of Genes with" ~ theta ~ "Aberrant Splicing")
x_lab <- expression(atop("Samples ordered by Number of Genes with",
                     "Significant" ~ theta ~ "Events"))
y_lab <- expression("Number of Genes with Significant" ~ theta ~ "Events")

#check which distribution this is considering
#give a z-score
theta_by_gene_plot <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
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

theta_by_gene_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/excess_junctions/by_gene_measures/theta_by_gene.pdf", plot=theta_by_gene_plot,  limitsize = FALSE, units = "in")

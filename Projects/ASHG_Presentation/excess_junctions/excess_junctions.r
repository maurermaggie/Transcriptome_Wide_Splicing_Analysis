#library(EnsDb.Hsapiens.v79)
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
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

################-----Filter All Junctions-----#################
all_filtered <- all_uncompiled %>% filter(padjust <=0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
all_count <- all_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_all <- setdiff(files$sampleID, all_filtered$sampleID)
missing_all_df <- data.frame(sampleID = missing_all, n = 0)
all_count <- bind_rows(all_count, missing_all_df)
all <- all_count %>% filter(! sampleID %in% low_RIN$sampleID)

psi5_filtered <- all_filtered %>% filter(type == "psi5") %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
psi5_count <- psi5_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi5 <- setdiff(files$sampleID, psi5_filtered$sampleID)
missing_psi5_df <- data.frame(sampleID = missing_psi5, n = 0)
psi5_count <- bind_rows(psi5_count, missing_psi5_df)
tally_psi5 <- psi5_count %>% filter(! sampleID %in% low_RIN$sampleID)

psi3_filtered <- all_filtered %>% filter(type == "psi3") %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
psi3_count <- psi3_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi3 <- setdiff(files$sampleID, psi3_filtered$sampleID)
missing_psi3_df <- data.frame(sampleID = missing_psi3, n = 0)
psi3_count <- bind_rows(psi3_count, missing_psi3_df)
tally_psi3 <- psi3_count %>% filter(! sampleID %in% low_RIN$sampleID)



tally_psi3 <- all_filtered %>% filter(type == "psi3") %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
tally_theta <- all_filtered %>% filter(type == "theta") %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
all <- all_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))

##### ADD IN PEOPLE WITH ZEROS!!!!!

################-----Psi5-----#################
joined <- tally_psi5
title_lab <- paste("Number of Psi5 events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

psi5 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point() +
  xlab("Samples Ordered by Number of Psi5 Events") +
  ylab("Number of Psi5 Events") +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= sampleID)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, n), y = n, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

psi5

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/psi5_excess.pdf", plot=psi5,  limitsize = FALSE, units = "in")

################-----Psi3-----#################
joined <- tally_psi3
title_lab <- paste("Number of Psi3 events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

psi3 <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point() +
  xlab("Samples Ordered by Number of Psi3 Events") +
  ylab("Number of Psi3 Events") +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= sampleID)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, n), y = n, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

psi3

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/psi3_excess.pdf", plot=psi3,  limitsize = FALSE, units = "in")

################-----Theta-----#################
joined <- tally_theta
title_lab <- paste("Number of Theta events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

theta <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point() +
  xlab("Samples Ordered by Number of Theta Events") +
  ylab("Number of Theta Events") +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= sampleID)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, n), y = n, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

theta

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/theta_excess.pdf", plot=theta,  limitsize = FALSE, units = "in")

################-----All-----#################
joined <- all
title_lab <- paste("Number of Significant events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

all_figure <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point() +
  xlab("Samples Ordered by Number of Significant Events") +
  ylab("Number of Significant Events") +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= sampleID)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, n), y = n, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

all_figure

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/all_excess.pdf", plot=all_figure,  limitsize = FALSE, units = "in")

################-----Theta Point-----#################
joined <- tally_theta
title_lab <- paste("Number of Theta events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")
g1 <- filter(joined, sampleID %in% c("RD268", "GSS225379"))
top_10 <- g1 %>% pull(sampleID)

theta_RNU4ATAC <- ggplot(joined, aes(x = reorder(sampleID, n), y = n)) +
  geom_point() +
  xlab("Samples Ordered by Number of Theta Events") +
  ylab("Number of Theta Events") +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= sampleID)) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(sampleID %in% top_10)), size=3, aes(x = reorder(sampleID, n), y = n, label = sampleID), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

theta_RNU4ATAC

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/ASHG_Presentation/theta_excess_RNU4ATAC.pdf", plot=theta_RNU4ATAC,  limitsize = FALSE, units = "in")


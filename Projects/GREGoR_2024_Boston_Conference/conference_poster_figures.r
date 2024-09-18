library(biomaRt)
library(tidyverse)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
BiocManager::install("compbio4all")
BiocManager::install("rentrez")
BiocManager::install("seqinr")
BiocManager::install("msa")

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
filepath_blood_RNU6ATAC <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/minor_intron_events_RNU6ATAC.pdf"
filepath_minor_introns_box_RNU6ATAC <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/minor_introns_box_RNU6ATAC.pdf"

################-----Number Minor Intron Retention Events-----#################
title_lab <- paste("Number of Minor Intron Retention Events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% filter(Sample_ID %in% c("RD268", "GSS225379"))
RNU6ATAC <- joined %>% filter(Sample_ID %in% c("RD380"))
normal <- joined %>% filter(! Sample_ID %in% c("RD380", "GSS225379", "RD268"))

RNU4ATAC$type <- "RNU4ATAC-opathy"
RNU6ATAC$type <- "Biallelic rare variants in RNU6ATAC"
normal$type <- NA

joined <- bind_rows(RNU4ATAC, RNU6ATAC, normal)

top_10 <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(Sample_ID, no_MIGs_theta_juncs, type) %>% head(3) %>% pull(type)
g1 <- subset(joined, type %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood <- ggplot(joined, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Samples Ordered by Number of Minor Intron Retention Events") +
  ylab("Number of Minor Intron Retention Events") +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= type)) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = type), position=position_dodge(width=0.9),
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
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_MIGs_theta_juncs)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 12*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) - sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_blood

ggsave(filename=filepath_blood_RNU6ATAC, plot=function_plot_blood,  limitsize = FALSE, units = "in")

################-----Number Minor Introns Retained-----#################
title_lab <- paste("Number of Minor Intron Retention Events", "\n", "Per Sample in the UDN and GREGoR Stanford Sites", "\n", sep=" ")
joined_MS <- joined %>% filter(Sample_ID %in% c("RD268", "GSS225379"))
joined_MS$MS <- "RNU4ATAC-opathies"

joined_RN6ATAC <- joined %>% filter(Sample_ID == "RD380")
joined_RN6ATAC$MS <- "RNU6ATAC-opathy"

joined_non_MS <- joined %>% filter(! Sample_ID %in% c("RD268", "GSS225379", "RD380"))
joined_non_MS$MS <- "No Minor Spliceosome-opathies"

joined <- bind_rows(joined_MS, joined_non_MS, joined_RN6ATAC)

function_counts_RNU4ATAC <- ggplot(joined, aes(x=MS, y=no_MIGs_theta_juncs)) + 
  geom_boxplot() +
  xlab("Mutation Type") +
  ylab("Number of Minor Intron Retention Events per Sample") +
  labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_grey(base_size = 17)+
  theme(plot.title = element_text(hjust = 0.5))

function_counts_RNU4ATAC
ggsave(filename=filepath_minor_introns_box_RNU6ATAC, plot=function_counts_RNU4ATAC,  limitsize = FALSE, units = "in")

joined_MS <- bind_rows(joined_MS, joined_RN6ATAC)
mean(joined_non_MS$no_MIGs_theta_juncs)
sd(joined_non_MS$no_MIGs_theta_juncs)
mean(joined_MS$no_MIGs_theta_juncs)
sd(joined_MS$no_MIGs_theta_juncs)
299/5.158416


library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions in",
                     "Minor Intron containing Genes (MIGs) per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs, numbers) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs, numbers) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
non_RNU4ATAC$type <- "NA"
RNU6ATAC <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs, numbers) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "RNU6ATAC-opathy"

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC, RNU6ATAC)

top_10 <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(sampleID, no_MIGs_theta_juncs, numbers, type) %>% head(5) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

top_5 <- joined %>% filter(sampleID == "RD380") %>% pull(type)
g2 <- subset(joined, type %in% top_5)


#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)
eight <- mean(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)
nine <- mean(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)
ten <- mean(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)
eleven <- mean(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)+sd(joined$no_MIGs_theta_juncs)


MIGs <- ggplot(joined, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, fill= type)) +
  geom_point(shape=21, alpha=.6, size=5) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend=TRUE) +
  #geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=8, aes(x = reorder(type, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = type), position=position_dodge(width=0.9),
                    #hjust=1.2, vjust= 0.5, colour="black", max.overlaps = Inf) +
  #geom_point(data=g2, aes(fill= type), shape=21, alpha=.6, size=5, fill="blue", show.legend=TRUE) +
  #geom_text_repel(data=(joined %>% filter(type %in% top_5)), size=8, aes(x = reorder(type, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = type), position=position_dodge(width=0.9),
  #                  hjust=1.2, vjust= 1.8, colour="black", max.overlaps = Inf) +
  labs(colour="type") +
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
  geom_hline(data=joined, aes(yintercept=seven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=eight, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=nine, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept= ten, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=eleven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_MIGs_theta_juncs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) - sd(joined$no_MIGs_theta_juncs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 22)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  theme(legend.position = "bottom")+
  scale_fill_manual(breaks = c("RNU4ATAC-opathy", "RNU6ATAC-opathy"),
    values = c("RNU6ATAC-opathy" = "red", "RNU4ATAC-opathy" = "blue", "NA" = "gray"))+
  guides(fill=guide_legend(title="Mutation Type:"))+
   coord_cartesian(clip = "off")

MIGs

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Poster/MIGs.pdf", plot=MIGs, limitsize = FALSE, units = "in")









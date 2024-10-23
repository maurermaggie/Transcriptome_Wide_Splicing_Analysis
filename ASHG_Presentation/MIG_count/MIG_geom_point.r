#library(EnsDb.Hsapiens.v79)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

################-----Number Minor Intron Retention Events UnLabeled-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)


MIGs_unlabeled <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=8, aes(x = reorder(type, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = type), position=position_dodge(width=0.9),
  #                 hjust=1.2,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
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
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

MIGs_unlabeled

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/MIG_count/MIGs_unlabeled.pdf", plot=MIGs_unlabeled,  limitsize = FALSE, units = "in")

################-----Number Minor Intron Retention Events Labeled-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) 
joined$numbers <- rownames(joined)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% head(5) %>% pull(numbers)
g1 <- subset(joined, numbers %in% top_10)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)

MIGs_labeled <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs)) +
  geom_point(shape=21, alpha=.6, size=5, fill=MyColour) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(fill= numbers), shape=21, alpha=.6, size=5, show.legend = FALSE, fill = "blue") +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=10, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, label = numbers), position=position_dodge(width=0.9),
                   hjust=3,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
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
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

MIGs_labeled

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/MIG_count/MIGs_labeled.pdf", plot=MIGs_labeled,  limitsize = FALSE, units = "in")

################-----Number Minor Intron Retention Events Labeled RNU4ATAC-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
not_RNU4ATAC <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
RD380 <- joined %>% filter(sampleID == "RD380")
RNU4ATAC$type <- "RNU4ATAC-opathy"
not_RNU4ATAC$type <- "None"
RD380$type <- "RNU6ATAC-opathy"
joined <- bind_rows(RNU4ATAC, not_RNU4ATAC, RD380)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, type) %>% head(4) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)

MyColour <- c("#EC936A",  "#FF0000", "blue")
names(MyColour) <- c("RNU4ATAC-opathy", "None", "RNU6ATAC-opathy")

MIGs_labeled_RNU4ATAC <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=5, show.legend=FALSE) +
  #scale_color_manual(values = MyColour)+
  scale_fill_manual(values = MyColour)+
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=6, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, label = type), position=position_dodge(),
  #                 hjust=1.25,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
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
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

MIGs_labeled_RNU4ATAC

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/MIG_count/MIGs_labeled_RNU4ATAC.pdf", plot=MIGs_labeled_RNU4ATAC,  limitsize = FALSE, units = "in")

################-----Number Minor Intron Retention Events Labeled RNU4ATAC-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

RNU6ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "Potential RNU6ATAC-opathy"

non_RNU6ATAC <- joined %>% filter(! sampleID %in% c("RD380"))  %>% select(sampleID, no_theta_juncs_in_MIGs, numbers)
non_RNU6ATAC$type <- NA

joined <- bind_rows(non_RNU6ATAC, RNU6ATAC)
top_10 <- RNU6ATAC %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)


MIGs_labeled_RNU6ATAC <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs)) +
  geom_point(shape=21, alpha=.6, size=5, fill="red") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(fill= numbers), shape=21, alpha=.6, size=5, show.legend = FALSE, fill = "blue") +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=6, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, label = type), position=position_dodge(),
                   hjust=1.25,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
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
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

MIGs_labeled_RNU6ATAC

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/MIG_count/MIGs_labeled_RNU6ATAC.pdf", plot=MIGs_labeled_RNU6ATAC,  limitsize = FALSE, units = "in")

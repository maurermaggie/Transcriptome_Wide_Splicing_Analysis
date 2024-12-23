library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
joined_fp <- args[1]
joined <- read_csv(joined_fp)

geom_point_RNU6ATAC <- args[2]

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
RNU6ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "RNU6ATAC-opathy"
non_RNU6ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD380", "RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
non_RNU6ATAC$type <- NA

joined <- bind_rows(RNU6ATAC, non_RNU6ATAC) %>% bind_rows(RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% filter(sampleID == "RD380") %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_MIGs_no_theta_juncs_in_MIGstheta_juncs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)

colours <- c("#D7DCD8", "#D7DCD8", "#D3B054")
names(colours) <- c(NA, "RNU4ATAC-opathy", "RNU6ATAC-opathy")

MIGs_RNU6ATAC <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU6ATAC-opathy")), size=8, aes(x = reorder(type, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, label = type), position=position_dodge(width=0.9),
  #                  hjust=1.2,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=seven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=eight, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=nine, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept= ten, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)

MIGs_RNU6ATAC

ggsave(filename=geom_point_RNU6ATAC, plot=MIGs_RNU6ATAC,  limitsize = FALSE, units = "in", height=12, width=10)

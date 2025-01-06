library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

args <- commandArgs(TRUE)
joined_fp <- args[1]
joined <- read_csv(joined_fp)

geom_point_RNU4ATAC <- args[2]
number_theta_in_MIGs <- args[3]
means <- args[4]

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
title_lab <- expression(atop("Number of" ~ theta ~ "Outlier Junctions",
                     "in MIGs per Sample"))
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("A1", "B1", "C1", "C2"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
RNU6ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("D1"))
RNU6ATAC$type <- "RNU6ATAC-opathy"
non_MS <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("D1", "A1", "B1", "C1", "C2"))
non_MS$type <- NA

joined <- bind_rows(RNU4ATAC, RNU6ATAC, non_MS)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(4) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eight <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
ten <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)

colours <- c("#D7DCD8", "#DA786C", "#D3B054")
names(colours) <- c(NA, "RNU4ATAC-opathy", "RNU6ATAC-opathy")

MIGs_RNU4ATAC <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(position = position_jitterdodge(0.1, dodge.width = .1), shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
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
  geom_hline(data=joined, aes(yintercept=nine, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept= ten, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_theta_juncs_in_MIGs)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 7*SD", x = 25, y=seven +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 8*SD", x = 25, y=eight +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 25, y=nine +10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=ten +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC

ggsave(filename=geom_point_RNU4ATAC, plot=MIGs_RNU4ATAC,  limitsize = FALSE, units = "in", height=12, width=10)

################-----Get Stats-----#################
RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("A1", "B1", "C1", "C2"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("D1", "B1", "C1", "C2"))
non_RNU4ATAC$type <- "non-RNU4ATAC"

A1 <- joined %>% filter(sampleID == "A1") %>% pull(no_theta_juncs_in_MIGs) - mean(joined$no_theta_juncs_in_MIGs)
A1_zscore <- A1/ sd(joined$no_theta_juncs_in_MIGs)
A1_no <- joined %>% filter(sampleID == "A1") %>% pull(no_theta_juncs_in_MIGs)

B1 <- joined %>% filter(sampleID == "B1") %>% pull(no_theta_juncs_in_MIGs) - mean(joined$no_theta_juncs_in_MIGs)
B1_zscore <- B1/ sd(joined$no_theta_juncs_in_MIGs)
B1_no <- joined %>% filter(sampleID == "B1") %>% pull(no_theta_juncs_in_MIGs)

C1 <- joined %>% filter(sampleID == "C1") %>% pull(no_theta_juncs_in_MIGs) - mean(joined$no_theta_juncs_in_MIGs)
C1_zscore <- C1/ sd(joined$no_theta_juncs_in_MIGs)
C1_no <- joined %>% filter(sampleID == "C1") %>% pull(no_theta_juncs_in_MIGs)

C2 <- joined %>% filter(sampleID == "C2") %>% pull(no_theta_juncs_in_MIGs) - mean(joined$no_theta_juncs_in_MIGs)
C2_zscore <- C2/ sd(joined$no_theta_juncs_in_MIGs)
C2_no <- joined %>% filter(sampleID == "C2") %>% pull(no_theta_juncs_in_MIGs)

D1 <- joined %>% filter(sampleID == "D1") %>% pull(no_theta_juncs_in_MIGs) - mean(joined$no_theta_juncs_in_MIGs)
D1_zscore <- D1/ sd(joined$no_theta_juncs_in_MIGs)
D1_no <- joined %>% filter(sampleID == "D1") %>% pull(no_theta_juncs_in_MIGs)

zscores <- data.frame(values=c("A1", "B1", "C1", "C2", "D1"),
                      zscores_theta_in_MIGs=c(A1_zscore, B1_zscore, C1_zscore, C2_zscore, D1_zscore),
                      number_theta_juncs_in_MIGs = c(A1_no, B1_no, C1_no, C1_no, D1_no))

write_csv(zscores, number_theta_in_MIGs)

RNU4ATAC_mean <- mean(RNU4ATAC$no_theta_juncs_in_MIGs)
RNU4ATAC_median <- median(RNU4ATAC$no_theta_juncs_in_MIGs)
RNU4ATA_sd <- sd(RNU4ATAC$no_theta_juncs_in_MIGs)
RNU4ATA_first <- quantile(RNU4ATAC$no_theta_juncs_in_MIGs, 0.25)
RNU4ATA_last <- quantile(RNU4ATAC$no_theta_juncs_in_MIGs, 0.75)

non_RNU4ATAC_mean <- mean(non_RNU4ATAC$no_theta_juncs_in_MIGs)
non_RNU4ATAC_median <- median(non_RNU4ATAC$no_theta_juncs_in_MIGs)
non_RNU4ATAC_sd <- sd(non_RNU4ATAC$no_theta_juncs_in_MIGs)
non_RNU4ATAC_first <- quantile(non_RNU4ATAC$no_theta_juncs_in_MIGs, 0.25)
non_RNU4ATAC_last <- quantile(non_RNU4ATAC$no_theta_juncs_in_MIGs, 0.75)

overall_mean <- mean(joined$no_theta_juncs_in_MIGs)
overall_median <- median(joined$no_theta_juncs_in_MIGs)
overall_sd <- sd(joined$no_theta_juncs_in_MIGs)
overall_first <- quantile(joined$no_theta_juncs_in_MIGs, 0.25)
overall_last <- quantile(joined$no_theta_juncs_in_MIGs, 0.75)

means_df <- data.frame(values= c("RNU4ATAC", "non-RNU4ATAC", "overall"), mean=c(RNU4ATAC_mean, non_RNU4ATAC_mean, overall_mean), sd=c(RNU4ATA_sd, non_RNU4ATAC_sd, overall_sd),
            first_quartile =c(RNU4ATA_first, non_RNU4ATAC_first, overall_first),
            last_quartile=c(RNU4ATA_last, non_RNU4ATAC_last, overall_last),
            median=c(RNU4ATAC_median, non_RNU4ATAC_median, overall_median))

write_csv(means_df, means)





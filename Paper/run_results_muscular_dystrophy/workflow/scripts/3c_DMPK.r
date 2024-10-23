library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv")
number_psi3_psi5_dmpk <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/psi3_psi5_dmpk.csv"
number_psi3_psi5_dmpk_plot <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/psi3_psi5_dmpk.pdf"

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

joined <- joined %>% select(sampleID, Psi3_Psi5_Junctions)
colnames(joined) <- c("sampleID", "n")
################-----Number Psi3 Psi5 Events DMPL Labeled-----#################
joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3 ~ "and" ~ psi ~ 5)
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant" ~ psi ~ 3 ~ "and" ~ psi ~ "5 Outliers"))
y_lab <- expression("Number of Junctions Significant" ~ psi ~ 3 ~ "and" ~ psi ~ "5 Outliers")
joined$numbers <- rownames(joined)

DMPK <- joined %>% filter(sampleID == "BEG_1230-1_T1227")
DMPK$type <- "DMPK Pathogenic Tandem Repeat"
non_DMPK <- joined %>% filter(sampleID != "BEG_1230-1_T1227")
non_DMPK$type <- NA

joined <- bind_rows(DMPK, non_DMPK)

top_10 <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score

DMPK_psi3_psi5 <- ggplot(joined, aes(x = reorder(numbers, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=8, aes(x = reorder(type, n), y = n, label = type), position=position_dodge(width=0.9),
                    hjust=1.2,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+10, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+10),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

DMPK_psi3_psi5

ggsave(filename=number_psi3_psi5_dmpk_plot, plot=DMPK_psi3_psi5,  limitsize = FALSE, units = "in")

################-----Get Stats-----#################
DMPK_pro <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n) - mean(joined$n)
DMPK_pro_zscore <- DMPK_pro/ sd(joined$n)
DMPK_pro_n <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n)

zscores <- data.frame(values=c("BEG_1230-1_T1227"),
                      zscores_theta_in_MIGs=c(DMPK_pro_zscore),
                      number_theta_juncs_in_MIGs = c(DMPK_pro_n))

write_csv(zscores, number_psi3_psi5_dmpk)





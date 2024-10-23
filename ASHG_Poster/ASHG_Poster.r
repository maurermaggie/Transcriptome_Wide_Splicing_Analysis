library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
RNU6ATAC <- joined %>% filter(sampleID == "RD380") %>% pull(no_MIGs_with_theta_juncs)
RD268 <- joined %>% filter(sampleID == "RD268") %>% pull(no_MIGs_with_theta_juncs)
GSS225379 <- joined %>% filter(sampleID == "GSS225379") %>% pull(no_MIGs_with_theta_juncs)
UDN550488 <- joined %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)
UDN238929 <- joined %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(no_MIGs_with_theta_juncs)

non_RNU4ATAC <- joined %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")) %>% pull(no_MIGs_with_theta_juncs)

non_RNU4ATAC_6ATAC <- joined %>% filter(! sampleID %in% c("RD380", "RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")) %>% pull(no_MIGs_with_theta_juncs)

################-----Number Minor Intron Retention Events RNU4ATAC Labeled Boxplot-----#################
RNU4ATAC <- joined %>% filter(no_MIGs_with_theta_juncs %in% c(153, 181, 174, 209))
RNU4ATAC$type <- "RNU4ATAC"
non_RNU4ATAC <- joined %>% filter(! no_MIGs_with_theta_juncs %in% c(153, 181, 174, 209, 160))
non_RNU4ATAC$type <- "None"
RNU6ATAC <- joined %>% filter(no_MIGs_with_theta_juncs == 160)
RNU6ATAC$type <- "RNU6ATAC"

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC, RNU6ATAC)
y_lab <- paste("Number of MIGs", "\n", "with Intron Retention")
title <- paste("Number of MIGs", "\n", "with Intron Retention per Sample")

dark_colors = c("gray", "blue", "red") %>% 
  col2rgb %>% #convert HEX colors to RGB Matrix
  "*"(0.7) %>% # make each component "darker"
  apply(2, lift_dv(rgb, maxColorValue = 255))
color <- c("gray", "blue", "red")
alpha_color <- 10
alpha_fill <- 0.5

MIG_boxplot <- ggplot(joined, aes(x=type, y=no_MIGs_with_theta_juncs, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  ylab(y_lab) +
  xlab("Minor Spliceosome Variants")+
  labs(title_lab=title)+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 75)+
  scale_fill_manual(values = alpha(color, alpha_fill)) + 
  scale_color_manual(values = alpha(color, alpha_color))+
  theme(panel.grid.minor = element_line(colour = "black", size = 0.15), 
        panel.grid.major = element_line(colour = "black", size = 0.15))+
   geom_boxplot(color = "black", fill = NA, fatten = NULL)+
  theme(plot.title = element_text(hjust = 0.5))



MIG_boxplot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Poster/MIGs_boxplot.pdf", plot=MIG_boxplot, limitsize = FALSE, units = "in", height=20, width=35)

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
title_lab <- paste("Number of Intron Retention Outliers in", "\n",
                     "Minor Intron containing Genes (MIGs) per Sample")
x_lab <- paste("Samples Ordered by Number of", "\n",
                       "Intron Retention Outliers in MIGs")
y_lab <- "Number of Intron Retention Outliers in MIGs"
joined$numbers <- rownames(joined)

RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam"))
RNU4ATAC$type <- "RNU4ATAC-opathy"
non_RNU4ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam", "RD380"))
non_RNU4ATAC$type <- "NA"
RNU6ATAC <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD380"))
RNU6ATAC$type <- "RNU6ATAC-opathy"

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC, RNU6ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(5) %>% pull(type)
g1 <- subset(joined, type %in% top_10)



#check which distribution this is considering
#give a z-score
seven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
nine <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)
eleven <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)


MIGs <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill= type)) +
  geom_point(shape=21, alpha=.6, size=8) +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  labs(colour="type") +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=seven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=nine, color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=eleven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  theme_classic(base_size = 35)+
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

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Poster/MIGs.pdf", plot=MIGs, limitsize = FALSE, units = "in",  width = 15, height = 15)





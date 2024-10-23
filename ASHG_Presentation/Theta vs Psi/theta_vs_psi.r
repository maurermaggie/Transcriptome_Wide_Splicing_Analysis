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

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

################-----Get Set-----#################
joined$percent <- joined$Psi3_Junctions/joined$Theta_Junctions

mean <- mean(joined$Theta_Junctions)
sd <- sd(joined$Theta_Junctions)
two_sd <- mean + sd + sd

joined %>% filter(Theta_Junctions >= two_sd) %>% arrange(percent) %>% print(n=36)
examples_theta_high <- c("RD142", "UDN281599", "GSS178184")
examples_same <- c("RD416", "RD384")

mean <- mean(joined$Psi3_Junctions)
sd <- sd(joined$Psi3_Junctions)
two_sd <- mean + sd + sd

joined %>% filter(Psi3_Junctions >= two_sd) %>% arrange(percent) %>% print(n=36)

examples_psi3_high <- c("RD235")

examples <- c(examples_psi3_high, examples_same, examples_theta_high)

joined_filter <- joined %>% filter(sampleID %in% examples) %>% select(sampleID, Psi3_Junctions, Theta_Junctions)
joined_filter$numbers <- c("E","F","D","A","B","C")

joined_filter_long <- melt(setDT(joined_filter), id.vars = c("sampleID", "numbers"), variable.name = "count")

joined_filter_long$Type <- ifelse(joined_filter_long$count == "Theta_Junctions", "Theta", "Psi3")

################-----Plot-----#################
Theta_psi3_comparison <- ggplot(joined_filter_long, aes(x = numbers, y = value)) +
  geom_bar(
    aes(color = Type, fill = Type),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7
    ) +
  ylab("Number of Aberrant Junctions")+
  scale_color_manual(values = c("Darkgreen", "Blue"))+
  scale_fill_manual(values = c("Darkgreen", "Blue")) +
  theme_gray(base_size=30)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=30))+
  xlab("Samples")

Theta_psi3_comparison
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/Theta vs Psi/theta_psi_comp.pdf", plot=Theta_psi3_comparison,  limitsize = FALSE, units = "in")

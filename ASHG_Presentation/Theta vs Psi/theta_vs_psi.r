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
colnames(all) <- c("sampleID", "all_count")

psi5_filtered <- all_filtered %>% filter(type == "psi5") %>% filter(padjust <=0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
psi5_count <- psi5_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi5 <- setdiff(files$sampleID, psi5_filtered$sampleID)
missing_psi5_df <- data.frame(sampleID = missing_psi5, n = 0)
psi5_count <- bind_rows(psi5_count, missing_psi5_df)
tally_psi5 <- psi5_count %>% filter(! sampleID %in% low_RIN$sampleID)
colnames(tally_psi5) <- c("sampleID", "psi5_count")

psi3_filtered <- all_filtered %>% filter(type == "psi3") %>% filter(padjust <=0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
psi3_count <- psi3_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_psi3 <- setdiff(files$sampleID, psi3_filtered$sampleID)
missing_psi3_df <- data.frame(sampleID = missing_psi3, n = 0)
psi3_count <- bind_rows(psi3_count, missing_psi3_df)
tally_psi3 <- psi3_count %>% filter(! sampleID %in% low_RIN$sampleID)
colnames(tally_psi3) <- c("sampleID", "psi3_count")

theta_filtered <- all_filtered %>% filter(type == "theta") %>% filter(padjust <=0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2)
theta_count <- theta_filtered %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
missing_theta <- setdiff(files$sampleID, theta_filtered$sampleID)
missing_theta_df <- data.frame(sampleID = missing_theta, n = 0)
theta_count <- bind_rows(theta_count, missing_theta_df)
tally_theta <- theta_count %>% filter(! sampleID %in% low_RIN$sampleID)
colnames(tally_theta) <- c("sampleID", "theta_count")

joined <- left_join(all, tally_psi3) %>% left_join(tally_psi5) %>% left_join(tally_theta)

################-----Get Set-----#################
joined$percent <- joined$psi3_count/joined$theta_count

mean <- mean(joined$theta_count)
sd <- sd(joined$theta_count)
two_sd <- mean + sd + sd

joined %>% filter(theta_count >= two_sd) %>% arrange(percent) %>% print(n=36)
examples_theta_high <- c("RD142", "UDN281599", "GSS178184")
examples_same <- c("RD416", "RD384")

mean <- mean(joined$psi3_count)
sd <- sd(joined$psi3_count)
two_sd <- mean + sd + sd

joined %>% filter(psi3_count >= two_sd) %>% arrange(percent) %>% print(n=36)

examples_psi3_high <- c("RD235")

examples <- c(examples_psi3_high, examples_same, examples_theta_high)

joined_filter <- joined %>% filter(sampleID %in% examples) %>% select(sampleID, theta_count, psi3_count)
joined_filter$numbers <- c("E","F","D","A","B","C")

joined_filter_long <- melt(setDT(joined_filter), id.vars = c("sampleID", "numbers"), variable.name = "count")

joined_filter_long$Type <- ifelse(joined_filter_long$count == "theta_count", "Theta", "Psi3")

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

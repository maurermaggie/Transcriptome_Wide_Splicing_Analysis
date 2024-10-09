require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)

display.brewer.all
pal <- brewer.pal(10, "Paired")

genes_junctions_joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_filtered_with_sibs.csv")
all_count_genes <- genes_junctions_joined %>% select(sampleID, contains("Genes"))
all_count_junctions <- genes_junctions_joined %>% select(sampleID, contains("Junctions"))

################-----Get Number by Gene per Person-----#################
mean(all_count_genes$All_Genes)
sd(all_count_genes$All_Genes)

mean(all_count_genes$Psi3_Genes)
sd(all_count_genes$Psi3_Genes)

mean(all_count_genes$Psi5_Genes)
sd(all_count_genes$Psi5_Genes)

mean(all_count_genes$Theta_Genes)
sd(all_count_genes$Theta_Genes)


################-----Get Number by Junction per Person-----#################
mean(all_count_junctions$All_Junctions)
sd(all_count_junctions$All_Junctions)

mean(all_count_junctions$Psi3_Junctions)
sd(all_count_junctions$Psi3_Junctions)

mean(all_count_junctions$Psi5_Junctions)
sd(all_count_junctions$Psi5_Junctions)

mean(all_count_junctions$Theta_Junctions)
sd(all_count_junctions$Theta_Junctions)

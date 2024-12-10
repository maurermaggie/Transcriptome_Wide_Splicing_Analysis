library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

########################################################################
###################-----All RNU4ATAC v RNU6ATAC-----####################
########################################################################
RNU4ATAC <- c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")
RNU6ATAC <- c("RD380")

all_RNU4ATAC <- all_uncompiled %>% filter(sampleID %in% RNU4ATAC) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
all_RNU6ATAC <- all_uncompiled %>% filter(sampleID %in% RNU6ATAC) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique

gene_list <- list(`RNU4ATAC Probands (n=4)`= all_RNU4ATAC, `RNU6ATAC Proband`=all_RNU6ATAC)

all <- ggVennDiagram(gene_list, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none",
                sets.bar.color = c("grey", "red"),
                top.bar.color = c("grey", "red" , "red"),
                intersection.matrix.color = c("grey", "red", "grey", "red")) 
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/Upset_plot/RNU6_v_all.pdf", plot=all,  limitsize = FALSE, units = "in")

########################################################################
###################-----All RNU4ATAC v RNU6ATAC-----####################
########################################################################
RD268 <- all_uncompiled %>% filter(sampleID =="RD268") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
all_RNU6ATAC <- all_uncompiled %>% filter(sampleID %in% RNU6ATAC) %>% %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
GSS225379 <- all_uncompiled %>% filter(sampleID =="GSS225379") %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
UDN550488 <- all_uncompiled %>% filter(sampleID =="UDN550488.Aligned.sortedByCoord.out.bam") %>%  filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
UDN238929 <- all_uncompiled %>% filter(sampleID =="UDN238929.Aligned.sortedByCoord.out.bam") %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
RD380 <- all_uncompiled %>% filter(sampleID =="RD380") %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta")%>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique

gene_list_split <- list(`RNU4ATAC Proband 1`= RD268, `RNU4ATAC Proband 2`=GSS225379, `RNU4ATAC Proband 3`=UDN550488, `RNU4ATAC Proband 4`=UDN238929, `RNU6ATAC Proband`=RD380)

all_split <- ggVennDiagram(gene_list_split, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none",   top.bar.show.numbers = TRUE,  top.bar.numbers.size = 3, sets.bar.show.numbers = TRUE,
                sets.bar.color = c("grey", "grey", "grey", "grey", "red"),
                top.bar.color = c("grey", "grey", "grey", "grey", "red", "grey", "grey", "grey", "red", "grey", "grey", "red", "grey", 'red', "red", "grey", "grey", "red", "grey", "red"),
                intersection.matrix.color = c("grey", "grey", "grey", "grey", "red", 
                                                "grey", "grey", 
                                                "grey", "grey", 
                                                "grey", "grey", 
                                                "grey", "red", 
                                                "grey", "grey", 
                                                "grey", "grey", 
                                                "grey", "red",
                                                "grey", "grey",
                                                "grey", "red",
                                                "grey", "red",
                                                "grey", "grey", "grey",
                                                "grey", "grey", "grey",
                                                "grey", "grey", "red",
                                                "grey", "grey", "grey",
                                                "grey", "grey", "red"
                                                ))

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/Upset_plot/RNU6_v_all_split.pdf", plot=all_split,  limitsize = FALSE, units = "in")

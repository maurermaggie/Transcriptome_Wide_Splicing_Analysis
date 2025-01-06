library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(VennDiagram)
library(tidyverse)

args <- commandArgs(TRUE)

all_uncompiled_fp <- args[1]
all_uncompiled <- read_csv(all_uncompiled_fp)

MIGs_fp <- args[2]
MIGs <- read_csv(MIGs_fp)

Venn_fp <- args[3]

MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

########################################################################
###################-----All RNU4ATAC v RNU6ATAC-----####################
########################################################################
RNU4ATAC <- c("A1", "B1", "C1", "C2")
RNU6ATAC <- c("D1")

all_RNU4ATAC <- all_uncompiled %>% filter(sampleID %in% RNU4ATAC) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
all_RNU6ATAC <- all_uncompiled %>% filter(sampleID %in% RNU6ATAC) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique

inter<- intersect(all_RNU4ATAC, all_RNU6ATAC) %>% length
RNU6 <- all_RNU6ATAC %>% length
RNU4 <- all_RNU4ATAC %>% length

plt <-draw.pairwise.venn(area1=RNU6, area2=RNU4,cross.area=inter, 
                   category=c("RNU6ATAC","RNU4ATAC"),fill=c("#88ACAA","#F2D177")) 
ggsave(plt, file=Venn_fp, device = "pdf", limitsize = FALSE, units = "in", height=20, width=20)

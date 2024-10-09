library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")

########################################################################
###################-----All RNU4ATAC v RNU6ATAC-----####################
########################################################################
RNU4ATAC <- c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")
RNU6ATAC <- c("RD380")

all_RNU4ATAC <- all_uncompiled %>% filter(sampleID %in% RNU4ATAC) %>% pull(hgncSymbol) %>% unique
all_RNU6ATAC <- all_uncompiled %>% filter(sampleID %in% RNU6ATAC)  %>% pull(hgncSymbol) %>% unique

inter<- intersect(all_RNU4ATAC, all_RNU6ATAC) 
diff_RNU6 <- setdiff(all_RNU6ATAC, all_RNU4ATAC) 
diff_RNU4 <- setdiff(all_RNU4ATAC, all_RNU6ATAC)

gene_list <- list(RNU4ATAC= all_RNU4ATAC, RNU6ATAC=all_RNU6ATAC)
venn_genes_all <- ggVennDiagram(gene_list) + scale_fill_gradient(low="grey90",high = "red")+
                    scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/VennDiagrams/RNU6_v_all.pdf", plot=venn_genes_all,  limitsize = FALSE, units = "in")

########################################################################
###################-----All RNU4ATAC v RNU6ATAC-----####################
########################################################################
RNU4ATAC <- c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")
RNU6ATAC <- c("RD380")

RD268 <- all_uncompiled %>% filter(sampleID =="RD268") %>% pull(hgncSymbol) %>% unique
all_RNU6ATAC <- all_uncompiled %>% filter(sampleID %in% RNU6ATAC)  %>% pull(hgncSymbol) %>% unique
GSS225379 <- all_uncompiled %>% filter(sampleID =="GSS225379") %>% pull(hgncSymbol) %>% unique
UDN550488 <- all_uncompiled %>% filter(sampleID =="UDN550488.Aligned.sortedByCoord.out.bam") %>% pull(hgncSymbol) %>% unique
UDN238929 <- all_uncompiled %>% filter(sampleID =="UDN238929.Aligned.sortedByCoord.out.bam") %>% pull(hgncSymbol) %>% unique

gene_list <- list(`RNU4ATAC Proband 1`= RD268, `RNU4ATAC Proband 2`=GSS225379, `RNU4ATAC Proband 3`=UDN550488, `RNU4ATAC Proband 4`=UDN238929, `RNU6ATAC Proband`=all_RNU6ATAC)
venn_genes_all_split <- ggVennDiagram(gene_list) + scale_fill_gradient(low="grey90",high = "red") +
            scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/VennDiagrams/RNU6_v_all_split.pdf", plot=venn_genes_all_split,  limitsize = FALSE, units = "in")

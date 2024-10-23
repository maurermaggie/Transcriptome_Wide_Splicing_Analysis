library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_AND_outlier_status_metadata_filtered_with_sibs.csv")
missing_metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/missing_metadata.csv")

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

RNU6ATAC <- all_uncompiled %>% filter(sampleID == "RD380") %>% filter(type == "theta") %>% filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique
RNU4ATAC <- all_uncompiled %>% filter(sampleID %in% c("RD268", "GSS225379", "UDN550488.Aligned.sortedByCoord.out.bam", "UDN238929.Aligned.sortedByCoord.out.bam")) %>% filter(type == "theta") %>% filter(hgncSymbol %in% MIG_table_select) %>% pull(hgncSymbol) %>% unique()

not_RNU6ATAC <- setdiff(MIG_table_select, RNU6ATAC)
not_RNU4ATAC <- setdiff(MIG_table_select, RNU4ATAC)

RNU4ATAC_RNU6ATAC <- intersect(RNU4ATAC, RNU6ATAC) %>% length
RNU4ATAC_notRNU6ATAC <- intersect(RNU4ATAC, not_RNU6ATAC) %>% length
notRNU4ATAC_RNU6ATAC <- intersect(not_RNU4ATAC, RNU6ATAC) %>% length
notRNU4ATAC_notRNU6ATAC <- intersect(not_RNU4ATAC, not_RNU6ATAC) %>% length

MS_ct <- data.frame(
    "RNU6ATAC" = c(RNU4ATAC_RNU6ATAC, RNU4ATAC_notRNU6ATAC),
    "RNU4ATAC" = c(notRNU4ATAC_RNU6ATAC, notRNU4ATAC_notRNU6ATAC),
    row.names = c("RNU6ATAC", "notRNU6ATAC"),
    stringsAsFactors = FALSE
)
            
colnames(MS_ct) <- c("RNU4ATAC", "notRNU4ATAC")

MS_f <- fisher.test(MS_ct)
MS_tidy <- tidy(MS_f)

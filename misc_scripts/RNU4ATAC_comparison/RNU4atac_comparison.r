library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(broom)
library(tidyverse)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/output/FRASER_Stanford_Filters/FRASER_output.csv")
MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

########################################################################
##########-----Get Significant MIG Intron Retention Events-----#########
########################################################################
#get MIG list of genes
MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG") %>% pull(gene_symbol)

########################################################################
###################-----All RNU4ATAC Upset Plot-----####################
########################################################################
RD268 <- all_uncompiled %>% filter(sampleID == "RD268") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
GSS225379 <- all_uncompiled %>% filter(sampleID == "GSS225379") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
UDN550488 <- all_uncompiled %>% filter(sampleID == "UDN550488.Aligned.sortedByCoord.out.bam") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
UDN238929 <- all_uncompiled %>% filter(sampleID == "UDN238929.Aligned.sortedByCoord.out.bam") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique

gene_list <- c(RD268, GSS225379, UDN550488, UDN238929) %>% unique 
genes <- data.frame(gene_list)
genes$A1 <- ifelse(genes$gene_list %in% RD268, TRUE, FALSE)
genes$B1 <- ifelse(genes$gene_list %in% GSS225379, TRUE, FALSE)
genes$C1 <- ifelse(genes$gene_list %in% UDN550488, TRUE, FALSE)
genes$C2 <- ifelse(genes$gene_list %in% UDN238929, TRUE, FALSE)
rownames(genes) <- genes$gene_list
#genes <- genes %>% select(-gene_list)
y_lab <- expression("Number of Shared MIGs with" ~ theta ~ "outliers")

upset <- upset(genes, colnames(genes)[2:5], name='gene_list', width_ratio=0.1,
    matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=3))
        + scale_color_manual(
            values=c('A1'='#DA786C', 'B1'='#88ACAA', 'C1'='#F2D177', 'C2'= '#9F5E64'),
            guide=guide_legend(override.aes=list(shape='circle'))
        )
    ),
    queries=list(
        upset_query(set='A1', fill='#DA786C'),
        upset_query(set='B1', fill='#88ACAA'),
        upset_query(set='C1', fill='#F2D177'),
        upset_query(set='C2', fill='#9F5E64')
    )
) +
theme_classic(base_size=30)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  xlab("RNU4ATAC-opathy")+
  ylab(y_lab)
upset
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/RNU4ATAC_comparison/RNU4ATAC_upset.pdf", plot=upset, limitsize = FALSE, units = "in", height=15, width=15)

########################################################################
################-----RNU4ATAC Contingency Tables-----###################
########################################################################
make_table <- function(list1, list2){
    all <- all_uncompiled %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
    not_list1 <- setdiff(all, list1)
    not_list2 <- setdiff(all, list2)

    list1_list2 <- intersect(list1, list2) %>% length
    list1_notlist2 <- intersect(list1, not_list2) %>% length
    notlist1_list2 <- intersect(not_list1, list2) %>% length
    notlist1_notlist2 <- intersect(not_list1, not_list2) %>% length

    MS_ct <- data.frame(
    "list1" = c(list1_list2, list1_notlist2),
    "not_list1" = c(notlist1_list2, notlist1_notlist2),
    row.names = c("list2", "not_list2"),
    stringsAsFactors = FALSE)
            
colnames(MS_ct) <- c("list1", "list2")

MS_f <- fisher.test(MS_ct)
MS_tidy <- tidy(MS_f)
MS_tidy
}

RD_GSS <- make_table(RD268, GSS225379)
RD_sib1 <- make_table(RD268, UDN550488)
RD_sib2 <- make_table(RD268, UDN238929)

GSS_sib1 <- make_table(GSS225379, UDN550488)
GSS_sib2 <- make_table(GSS225379, UDN238929)

sib1_sib2 <- make_table(UDN550488, UDN238929)

df <- data.frame(value=c("RD268_GSS225379", "RD268_UDN550488", "RD268_UDN238929", "GSS225379_UDN550488", "GSS225379_UDN238929", "UDN550488_UDN238929"),
                    pvalue=c(RD_GSS$p.value, RD_sib1$p.value, RD_sib2$p.value, GSS_sib1$p.value, GSS_sib2$p.value, sib1_sib2$p.value),
                    estimate=c(RD_GSS$estimate, RD_sib1$estimate, RD_sib2$estimate, GSS_sib1$estimate, GSS_sib2$estimate, sib1_sib2$estimate))
df$q <- p.adjust(df$pvalue, method="fdr", n=length(df$pvalue))

########################################################################
##########-----RNU4ATAC vs RNU6ATAC Contingency Tables-----#############
########################################################################
RNU4ATAC <- c(RD268, GSS225379, UDN550488, UDN238929) %>% unique
RNU6ATAC <- all_uncompiled %>% filter(sampleID == "RD380") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
RNU4ATAC_RNU6ATAC <- intersect(RNU4ATAC, RNU6ATAC) %>% length

all <- all_uncompiled %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
not_RNU4ATAC <- setdiff(all, RNU4ATAC)
not_RNU6ATAC <- setdiff(all, RNU6ATAC)

RNU4ATAC_notRNU6ATAC <- intersect(RNU4ATAC, not_RNU6ATAC) %>% length
notRNU4ATAC_RNU6ATAC <- intersect(not_RNU4ATAC, RNU6ATAC) %>% length
notRNU4ATAC_notRNU6ATAC <- intersect(not_RNU4ATAC, not_RNU6ATAC) %>% length


MS_ct <- data.frame(
    "RNU4ATAC" = c(RNU4ATAC_RNU6ATAC, RNU4ATAC_notRNU6ATAC),
    "notRNU4ATAC" = c(notRNU4ATAC_RNU6ATAC, notRNU4ATAC_notRNU6ATAC),
    row.names = c("RNU6ATAC", "notRNU6ATAC"),
    stringsAsFactors = FALSE
)
            
colnames(MS_ct) <- c("RNU4ATAC", "notRNU4ATAC")

MS_f <- fisher.test(MS_ct)
MS_tidy <- tidy(MS_f)

########################################################################
##########-----RNU4ATAC vs RNU6ATAC Contingency Tables-----#############
########################################################################
RNU4ATAC <- c(RD268, GSS225379, UDN550488, UDN238929) %>% unique
RNU6ATAC <- all_uncompiled %>% filter(sampleID == "RD380") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
RNU4ATAC_RNU6ATAC <- intersect(RNU4ATAC, RNU6ATAC) %>% length

all <- all_uncompiled %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi)>=0.3) %>% filter(type == "theta") %>% 
                    filter(hgncSymbol %in% MIG_table_select) %>%
                     pull(hgncSymbol) %>% unique
not_RNU4ATAC <- setdiff(all, RNU4ATAC)
not_RNU6ATAC <- setdiff(all, RNU6ATAC)

RNU4ATAC_notRNU6ATAC <- intersect(RNU4ATAC, not_RNU6ATAC) %>% length
notRNU4ATAC_RNU6ATAC <- intersect(not_RNU4ATAC, RNU6ATAC) %>% length
notRNU4ATAC_notRNU6ATAC <- intersect(not_RNU4ATAC, not_RNU6ATAC) %>% length


MS_ct <- data.frame(
    "RNU4ATAC" = c(RNU4ATAC_RNU6ATAC, RNU4ATAC_notRNU6ATAC),
    "notRNU4ATAC" = c(notRNU4ATAC_RNU6ATAC, notRNU4ATAC_notRNU6ATAC),
    row.names = c("RNU6ATAC", "notRNU6ATAC"),
    stringsAsFactors = FALSE
)
            
colnames(MS_ct) <- c("RNU4ATAC", "notRNU4ATAC")

MS_f <- fisher.test(MS_ct)
MS_tidy <- tidy(MS_f)

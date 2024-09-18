#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(tidyverse)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(GO.db)

all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
joined <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/compiled_counts/FRASER_results_raw.csv")
filepath_blood <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/minor_intron_events.pdf"
filepath_minor_introns_box <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/minor_introns_box.pdf"

subset_uncompiled <- read_csv("/home/maurertm/smontgom/maurertm/FRASER_analysis/outputs_old/FRASER_Output_Apr_1/FRASER_output_2024-04-01_14:51.csv")
subset <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/NIH_talk/FRASER_results_raw.csv")
filepath_juncs_no_comparison <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/juncs_no_comparison.pdf"
filepath_subset_MIGs <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/subset_minor_introns.pdf"
filepath_compare_RD268 <- "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/filepath_compare_RD268.pdf"

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
################-----Number Minor Intron Retention Events-----#################
title_lab <- paste("Number of Minor Intron Retention Events", "\n", "per Sample in the UDN and GREGoR Stanford Sites", sep=" ")

joined$numbers <- rownames(joined)

top_10 <- joined %>% arrange(desc(no_MIGs_theta_juncs)) %>% select(Sample_ID, no_MIGs_theta_juncs, numbers) %>% head(3) %>% pull(numbers)
g1 <- subset(joined, numbers %in% top_10)


#check which distribution this is considering
#give a z-score
function_plot_blood <- ggplot(joined, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Samples Ordered by Number of Minor Intron Retention Events") +
  ylab("Number of Minor Intron Retention Events") +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)-sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs)+sd(no_MIGs_theta_juncs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$no_MIGs_theta_juncs)+8),  size=2) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 10*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean + 12*SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs) + sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_MIGs_theta_juncs) - sd(joined$no_MIGs_theta_juncs)+8),  size=2) +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

function_plot_blood

ggsave(filename=filepath_blood, plot=function_plot_blood,  limitsize = FALSE, units = "in")

################-----Compare groups-----#################
title_lab <- paste("Comparison of the Number of Aberrant Splicing Events", "\n", "Per Sample in Full (n=407) vs Subsetted (n=32) Datasets", sep=" ")

total_subset <- subset_uncompiled %>% group_by(sampleID) %>% tally()
total_subset$Type <- "32"
total_full <- all_uncompiled %>% group_by(sampleID) %>% tally()
total_full$Type <- "407"

combined <- bind_rows(total_subset, total_full)
colnames(combined) <- c("sampleID", "Number_of_Aberrant_Junctions", "Type")

function_counts_per <- ggplot(combined, aes(x=Type, y=Number_of_Aberrant_Junctions)) + 
  geom_boxplot() +
  xlab("Populations Sizes") +
  ylab("Number of Total Aberrant Junctions Per Sample") +
  labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_grey(base_size = 20)

function_counts_per
ggsave(filename=filepath_juncs_no_comparison, plot=function_counts_per,  limitsize = FALSE, units = "in")

################-----Number Subset Blood-----#################
y_lab <- "Number of Minor Intron Retention Rvents"
title_lab <- paste("Number of Minor Intron Retention Rvents in a Subset (n=32)", "\n", " of the UDN and GREGoR Stanford Sites", "\n", "with a RNU4ATAC-opathy proband (29)", sep=" ")

subset$numbers <- rownames(subset)

g1 <- subset(subset, Sample_ID %in% c("RD268"))

#check which distribution this is considering
#give a z-score

function_plot_subset <- ggplot(subset, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs)) +
  geom_point() +
  xlab("Samples Ordered by Number of Minor Intron Retention Events") +
  ylab("Number of Minor Intron Retention Events") +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(colour= Sample_ID)) +
  labs(colour="Sample ID") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5))

function_plot_subset

ggsave(filename=filepath_subset_MIGs, plot=function_plot_subset,  limitsize = FALSE, units = "in")

################-----Number Minor Introns Retained-----#################
title_lab <- paste("Comparison of the Number of Minor Intron Retention Events", "\n", "Per Sample between Probands with and without RNU4ATAC-opathies", sep=" ")
joined_MS <- joined %>% filter(Sample_ID %in% c("RD268", "GSS225379"))
joined_MS$MS <- "RNU4ATAC-opathy Probands"

joined_non_MS <- joined %>% filter(! Sample_ID %in% c("RD268", "GSS225379"))
joined_non_MS$MS <- "Non RNU4ATAC-opathy Probands"

joined <- bind_rows(joined_MS, joined_non_MS)

function_counts_RNU4ATAC <- ggplot(joined, aes(x=MS, y=no_MIGs_theta_juncs)) + 
  geom_boxplot() +
  xlab("Mutation Type") +
  ylab("Number of Minor Intron Retention Events per Sample") +
  labs(title= title_lab) +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_grey(base_size = 17)

function_counts_RNU4ATAC
ggsave(filename=filepath_minor_introns_box, plot=function_counts_RNU4ATAC,  limitsize = FALSE, units = "in")

mean(joined_non_MS$no_MIGs_theta_juncs)
sd(joined_non_MS$no_MIGs_theta_juncs)
mean(joined_MS$no_MIGs_theta_juncs)
sd(joined_MS$no_MIGs_theta_juncs)
309.5/5.832099

################-----Compare Subset vs Full RD268-----#################
total_subset <- subset %>% filter(Sample_ID == "RD268") %>% pull(no_MIGs_theta_juncs)
total_full <- joined %>% filter(Sample_ID == "RD268") %>% pull(no_MIGs_theta_juncs)

data_RD268 <- data.frame(Number_of_Minor_Intron_Retention_Events = c(as.numeric(total_subset), as.numeric(total_full)), Sample_Size = c("32", "407"))

title_lab <- paste("Comparing the Number of Minor Intron Retention Events", "\n", "for Proband 44 when Analysed in Full (n=407)", "\n", "and Subsetted (n=32) Datasets", sep= " ")
function_RD268 <- ggplot(data_RD268, aes(Sample_Size, Number_of_Minor_Intron_Retention_Events, fill = Sample_Size)) +     
                    geom_col(position = 'dodge') +
                    geom_line()+
                    geom_point()+
                    xlab("Sample Size") +
                    ylab("Number of Minor Intron Retention Events Detected in Proband 44") +
                    labs(title= title_lab) +
                    theme_update(plot.title = element_text(hjust = 0.5))+
                    theme_grey(base_size = 20)+
                    theme(plot.title = element_text(hjust = 0.5))
                    
function_RD268
ggsave(filename=filepath_compare_RD268, plot=function_RD268,  limitsize = FALSE, units = "in")

################-----MIGs GO-----#################
MIGs <- MIGs %>% filter(gene_class == "MIG") %>% dplyr::select("gene_symbol","ensembl_gene_id")

ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

go_terms <- getBM(
  mart=ensembl, attributes=c("hgnc_symbol", "uniprot_gn_id", "uniprot_gn_symbol", "go_id", "namespace_1003", "name_1006"),
  filters="hgnc_symbol", values=MIGs$gene_symbol)

write_csv(go_terms, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/go_terms.csv")
go_terms <- go_terms[!(!is.na(go_terms$go_id) & go_terms$go_id==""), ]

go_meaning <- select(GO.db, columns=c("GOID","TERM"), key = go, keytype="GOID")
colnames(go_meaning) <- c("go_id", "go_meaning")
go_terms <- left_join(go_terms, go_meaning)

write_csv(go_terms, "/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/go_terms_with_meaning.csv")
go_terms_joined <- go_terms %>% group_by(go_meaning) %>% tally() %>% arrange(desc(n))

go_terms_top_10 <- go_terms_joined[1:10,]

title_lab <- paste("Examination of the Distribution of", "\n", "Gene Ontology IDs assocaited with MIGs", sep= " ")

p<-ggplot(data=go_terms_joined, aes(x=go_meaning, y=n)) +
  geom_bar(stat="identity") +
  geom_point(data=go_terms_top_10, aes(colour= go_meaning))+
  xlab("Gene Ontology ID") +
  ylab("Number of Genes with Gene Ontology ID") +
  labs(title= title_lab) +
  guides(color = guide_legend(title = "Gene Ontology ID Meaning"))+
  theme_update(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/projects/GREGoR_2024_Boston_Conference/go_terms.pdf", plot=p,  limitsize = FALSE, units = "in")

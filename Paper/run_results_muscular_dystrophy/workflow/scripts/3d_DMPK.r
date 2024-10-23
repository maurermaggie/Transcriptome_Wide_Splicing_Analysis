library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(broom)
library(FRASER)
library("org.Hs.eg.db")
library(tidyverse)

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

psi3 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi3_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
psi5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi5_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
theta <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_theta_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_muscle_jaccard_results_jaccard.csv")
output_file <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/GSEA_analysis.pdf"
all_count <- bind_rows(psi3, psi5, theta, jaccard)

count <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/dataframes/counts.csv")
dmpk_dmpk_genes <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/3d_dmpk_misspliced.pdf" 
dmpk_dmpk_genes_stats <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/3d_dmpk_misspliced_zscores.csv" 
dmpk_dmpk_jaccard <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/3d_dmpk_misspliced_jaccard.pdf" 
dmpk_dmpk_genes_stats_jaccard <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/3d_dmpk_misspliced_zscores_jaccard.csv" 
dmpk_dmpk_genes_means <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/stats/3d_dmpk_misspliced_means.csv" 

################-----Get ENSG IDs-----#################
Genes <- c("CLCN1",	"TNNT2"	, "DMD", "RYR1", "DTNA", "ATP2A1", "TTN", "ATP2A2", "LDB3", "BIN1", "PDLIM3",
"CACNA1S", "CAPZB", "CACNA2D1", "KCNAB1", "DMD", "SCN5A", "DTNA", "TTN", "MTMR1", "LDB3",	"MAPT", "CAPZB", "MBNL1", "TNNT3", "MBNL2",	
"PDLIM3", "DGKH", "PDLIM7", "FHOD1", "APP", "NRAP", "GRIN1", "ABLIM2", "MAPT", "MYBPC1", "CAMK2D", "MYOM1", "SORBS1",
"ANK2", "DCLK1", "MPRIP", "INSR", "TANC2", "SOS1", "KCNMA1", "ALPK3", "CSNK1D", "CAMK2B", "CACNA1D", "MAP4K4", "LIMCH1",
"ADD1", "MBNL1", "ADD3", "MBNL2", "PPP1R12A", "NFIX", "CLASP2", "MEF2C", "RYR2", "NCOR2", "DLGAP1", "FXR1", "NTM", "CAMKK2",	
"GFPT1", "WNK1", "IMPDH2", "DGKH", "MTMR1", "NDUFV3", "CAPN3", "GRIP1", "PHKA1", "PHKA2", "KIF13A","VPS39",	"ARFGAP2",	
"COPZ2",	"TBC1D15","ANK2","OPA1",	"UBE2D3",	"USP25",	"TXNL4A",	"VEGFA",	"MAPT", "MLF1",	"CLASP1",	"NDUFV3",	"CLTC")	

ensg = mapIds(org.Hs.eg.db,
                    keys=Genes, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")

number_DMPK_junctions_psi3_psi5 <- all_count %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
joined <- number_DMPK_junctions_psi3_psi5
number_DMPK_genes_psi3_psi5 <- all_count %>% filter(type %in% c("psi3", "psi5")) %>% filter(hgncSymbol %in% ensg) %>% select(sampleID, hgncSymbol) %>% unique() %>% group_by(sampleID)  %>% tally() %>% arrange(desc(n))

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3 ~ "and" ~ psi ~ "5 in DMPK Mis-Spliced Genes")
x_lab <- expression(atop("Samples ordered by Number of Junctions with Significant",
                     psi ~ 3 ~ "and" ~ psi ~ "5 Outliers in DMPK Mis-Spliced Genes"))
y_lab <- expression(atop("Number of Junctions with Significant" ~ psi ~ 3 ~ "and" ~ psi ~ "5 Outliers",
                        "in DMPK Mis-Spliced Genes"))
joined$numbers <- rownames(joined)

DMPK <- joined %>% filter(sampleID == "BEG_1230-1_T1227")
DMPK$type <- "DMPK Pathogenic Tandem Repeat"
non_DMPK <- joined %>% filter(sampleID != "BEG_1230-1_T1227")
non_DMPK$type <- NA

joined <- bind_rows(DMPK, non_DMPK)

top_10 <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score

DMPK_psi3_psi5 <- ggplot(joined, aes(x = reorder(numbers, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=8, aes(x = reorder(type, n), y = n, label = type), position=position_dodge(width=0.9),
                    hjust=1.2,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+2, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+2),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

DMPK_psi3_psi5

ggsave(filename=dmpk_dmpk_genes, plot=DMPK_psi3_psi5,  limitsize = FALSE, units = "in")

################-----Get Stats-----#################
DMPK_pro <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n) - mean(joined$n)
DMPK_pro_zscore <- DMPK_pro/ sd(joined$n)
DMPK_pro_n <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n)

not_DMPK <- joined %>% filter(sampleID != "BEG_1230-1_T1227") 
mean_not_DMPK <- mean(not_DMPK$n)
sd_not_DMPK <- sd(not_DMPK$n)

DMPK_pro_zscore <- DMPK_pro/ sd(joined$n)
DMPK_pro_n <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n)

zscores <- data.frame(values=c("BEG_1230-1_T1227"),
                      zscores_DMPK_genes_in_jaccard=c(DMPK_pro_zscore),
                      number_DMPK_genes_in_jaccard = c(DMPK_pro_n))

write_csv(zscores, dmpk_dmpk_genes_stats)

################-----Number Minor Intron Retention Events RNU4ATAC Labeled-----#################
number_DMPK_junctions_jaccard <- all_count %>% filter(type %in% c("jaccard")) %>% filter(hgncSymbol %in% ensg) %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))
joined <- number_DMPK_junctions_jaccard
number_DMPK_junctions_jaccard_genes <- all_count %>% filter(type %in% c("jaccard")) %>% filter(hgncSymbol %in% ensg) %>% select(sampleID, hgncSymbol) %>% unique %>% group_by(sampleID) %>% tally() %>% arrange(desc(n))

joined$numbers <- rownames(joined)
title_lab <- expression(psi ~ 3 ~ "and" ~ psi ~ "5 in DMPK Mis-Spliced Genes")
x_lab <- expression(atop("Samples ordered by Number of Junctions with Significant",
                     psi ~ 3 ~ "and" ~ psi ~ "5 Outliers in DMPK Mis-Spliced Genes"))
y_lab <- expression(atop("Number of Junctions with Significant" ~ psi ~ 3 ~ "and" ~ psi ~ "5 Outliers",
                        "in DMPK Mis-Spliced Genes"))
joined$numbers <- rownames(joined)

DMPK <- joined %>% filter(sampleID == "BEG_1230-1_T1227")
DMPK$type <- "DMPK Pathogenic Tandem Repeat"
non_DMPK <- joined %>% filter(sampleID != "BEG_1230-1_T1227")
non_DMPK$type <- NA

joined <- bind_rows(DMPK, non_DMPK)

top_10 <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score

DMPK_jaccard <- ggplot(joined, aes(x = reorder(numbers, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="limegreen") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  geom_text_repel(data=(joined %>% filter(type %in% top_10)), size=8, aes(x = reorder(type, n), y = n, label = type), position=position_dodge(width=0.9),
                    hjust=1.2,colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(n), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  geom_hline(data=joined, aes(yintercept=mean(n)-sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(n)+sd(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$n)+2, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$n) + sd(joined$n) + sd(joined$n) + sd(joined$n)+2),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$n) - sd(joined$n)+2),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))

DMPK_jaccard

ggsave(filename=dmpk_dmpk_jaccard, plot=DMPK_jaccard,  limitsize = FALSE, units = "in")

################-----Get Stats-----#################
DMPK_pro_jac <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n) - mean(joined$n)
DMPK_pro_zscore_jac <- DMPK_pro_jac/ sd(joined$n)
DMPK_pro_n_jac <- joined %>% filter(sampleID == "BEG_1230-1_T1227") %>% pull(n)

not_DMPK_jac <- joined %>% filter(sampleID != "BEG_1230-1_T1227") 
mean_not_DMPK_jac <- mean(not_DMPK_jac$n)
sd_not_DMPK_jac <- sd(not_DMPK_jac$n)

zscores <- data.frame(values=c("BEG_1230-1_T1227"),
                      zscores_DMPK_genes_in_jaccard=c(DMPK_pro_zscore_jac),
                      number_DMPK_genes_in_jaccard = c(DMPK_pro_n_jac))

write_csv(zscores, dmpk_dmpk_genes_stats_jaccard)


means <- data.frame(values=c("Jaccard_non_DMPK", "PSI3_PSI5_DMPK"),
                    mean=c(mean_not_DMPK_jac, mean_not_DMPK),
                    sd=c(sd_not_DMPK_jac, sd_not_DMPK))

write_csv(means, dmpk_dmpk_genes_stats_jaccard)







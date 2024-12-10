library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(tidyverse)

twenty5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_25_filtered_stanford/FRASER_output.csv")
twenty5_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_25_filtered_stanford/FRASER_input.csv")
fifty <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_50_filtered_stanford/FRASER_output.csv")
fifty_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_50_filtered_stanford/FRASER_input.csv")
one00  <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_100_filtered_stanford/FRASER_output.csv")
one00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_100_filtered_stanford/FRASER_input.csv")
one50 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_150_filtered_stanford/FRASER_output.csv")
one50_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_150_filtered_stanford/FRASER_input.csv")
two00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_200_filtered_stanford/FRASER_output.csv")
two00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_200_filtered_stanford/FRASER_input.csv")
three00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_300_filtered_stanford/FRASER_output.csv")
three00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_300_filtered_stanford/FRASER_input.csv")
four00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_400_filtered_stanford/FRASER_output.csv")
four00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER_snakemake_old_filters/output/Subset_400_filtered_stanford/FRASER_input.csv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

########################################################################
######################-----Make Dataframe-----##########################
########################################################################

######################-----25-----##########################
get_number_MIGs_theta <- function(dataframe, files, col_name){
        MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
        colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
        MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG")

        #get theta dataframe
        results_filtered_theta <- dataframe %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta")
        results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol")
        colnames(results_filtered_theta) <- c("sampleID", "gene_symbol")

        #combine dataframes
        sig_genes_type_all_samples <- left_join(MIG_table_select, results_filtered_theta) %>% filter(!is.na(sampleID))

        joined <- sig_genes_type_all_samples %>% select(sampleID, gene_symbol) %>% 
            group_by(sampleID) %>% tally() %>% arrange(desc(n))

        missing_all <- setdiff(files$sampleID, joined$sampleID)
        missing_all_df <- data.frame(sampleID = missing_all, n = 0)
        MIG_count <- bind_rows(joined, missing_all_df)

        MIG_colname <- paste("no_MIGs_theta_juncs", col_name, sep="_")
        colnames(MIG_count) <- c("sampleID", MIG_colname)
        MIG_count

}

twenty5_MIGs <- get_number_MIGs_theta(twenty5, twenty5_files, "25") 
colnames(twenty5_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
fifty_MIGs <- get_number_MIGs_theta(fifty, fifty_files, "50") 
colnames(fifty_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
one00_MIGs <- get_number_MIGs_theta(one00, one00_files, "100") 
colnames(one00_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
one50_MIGs <- get_number_MIGs_theta(one50, one50_files, "150")
colnames(one50_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
two00_MIGs <- get_number_MIGs_theta(two00, two00_files, "200") 
colnames(two00_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
three00_MIGs <- get_number_MIGs_theta(three00, three00_files, "300") 
colnames(three00_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")
four00_MIGs <- get_number_MIGs_theta(four00, four00_files, "400") 
colnames(four00_MIGs) <- c("sampleID", "no_theta_juncs_in_MIGs")

########################################################################
######################-----Plot-----##########################
########################################################################

######################-----25-----##########################
x_lab <- expression(atop("Samples Ordered by Number of",
                       ~ theta ~ "Outlier Junctions in MIGs"))
y_lab <- expression("Number of" ~ theta~  "Outlier Junctions in MIGs")
twenty5_MIGs$numbers <- rownames(twenty5_MIGs)

RNU4ATAC <- twenty5_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- twenty5_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
four <- mean(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)+sd(joined$no_theta_juncs_in_MIGs)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_25 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=four, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 4*SD", x = 1, y=four +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_25

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/25.pdf", plot=MIGs_RNU4ATAC_25,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----50-----##########################
fifty_MIGs$numbers <- rownames(fifty_MIGs)

RNU4ATAC <- fifty_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- fifty_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
six <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*6)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_50 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=six, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 6*SD", x = 1, y=six +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_50

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/50.pdf", plot=MIGs_RNU4ATAC_50,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----100-----##########################
one00_MIGs$numbers <- rownames(one00_MIGs)

RNU4ATAC <- one00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- one00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
nine <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*9)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_100 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=nine, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 9*SD", x = 1, y=nine +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_100

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/100.pdf", plot=MIGs_RNU4ATAC_100,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----150-----##########################
one50_MIGs$numbers <- rownames(one50_MIGs)

RNU4ATAC <- one50_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- one50_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
eleven <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*11)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_150 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=eleven, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 11*SD", x = 1, y=eleven +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_150

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/150.pdf", plot=MIGs_RNU4ATAC_150,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----200-----##########################
two00_MIGs$numbers <- rownames(two00_MIGs)

RNU4ATAC <- two00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- two00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
thir <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*13)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_200 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=thir, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 13*SD", x = 1, y=thir +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_200

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/200.pdf", plot=MIGs_RNU4ATAC_200,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----200-----##########################
three00_MIGs$numbers <- rownames(three00_MIGs)

RNU4ATAC <- three00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- three00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
sixt <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*16)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_300 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=sixt, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 16*SD", x = 1, y=sixt +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_300

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/300.pdf", plot=MIGs_RNU4ATAC_300,  limitsize = FALSE, units = "in", height=12, width=10)

######################-----200-----##########################
four00_MIGs$numbers <- rownames(four00_MIGs)

RNU4ATAC <- four00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(sampleID %in% c("RD268"))
RNU4ATAC$type <- "A1"
non_RNU4ATAC <- four00_MIGs %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers) %>% filter(! sampleID %in% c("RD268"))
non_RNU4ATAC$type <- NA

joined <- bind_rows(RNU4ATAC, non_RNU4ATAC)

top_10 <- joined %>% arrange(desc(no_theta_juncs_in_MIGs)) %>% select(sampleID, no_theta_juncs_in_MIGs, numbers, type) %>% head(1) %>% pull(type)
g1 <- subset(joined, type %in% top_10)

#check which distribution this is considering
#give a z-score
ninet <- mean(joined$no_theta_juncs_in_MIGs)+(sd(joined$no_theta_juncs_in_MIGs)*19)

colours <- c("#D7DCD8", "#DA786C")
names(colours) <- c(NA, "A1")

MIGs_RNU4ATAC_400 <- ggplot(joined, aes(x = reorder(numbers, no_theta_juncs_in_MIGs), y = no_theta_juncs_in_MIGs, fill=type)) +
  geom_point(shape=21, alpha=.6, size=10, show.legend=FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  #labs(title=title_lab)  +
  #geom_point(data=g1, aes(fill= type), shape=21, alpha=.6, size=5, show.legend = FALSE) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(type == "RNU4ATAC-opathy")), size=8, aes(label = type), #position=position_dodge(width=0.9),
  #                  hjust=1.2, colour="black", max.overlaps = Inf) +
  geom_hline(data= joined, aes(yintercept=mean(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="solid", size=0.5)+
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5)+        
  #geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)-sd(no_theta_juncs_in_MIGs), color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=mean(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs)+sd(no_theta_juncs_in_MIGs), color = group),
              color="blue", linetype="dashed", size=0.5) +
  geom_hline(data=joined, aes(yintercept=ninet, color = group),
              color="blue", linetype="dashed", size=0.5) +
  #geom_hline(data=joined, aes(yintercept=eleven, color = group),
  #            color="blue", linetype="dashed", size=0.5) +
  geom_text(data = joined, aes(label = "Mean", x = 1, y=mean(joined$no_theta_juncs_in_MIGs)+4, hjust=.00025),  size=8) + 
  geom_text(data = joined, aes(label = "Mean + SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 2*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 3*SD", x = 1, y=mean(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs) + sd(joined$no_theta_juncs_in_MIGs)+4),  size=8, hjust=0.000025) +
  geom_text(data = joined, aes(label = "Mean + 19*SD", x = 1, y=ninet +4),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean + 11*SD", x = 25, y=eleven +10),  size=8, hjust=0.000025) +
  #geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$no_theta_juncs_in_MIGs) - sd(joined$no_theta_juncs_in_MIGs)+10),  size=8, hjust=0.000025) +
  theme_classic(base_size = 25)+
  coord_cartesian(clip="off")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  scale_fill_manual(values=colours)+
  expand_limits(y=0)


MIGs_RNU4ATAC_400

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/sample_size_MIGs/400.pdf", plot=MIGs_RNU4ATAC_400,  limitsize = FALSE, units = "in", height=12, width=10)

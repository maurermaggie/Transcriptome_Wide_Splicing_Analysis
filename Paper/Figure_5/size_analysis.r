library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

twenty5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_25_filtered/FRASER_output.csv")
twenty5_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_25_filtered/FRASER_input.csv")
fifty <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_50_filtered/FRASER_output.csv")
fifty_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_50_filtered/FRASER_input.csv")
one00  <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_100_filtered/FRASER_output.csv")
one00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_100_filtered/FRASER_input.csv")
one50 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_150_filtered/FRASER_output.csv")
one50_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_150_filtered/FRASER_input.csv")
two00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_200/FRASER_output.csv")
two00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_200/FRASER_input.csv")
three00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_300_filtered/FRASER_output.csv")
three00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_300_filtered/FRASER_input.csv")
four00 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_400_filtered/FRASER_output.csv")
four00_files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/Subset_400_filtered/FRASER_input.csv")

MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")

low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

########################################################################
######################-----Make Dataframe-----##########################
########################################################################

######################-----MIGs with IR in RD268-----##########################
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
        MIG_count <- MIG_count %>% filter(! sampleID %in% low_RIN$sampleID)

        MIG_colname <- paste("no_MIGs_theta_juncs", col_name, sep="_")
        colnames(MIG_count) <- c("sampleID", MIG_colname)
        MIG_count

}

twenty5_MIGs <- get_number_MIGs_theta(twenty5, twenty5_files, "25") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_25)
fifty_MIGs <- get_number_MIGs_theta(fifty, fifty_files, "50") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_50)
one00_MIGs <- get_number_MIGs_theta(one00, one00_files, "100") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_100)
one50_MIGs <- get_number_MIGs_theta(one50, one50_files, "150") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_150)
two00_MIGs <- get_number_MIGs_theta(two00, two00_files, "200") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_200)
three00_MIGs <- get_number_MIGs_theta(three00, three00_files, "300") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_300)
four00_MIGs <- get_number_MIGs_theta(four00, four00_files, "400") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_400)

number_MIGs_with_IR <- data.frame(sample_size = c(25, 50, 100, 150, 200, 300, 400), 
                            number_MIGs_theta=c(twenty5_MIGs, fifty_MIGs, one00_MIGs, one50_MIGs, two00_MIGs, three00_MIGs, four00_MIGs))

######################-----MIGs with IR mean/ person-----##########################
twenty5_MIGs_mean <- get_number_MIGs_theta(twenty5, twenty5_files, "25") 
mean_25 <- mean(twenty5_MIGs_mean$no_MIGs_theta_juncs_25)
sd_25 <- sd(twenty5_MIGs_mean$no_MIGs_theta_juncs_25)

fifty_MIGs_mean <- get_number_MIGs_theta(fifty, fifty_files, "50") 
mean_50 <- mean(fifty_MIGs_mean$no_MIGs_theta_juncs_50)
sd_50 <- sd(fifty_MIGs_mean$no_MIGs_theta_juncs_50)

one00_MIGs_mean <- get_number_MIGs_theta(one00, one00_files, "100")
mean_100 <-  mean(one00_MIGs_mean$no_MIGs_theta_juncs_100)
sd_100 <-  sd(one00_MIGs_mean$no_MIGs_theta_juncs_100)

one50_MIGs_mean <- get_number_MIGs_theta(one50, one50_files, "150") 
mean_150 <- mean(one50_MIGs_mean$no_MIGs_theta_juncs_150)
sd_150 <- sd(one50_MIGs_mean$no_MIGs_theta_juncs_150)

two00_MIGs_mean <- get_number_MIGs_theta(two00, two00_files, "200")
mean_200 <-  mean(two00_MIGs_mean$no_MIGs_theta_juncs_200)
sd_200 <-  sd(two00_MIGs_mean$no_MIGs_theta_juncs_200)

three00_MIGs_mean <- get_number_MIGs_theta(three00, three00_files, "300")
mean_300 <- mean(three00_MIGs_mean$no_MIGs_theta_juncs_300) 
sd_300 <- sd(three00_MIGs_mean$no_MIGs_theta_juncs_300) 

four00_MIGs_mean <- get_number_MIGs_theta(four00, four00_files, "400")
mean_400<-  mean(four00_MIGs_mean$no_MIGs_theta_juncs_400)
sd_400<-  sd(four00_MIGs_mean$no_MIGs_theta_juncs_400)

number_MIGs_with_IR_all <- data.frame(sample_size = c(25, 50, 100, 150, 200, 300, 400), 
                            mean_number_MIGs_theta=c(mean_25, mean_50, mean_100, mean_150, mean_200, mean_300, mean_400),
                            sd_number_MIGs_theta=c(sd_25, sd_50, sd_100, sd_150, sd_200, sd_300, sd_400))

all <- left_join(four00_MIGs_mean, three00_MIGs_mean) %>% left_join(two00_MIGs_mean) %>% left_join(one50_MIGs_mean) %>%
        left_join(one00_MIGs_mean) %>% left_join(fifty_MIGs_mean) %>% left_join(twenty5_MIGs_mean)

########################################################################
######################-----Plot-----##########################
########################################################################

######################-----MIGs with IR in RD268-----##########################
MIG_plot <- ggplot(number_MIGs_with_IR,aes(x=sample_size, y=number_MIGs_theta)) + 
  geom_point(method='lm', formula= sample_size~number_MIGs_theta) +
  xlab("Sample Size") +
  ylab("Number of MIGs with Theta Outliers in Proband 1") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5)) 

MIG_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_5/Size Analysis/MIG_analysis_RD268.pdf", plot=MIG_plot,  limitsize = FALSE, units = "in")

get_lm <- lm(sample_size ~ number_MIGs_theta, data = number_MIGs_with_IR)
print(summary(get_lm))

######################-----MIGs with IR in all-----##########################
MIG_plot_all <- ggplot(number_MIGs_with_IR_all,aes(x=sample_size, y=mean_number_MIGs_theta)) + 
  geom_point(method='lm', formula= sample_size~mean_number_MIGs_theta) +
  geom_errorbar(aes(ymin = mean_number_MIGs_theta-sd_number_MIGs_theta, ymax = mean_number_MIGs_theta+sd_number_MIGs_theta),width = 0.2) +
  xlab("Sample Size") +
  ylab("Number of MIGs with Theta Outliers in Proband 1") +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5))



MIG_plot_all

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_5/Size Analysis/MIG_analysis_all.pdf", plot=MIG_plot_all,  limitsize = FALSE, units = "in")

get_lm <- lm(sample_size ~ mean_number_MIGs_theta, data = number_MIGs_with_IR_all)
print(summary(get_lm))

######################-----MIGs with IR in all-----##########################
colnames(all) <- c("sampleID", "400", "300", "200", "150", "100", "50", "25")
long <- melt(setDT(all), id.vars = c("sampleID"), variable.name = "count")

top_10 <- long %>% filter(sampleID == "RD268") %>% pull(sampleID)
g1 <- subset(long, sampleID %in% top_10)

MIG_violin_plot<- ggplot(long,aes(x=count, y=value)) + 
  geom_violin() +
  geom_point(data=g1, aes(fill= sampleID), shape=21, alpha=.6, size=10) +
  labs(colour="sampleID") +
  xlab("Sample Size") +
  ylab("Number of MIGs with Theta Outliers in Proband 1") +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  theme_gray(base_size = 30)+
  theme(axis.text.x = element_text( vjust = 1))+
  stat_summary(fun.y=mean, geom="point", size=2, color="red")


MIG_violin_plot

ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Figure_5/Size Analysis/MIG_violin_plot.pdf", plot=MIG_violin_plot,  limitsize = FALSE, units = "in")


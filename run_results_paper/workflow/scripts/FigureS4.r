library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(scales)
library(tidyverse)

args <- commandArgs(TRUE)
size_dir <- args[1]
MIGs_fp <- args[2]
output_file <- args[3]

twenty5_fp <- paste(size_dir, "Subset_25_filtered_stanford/FRASER_output.csv", sep="/")
twenty5 <- read_csv(twenty5_fp)
twenty_f_fp <- paste(size_dir, "Subset_25_filtered_stanford/FRASER_input.csv", sep="/")
twenty5_files <- read_csv(twenty_f_fp)

fifty_fp <- paste(size_dir, "Subset_50_filtered_stanford/FRASER_output.csv", sep="/")
fifty <- read_csv(fifty_fp)
fifty_f_fp <- paste(size_dir, "Subset_50_filtered_stanford/FRASER_input.csv", sep="/")
fifty_files <- read_csv(fifty_f_fp)

one00_fp <- paste(size_dir, "Subset_100_filtered_stanford/FRASER_output.csv", sep="/")
one00  <- read_csv(one00_fp)
one00_f_fp <- paste(size_dir, "Subset_100_filtered_stanford/FRASER_input.csv", sep="/")
one00_files <- read_csv(one00_f_fp)

one50_fp <- paste(size_dir, "Subset_150_filtered_stanford/FRASER_output.csv", sep="/")
one50  <- read_csv(one50_fp)
one50_f_fp <- paste(size_dir, "Subset_150_filtered_stanford/FRASER_input.csv", sep="/")
one50_files <- read_csv(one50_f_fp)

two00_fp <- paste(size_dir, "Subset_200_filtered_stanford/FRASER_output.csv", sep="/")
two00  <- read_csv(two00_fp)
two00_f_fp <- paste(size_dir, "Subset_200_filtered_stanford/FRASER_input.csv", sep="/")
two00_files <- read_csv(two00_f_fp)

three00_fp <- paste(size_dir, "Subset_300_filtered_stanford/FRASER_output.csv", sep="/")
three00  <- read_csv(three00_fp)
three00_f_fp <- paste(size_dir, "Subset_300_filtered_stanford/FRASER_input.csv", sep="/")
three00_files <- read_csv(three00_f_fp)

three00_fp <- paste(size_dir, "Subset_300_filtered_stanford/FRASER_output.csv", sep="/")
three00  <- read_csv(three00_fp)
three00_f_fp <- paste(size_dir, "Subset_300_filtered_stanford/FRASER_input.csv", sep="/")
three00_files <- read_csv(three00_f_fp)

four00_fp <- paste(size_dir, "Subset_400_filtered_stanford/FRASER_output.csv", sep="/")
four00  <- read_csv(four00_fp)
four00_f_fp <- paste(size_dir, "Subset_400_filtered_stanford/FRASER_input.csv", sep="/")
four00_files <- read_csv(four00_f_fp)

MIGs <- read_csv(MIGs_fp)

lm_MIGs_fp <- args[4]
lm_MIGs_z_fp <- args[5]

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

######################-----MIGs with IR in RD268-----##########################
get_number_MIGs_genes_theta <- function(dataframe, files, col_name){
        MIG_table_select <- MIGs %>% select(gene_symbol, gene_class) %>% unique
        colnames(MIG_table_select) <- c("gene_symbol", "gene_class")
        MIG_table_select <- MIG_table_select %>% filter(gene_class == "MIG")

        #get theta dataframe
        results_filtered_theta <- dataframe %>% filter(!is.na(hgncSymbol)) %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(type == "theta")
        results_filtered_theta <- results_filtered_theta %>% select("sampleID", "hgncSymbol") %>% unique()
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

twenty5_MIGs_genes <- get_number_MIGs_genes_theta(twenty5, twenty5_files, "25") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_25)
fifty_MIGs_genes <- get_number_MIGs_genes_theta(fifty, fifty_files, "50") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_50)
one00_MIGs_genes <- get_number_MIGs_genes_theta(one00, one00_files, "100") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_100)
one50_MIGs_genes <- get_number_MIGs_genes_theta(one50, one50_files, "150") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_150)
two00_MIGs_genes <- get_number_MIGs_genes_theta(two00, two00_files, "200") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_200)
three00_MIGs_genes <- get_number_MIGs_genes_theta(three00, three00_files, "300") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_300)
four00_MIGs_genes <- get_number_MIGs_genes_theta(four00, four00_files, "400") %>% filter(sampleID=="RD268") %>% pull (no_MIGs_theta_juncs_400)

number_MIGs_with_IR_genes <- data.frame(sample_size = c(25, 50, 100, 150, 200, 300, 400), 
                            number_MIGs_theta=c(twenty5_MIGs_genes, fifty_MIGs_genes, one00_MIGs_genes, one50_MIGs_genes, two00_MIGs_genes, three00_MIGs_genes, four00_MIGs_genes))

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

all$mean_400 <- mean(all$`no_MIGs_theta_juncs_400`)
all$mean_300 <- mean(all$`no_MIGs_theta_juncs_300`, na.rm = T)
all$mean_200 <- mean(all$`no_MIGs_theta_juncs_200`, na.rm = T)
all$mean_150 <- mean(all$`no_MIGs_theta_juncs_150`, na.rm = T)
all$mean_100 <- mean(all$`no_MIGs_theta_juncs_100`, na.rm = T)
all$mean_50 <- mean(all$`no_MIGs_theta_juncs_50`, na.rm = T)
all$mean_25 <- mean(all$`no_MIGs_theta_juncs_25`, na.rm = T)

all$sd_400 <- sd(all$`no_MIGs_theta_juncs_400`)
all$sd_300 <- sd(all$`no_MIGs_theta_juncs_300`, na.rm = T)
all$sd_200 <- sd(all$`no_MIGs_theta_juncs_200`, na.rm = T)
all$sd_150 <- sd(all$`no_MIGs_theta_juncs_150`, na.rm = T)
all$sd_100 <- sd(all$`no_MIGs_theta_juncs_100`, na.rm = T)
all$sd_50 <- sd(all$`no_MIGs_theta_juncs_50`, na.rm = T)
all$sd_25 <- sd(all$`no_MIGs_theta_juncs_25`, na.rm = T)

number_MIGs_with_IR$mean <- c(unique(all$mean_25), unique(all$mean_50),unique(all$mean_100), unique(all$mean_150),
                              unique(all$mean_200), unique(all$mean_300),unique(all$mean_400))

number_MIGs_with_IR$sd <- c(unique(all$sd_25), unique(all$sd_50),unique(all$sd_100), unique(all$sd_150),
                              unique(all$sd_200), unique(all$sd_300),unique(all$sd_400))

y_lab <- paste("Individual A1:", "\n", "Number of Significant Intron Retention Events in MIGs")

number_MIGs_with_IR$zscore <- (number_MIGs_with_IR$number_MIGs_theta - number_MIGs_with_IR$mean)/number_MIGs_with_IR$sd
MIG_plot <- ggplot(number_MIGs_with_IR, aes(x=sample_size, y=number_MIGs_theta, fill=zscore)) + 
  geom_point(shape=21, size=20, colour="black", aes(y=number_MIGs_theta, x = reorder(sample_size, number_MIGs_theta))) +
  #geom_point(shape=21, size=10, fill="#B8CFCB", colour="black", aes(y=mean, x = reorder(sample_size, number_MIGs_theta)))+
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, x = reorder(sample_size, number_MIGs_theta)),width = 0.2) +
  xlab("Sample Size") +
  ylab(y_lab) +
  #geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 27)+
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(x = 0, y = 0)+
  scale_fill_gradient2(position="bottom" , low = "#3A615D", mid = "#B4B689", high = "#F3D177", 
                       midpoint = median(number_MIGs_with_IR$zscore)) 

MIG_plot

ggsave(filename=output_file, plot=MIG_plot,  limitsize = FALSE, units = "in")

sink(file=lm_MIGs_fp)
get_lm <- lm(sample_size ~ number_MIGs_theta, data = number_MIGs_with_IR)
summary(get_lm)
sink()

sink(file=lm_MIGs_z_fp)
get_lm_z <- lm(zscore ~ number_MIGs_theta, data = number_MIGs_with_IR)
summary(get_lm_z)
sink()
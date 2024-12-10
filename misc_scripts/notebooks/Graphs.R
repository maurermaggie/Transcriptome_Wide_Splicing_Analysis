#library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)

joined <- read_csv("/home/maurertm/smontgom/shared/UDN/FRASER_analysis/scripts/FRASER_snakemake/output/FRASER_all/FRASER_results_raw.csv")
results <- read_csv("/home/maurertm/smontgom/shared/UDN/FRASER_analysis/scripts/FRASER_snakemake/output/FRASER_all/FRASER_output.csv")

########################################################
################-----Plot Function-----#################
########################################################
plot_number_junctions <- function(sample, column_name, variable, type, output_directory) {

  y_lab <- paste("Number of", type,  "Junctions with", variable, sep=" ")
  title_lab <- paste("Number of", type, "Junctions with", variable, "\n", "in the UDN and GREGoR Stanford Sites", sep=" ")
  file_name <- paste(type, "_", variable, ".pdf", sep="")
  filepath <- paste(output_directory, file_name, sep="/")

  joined$variable_of_interest <- column_name
  g1 <- subset(joined, Sample_ID == sample)

  top_10 <- joined %>% arrange(desc(variable_of_interest)) %>% select(Sample_ID, variable_of_interest) %>% head(10) %>% pull(Sample_ID)
  
  #check which distribution this is considering
  #give a z-score
  function_plot <- ggplot(joined, aes(x = reorder(Sample_ID, variable_of_interest), y = variable_of_interest)) +
    geom_point() +
    xlab("Sample_ID") +
    ylab(y_lab) +
    labs(title=title_lab)  +
    geom_point(data=g1, aes(colour= Sample_ID)) +
    labs(colour="Sample_ID") +
    geom_text_repel(data=(joined %>% filter(Sample_ID %in% top_10)), size=3, aes(x = reorder(Sample_ID, variable_of_interest), y = variable_of_interest, label = Sample_ID), position=position_dodge(width=0.9),
                    hjust=1.5,colour="black", max.overlaps = Inf) +
    geom_hline(data= joined, aes(yintercept=mean(variable_of_interest), color = group),
              color="blue", linetype="solid", size=0.5)+
    geom_hline(data=joined, aes(yintercept=mean(variable_of_interest)+sd(variable_of_interest), color = group),
              color="blue", linetype="dashed", size=0.5)+ 
    geom_hline(data=joined, aes(yintercept=mean(variable_of_interest)+sd(variable_of_interest)+sd(variable_of_interest), color = group),
              color="blue", linetype="dashed", size=0.5)+        
    geom_hline(data=joined, aes(yintercept=mean(variable_of_interest)-sd(variable_of_interest), color = group),
              color="blue", linetype="dashed", size=0.5) +
    geom_hline(data=joined, aes(yintercept=mean(variable_of_interest)+sd(variable_of_interest)+sd(variable_of_interest)+sd(variable_of_interest), color = group),
              color="blue", linetype="dashed", size=0.5) +
    geom_text(data = joined, aes(label = "Mean", x = 25, y=mean(joined$variable_of_interest)+8),  size=2) + 
    geom_text(data = joined, aes(label = "Mean + SD", x = 25, y=mean(joined$variable_of_interest) + sd(joined$variable_of_interest)+8),  size=2) +
    geom_text(data = joined, aes(label = "Mean + 2*SD", x = 25, y=mean(joined$variable_of_interest) + sd(joined$variable_of_interest) + sd(joined$variable_of_interest)+8),  size=2) +
    geom_text(data = joined, aes(label = "Mean + 3*SD", x = 25, y=mean(joined$variable_of_interest) + sd(joined$variable_of_interest) + sd(joined$variable_of_interest) + sd(joined$variable_of_interest)+8),  size=2) +
    geom_text(data = joined, aes(label = "Mean - SD", x = 25, y=mean(joined$variable_of_interest) - sd(joined$variable_of_interest)+8),  size=2) +
    theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4))

  function_plot

  #joined <- joined %>% select(-variable_of_interest)
  function_plot
  ggsave(filename=filepath, plot=function_plot,  limitsize = FALSE, units = "in")
}

########################################################
##################-----Plot Graphs-----#################
########################################################
joined_values <- joined %>% select(-c(Sample_ID, Individual_ID))
joined_values_colnames <- colnames(joined_values)

for (i in joined_values_colnames) {
  if (grepl(i, "^psi3_psi5$"))
  {
    Type <- "Psi3_Psi5"
  } else if (grepl("theta", i) == TRUE) {
    Type <- "Theta"
  } else if (grepl("psi3", i) == TRUE) {
    Type <- "Psi3"
  } else if (grepl("psi5", i) == TRUE) {
    Type <- "Psi5"
  } else {
    print(i)
  }

  if (grepl("juncs$", i) == TRUE) 
  {
    Variable <- "padjust < 0.05"
  } else if (grepl("abs_z_score", i) == TRUE) {
    Variable <- "Abs(Z Score) > 2" 
  } else if (grepl("z_score_under", i) == TRUE) {
    Variable <- "Z Score < 2"
  } else if (grepl("z_score_over", i) == TRUE) {
    Variable <- "Z Score > 2" 
  } else if (grepl("pval_1E_neg6$", i) == TRUE) {
    Variable <- "padjust < 1E-6"
  } else if (grepl("pval_1E_neg6_abs_deltaPsi", i) == TRUE) {
    Variable <- "padjust < 1E-6 and abs(deltaPsi) > 0.3"
  } else if (grepl("abs_deltaPsi", i) == TRUE) {
    Variable <- "abs(deltaPsi) > 0.3"
  } else if (grepl("deltaPsi_over", i) == TRUE) {
    Variable <- "deltaPsi > 0.3" 
  } else if (grepl("deltaPsi_under", i) == TRUE) {
    Variable <- "deltaPsi < 0.3"
  } else if (grepl("^psi3_psi5$", i) == TRUE) {
    Variable <- "padjust < 0.05"
  } else {
    print(i)
  }

  if (grepl("MIGs", i) == TRUE) 
  {
    Variable <- paste(Variable, "\n", " in Minor Intron Containing Genes", sep = "")
  }

  column <- joined %>% select(all_of(i))
  print("column")
  colnames(column) <- "filler"
  print(column)
  plot_number_junctions("GSS225379", column$filler, Variable, Type, "/oak/stanford/groups/smontgom/shared/UDN/FRASER_analysis/scripts/FRASER_snakemake/images/")
}

 
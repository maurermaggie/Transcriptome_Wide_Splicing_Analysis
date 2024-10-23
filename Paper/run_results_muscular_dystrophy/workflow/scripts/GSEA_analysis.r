library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

psi3 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi3_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
psi5 <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_psi5_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
theta <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_p_value_0.3_deltapsi_0.1_muscle_theta_185_samples_pdj_0.3_deltapsi_0.1_results.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Muscle_results/filtered_muscle_jaccard_results_jaccard.csv")
output_file <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/run_results_muscular_dystrophy/output/plots/GSEA_analysis.pdf"

all <- bind_rows(psi3, psi5, theta, jaccard)

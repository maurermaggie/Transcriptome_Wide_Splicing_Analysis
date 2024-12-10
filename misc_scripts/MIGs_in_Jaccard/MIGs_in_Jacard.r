library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_no_sibs/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")
MIGs <- read_csv("/home/maurertm/smontgom/shared/UDN/ReferenceFiles/Minor_Intron_Database/Homo_sapiens_gene.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_no_sibs/FRASER2/FRASER_ouput_filtered_padjust_deltaPsi.csv")


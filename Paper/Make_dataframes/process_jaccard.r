library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(sqldf)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")
jaccard <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_no_sibs/FRASER2/FRASER_ouput2_filtered.csv", col_names=FALSE)
colnames(jaccard) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi","counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")
jaccard_filtered <- jaccard %>% filter(padjust <=0.05) %>% filter(abs(deltaPsi) >= 0.3)
write_csv(jaccard_filtered, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_no_sibs/FRASER2/FRASER_ouput_filtered_padjust_deltaPsi.csv")

#jaccard <- fread("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_Sept_13/FRASER2/FRASER_output.csv")
jaccard_fp <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_Sept_13/FRASER2/FRASER_output.csv"
jaccard <- read.csv.sql(
    file = jaccard_fp,
    sql = "
       select
        V1, V2, V3, V6, V7, V8, V10, V11, V12
       from
           file
       where V10 < 0.051",
    header = FALSE,
    sep = ",",
    eol = "\n",
    `field.types` = list(
        col1 = c("V1"),
        col2 = c("V2"),
        col3 = c("V3"),
        col4 = c("V4"),
        col5 = c("V5"),
        col6 = c("V6"),
        col7 = c("V7"),
        col8 = c("V8"),
        col9 = c("V9"),
        col10 = c("V10"),
        col11 = c("V11"),
        col12 = c("V12"),
        col13 = c("V13"),
        col14 = c("V14"),
        col15 = c("V15"),
        col16 = c("V16"),
        col17 = c("V17"),
        col18 = c("V18"),
        col19 = c("V19")
    ),
    dbname = tempfile(),
    drv = "SQLite")

file_in <- file("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake/output/FRASER2_Sept_13/FRASER2/FRASER_output.csv","r")
chunk_size <- 2000000 # choose the best size for you
first <- readLines(file_in, n=chunk_size)
first <- strsplit(first[], ",")
first <- do.call(rbind, first)
first <- as.data.frame(first)
colnames(first) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi","counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")
first$deltaPsi <- as.numeric(first$deltaPsi)
first <- first %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
write_csv(first, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Make_dataframes/jaccard.csv")

second <- readLines(file_in,skip=2000000,n=2000000)
second <- strsplit(second[], ",")
second <- do.call(rbind, second)
second <- as.data.frame(second)
colnames(second) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi","counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")
second$deltaPsi <- as.numeric(second$deltaPsi)
second <- second %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
first_second <- bind_rows(first, second)
write_csv(first_second, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Make_dataframes/jaccard.csv")



third <- readLines(file_in,skip=4000000,n=2000000)
fourth <- readLines(file_in,skip=6000000,n=2000000)
fifth <- readLines(file_in,skip=8000000,n=2000000)
sixth <- readLines(file_in,skip=10000000,n=2000000)
seventh <- readLines(file_in,skip=12000000,n=2000000)
eighth <- readLines(file_in,skip=14000000,n=2000000)
ninth <- readLines(file_in,skip=16000000,n=2000000)
tenth <- readLines(file_in,skip=18000000,n=2000000)
eleventh <- readLines(file_in,skip=20000000,n=2000000)
twelvth <- readLines(file_in,skip=22000000,n=2000000)
thirtenth <- readLines(file_in,skip=24000000,n=2000000)
fourtenth <- readLines(file_in,skip=26000000,n=2000000)
fiftenth <- readLines(file_in,skip=28000000,n=2000000)
sixtenth <- readLines(file_in,skip=30000000,n=2000000)
sevententh <- readLines(file_in,skip=32000000,n=2000000)
eightenth <- readLines(file_in,skip=34000000,n=2000000)
ninetenth <- readLines(file_in,skip=36000000,n=2000000)






twenty <- readLines(file_in,skip=38000000,n=2000000)
twenty1 <- readLines(file_in,skip=40000000,n=2000000)
twenty2 <- readLines(file_in,skip=42000000,n=2000000)
twenty3 <- readLines(file_in,skip=44000000,n=2000000)
twenty4 <- readLines(file_in,skip=46000000,n=2000000)
twenty5 <- readLines(file_in,skip=48000000,n=2000000)
twenty6 <- readLines(file_in,skip=50000000,n=2000000)
twenty7 <- readLines(file_in,skip=52000000,n=2000000)


df1 <- strsplit(x[], ",")
df1 <- do.call(rbind, df1)
df1 <- as.data.frame(df1)

colnames(df1) <- c("seqnames","start","end","width","strand","sampleID","hgncSymbol","type","pValue","padjust","psiValue","deltaPsi","counts","totalCounts","meanCounts","meanTotalCounts","nonsplitCounts","nonsplitProportion","nonsplitProportion_99quantile")

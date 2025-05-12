library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db

FRASER2 <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER_snakemake_old_filters/output/FRASER_no_missing_rm_seqbatch_8/FRASER2/FRASER_filtered.rds")
#FRASER1 <- read_rds("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_no_missing_rm_seqbatch_8/FRASER_filtered.rds")


fds_annotated <- annotateRangesWithTxDb(FRASER2, txdb=txdb, orgDb=orgDb)

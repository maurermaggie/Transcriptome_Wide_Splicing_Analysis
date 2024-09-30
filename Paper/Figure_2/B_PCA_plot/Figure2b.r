#library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
#library(ggforce)
library(GenomicRanges)
library(patchwork)
library(reshape)
library(GO.db)
library(tidyverse)


all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")

################-----Make PCA Matrix-----#################
all_filtered <- all_uncompiled %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(abs(zScore) >= 2) %>% filter(!is.na(hgncSymbol))
all_select <- all_filtered %>% select(sampleID, hgncSymbol, counts, start, end)

all_select$junction <- paste(all_select$hgncSymbol,":",all_select$start,"-",all_select$end, sep="")

all_select <- all_select %>% select(-hgncSymbol, -start, -end) %>% unique

all_select_wide <- dcast(all_select, sampleID ~ junction, value.var = "counts")
rownames(all_select_wide) <- all_select_wide$sampleID
all_select_wide <- all_select_wide %>% select(-sampleID)
all_select_wide[is.na(all_select_wide)] <- 0

all_normalized <- scale(all_select_wide)

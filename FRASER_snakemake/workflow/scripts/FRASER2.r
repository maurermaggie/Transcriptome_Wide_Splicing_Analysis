library(FRASER)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyverse)

args <- commandArgs(TRUE)
input <- args[1]
output_directory <- args[2]
print(input)

#you HAVE to include as.data.table here or the "settings <-" command will NOT work
sampleTable <- read_csv(input) %>% as.data.table
sampleTable <- sampleTable
bamFiles <- sampleTable[,2]
sampleTable[,bamFile:=bamFiles]

settings <- FraserDataSet(colData=sampleTable, workingDir=output_directory)

#if(.Platform$OS.type == "unix") {
#register(MulticoreParam(workers=min(30, multicoreWorkers())))
#} else {
#register(SnowParam(workers=min(30, multicoreWorkers())))
#}

#this next line is the only part of the script that could be parallelized
#but it would be complicated because you would need to grab the n number from the sample file and makes sure thay all have the same
#output directory
fds <- countRNAData(settings, recount = TRUE)

fds <- calculatePSIValues(fds)

fds_filtered <- filterExpressionAndVariability(fds, minExpressionInOneSample=20, minDeltaPsi=0.0, filter=FALSE)

#you would think that these 6 commands could be parallelized, but if you do, you will have to save the data (which takes hours) and reload it after the
#6 lines (which also takes hours), so it's actually faster not to parallelize this. 
jaccard <- optimHyperParams(fds_filtered, type="jaccard", plot=FALSE)
jaccard_q <- bestQ(jaccard, type="jaccard")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_annotated <- annotateRangesWithTxDb(fds_filtered, txdb=txdb, orgDb=orgDb)

# retrieve the estimated optimal dimension of the latent space
fds_fraser <- FRASER(fds_annotated, q=c(jaccard=jaccard_q))

results_fraser2 <- results(fds_fraser, all=TRUE)
results_df <- as.data.frame(results_fraser2)

filename <- paste0(output_directory, "/FRASER_output.csv", sep="")
write_csv(results_df, filename)
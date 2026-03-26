library(GenomicFeatures)  # TxDb, intronsByTranscript
library(dplyr)
library(tibble)
library(readxl)
library(magrittr)
library(biomaRt)
require(data.table)
library(tidyverse)

# 1) Download a matching GTF for your genome build
#    e.g. Ensembl GRCh38 or GENCODE; make sure it matches MIDB's assembly. [web:70][web:76]
gtf_file <- "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/hg38.ensGene.gtf"  # path to your GTF

# 2) Your MIG list from MIDB (as a vector of Ensembl gene IDs, e.g. ENSG...)
#    Suppose you already have something like:
mig_genes <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/Homo_sapiens_intron.csv") %>% 
             filter(intron_class == "minor") %>%
             pull(ensembl_gene_id) %>% unique

#################################################################################
################----Build TxDb and get introns per transcript----################
#################################################################################
# Build TxDb from the GTF (only needs to be done once per GTF)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract introns grouped by transcript [web:71][web:80]
itr_by_tx <- intronsByTranscript(txdb, use.names = TRUE)

# Convert to a data frame: one row per intron
itr_df <- as.data.frame(itr_by_tx)

# itr_df will have columns like: seqnames, start, end, width, strand, tx_id, tx_name
# tx_name is usually the transcript ID (ENST...), tx_id is an internal numeric ID. [web:68]

introns_per_tx <- itr_df %>%
  as_tibble() %>%
  count(group_name, name = "n_introns_tx")

#################################################################################
#########################----Map transcripts to genes----########################
#################################################################################
# Get transcript-to-gene mapping from TxDb [web:71]
tx2gene <- transcripts(txdb, columns = c("tx_name", "gene_id")) %>%
  as_tibble() %>%
  dplyr::select(tx_name, gene_id)

# Join intron counts to their genes
colnames(introns_per_tx) <- c("tx_name", "n_introns_tx")

tx_introns_gene <- introns_per_tx %>%
  inner_join(tx2gene, by = "tx_name")

#################################################################################
##################----Collapse to one intron count per gene----##################
#################################################################################

# 1) Max introns across transcripts per gene (often robust)
tx_introns_gene$gene_id <- unlist(tx_introns_gene$gene_id)
introns_per_gene_max <- tx_introns_gene %>%
  group_by(gene_id) %>%
  summarise(n_introns_gene_max = max(n_introns_tx), .groups = "drop")

# 2) Mean introns per transcript per gene
introns_per_gene_mean <- tx_introns_gene %>%
  group_by(gene_id) %>%
  summarise(n_introns_gene_mean = mean(n_introns_tx), .groups = "drop")

# 3) Canonical transcript only – if you prefer Ensembl’s canonical, you’d first
#    get that list via BioMart or Ensembl REST and filter tx_introns_gene to those. [web:42][web:57]

#################################################################################
##############----Restrict to MIGs and compute the MIG average----###############
#################################################################################
# Keep only MIG genes
migs_introns <- introns_per_gene_max %>%
  filter(gene_id %in% mig_genes)

# Average total introns across all MIGs
mean_introns_migs <- mean(migs_introns$n_introns_gene_max, na.rm = TRUE)

mean_introns_migs_df <- data.frame(mean="number_of_introns_in_MIGs", value=mean_introns_migs)
write_csv(migs_introns, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/MIGs_introns_summary.csv")

write_csv(mean_introns_migs_df, "/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Fibroblast_RNU6ATAC/mean_no_introns_in_MIGs.csv")

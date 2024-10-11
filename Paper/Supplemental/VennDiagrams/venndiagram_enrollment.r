library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ggVennDiagram)
library(tidyverse)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata_updated_batches.csv")

########################################################################
###################-----Venn Diagram Enrollment-----####################
########################################################################
GSS <- metadata %>% filter(Enrollment_Type %in% c("GSS", "GSS_UDN")) %>% pull(ID) %>% unique
UDN <- metadata %>% filter(Enrollment_Type %in% c("UDN", "GSS_UDN"))  %>% pull(ID) %>% unique

enrollment_list <- list(`GREGoR`= GSS, `UDN`=UDN)
venn_enrollment_all <- ggVennDiagram(enrollment_list) + scale_fill_gradient(low="grey90",high = "red")+
                    scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Paper/Supplemental/VennDiagrams/Enrollment_VennDiagram.pdf", plot=venn_enrollment_all,  limitsize = FALSE, units = "in")


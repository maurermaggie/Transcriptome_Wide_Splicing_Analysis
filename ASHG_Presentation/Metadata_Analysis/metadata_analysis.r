require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_output.csv")
#files <- read_csv("/home/maurertm/smontgom/shared/UDN/Output/hg38/Splicing/FRASER/blood/output_5_2_24/FRASER_input.csv")
files <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_input.csv")
metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/FRASER_metadata.csv")
low_RIN <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/low_RIN.csv")

#all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13/FRASER_output.csv")
all_uncompiled <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/pipelines/FRASER_snakemake/output/FRASER_Sep_13 copy/FRASER_output.csv")
display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Age Distribution-----#################
metadata_with_age <- metadata %>% filter(!is.na(age))
metadata$ID <- paste(metadata$RDID, metadata$GSS_ID, metadata$UDN_ID, sep="_")
age_density <- metadata %>% group_by(age, affected_status) %>% tally()
colnames(age_density) <- c("age", "Affected Status", "n")

Age_plot <- ggplot(age_density, aes(y=n, x=age, group=`Affected Status`, color=`Affected Status`)) + 
  geom_smooth()+
  labs(title_lab = "Age Distribution")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  xlab("Age") +
  ylab("Number of Samples") +
  ylim(0, 15)+
  theme_classic(base_size = 27)

Age_plot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/Metadata_Analysis/age.pdf", plot=Age_plot,  limitsize = FALSE, units = "in")

################-----Sex Distribution-----#################
metadata_with_sex <- metadata %>% filter(!is.na(sex)) %>% select(affected_status, sex, ID)
male_case <- metadata_with_sex %>% filter(affected_status == "Case") %>% filter(sex == "M") %>% nrow
male_control <- metadata_with_sex %>% filter(affected_status == "Control") %>% filter(sex == "M") %>% nrow
female_case <- metadata_with_sex %>% filter(affected_status == "Case") %>% filter(sex == "F") %>% nrow
female_control <- metadata_with_sex %>% filter(affected_status == "Control") %>% filter(sex == "F") %>% nrow

sex <- data.frame(Sex=c("Male", "Male", "Female", "Female"), Affected_status=c("Case", "Control", "Case", "Control"), Values= c(male_case, male_control, female_case, female_control))
sex_plot <- ggplot(sex,aes(x = Affected_status, y = Values, fill = Sex)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title_lab="Sex Distribution")+
  theme_classic(base_size=30)+
  xlab("Sex and Affected Status") +
  ylab("Number of Individuals")

sex_plot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/Metadata_Analysis/sex.pdf", plot=sex_plot,  limitsize = FALSE, units = "in")

require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(tidyverse)

metadata <- read_csv("/home/maurertm/smontgom/shared/UDN/Analysis/FRASER/inputs/metadata/counts_metadata_filtered_with_sibs.csv")

display.brewer.all
pal <- c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"))

################-----Age Distribution-----#################
metadata$ID <- paste(metadata$RDID, metadata$GSS_ID, metadata$UDN_ID, sep="_")
age_density <- metadata %>% group_by(age, affected_status) %>% tally()
colnames(age_density) <- c("age", "Affected Status", "n")

Age_plot <- ggplot(age_density, aes(y=n, x=age, group=`Affected Status`, color=`Affected Status`)) + 
  geom_smooth()+
  labs(title_lab = "Age Distribution")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  xlab("Age") +
  ylab("Number of Samples") +
  ylim(0, 12)+
  theme_classic(base_size = 27)

Age_plot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/Metadata_Analysis/age.pdf", plot=Age_plot,  limitsize = FALSE, units = "in")

cases <- metadata %>% filter(affected_status == "Case")
mean(cases$age)
sd(cases$age)

controls <- metadata %>% filter(affected_status == "Control")
mean(controls$age)
sd(controls$age)

################-----Sex Distribution-----#################
male_case <- metadata %>% filter(affected_status == "Case") %>% filter(sex == "M") %>% nrow
male_control <- metadata %>% filter(affected_status == "Control") %>% filter(sex == "M") %>% nrow
female_case <- metadata%>% filter(affected_status == "Case") %>% filter(sex == "F") %>% nrow
female_control <- metadata %>% filter(affected_status == "Control") %>% filter(sex == "F") %>% nrow

sex <- data.frame(Sex=c("Male", "Male", "Female", "Female"), Affected_status=c("Case", "Control", "Case", "Control"), Values= c(male_case, male_control, female_case, female_control))
sex_plot <- ggplot(sex,aes(x = Affected_status, y = Values, fill = Sex)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title_lab="Sex Distribution")+
  theme_classic(base_size=30)+
  xlab("Sex and Affected Status") +
  ylab("Number of Individuals")

sex_plot
ggsave(filename="/home/maurertm/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/ASHG_Presentation/Metadata_Analysis/sex.pdf", plot=sex_plot,  limitsize = FALSE, units = "in")

################-----Affected Status Distribution-----#################
metadata %>% group_by(affected_status) %>% tally()

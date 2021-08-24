
#############################
#  load packages
#############################


library(tidyverse)
require(ggpubr)
require(ggsci)
require(readxl)

########################################################################################################
# load data
########################################################################################################

single_peptides <- read_csv("data/single_peptides_per_vaccinee.csv")

# load data
manuscript_data <- read_tsv(file.path("results/fcm/df_1.tsv"))

# order Status_2
manuscript_data$Status_2 <- factor(manuscript_data$Status_2, levels = c("Early-converter", "Late-converter", "Non-converter"))

# Create folders
dir.create("tables")

#---------
#
# clean data
#
#---------
# create metadata
metadata <- filter(manuscript_data, Time_Point == 0) %>% select(Vaccinee, Gender, Age, Status_2)

# inner join
single_peptides_per_vaccinee <- inner_join(metadata, single_peptides, by = "Vaccinee")

# paste peptides while removing NA
single_peptides_per_vaccinee$peptides <- apply(single_peptides_per_vaccinee[ , c('peptide_1', 'peptide_2', 'peptide_3', 'peptide_4', 'peptide_5', 'peptide_6')], 1, 
                                               function(x) paste(x[!is.na(x)], collapse = ", "))
# Check
table(epitope_mapping_ex_vivo=single_peptides_per_vaccinee$epitope_mapping_ex_vivo, epitope_mapping_in_vitro=single_peptides_per_vaccinee$epitope_mapping_in_vitro)

df <- 
  data.frame(table(filter(CFSE_assay, Stimulation_Condition == "Single Peptide")$Vaccinee), 
             single_peptides_per_vaccinee$peptide_number,
             table(filter(CFSE_assay, Stimulation_Condition == "Single Peptide")$Vaccinee) == single_peptides_per_vaccinee$peptide_number)
names(df) <- c("Vaccinee", "CFSE_assay", "single_peptides", "Status")
table(filter(CFSE_assay, Stimulation_Condition == "Single Peptide")$Vaccinee) == single_peptides_per_vaccinee$peptide_number

# export
single_peptides_export <- single_peptides_per_vaccinee[ , c('Vaccinee', 'Gender', 'Age', 'Status_2', 'peptide_number', 'peptides', 'sorting_of_single_peptide_specific_CD4')]
names(single_peptides_export) <- c('Vaccinee', 'Gender', 'Age', 'Status', 'peptide_number', 'single peptides', 'sorting_of_single_peptide_specific_CD4')
write_csv(single_peptides_export, "tables/single_peptides_export.csv")




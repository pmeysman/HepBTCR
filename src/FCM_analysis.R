
########################################################################################################
# load packages
########################################################################################################

library(dplyr)
library(tidyverse)
library(ggplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggsci)
library(readxl)
library(ggrepel)
library(rstatix)

###################################################################################################
# load metadata
###################################################################################################

# load
df_metadata <- read_csv("data/metadata.csv")

# Select columns
df_metadata <- df_metadata %>% dplyr::select(Vaccinee, Gender, Age)

# Vaccinee as character
df_metadata[ , c(1)] <- map(df_metadata[ , c(1)], as.character)
# Gender as factor
df_metadata[ , c(2)] <- map(df_metadata[ , c(2)], as.factor)

########################################################################################################
# load groups with serology
########################################################################################################

# load
df_Status <- read_csv("data/groups.csv")

# as factor
df_Status[ , c(4:10)] <- map(df_Status[ , c(4:10)], as.factor)
# Vaccinee as character
df_Status[ , c(1)] <- map(df_Status[ , c(1)], as.character)
# Time_Point as factor
df_Status[ , c(2)] <- map(df_Status[ , c(2)], as.factor)
# order Time_Point
df_Status$Time_Point <- factor(df_Status$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
# calculation of difference in Antibody_titre
########################################################################################################
# recode data to enable usage of mutate

temp <- df_Status

temp$Time_Point <- recode(temp$Time_Point, "0" = "Time_Point_0", "60" = "Time_Point_60", "180" = "Time_Point_180", "365" = "Time_Point_365")
# select data
FC <- dplyr::select(temp, Vaccinee, Time_Point, Antibody_titre) %>%
  # wide formate
  tidyr::spread(key = Time_Point, value = Antibody_titre) %>%
  # mutate
  mutate(Antibody_titre_FC_60_0 = log2((Time_Point_60-Time_Point_0)+1),
         Antibody_titre_FC_180_0 = log2((Time_Point_180-Time_Point_0)+1),
         Antibody_titre_FC_365_0 = log2((Time_Point_365-Time_Point_0)+1))

# bind df_Status and FC
df_Status_with_FC <- left_join(df_Status, FC[c(1, 6:8)], by = c("Vaccinee"))

# Plot Status_2
temp <- df_Status_with_FC %>%
       dplyr::filter(Time_Point == 0) %>%
       group_by(Status_2) %>%
       summarise(counts = n()) 

temp
# Status_2        counts
# 1 Early-converter     21
# 2 Late-converter       9
# 3 Non-converter        4

ggplot(temp, aes(x = Status_2, y = counts)) +
  geom_bar(
    aes(color = Status_2, fill = Status_2),
    stat = "identity", position = position_stack(), width = 0.4)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#CD534CFF"))+
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#CD534CFF"))

########################################################################################################
# load Cell count data
########################################################################################################

# load
df_counts <- read_csv("data/Cell Count.csv")

# remove columnes with "freq"
df_counts <- dplyr::select(df_counts, -contains("freq"))

# Vaccinee as character
df_counts[ , c(1)] <- map(df_counts[ , c(1)], as.character)
# Time_Point as factor
df_counts[ , c(2)] <- map(df_counts[ , c(2)], as.factor)
# order Time_Point
df_counts$Time_Point <- factor(df_counts$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
# load Antibody and Cytokine Data
########################################################################################################

# load
antibody_cytokine <- read_csv("data/Antibody and Cytokine Data.csv")

# Convert Vaccine
antibody_cytokine$Vaccinee <- as.character(antibody_cytokine$Vaccinee)

# Time_Point as factor
antibody_cytokine[ , c(2:5)] <- map(antibody_cytokine[ ,c(2:5)], as.factor)

# Renaming levels of factors
antibody_cytokine$CMV_Status <- recode(antibody_cytokine$CMV_Status, "Negative"="CMV(-)", "Positive"="CMV(+)")
antibody_cytokine$EBV_Status <- recode(antibody_cytokine$EBV_Status, "Negative"="EBV(-)", "Positive"="EBV(+)")
antibody_cytokine$HSV1_2_Status <- recode(antibody_cytokine$HSV1_2_Status, "Negative"="HSV(-)", "Positive"="HSV(+)")
antibody_cytokine$HHV6_Status <- recode(antibody_cytokine$HHV6_Status, "Negative"="HHV6(-)", "Positive"="HHV6(+)")

# Check for dimensions and names of variables
dim(antibody_cytokine)
names(antibody_cytokine)

# Check class of variables
map(antibody_cytokine, class)

# Check for sum of missing values 
sum(is.na(antibody_cytokine))

########################################################################################################
# left join df_metadata, df_counts,  df_Status_with_FC and antibody_cytokine
########################################################################################################

df_1 <- df_metadata %>% left_join(df_counts, by = c("Vaccinee")) %>% 
  left_join(df_Status_with_FC, by = c("Vaccinee", "Time_Point")) %>% 
  left_join(antibody_cytokine, by = c("Vaccinee"))

# Vaccinee as.charchter
df_1$Vaccinee <- as.character(df_1$Vaccinee)

# Gender as.factor
df_1$Gender <- as.factor(df_1$Gender)

# order Time_Point
df_1$Time_Point <- factor(df_1$Time_Point, levels = c("0", "60", "180", "365"))

# export as tsv
write_tsv(df_1, file.path("Processed data/df_1.tsv"))

########################################################################################################
# Load lymphocytes from FCM data
########################################################################################################

# load
df_lymphocytes <- read_csv("data/Lymphocytes.csv")
# rename first column
colnames(df_lymphocytes)[1] <- "Sample"
# Remove Mean and SD rows
df_lymphocytes <- df_lymphocytes %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
df_lymphocytes <- df_lymphocytes[ grep("Compensation", df_lymphocytes$Sample, invert = TRUE) , ]

# remove ".fcs" from end of Sample
# extract Vaccinee, Condition and Time_Point from Sample
df_lymphocytes <- df_lymphocytes %>% mutate(Sample = stringr::str_remove(df_lymphocytes$Sample, pattern = ".fcs"),
                                            Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                                            Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                                            Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())
# check
table(df_lymphocytes$Condition, df_lymphocytes$Time_Point)

# Vaccinee as character
df_lymphocytes[ , c(2)] <- map(df_lymphocytes[ , c(2)], as.character)
# Time_Point and Condition as factor
df_lymphocytes[ , c(3:4)] <- map(df_lymphocytes[ , c(3:4)], as.factor)

# order Time_Point
df_lymphocytes$Time_Point <- factor(df_lymphocytes$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
# Load FCM data
########################################################################################################

# Load FCM data
Table_1 <- read_csv("data/Table_1.csv")
Table_1$...18 <- NULL
colnames(Table_1)[1] <- "Sample"
Table_2 <- read_csv("data/Table_2.csv")
Table_2$...14 <- NULL
colnames(Table_2)[1] <- "Sample"
Table_3 <- read_csv("data/Table_3.csv")
Table_3$...14 <- NULL
colnames(Table_3)[1] <- "Sample"
Table_4 <- read_csv("data/Table_4.csv")
Table_4$...14 <- NULL
colnames(Table_4)[1] <- "Sample"
Table_5 <- read_csv("data/Table_5.csv")
Table_5$...14 <- NULL
colnames(Table_5)[1] <- "Sample"

# left join Table_1 to 5
df_FCM_basic <- Table_1 %>% dplyr::left_join(Table_2, by = c("Sample")) %>% 
                    dplyr::left_join(Table_3, by = c("Sample")) %>% 
                    dplyr::left_join(Table_4, by = c("Sample")) %>% 
                    dplyr::left_join(Table_5, by = c("Sample")) 
  
# replace "[ \t]\\|[ \t]Count" with ""
colnames(df_FCM_basic) <- str_replace(colnames(df_FCM_basic), pattern = "[ \t]\\|[ \t]Count", replacement = "")

# Remove Mean and SD rows
df_FCM_basic <- df_FCM_basic %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
df_FCM_basic <- df_FCM_basic[ grep("Compensation", df_FCM_basic$Sample, invert = TRUE) , ]

# remove "_CD4 T cells.fcs" from end of Sample
# remove _R at the end of Sample file
df_FCM_basic <- df_FCM_basic %>% mutate(Sample = stringr::str_remove(df_FCM_basic$Sample, pattern = "_CD4 T cells.fcs"))
df_FCM_basic <- df_FCM_basic %>% mutate(Sample = stringr::str_remove(df_FCM_basic$Sample, pattern = "_R"))

# extract Vaccinee, Condition and Time_Point from Sample and reorder
df_FCM_basic <- df_FCM_basic %>% mutate(Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                        Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                        Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())

# Check
table(df_FCM_basic$Condition, df_FCM_basic$Time_Point)
#        0 180 365 60
# HBsAg 34  34  34 34
# NC    34  34  34 34
# PC     4   3   3  2

# Vaccinee as character
df_FCM_basic[ , c(2)] <- map(df_FCM_basic[ , c(2)], as.character)

# Time_Point and Condition as factor
df_FCM_basic[ , c(3:4)] <- map(df_FCM_basic[ , c(3:4)], as.factor)

# order Time_Point
df_FCM_basic$Time_Point <- factor(df_FCM_basic$Time_Point, levels = c("0", "60", "180", "365"))

###################################################################################################
# add more variables 
###################################################################################################

df_FCM_basic <- df_FCM_basic %>% 
  
  # TH
  dplyr::mutate(TH = memory_TH + naive_TH) %>%
  # TFH
  dplyr::mutate(TFH = memory_TFH + naive_TFH) %>%
  # Treg
  dplyr::mutate(Treg = CD45RAneg_Treg + CD45RApos_Treg) %>%
  # TFR
  dplyr::mutate(TFR = CD45RAneg_TFR + CD45RApos_TFR) %>%
  
  # CD137neg_CD154pos_TH
  dplyr::mutate(CD137neg_CD154pos_TH   =   CD137neg_CD154pos_memory_TH + CD137neg_CD154pos_naive_TH) %>%
  # CD137neg_CD154pos_TFH
  dplyr::mutate(CD137neg_CD154pos_TFH  =  CD137neg_CD154pos_memory_TFH + CD137neg_CD154pos_naive_TFH) %>%
  # CD137neg_CD154pos_Treg
  dplyr::mutate(CD137neg_CD154pos_Treg = CD137neg_CD154pos_CD45RAneg_Treg + CD137neg_CD154pos_CD45RApos_Treg) %>%
  # CD137neg_CD154pos_TFR
  dplyr::mutate(CD137neg_CD154pos_TFR  =  CD137neg_CD154pos_CD45RAneg_TFR + CD137neg_CD154pos_CD45RApos_TFR) %>%
  
  # CD137pos_CD154neg_TH
  dplyr::mutate(CD137pos_CD154neg_TH =   CD137pos_CD154neg_memory_TH + CD137pos_CD154neg_naive_TH) %>%
  # CD137pos_CD154neg_TFH
  dplyr::mutate(CD137pos_CD154neg_TFH =  CD137pos_CD154neg_memory_TFH + CD137pos_CD154neg_naive_TFH) %>%
  # CD137pos_CD154neg_Treg
  dplyr::mutate(CD137pos_CD154neg_Treg = CD137pos_CD154neg_CD45RAneg_Treg + CD137pos_CD154neg_CD45RApos_Treg) %>%
  # CD137pos_CD154neg_TFR
  dplyr::mutate(CD137pos_CD154neg_TFR =  CD137pos_CD154neg_CD45RAneg_TFR + CD137pos_CD154neg_CD45RApos_TFR)

###################################################################################################
# left join df_lymphocytes and df_FCM_basic
###################################################################################################

df_FCM_basic_with_lymphocytes <- df_lymphocytes %>% left_join(df_FCM_basic[ ,-1], by = c("Vaccinee", "Time_Point", "Condition"))

########################################################################################################
# Select for numeric variables and calculate freq out of CD4
########################################################################################################

Out_CD4_freq <- subset(df_FCM_basic_with_lymphocytes, select=-c(Sample:Lymphocytes)) %>%
  # out of CD4 T cells
  dplyr::transmute_all(funs(out_CD4_freq = ./CD4*1000000)) 

# reattach
Out_CD4_freq <- bind_cols(subset(df_FCM_basic_with_lymphocytes, select=c(Sample:Condition)) , Out_CD4_freq)

###################################################################################################
# Select for numeric variables and calculate freq out of specific subsets
###################################################################################################

Out_subset_freq <- subset(df_FCM_basic_with_lymphocytes) %>% 
  
  # out of CD45RAneg
  dplyr::mutate(memory_TH_out_memory_CD4_freq = memory_TH/CD45RAneg*100) %>%
  dplyr::mutate(memory_TFH_out_memory_CD4_freq = memory_TFH/CD45RAneg*100) %>%
  dplyr::mutate(CD45RAneg_Treg_out_memory_CD4_freq = CD45RAneg_Treg/CD45RAneg*100) %>%
  dplyr::mutate(CD45RAneg_TFR_out_memory_CD4_freq = CD45RAneg_TFR/CD45RAneg*100) %>%
  
  dplyr::mutate(CD45RAneg_CD127neg_pos_CD25neg_out_memory_CD4_freq = CD45RAneg_CD127neg_pos_CD25neg/CD45RAneg*100) %>%
  dplyr::mutate(CD45RAneg_CD127neg_CD25pos_out_memory_CD4_freq = CD45RAneg_CD127neg_CD25pos/CD45RAneg*100) %>%
  
  dplyr::mutate(CD127neg_CD25neg_memory_TH_out_memory_CD4_freq = CD127neg_CD25neg_memory_TH/CD45RAneg*100) %>%
  dplyr::mutate(CD127pos_CD25neg_memory_TH_out_memory_CD4_freq = CD127pos_CD25neg_memory_TH/CD45RAneg*100) %>%
  
  # out of Treg, TFR, TH, TFH 
  dplyr::mutate(naive_TH_out_TH_freq = naive_TH/TH*100) %>%
  dplyr::mutate(memory_TH_out_TH_freq = memory_TH/TH*100) %>%
  
  dplyr::mutate(naive_TFH_out_TFH_freq = naive_TFH/TFH*100) %>%
  dplyr::mutate(memory_TFH_out_TFH_freq = memory_TFH/TFH*100) %>%
  
  dplyr::mutate(CD45RApos_Treg_out_Treg_freq = CD45RApos_Treg/Treg*100) %>%
  dplyr::mutate(CD45RAneg_Treg_out_Treg_freq = CD45RAneg_Treg/Treg*100) %>%
  
  dplyr::mutate(CD45RApos_TFR_out_TFR_freq = CD45RApos_TFR/TFR*100) %>%
  dplyr::mutate(CD45RAneg_TFR_out_TFR_freq = CD45RAneg_TFR/TFR*100) %>%
  
  # out of Treg, TFR, TH, TFH 
  dplyr::mutate(CD137pos_CD154neg_Treg_out_Treg_freq = CD137pos_CD154neg_Treg/Treg*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_Treg_out_Treg_freq = CD137neg_CD154pos_Treg/Treg*1000000) %>%
  
  dplyr::mutate(CD137pos_CD154neg_TFR_out_TFR_freq = CD137pos_CD154neg_TFR/TFR*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_TFR_out_TFR_freq = CD137neg_CD154pos_TFR/TFR*1000000) %>%
  
  dplyr::mutate(CD137pos_CD154neg_TH_out_TH_freq = CD137pos_CD154neg_TH/TH*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_TH_out_TH_freq = CD137neg_CD154pos_TH/TH*1000000) %>%  
  
  dplyr::mutate(CD137pos_CD154neg_TFH_out_TFH_freq = CD137pos_CD154neg_TFH/TFH*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_TFH_out_TFH_freq = CD137neg_CD154pos_TFH/TFH*1000000) %>%  
  
  # out of CD127neg_CD25neg_memory_TH, CD127pos_CD25neg_memory_TH, CD127neg_CD25neg_naive_TH, CD127pos_CD25neg_naive_TH
  
  dplyr::mutate(CD137neg_CD154pos_CD127neg_CD25neg_memory_TH_out_CD127neg_CD25neg_memory_TH_freq = CD137neg_CD154pos_CD127neg_CD25neg_memory_TH/CD127neg_CD25neg_memory_TH*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_CD127pos_CD25neg_memory_TH_out_CD127pos_CD25neg_memory_TH_freq = CD137neg_CD154pos_CD127pos_CD25neg_memory_TH/CD127pos_CD25neg_memory_TH*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_CD127neg_CD25neg_naive_TH_out_CD127neg_CD25neg_naive_TH_freq = CD137neg_CD154pos_CD127neg_CD25neg_naive_TH/CD127neg_CD25neg_naive_TH*1000000) %>%
  dplyr::mutate(CD137neg_CD154pos_CD127pos_CD25neg_naive_TH_out_CD127pos_CD25neg_naive_TH_freq = CD137neg_CD154pos_CD127pos_CD25neg_naive_TH/CD127pos_CD25neg_naive_TH*1000000) %>%
  
  # out CD137neg_CD154pos
  
  dplyr::mutate(TH_out_CD137neg_CD154pos_freq   = CD137neg_CD154pos_TH/CD137neg_CD154pos*100) %>%
  dplyr::mutate(TFH_out_CD137neg_CD154pos_freq  = CD137neg_CD154pos_TFH/CD137neg_CD154pos*100) %>%
  dplyr::mutate(Treg_out_CD137neg_CD154pos_freq = CD137neg_CD154pos_Treg/CD137neg_CD154pos*100) %>%
  dplyr::mutate(TFR_out_CD137neg_CD154pos_freq  = CD137neg_CD154pos_TFR/CD137neg_CD154pos*100) %>%
  
  # out CD137pos_CD154neg
  
  dplyr::mutate(TH_out_CD137pos_CD154neg_freq   = CD137pos_CD154neg_TH/CD137pos_CD154neg*100) %>%
  dplyr::mutate(TFH_out_CD137pos_CD154neg_freq  = CD137pos_CD154neg_TFH/CD137pos_CD154neg*100) %>%
  dplyr::mutate(Treg_out_CD137pos_CD154neg_freq = CD137pos_CD154neg_Treg/CD137pos_CD154neg*100) %>%
  dplyr::mutate(TFR_out_CD137pos_CD154neg_freq  = CD137pos_CD154neg_TFR/CD137pos_CD154neg*100) 

# keep only freq
Out_subset_freq <- subset(Out_subset_freq, select = -c(Lymphocytes:CD137pos_CD154neg_TFR))

########################################################################################################
# Load MFI data
########################################################################################################

Table_6_CD137 <- read_csv("data/Table (CD137 MFI).csv")
Table_6_CD137$...15 <- NULL
colnames(Table_6_CD137)[1] <- "Sample"
Table_6_CD154 <- read_csv("data/Table (CD154 MFI).csv")
colnames(Table_6_CD154)[1] <- "Sample"
Table_6_CD25 <- read_csv("data/Table (CD25 MFI).csv")
colnames(Table_6_CD25)[1] <- "Sample"
Table_6_CD127 <- read_csv("data/Table (CD127 MFI).csv")
Table_6_CD127$X15 <- NULL
colnames(Table_6_CD127)[1] <- "Sample"

# left join Table_6_CD137 to Table_6_CD154
df_MFI <- Table_6_CD137 %>% dplyr::left_join(Table_6_CD154, by = c("Sample")) %>%
                            dplyr::left_join(Table_6_CD25, by = c("Sample"))  %>%
                            dplyr::left_join(Table_6_CD127, by = c("Sample"))
  
# replace "[ \t]\\|[ \t]Count" with ""
colnames(df_MFI) <- str_replace(colnames(df_MFI), pattern = "[ \t]\\|[ \t]Median", replacement = "")

# replace "..Comp-APC-A." with "_CD154_MFI"
colnames(df_MFI) <- str_replace(colnames(df_MFI), pattern = "..Comp-APC-A.", replacement = "_CD154_MFI")

# replace "..Comp-PE-A." with "_CD137_MFI"
colnames(df_MFI) <- str_replace(colnames(df_MFI), pattern = "..Comp-PE-A.", replacement = "_CD137_MFI")

# replace "..Comp-Pacific Blue-A." with "_CD25_MFI"
colnames(df_MFI) <- str_replace(colnames(df_MFI), pattern = "..Comp-Pacific Blue-A.", replacement = "_CD25_MFI")

# replace "..Comp-BV786-A." with "_CD127_MFI"
colnames(df_MFI) <- str_replace(colnames(df_MFI), pattern = "..Comp-BV786-A.", replacement = "_CD127_MFI")

# Remove Mean and SD rows
df_MFI <- df_MFI %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
df_MFI <- df_MFI[ grep("Compensation", df_MFI$Sample, invert = TRUE) , ]

# remove "_CD4 T cells.fcs" from end of Sample
# remove _R at the end of Sample file
df_MFI <- df_MFI %>% mutate(Sample = stringr::str_remove(df_MFI$Sample, pattern = "_CD4 T cells.fcs"))
df_MFI <- df_MFI %>% mutate(Sample = stringr::str_remove(df_MFI$Sample, pattern = "_R"))

# extract Vaccinee, Condition and Time_Point from Sample and reorder
df_MFI <- df_MFI %>% mutate(Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                        Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                        Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())

# Check
table(df_MFI$Condition, df_MFI$Time_Point)
#        0 180 365 60
# HBsAg 34  34  34 34
# NC    34  34  34 34
# PC     4   3   3  2

# Vaccinee as character
df_MFI[ , c(2)] <- map(df_MFI[ , c(2)], as.character)

# Time_Point and Condition as factor
df_MFI[ , c(3:4)] <- map(df_MFI[ , c(3:4)], as.factor)

# order Time_Point
df_MFI$Time_Point <- factor(df_MFI$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
# Load Table_Treg_profile data
########################################################################################################

Table_Treg_profile <- read_csv("data/Table_Treg_profile.csv")
Table_Treg_profile$...18 <- NULL
colnames(Table_Treg_profile)[1] <- "Sample"

# Remove Mean and SD rows
Table_Treg_profile <- Table_Treg_profile %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
Table_Treg_profile <- Table_Treg_profile[ grep("Compensation", Table_Treg_profile$Sample, invert = TRUE) , ]

# remove "_CD4 T cells.fcs" from end of Sample
# remove _R at the end of Sample file
Table_Treg_profile <- Table_Treg_profile %>% mutate(Sample = stringr::str_remove(Table_Treg_profile$Sample, pattern = "_CD4 T cells.fcs"))
Table_Treg_profile <- Table_Treg_profile %>% mutate(Sample = stringr::str_remove(Table_Treg_profile$Sample, pattern = "_R"))

# extract Vaccinee, Condition and Time_Point from Sample and reorder
Table_Treg_profile <- Table_Treg_profile %>% mutate(Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                            Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                            Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())

# Check
table(Table_Treg_profile$Condition, Table_Treg_profile$Time_Point)
#        0 180 365 60
# HBsAg 34  34  34 34
# NC    34  34  34 34
# PC     4   3   3  2

# Vaccinee as character
Table_Treg_profile[ , c(2)] <- map(Table_Treg_profile[ , c(2)], as.character)

# Time_Point and Condition as factor
Table_Treg_profile[ , c(3:4)] <- map(Table_Treg_profile[ , c(3:4)], as.factor)

# order Time_Point
Table_Treg_profile$Time_Point <- factor(Table_Treg_profile$Time_Point, levels = c("0", "60", "180", "365"))

# Select for numeric variables and calculate freq out of specific subsets
Table_Treg_profile <- subset(Table_Treg_profile) %>% 
  
  # out of CD4
  dplyr::mutate(CD137neg_CD45RAneg_Treg_out_CD4_freq = CD137neg_CD45RAneg_Treg/CD4*1000000) %>%
  dplyr::mutate(CD137pos_CD45RAneg_Treg_out_CD4_freq = CD137pos_CD45RAneg_Treg/CD4*1000000) %>%
  
  dplyr::mutate(CD137neg_CD45RApos_Treg_out_CD4_freq = CD137neg_CD45RApos_Treg/CD4*1000000) %>%
  dplyr::mutate(CD137pos_CD45RApos_Treg_out_CD4_freq = CD137pos_CD45RApos_Treg/CD4*1000000) %>%
  
  # out of Treg
  dplyr::mutate(CD137neg_CD45RAneg_Treg_out_Treg_freq = CD137neg_CD45RAneg_Treg/Treg*100) %>%
  dplyr::mutate(CD137pos_CD45RAneg_Treg_out_Treg_freq = CD137pos_CD45RAneg_Treg/Treg*100) %>%
  
  dplyr::mutate(CD137neg_CD45RApos_Treg_out_Treg_freq = CD137neg_CD45RApos_Treg/Treg*100) %>%
  dplyr::mutate(CD137pos_CD45RApos_Treg_out_Treg_freq = CD137pos_CD45RApos_Treg/Treg*100) %>%
  
  # out of CD45RAneg_Treg
  dplyr::mutate(CD137neg_CD45RAneg_Treg_out_CD45RAneg_Treg_freq = CD137neg_CD45RAneg_Treg/CD45RAneg_Treg*100) %>%
  dplyr::mutate(CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq = CD137pos_CD45RAneg_Treg/CD45RAneg_Treg*100) %>%
  
  # out of CD45RApos_Treg
  dplyr::mutate(CD137neg_CD45RApos_Treg_out_CD45RApos_Treg_freq = CD137neg_CD45RApos_Treg/CD45RApos_Treg*100) %>%
  dplyr::mutate(CD137pos_CD45RApos_Treg_out_CD45RApos_Treg_freq = CD137pos_CD45RApos_Treg/CD45RApos_Treg*100) %>%
  
  
  subset(select = -c(CD4:CD137pos_CD45RAneg_Treg))

########################################################################################################
## calculation of significance of CD4 T cell responses
########################################################################################################

###########################
###########################

## calculate percentage and subtract percentage in the negative control, return only if higher than 0
## return the response_freq as.numeric and out of 1000000 CD4 T cells
response_freq <- function(mtx) {
  freq <- (mtx[2,1]/mtx[2,2] - mtx[1,1]/mtx[1,2])*1000000
  if (freq <= 0){
    return(as.numeric(0))} else {as.numeric(freq)}}

###########################
###########################

# create vector of Vaccinees to loop over
Vaccinees <- unique(df_FCM_basic_with_lymphocytes$Vaccinee)

# create vector of target_pop to loop over
target_pops <- c(names(dplyr::select(df_FCM_basic_with_lymphocytes, contains("CD137pos_CD154neg"))),
                 names(dplyr::select(df_FCM_basic_with_lymphocytes, contains("CD137neg_CD154pos"))),
                 names(dplyr::select(df_FCM_basic_with_lymphocytes, contains("CD137pos_CD154pos"))))

# Create an empty list (higher level list)
data_h = list()

# loop through Vaccinees to fill data
for (i in Vaccinees){
  
  print(paste(c("Vaccinee ", i), collapse = ""))
  flush.console()
  # Select Vaccine data
  Vaccinee_data <- dplyr::filter(df_FCM_basic_with_lymphocytes, Vaccinee == i)
  
  # Create an empty list (lower level list)
  data_l = list()
  
  # loop through target_pops to fill data
  for (target_pop in target_pops){

    ### retrieve counts of CD4 and target_pop at Time_Point 0
    # of control condition
    control_0     <- Vaccinee_data %>% dplyr::filter(Time_Point == "0"  & Condition == "NC") %>% dplyr::select(target_pop, CD4)
    # of stimulation condition
    stimulation_0     <- Vaccinee_data %>% dplyr::filter(Time_Point == "0"  & Condition == "HBsAg") %>% dplyr::select(target_pop, CD4)
    # generate matrix
    matrix_0     <- bind_rows(control_0, stimulation_0)
    # calculate freq before subtracting target_pop from CD4
    freq_0     <- response_freq(matrix_0)
    # subtract target_pop from CD4
    matrix_0[1,2] <- matrix_0[1,2]-matrix_0[1,1]
    matrix_0[2,2] <- matrix_0[2,2]-matrix_0[2,1]
    # calculate fisher exaxt t test and retrieve p value after subtracting target_pop from CD4
    fisher_test_0     <- fisher.test(matrix_0, alternative = "less")$p.value
    
    ### retrieve counts of CD4 and target_pop at Time_Point 60
    # of control condition
    control_60     <- Vaccinee_data %>% dplyr::filter(Time_Point == "60"  & Condition == "NC") %>% dplyr::select(target_pop, CD4)
    # of stimulation condition
    stimulation_60     <- Vaccinee_data %>% dplyr::filter(Time_Point == "60"  & Condition == "HBsAg") %>% dplyr::select(target_pop, CD4)
    # generate matrix
    matrix_60     <- bind_rows(control_60, stimulation_60)
    # calculate freq before subtracting target_pop from CD4
    freq_60     <- response_freq(matrix_60)
    # subtract target_pop from CD4
    matrix_60[1,2] <- matrix_60[1,2]-matrix_60[1,1]
    matrix_60[2,2] <- matrix_60[2,2]-matrix_60[2,1]
    # calculate fisher exaxt t test and retrieve p value after subtracting target_pop from CD4
    fisher_test_60     <- fisher.test(matrix_60, alternative = "less")$p.value
    
    ### retrieve counts of CD4 and target_pop at Time_Point 180
    # of control condition
    control_180     <- Vaccinee_data %>% dplyr::filter(Time_Point == "180"  & Condition == "NC") %>% dplyr::select(target_pop, CD4)
    # of stimulation condition
    stimulation_180     <- Vaccinee_data %>% dplyr::filter(Time_Point == "180"  & Condition == "HBsAg") %>% dplyr::select(target_pop, CD4)
    # generate matrix
    matrix_180     <- bind_rows(control_180, stimulation_180)
    # calculate freq before subtracting target_pop from CD4
    freq_180     <- response_freq(matrix_180)
    # subtract target_pop from CD4
    matrix_180[1,2] <- matrix_180[1,2]-matrix_180[1,1]
    matrix_180[2,2] <- matrix_180[2,2]-matrix_180[2,1]
    # calculate fisher exaxt t test and retrieve p value after subtracting target_pop from CD4
    fisher_test_180     <- fisher.test(matrix_180, alternative = "less")$p.value
    
    
    ### retrieve counts of CD4 and target_pop at Time_Point 365
    # of control condition
    control_365     <- Vaccinee_data %>% dplyr::filter(Time_Point == "365"  & Condition == "NC") %>% dplyr::select(target_pop, CD4)
    # of stimulation condition
    stimulation_365     <- Vaccinee_data %>% dplyr::filter(Time_Point == "365"  & Condition == "HBsAg") %>% dplyr::select(target_pop, CD4)
    # generate matrix
    matrix_365     <- bind_rows(control_365, stimulation_365)
    # calculate freq before subtracting target_pop from CD4
    freq_365     <- response_freq(matrix_365)
    # subtract target_pop from CD4
    matrix_365[1,2] <- matrix_365[1,2]-matrix_365[1,1]
    matrix_365[2,2] <- matrix_365[2,2]-matrix_365[2,1]
    # calculate fisher exaxt t test and retrieve p value after subtracting target_pop from CD4
    fisher_test_365     <- fisher.test(matrix_365, alternative = "less")$p.value
    
    # combine results and add to data_l using tibble
    data_l[[target_pop]] <- tibble(
      Time_Point = c("0", "60", "180", "365"),
      response_freq       = c(freq_0, freq_60, freq_180, freq_365),
      response_p_value    = c(fisher_test_0, fisher_test_60, fisher_test_180, fisher_test_365)
    )
  }
  
  # combine the list outside the loop to a tbl
  # add it to data_h
  data_h[[i]] <- bind_rows(data_l, .id = "target_pop")
}

# combine the list outside the loop to a tbl
FCM_Analysis_NC_Subtracted <- bind_rows(data_h, .id = "Vaccinee")

########################################################################################################
# Check for a seq of thresholds by a loop
########################################################################################################

# create an empty list
freq_positive <- list()

# set loop value to zero
loop <- 0

for (n in c(5 %o% 10^(-(2:10)))) {
  # zero the value of response_freq if response_p_value is larger than the threshold
  FCM_Analysis_NC_Subtracted_threshold <- FCM_Analysis_NC_Subtracted %>% mutate(response_freq = ifelse(response_p_value >= n, 0, response_freq))
  
  # what is the precentage of of response_freq more than zero
  precentage <- count(dplyr::filter(FCM_Analysis_NC_Subtracted_threshold, response_freq > 0))/count(dplyr::filter(FCM_Analysis_NC_Subtracted_threshold))*100
  
  loop <- loop +1
  freq_positive[[loop]] <- tibble(
    threshold = c(n),
    freq       = c(as.numeric(precentage)))
}

########################################################################################################
########################################################################################################

# set a threshold
threshold <- 0.05

# zero the value of response_freq if response_p_value is larger than the threshold
FCM_Analysis_NC_Subtracted_threshold <- FCM_Analysis_NC_Subtracted %>% mutate(response_freq = ifelse(response_p_value >= threshold, 0, response_freq))

# what is the precentage of of response_freq more than zero
count(dplyr::filter(FCM_Analysis_NC_Subtracted_threshold, response_freq > 0))/count(dplyr::filter(FCM_Analysis_NC_Subtracted_threshold))*100

# plot a geom_histogram of response_freq faceted by Time_Point
ggplot(data=dplyr::filter(FCM_Analysis_NC_Subtracted_threshold, response_freq > 0), aes(response_freq)) + 
  geom_histogram(breaks=seq(0, 2000, by = 10), 
                 col="red", 
                 fill="green", 
                 alpha = .2)+
  facet_grid(Time_Point~.)

# select for certain columns
FCM_Analysis_NC_Subtracted_threshold_wide <- FCM_Analysis_NC_Subtracted_threshold %>%dplyr::select(Vaccinee, Time_Point, target_pop, response_freq) %>% 
  # change to a wide format
  tidyr::spread(key = target_pop, value = response_freq)

# Time_Point as factor
FCM_Analysis_NC_Subtracted_threshold_wide[ , c(2)] <- map(FCM_Analysis_NC_Subtracted_threshold_wide[ , c(2)], as.factor)
# order Time_Point
FCM_Analysis_NC_Subtracted_threshold_wide$Time_Point <- factor(FCM_Analysis_NC_Subtracted_threshold_wide$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
# Load FCM data
########################################################################################################

Subset_specific <- read_csv("data/Subset_specific_sign.csv")
colnames(Subset_specific)[1] <- "Sample"

# replace "[ \t]\\|[ \t]Count" with ""
colnames(Subset_specific) <- str_replace(colnames(Subset_specific), pattern = "[ \t]\\|[ \t]Count", replacement = "")

# Remove Mean and SD rows
Subset_specific <- Subset_specific %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
Subset_specific <- Subset_specific[ grep("Compensation", Subset_specific$Sample, invert = TRUE) , ]

# remove "_CD4 T cells.fcs" from end of Sample
# remove _R at the end of Sample file
Subset_specific <- Subset_specific %>% mutate(Sample = stringr::str_remove(Subset_specific$Sample, pattern = "_CD4 T cells.fcs"))
Subset_specific <- Subset_specific %>% mutate(Sample = stringr::str_remove(Subset_specific$Sample, pattern = "_R"))

# extract Vaccinee, Condition and Time_Point from Sample and reorder
Subset_specific <- Subset_specific %>% mutate(Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                                              Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                                              Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())

# Check
table(Subset_specific$Condition, Subset_specific$Time_Point)
#        0 180 365 60
# HBsAg 34  34  34 34
# NC    34  34  34 34
# PC     4   3   3  2

# Vaccinee as character
Subset_specific[ , c(2)] <- map(Subset_specific[ , c(2)], as.character)

# Time_Point and Condition as factor
Subset_specific[ , c(3:4)] <- map(Subset_specific[ , c(3:4)], as.factor)

# order Time_Point
Subset_specific$Time_Point <- factor(Subset_specific$Time_Point, levels = c("0", "60", "180", "365"))

########################################################################################################
########################################################################################################

Subset_specific$naive_CD4 <- Subset_specific$naive_TH + Subset_specific$naive_TFH + Subset_specific$CD45RApos_Treg + Subset_specific$CD45RApos_TFR
Subset_specific$CD137pos_CD154neg_naive_CD4 <- Subset_specific$CD137pos_CD154neg_naive_TH + Subset_specific$CD137pos_CD154neg_naive_TFH + Subset_specific$CD137pos_CD154neg_CD45RApos_Treg + Subset_specific$CD137pos_CD154neg_CD45RApos_TFR
Subset_specific$CD137neg_CD154pos_naive_CD4 <- Subset_specific$CD137neg_CD154pos_naive_TH + Subset_specific$CD137neg_CD154pos_naive_TFH + Subset_specific$CD137neg_CD154pos_CD45RApos_Treg + Subset_specific$CD137neg_CD154pos_CD45RApos_TFR
Subset_specific$CD137pos_CD154pos_naive_CD4 <- Subset_specific$CD137pos_CD154pos_naive_TH + Subset_specific$CD137pos_CD154pos_naive_TFH + Subset_specific$CD137pos_CD154pos_CD45RApos_Treg + Subset_specific$CD137pos_CD154pos_CD45RApos_TFR

Subset_specific$memory_CD4 <- Subset_specific$memory_TH + Subset_specific$memory_TFH + Subset_specific$CD45RAneg_Treg + Subset_specific$CD45RAneg_TFR
Subset_specific$CD137pos_CD154neg_memory_CD4 <- Subset_specific$CD137pos_CD154neg_memory_TH + Subset_specific$CD137pos_CD154neg_memory_TFH + Subset_specific$CD137pos_CD154neg_CD45RAneg_Treg + Subset_specific$CD137pos_CD154neg_CD45RAneg_TFR
Subset_specific$CD137neg_CD154pos_memory_CD4 <- Subset_specific$CD137neg_CD154pos_memory_TH + Subset_specific$CD137neg_CD154pos_memory_TFH + Subset_specific$CD137neg_CD154pos_CD45RAneg_Treg + Subset_specific$CD137neg_CD154pos_CD45RAneg_TFR
Subset_specific$CD137pos_CD154pos_memory_CD4 <- Subset_specific$CD137pos_CD154pos_memory_TH + Subset_specific$CD137pos_CD154pos_memory_TFH + Subset_specific$CD137pos_CD154pos_CD45RAneg_Treg + Subset_specific$CD137pos_CD154pos_CD45RAneg_TFR

Subset_specific$CD4 <- Subset_specific$naive_CD4 + Subset_specific$memory_CD4
Subset_specific$CD137pos_CD154neg_CD4 <- Subset_specific$CD137pos_CD154neg_naive_CD4 + Subset_specific$CD137pos_CD154neg_memory_CD4
Subset_specific$CD137neg_CD154pos_CD4 <- Subset_specific$CD137neg_CD154pos_naive_CD4 + Subset_specific$CD137neg_CD154pos_memory_CD4
Subset_specific$CD137pos_CD154pos_CD4 <- Subset_specific$CD137pos_CD154pos_naive_CD4 + Subset_specific$CD137pos_CD154pos_memory_CD4

# check dim
dim(Subset_specific) # 284  48

########################################################################################################
## calculation of significance of CD4 T cell responses
########################################################################################################

# create vector of Vaccinees to loop over
Vaccinees <- unique(Subset_specific$Vaccinee)

# create vector of mother_pops to loop over
mother_pops <- c("naive_TH", "memory_TH", "naive_TFH", "memory_TFH", "CD45RApos_Treg", "CD45RAneg_Treg", "CD45RApos_TFR", "CD45RAneg_TFR", "naive_CD4", "memory_CD4", "CD4")

# create vector of target_pop to loop over
target_pops <- c(names(dplyr::select(Subset_specific, contains("CD137pos_CD154neg"))),
                 names(dplyr::select(Subset_specific, contains("CD137neg_CD154pos"))),
                 names(dplyr::select(Subset_specific, contains("CD137pos_CD154pos"))))

## calculate percentage and subtract percentage in the negative control, return only if higher than 0
## return the response_freq as.numeric and out of 1000000 cells of the mother_pop
response_freq <- function(mtx) {
  freq <- (mtx[2,1]/mtx[2,2] - mtx[1,1]/mtx[1,2])*1000000
  if (freq <= 0){
    return(as.numeric(0))} else {as.numeric(freq)}}

# Create an empty list (higher level list)
data_h = list()

# loop through Vaccinees to fill data
for (Vaccinee_i in Vaccinees){
  
  print(paste(c("Vaccinee ", Vaccinee_i), collapse = ""))
  flush.console()
  # Select Vaccine data
  Vaccinee_data <- dplyr::filter(Subset_specific, Vaccinee == Vaccinee_i)
  
  # Create an empty list (lower level list)
  data_l = list()
  
  # loop through mother_pops to fill data
  for (mother_pop in mother_pops){
    
    print(paste(c("Mother Population ", mother_pop), collapse = ""))
    flush.console()
    
    # Create an empty list (lower level list)
    data_w = list()
    
    # create vector of daughter_pops to loop over
    daughter_pops <- names(dplyr::select(Subset_specific, contains(mother_pop)) %>% dplyr::select(contains("CD137")))
    
    # loop through daughter_pops to fill data
    for (daughter_pop in daughter_pops){
      
      #print(paste(c("Daughter Population ", daughter_pop), collapse = ""))
      #flush.console()
      
      ###########################
      ###########################
      ### retrieve counts of mother_pop and daughter_pop at Time_Point 0
      # of control condition
      control_0     <- Vaccinee_data %>% dplyr::filter(Time_Point == "0"  & Condition == "NC") %>% dplyr::select(daughter_pop, mother_pop)
      # of stimulation condition
      stimulation_0     <- Vaccinee_data %>% dplyr::filter(Time_Point == "0"  & Condition == "HBsAg") %>% dplyr::select(daughter_pop, mother_pop)
      # generate matrix
      matrix_0     <- bind_rows(control_0, stimulation_0)
      # calculate freq before subtracting daughter_pop from mother_pop
      freq_0     <- response_freq(matrix_0)
      # subtract daughter_pop from mother_pop
      matrix_0[1,2] <- matrix_0[1,2]-matrix_0[1,1]
      matrix_0[2,2] <- matrix_0[2,2]-matrix_0[2,1]
      # calculate fisher exaxt t test and retrieve p value after subtracting daughter_pop from mother_pop
      fisher_test_0     <- fisher.test(matrix_0, alternative = "less")$p.value
      
      ### retrieve counts of mother_pop and daughter_pop at Time_Point 60
      # of control condition
      control_60     <- Vaccinee_data %>% dplyr::filter(Time_Point == "60"  & Condition == "NC") %>% dplyr::select(daughter_pop, mother_pop)
      # of stimulation condition
      stimulation_60     <- Vaccinee_data %>% dplyr::filter(Time_Point == "60"  & Condition == "HBsAg") %>% dplyr::select(daughter_pop, mother_pop)
      # generate matrix
      matrix_60     <- bind_rows(control_60, stimulation_60)
      # calculate freq before subtracting daughter_pop from mother_pop
      freq_60     <- response_freq(matrix_60)
      # subtract daughter_pop from mother_pop
      matrix_60[1,2] <- matrix_60[1,2]-matrix_60[1,1]
      matrix_60[2,2] <- matrix_60[2,2]-matrix_60[2,1]
      # calculate fisher exaxt t test and retrieve p value after subtracting daughter_pop from mother_pop
      fisher_test_60     <- fisher.test(matrix_60, alternative = "less")$p.value
      
      ### retrieve counts of mother_pop and daughter_pop at Time_Point 180
      # of control condition
      control_180     <- Vaccinee_data %>% dplyr::filter(Time_Point == "180"  & Condition == "NC") %>% dplyr::select(daughter_pop, mother_pop)
      # of stimulation condition
      stimulation_180     <- Vaccinee_data %>% dplyr::filter(Time_Point == "180"  & Condition == "HBsAg") %>% dplyr::select(daughter_pop, mother_pop)
      # generate matrix
      matrix_180     <- bind_rows(control_180, stimulation_180)
      # calculate freq before subtracting daughter_pop from mother_pop
      freq_180     <- response_freq(matrix_180)
      # subtract daughter_pop from mother_pop
      matrix_180[1,2] <- matrix_180[1,2]-matrix_180[1,1]
      matrix_180[2,2] <- matrix_180[2,2]-matrix_180[2,1]
      # calculate fisher exaxt t test and retrieve p value after subtracting daughter_pop from mother_pop
      fisher_test_180     <- fisher.test(matrix_180, alternative = "less")$p.value
      
      
      ### retrieve counts of mother_pop and daughter_pop at Time_Point 365
      # of control condition
      control_365     <- Vaccinee_data %>% dplyr::filter(Time_Point == "365"  & Condition == "NC") %>% dplyr::select(daughter_pop, mother_pop)
      # of stimulation condition
      stimulation_365     <- Vaccinee_data %>% dplyr::filter(Time_Point == "365"  & Condition == "HBsAg") %>% dplyr::select(daughter_pop, mother_pop)
      # generate matrix
      matrix_365     <- bind_rows(control_365, stimulation_365)
      # calculate freq before subtracting daughter_pop from mother_pop
      freq_365     <- response_freq(matrix_365)
      # subtract daughter_pop from mother_pop
      matrix_365[1,2] <- matrix_365[1,2]-matrix_365[1,1]
      matrix_365[2,2] <- matrix_365[2,2]-matrix_365[2,1]
      # calculate fisher exaxt t test and retrieve p value after subtracting daughter_pop from mother_pop
      fisher_test_365     <- fisher.test(matrix_365, alternative = "less")$p.value
      
      
      # combine results and add to data_w using tibble
      data_w[[daughter_pop]] <- tibble(
        Time_Point = c("0", "60", "180", "365"),
        response_freq       = c(freq_0, freq_60, freq_180, freq_365),
        response_p_value    = c(fisher_test_0, fisher_test_60, fisher_test_180, fisher_test_365))
    }
    
    # combine the list outside the loop to a tbl
    # add it to data_l
    data_l[[mother_pop]] <- bind_rows(data_w, .id = "Daughter_Population")
  }
  
  # combine the list outside the loop to a tbl
  # add it to data_h
  data_h[[Vaccinee_i]] <- bind_rows(data_l, .id = "Mother_Population")
  
}

# combine the list outside the loop to a tbl
Subset_specific_NC_Subtracted <- bind_rows(data_h, .id = "Vaccinee")

########################################################################################################
# Check for a seq of thresholds by a loop
########################################################################################################

# create an empty list
freq_positive <- list()

# set loop value to zero
loop <- 0

for (n in c(1 %o% 10^(-(2:10)))) {
  # zero the value of response_freq if response_p_value is larger than the threshold
  Subset_specific_NC_Subtracted_threshold <- Subset_specific_NC_Subtracted %>% mutate(response_freq = ifelse(response_p_value >= n, 0, response_freq))
  
  # what is the precentage of of response_freq more than zero
  precentage <- count(dplyr::filter(Subset_specific_NC_Subtracted_threshold, response_freq > 0))/count(dplyr::filter(Subset_specific_NC_Subtracted_threshold))*100
  
  loop <- loop +1
  freq_positive[[loop]] <- tibble(
    threshold = c(n),
    freq       = c(as.numeric(precentage)))
}

########################################################################################################
#
########################################################################################################

# set a threshold
threshold <- 0.05

# zero the value of response_freq if response_p_value is larger than the threshold
Subset_specific_NC_Subtracted_threshold <- Subset_specific_NC_Subtracted %>% mutate(response_freq = ifelse(response_p_value >= threshold, 0, response_freq))

# what is the precentage of of response_freq more than zero
count(dplyr::filter(Subset_specific_NC_Subtracted_threshold, response_freq > 0))/count(dplyr::filter(Subset_specific_NC_Subtracted_threshold))*100

# plot a geom_histogram of response_freq more than 0 faceted by Time_Point
ggplot(data=dplyr::filter(Subset_specific_NC_Subtracted_threshold, response_freq > 0), aes(response_freq)) + 
  geom_histogram(breaks=seq(0, 2000, by = 10), 
                 col="red", 
                 fill="green", 
                 alpha = .2)+
  facet_grid(Time_Point~.)

########################################################################################################
#
########################################################################################################

# unite Daughter Population and Mother Population
Subset_specific_NC_Subtracted_threshold_united <- Subset_specific_NC_Subtracted_threshold %>% unite("Population", c("Daughter_Population","Mother_Population"), sep = "_out_")

# order Time_Point
Subset_specific_NC_Subtracted_threshold_united$Time_Point <- factor(Subset_specific_NC_Subtracted_threshold_united$Time_Point, levels = c("0", "60", "180", "365"))

# group_by Mother_Population and Time_Point and calculate the average of response_freq
# Subset_specific_NC_Subtracted_threshold_united %>% group_by(Population, Time_Point) %>% summarise(n = n(), response_freq_avg = mean(response_freq)) %>% View()

# select for certain columns
Subset_specific_NC_Subtracted_threshold_wide <- Subset_specific_NC_Subtracted_threshold_united %>%dplyr::select(Vaccinee, Time_Point, Population, response_freq) %>% 
  # change to a wide format
  tidyr::spread(key = Population, value = response_freq)


########################################################################################################
#
########################################################################################################

## Clean Out_CD4_freq
# Remove columns related to CD137 and CD154 signature
Out_CD4_freq_temp <- Out_CD4_freq %>% subset(select=-c(CD137neg_CD154pos_out_CD4_freq:CD137pos_CD154neg_out_CD4_freq, CD137pos_CD154neg_memory_TH_out_CD4_freq:CD137pos_CD154pos_CD45RApos_TFR_out_CD4_freq, CD137neg_CD154pos_TH_out_CD4_freq:CD137pos_CD154neg_TFR_out_CD4_freq))
# Keep only NC data
Out_CD4_freq_temp <- filter(Out_CD4_freq_temp, Condition == "NC")
# Remove Column
Out_CD4_freq_temp$Condition <- NULL

########################################################################################################
#
########################################################################################################

## Clean Out_subset_freq
# Remove columns related to CD137 and CD154 signature
Out_subset_freq_temp <- Out_subset_freq %>% subset(select=-c(CD137pos_CD154neg_Treg_out_Treg_freq:TFR_out_CD137pos_CD154neg_freq))
# Keep only NC data
Out_subset_freq_temp <- filter(Out_subset_freq_temp, Condition == "NC")
# Remove Column
Out_subset_freq_temp$Condition <- NULL

########################################################################################################
#
########################################################################################################

## Clean df_MFI
# Keep only NC data
df_MFI_temp <- filter(df_MFI, Condition == "NC")
# Remove Column
df_MFI_temp$Condition <- NULL

########################################################################################################
#
########################################################################################################

## Clean Table_Treg_profile
# Keep only NC data
Table_Treg_profile_temp <- filter(Table_Treg_profile, Condition == "NC")
# Remove Column
Table_Treg_profile_temp$Condition <- NULL

########################################################################################################
#
########################################################################################################

# left join df_1, Out_CD4_freq, Out_subset_freq, df_MFI, Table_Treg_profile, Subset_specific_end
df_all <- df_1 %>% 
  left_join(Out_CD4_freq_temp, by = c("Vaccinee", "Time_Point"))                                 %>% 
  left_join(Out_subset_freq_temp, by = c("Sample", "Vaccinee", "Time_Point"))                    %>%
  left_join(df_MFI_temp, by = c("Sample", "Vaccinee", "Time_Point"))                             %>% 
  left_join(Table_Treg_profile_temp, by = c("Sample", "Vaccinee", "Time_Point"))       %>%
  left_join(FCM_Analysis_NC_Subtracted_threshold_wide, by = c("Vaccinee", "Time_Point"))         %>%
  left_join(Subset_specific_NC_Subtracted_threshold_wide, by = c("Vaccinee", "Time_Point"))

###################################################################################################
#
###################################################################################################

# order Time_Point
df_all$Time_Point <- factor(df_all$Time_Point, levels = c("0", "60", "180", "365"))

# order Set
df_all$Set <- factor(df_all$Set, levels = c("1", "2", "3", "4", "5", "6", "7"))

# order Status
df_all$Status_1 <- factor(df_all$Status_1, levels = c("Non-converter", "Early-converter", "Mid-converter", "Late-converter"))

# order Status_2
df_all$Status_2 <- factor(df_all$Status_2, levels = c("Non-converter", "Early-converter", "Late-converter"))

# order Status_3
df_all$Status_3 <- factor(df_all$Status_3, levels = c("Non-converter", "Low early-converter", "High early-converter", "Late-converter"))

# order Status_4
df_all$Status_4 <- factor(df_all$Status_4, levels = c("Non-converter", "Low early-converter", "High early-converter", "Mid-converter", "Late-converter"))

# order variables 
df_all <- df_all %>% dplyr::select(Set, Sample, Vaccinee, Gender, Age, Time_Point, 
                                         Status_1, Status_2, Status_3, Status_4,
                                         Antibody_titre, Antibody_titre_FC_60_0, Antibody_titre_FC_180_0, Antibody_titre_FC_365_0, 
                                         Status_60_0, Status_180_0, Status_365_0, 
                                         CMV_Status, EBV_Status, HSV1_2_Status, HHV6_Status,
                                         WBC, RBC, HGB, HCT, PLT, LYM_absolute, MON_absolute, GRA_absolute, everything())

# Check for dimensions and names of variables
dim(df_all)
names(df_all)

# Check class of variables
map(df_all, class)

# Check for sum of missing values 
sum(is.na(df_all))

# Check for missing values 
#sapply(df_all, function(x) sum(is.na(x))) %>% View()

# export as tsv
write_tsv(df_all, file.path("results/fcm/df_all.tsv"))

###################################################################################################
# load memory CD4 T cell repertoires
###################################################################################################

# load
all_diversity <- read_tsv(file = "data/all_diversity.tsv")
# Check
table(all_diversity$Vaccinee, all_diversity$Time_Point)
# check nrow
nrow(all_diversity)/2 # 35

# Both 8 and 23 are duplicated
# Keep unique rows
all_diversity <- all_diversity[!duplicated(all_diversity[1:2]), ]
# Check
table(all_diversity$Vaccinee, all_diversity$Time_Point)
# check nrow
nrow(all_diversity)/2 # 33
# Vaccinee as.character
all_diversity[ , c(1)] <- map(all_diversity[ , c(1)], as.character)
# Time_Point as.factor
all_diversity[ , c(2)] <- map(all_diversity[ , c(2)], as.factor)

# Remove last two columns
all_diversity[ , c(5,6)] <- NULL

# left join df_all, and all_diversity
df_all_with_diversity <- filter(df_all, Time_Point == "0" | Time_Point == "60") %>% 
  inner_join(all_diversity, by = c("Vaccinee", "Time_Point")) 

# Check
table(df_all_with_diversity$Vaccinee, df_all_with_diversity$Time_Point)

# Time_Point as.factor
df_all_with_diversity[ , c(4)] <- map(df_all_with_diversity[ , c(4)], as.factor)

###################################################################################################
# end
###################################################################################################



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

#############################
# create directory to save figures
#############################

fig_dir <- 'figures'
dir.create(fig_dir)

supp_fig <- 'figures/supplementary figures'
dir.create(supp_fig)


########################################################################################################
# load data
########################################################################################################

manuscript_data <- read_tsv(file.path("results/fcm/df_1.tsv"))
# order Status_2
manuscript_data$Status_2 <- factor(manuscript_data$Status_2, levels = c("Early-converter", "Late-converter", "Non-converter"))

########################################################################################################
# load, clean and prep data 
########################################################################################################

# Load data
CFSE_assay_1 <- read_csv("Raw data/Sorting of expanded CD4 T cells (after ex vivo epitope mapping).csv" , col_select = c(1:3))
CFSE_assay_2 <- read_csv("Raw data/Sorting of expanded CD4 T cells (after in vitro epitope mapping).csv", col_select = c(1:3))

# change colnames
colnames(CFSE_assay_1) <- c("sample", "CD4_count", "CFSElow_CD4_count")
colnames(CFSE_assay_2) <- c("sample", "CD4_count", "CFSElow_CD4_count")

## Expansion of vaccine-specific CD4 T cells
# clean dataset
CFSE_assay_1 <- CFSE_assay_1 %>% dplyr::filter(sample != "Mean", sample != "SD") %>% # Remove Mean and SD
  dplyr::filter(!grepl("Compensation Controls",sample)) # Remove Compensation Controls

CFSE_assay_2 <- CFSE_assay_2 %>% filter(sample != "Mean", sample != "SD") %>% # Remove Mean and SD
  filter(!grepl("Compensation Controls",sample)) # Remove Compensation Controls


# extract Vaccinee, Antigen from sample
CFSE_assay_1 <- CFSE_assay_1 %>% mutate(Vaccinee = as.numeric(stringr::str_extract(CFSE_assay_1$sample, pattern = "\\d\\d")),
                                        Antigen = stringr::str_extract(CFSE_assay_1$sample, pattern = "[V|P|N][[:alnum:]]+"))
CFSE_assay_2 <- CFSE_assay_2 %>% mutate(Vaccinee = as.numeric(stringr::str_extract(CFSE_assay_2$sample, pattern = "\\d\\d")),
                                        Antigen = stringr::str_extract(CFSE_assay_2$sample, pattern = "[V|P|N][[:alnum:]]+"))
# Check
table(unique(CFSE_assay_1$Vaccinee) == unique(manuscript_data$Vaccinee)) # 34

# what is common between CFSE_assay_1 and CFSE_assay_2
intersect(unique(CFSE_assay_1$Vaccinee), unique(CFSE_assay_2$Vaccinee))
# what is different between CFSE_assay_1 and CFSE_assay_2
vaccinees_vector <- setdiff(unique(CFSE_assay_1$Vaccinee), unique(CFSE_assay_2$Vaccinee))

# keep only CFSE_assay_1 that is not in CFSE_assay_2
nrow(CFSE_assay_1) # 151
CFSE_assay_1 <- filter(CFSE_assay_1, Vaccinee %in% vaccinees_vector)
nrow(CFSE_assay_1) # 124

CFSE_assay <- bind_rows(list(ex_vivo = CFSE_assay_1, in_vitro = CFSE_assay_2), .id = "epitope_mapping")
nrow(CFSE_assay)
table(CFSE_assay$Vaccinee, CFSE_assay$Antigen)

# mutate
CFSE_assay <- CFSE_assay %>% mutate(CFSElow_CD4_freq = 100*CFSElow_CD4_count/CD4_count)

manuscript_data_day_60 <- filter(manuscript_data, Time_Point == 60)

# inner join
CFSE_assay <- inner_join(CFSE_assay, manuscript_data_day_60, by = "Vaccinee")


CFSE_assay$Stimulation_Condition <- CFSE_assay$Antigen
CFSE_assay$Stimulation_Condition <- if_else(CFSE_assay$Antigen != "NC" & CFSE_assay$Antigen != "PP", "SP", CFSE_assay$Antigen)

# recode 
CFSE_assay$Stimulation_Condition <- recode(CFSE_assay$Stimulation_Condition, 
                                           NC = "Negative Control", PP = "Peptide Pool", SP = "Single Peptide")

# arrange by Vaccinee
CFSE_assay <- CFSE_assay %>% dplyr::arrange(Vaccinee)
CFSE_assay$Vaccinee <- as.factor(CFSE_assay$Vaccinee)

########################################################################################################
# Fig_S2a: Antibody_titre
########################################################################################################

my_comparisons <- list(c("0", "60"),
                       c("60", "180"),
                       c("180", "365"),
                       c("0", "180"),
                       c("0", "365"))

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  dplyr::group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(
    Antibody_titre ~ Time_Point, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_x_position(x = "Time_Point") %>%
  rstatix::add_y_position(y.trans = log2)


# Antibody_titre
Fig_S2a <-
  ggboxplot(data = dplyr::filter(manuscript_data), x = "Time_Point", y = "Antibody_titre",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", step.increase = 0.035, tip.length = 0.02)  + # Add statistical test p-values
  geom_hline(yintercept = 10, color = "#d02310", alpha = 0.2)+
  geom_text(aes(x=1.6,y = 6,label = "protective anti-HBs titer ", vjust = -0.5), size = 2, alpha = 0.2, color = "#d02310") +
  xlab("Time Point in Days") +
  ylab(expression("Antibody titer")) +
  scale_y_continuous(
    trans = "log2",
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = c("0", "10", "100", "1000", "10000", "100000")) +
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 10),
        axis.title  = element_text(size = 10),
        legend.position = "none",
        aspect.ratio = 2.2)
ggsave(plot = Fig_S2a, file=paste(c(supp_fig, "Fig_S2a_Antibody_titer_over_time", ".png"), collapse = ""), height = 12, width = 15, dpi = 600, units = "cm")

########################################################################################################
# Fig_S2b: CMV_EBV_HSV
########################################################################################################

# Tabulate for CMV_EBV_HSV
CMV_EBV_HSV <- filter(manuscript_data, Time_Point == 0) %>% select(Status_2, CMV_Status, EBV_Status, HSV1_2_Status)
table(CMV = CMV_EBV_HSV$CMV_Status, HBV = CMV_EBV_HSV$Status_2)
chisq.test(CMV_EBV_HSV$CMV_Status, CMV_EBV_HSV$Status_2) 
fisher.test(CMV_EBV_HSV$CMV_Status, CMV_EBV_HSV$Status_2)

table(EBV = CMV_EBV_HSV$EBV_Status, HBV = CMV_EBV_HSV$Status_2)
chisq.test(CMV_EBV_HSV$EBV_Status, CMV_EBV_HSV$Status_2) 
fisher.test(CMV_EBV_HSV$EBV_Status, CMV_EBV_HSV$Status_2)

table(HSV = CMV_EBV_HSV$HSV1_2_Status, HBV = CMV_EBV_HSV$Status_2)
chisq.test(CMV_EBV_HSV$HSV1_2_Status, CMV_EBV_HSV$Status_2) 
fisher.test(CMV_EBV_HSV$HSV1_2_Status, CMV_EBV_HSV$Status_2)

######################################################################################################### 
converter_group_CMV_Status <- manuscript_data %>%
  dplyr::filter(Time_Point == 0) %>%
  dplyr::group_by(Status_2, CMV_Status) %>%
  dplyr::summarise(counts = n()) 

plot_1 <- 
  ggplot(converter_group_CMV_Status, aes(x = Status_2, y = counts)) +
  geom_bar(aes(fill = CMV_Status),
           stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values =  c("#0079bb","#bb4200")) +
  ylim(0,25) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "top", aspect.ratio = 2)
######################################################################################################### 
converter_group_EBV_Status <- manuscript_data %>%
  dplyr::filter(Time_Point == 0) %>%
  group_by(Status_2, EBV_Status) %>%
  summarise(counts = n()) 

plot_2 <- 
  ggplot(converter_group_EBV_Status, aes(x = Status_2, y = counts)) +
  geom_bar(aes(fill = EBV_Status),
           stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values =  c("#0079bb","#bb4200")) +
  ylim(0,25) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "top", aspect.ratio = 2)
######################################################################################################### 
converter_group_HSV_Status <- manuscript_data %>%
  dplyr::filter(Time_Point == 0) %>%
  group_by(Status_2, HSV1_2_Status) %>%
  summarise(counts = n()) 

plot_3 <- 
  ggplot(converter_group_HSV_Status, aes(x = Status_2, y = counts)) +
  geom_bar(aes(fill = HSV1_2_Status),
           stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values =  c("#0079bb","#bb4200")) +
  ylim(0,25) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "top", aspect.ratio = 2)

Fig_S2b <- gridExtra::grid.arrange(plot_1, plot_2, plot_3, ncol = 3, top = textGrob("CMV, EBV and HSV seropositivity in the cohort ",gp=gpar(fontsize=18,fontface="bold")))
ggsave(plot = Fig_S2b, file=paste(c(supp_fig, "Fig_S2b_CMV_EBV_HSV", ".png"), collapse = ""), height = 12, width = 16, dpi = 600, units = "cm")

########################################################################################################
########################################################################################################

my_comparisons <- list(c("Early-converter", "Late-converter"),
                       c("Late-converter", "Non-converter"),
                       c("Early-converter", "Non-converter"))

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Time_Point == 0) %>%
  rstatix::wilcox_test(Age ~ Status_2, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

plot_1 <-
  ggboxplot(data = dplyr::filter(manuscript_data, Time_Point == 0), x = "Status_2", y = "Age",
            color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression("Age")) +
  theme(axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 10),
        axis.title.y  = element_text(size = 11),
        axis.title.x  = element_blank(),
        legend.position = "none",
        aspect.ratio = 2)

Fig_S2c <- gridExtra::grid.arrange(plot_1, ncol = 1, top = textGrob("Age per group",gp=gpar(fontsize=18,fontface="bold")))

ggsave(plot = Fig_S2c, file=paste(c(supp_fig, "Fig_S2c_Age_per_group", ".png"), collapse = ""), height = 12, width = 6, dpi = 600, units = "cm")

########################################################################################################
# Fig_S4: plot CFSE_assay for all vaccinees and per vaccinee
########################################################################################################

my_comparisons <- list(c("Single Peptide", "Peptide Pool"),
                       c("Negative Control", "Peptide Pool"),
                       c("Negative Control", "Single Peptide"))

# wilcox_test
stat.test <- 
  dplyr::filter(CFSE_assay) %>%
  rstatix::wilcox_test(CFSElow_CD4_freq ~ Stimulation_Condition, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Stimulation_Condition")

plot_1 <-
  ggplot(data = dplyr::filter(CFSE_assay), aes(x = Stimulation_Condition, y = CFSElow_CD4_freq, color = Stimulation_Condition)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.4) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.001)  + # Add statistical test p-values
  guides(color=guide_legend(title="Stimulation Condition")) + 
  ylab(expression(paste("% CFSE"^"low", " CD4 T cells / ", "CD4 T cells"))) +
  scale_y_continuous(breaks=seq(0, 30, by = 10), limits = c(0, 55)) + 
  theme_minimal() + 
  theme(axis.text.x   = element_blank(),
        axis.text.y   = element_text(size = 8, face = "bold"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 12, face = "bold"), 
        legend.position = "none",
        aspect.ratio = 2.5) + scale_color_npg()

# wilcox_test
stat.test <- 
  dplyr::filter(CFSE_assay) %>%
  dplyr::group_by(Status_2) %>%
  rstatix::wilcox_test(CFSElow_CD4_freq ~ Stimulation_Condition, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Stimulation_Condition")

plot_2 <- 
  ggplot(data = dplyr::filter(CFSE_assay), aes(x = Stimulation_Condition, y = CFSElow_CD4_freq, color = Stimulation_Condition)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25, height = 0, size = 1, alpha = 0.4) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.001)  + # Add statistical test p-values
  guides(color=guide_legend(title="Stimulation Condition")) + 
  ylab("") +
  scale_y_continuous(breaks=seq(0, 30, by = 10), limits = c(0, 55)) + 
  theme_minimal() + 
  theme(axis.text.x   = element_blank(),
        axis.text.y   = element_text(size = 8, face = "bold"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 12, face = "bold"), 
        legend.position = "none",
        aspect.ratio = 2.5) + scale_color_npg() + facet_wrap(~Status_2)

plot_3 <-
  ggplot(data = filter(CFSE_assay), aes(x = Stimulation_Condition, y = CFSElow_CD4_freq, color = Stimulation_Condition)) + 
  geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.6) +
  ylim(0, 38) +
  ylab("") +
  guides(color=guide_legend(title="Stimulation Condition")) + 
  theme_minimal() + 
  theme(axis.text.x   = element_blank(),
        axis.text.y   = element_text(size = 8, face = "bold"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 12, face = "bold"), 
        legend.position = "none",
        aspect.ratio = 1) +
  scale_color_npg() + facet_wrap(~Vaccinee, ncol=5) 

plot_legend <-
  get_legend(p = ggplot(data = filter(CFSE_assay), aes(x = Stimulation_Condition, y = CFSElow_CD4_freq, color = Stimulation_Condition)) + 
               geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.6) +
               ylim(0, 38) +
               ylab(expression(paste("% CFSE"^"low", " CD4 T cells / ", "CD4 T cells"))) +
               guides(color=guide_legend(title="Stimulation Condition")) + 
               theme_minimal() + 
               theme(axis.text.x   = element_blank(),
                     axis.text.y   = element_text(size = 8, face = "bold"),
                     axis.title.x  = element_blank(),
                     axis.title.y  = element_text(size = 12, face = "bold"), 
                     legend.position = "bottom",
                     aspect.ratio = 1) +
               scale_color_npg() + facet_wrap(~Vaccinee))

grid.arrange_matrix <- matrix(data = c(rep(1,12), rep(c(2,2,3,3,3,3),4)), nrow = 6, ncol = 6, byrow = FALSE)
grid.arrange_matrix <- rbind(grid.arrange_matrix, rep(4,6))

Fig_S4 <- gridExtra::grid.arrange(plot_1, plot_2, plot_3, as_ggplot(plot_legend), layout_matrix = grid.arrange_matrix,
                                 top = textGrob("In vitro expansion of CD4 T cells (7 days)",gp=gpar(fontsize=18,fontface="bold")))

ggsave(plot = Fig_S4, file=paste(c(supp_fig, "Fig_S4_CFSE_assay", ".png"), collapse = ""), height = 27.7, width = 17, dpi = 600, units = "cm")

########################################################################################################
# Fig_S5b: CD25 and CD127 expression on vaccine-specific CD4 T cells  
########################################################################################################

# dissected by converse expression of CD40L and 4-1BB

# Load metadata
Table_extra <- read_csv("Raw data/CD25_MFI_CD127_MFI_Tcon_Treg.csv")

colnames(Table_extra)[1] <- "Sample"

# replace "[ \t]\\|[ \t]Count" with ""
colnames(Table_extra) <- str_replace(colnames(Table_extra), pattern = "[ \t]\\|[ \t]Median", replacement = "")

# replace "..Comp-Pacific Blue-A." with "_CD25_MFI"
colnames(Table_extra) <- str_replace(colnames(Table_extra), pattern = "..Comp-Pacific Blue-A.", replacement = "_CD25_MFI")

# replace "..Comp-BV786-A." with "_CD127_MFI"
colnames(Table_extra) <- str_replace(colnames(Table_extra), pattern = "..Comp-BV786-A.", replacement = "_CD127_MFI")

# replace "..Comp-FITC-A." with "_CD45RA_MFI"
colnames(Table_extra) <- str_replace(colnames(Table_extra), pattern = "..Comp-FITC-A.", replacement = "_CD45RA_MFI")

# Remove Mean and SD rows
Table_extra <- Table_extra %>% dplyr::filter(Sample != "Mean" & Sample != "SD")
# Remove Compensation
Table_extra <- Table_extra[ grep("Compensation", Table_extra$Sample, invert = TRUE) , ]

# remove "_CD4 T cells.fcs" from end of Sample
# remove _R at the end of Sample file
Table_extra <- Table_extra %>% mutate(Sample = stringr::str_remove(Table_extra$Sample, pattern = "_CD4 T cells.fcs"))
Table_extra <- Table_extra %>% mutate(Sample = stringr::str_remove(Table_extra$Sample, pattern = "_R"))

# extract Vaccinee, Condition and Time_Point from Sample and reorder
Table_extra <- Table_extra %>% mutate(Vaccinee = as.numeric(stringr::str_extract(Sample, pattern = "\\d\\d")),
                                      Condition =            stringr::str_extract(Sample, pattern = "(NC|HBsAg|PC)\\d*"),
                                      Time_Point =            stringr::str_extract(Sample, pattern = "[:digit:]+$")) %>%
  dplyr::select(Sample, Vaccinee, Time_Point, Condition, everything())

# Check
table(Table_extra$Condition, Table_extra$Time_Point)

# Vaccinee as character
Table_extra[ , c(2)] <- map(Table_extra[ , c(2)], as.character)

# Time_Point and Condition as factor
Table_extra[ , c(3:4)] <- map(Table_extra[ , c(3:4)], as.factor)

# order Time_Point
Table_extra$Time_Point <- factor(Table_extra$Time_Point, levels = c("0", "60", "180", "365"))

# fitler for Condition HBsAg and Time_Point 60
Table_extra_HBsAg_60 <- filter(Table_extra, Condition == "HBsAg" & Time_Point == 60)

# Remove Condition, Time_Point columns
Table_extra_HBsAg_60 <- dplyr::select(Table_extra_HBsAg_60, -c(Sample, Vaccinee, Condition, Time_Point))

# MFI_CD25_CD127_Tcon_Treg
# subset
MFI_CD25_CD127_Tcon_Treg <- subset(Table_extra_HBsAg_60, select = c(CD137neg_CD154pos_CD25_MFI:CD137pos_CD154pos_CD45RA_MFI))
# convert to long formate
MFI_CD25_CD127_Tcon_Treg <- tidyr::gather(MFI_CD25_CD127_Tcon_Treg, key = "Variable", value = "MFI")
# separate Variable by _
MFI_CD25_CD127_Tcon_Treg <- MFI_CD25_CD127_Tcon_Treg %>% separate(Variable, into = c("Subset1", "Subset2", "Marker", "Extra"), sep = "_", extra = "merge", fill = "left") %>%
  # unite Subset1 and Subset2
  unite("Subset", c("Subset1","Subset2"), sep = "_") %>% 
  # remove Extra column
  dplyr::select(-Extra)
# order Subset
MFI_CD25_CD127_Tcon_Treg$Subset <- factor(MFI_CD25_CD127_Tcon_Treg$Subset, levels = c("CD137neg_CD154pos","CD137pos_CD154neg","CD137pos_CD154pos"))

my_comparisons <- list(c("CD137neg_CD154pos", "CD137pos_CD154neg"))

plot_1 <-
  ggboxplot(data = filter(MFI_CD25_CD127_Tcon_Treg, Marker == "CD25" & Subset != "CD137pos_CD154pos"), x = "Subset", y = "MFI",
            color = "Subset", add = "jitter", add.params = list(alpha = 0.5))+
  scale_x_discrete(labels=c("CD137neg_CD154pos" = expression(paste("CD40L"^"+", "4-1BB"^"−", " CD4 T cells")), 
                            "CD137pos_CD154neg" = expression(paste("CD40L"^"−", "4-1BB"^"+", " CD4 T cells")))) +
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0, 0.1))) +
  ylab(expression(paste("CD25 MFI"))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme(axis.title.x =  element_blank(),
        axis.text.x = element_text(size = 14, face = "bold",  angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "none",
        aspect.ratio = 2) +
  ggsci::scale_color_aaas()

plot_2 <-
  ggboxplot(data = filter(MFI_CD25_CD127_Tcon_Treg, Marker == "CD127" & Subset != "CD137pos_CD154pos"), x = "Subset", y = "MFI",
            color = "Subset", add = "jitter", add.params = list(alpha = 0.5))+
  scale_x_discrete(labels=c("CD137neg_CD154pos" = expression(paste("CD40L"^"+", "4-1BB"^"−", " CD4 T cells")), 
                            "CD137pos_CD154neg" = expression(paste("CD40L"^"−", "4-1BB"^"+", " CD4 T cells")))) +
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0, 0.1))) +
  ylab(expression(paste("CD127 MFI"))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme(axis.title.x =  element_blank(),
        axis.text.x = element_text(size = 14, face = "bold",  angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "none",
        aspect.ratio = 2) +
  ggsci::scale_color_aaas()

Fig_S5b <-
  gridExtra::grid.arrange(plot_1, plot_2, ncol = 2,
                          top = textGrob("CD25 and CD127 expression on vaccine-specific CD4 T cells \n dissected by converse expression of CD40L and 4-1BB",gp=gpar(fontface="bold", fontsize=14)))

ggsave(plot = Fig_S5b, file=paste(c(supp_fig, "Fig_S5b_CD25_CD137_CD40L_41BB", ".png"), collapse = ""), height = 15, width = 18, dpi = 600, units = "cm")

########################################################################################################
# Fig_S6
########################################################################################################

correlation_plot_1 <-
  ggscatter(data = dplyr::filter(manuscript_data, Time_Point == 60), 
            x = "CD137neg_CD154pos_memory_CD4_out_memory_CD4", 
            y = "Antibody_titre_FC_365_0",
            fill = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"), 
            shape = 21, size = 4, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "#E7B800", fill = "lightgray", size = 1), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "spearman",label.sep = "\n", size = 5)) +
  scale_y_continuous(breaks=seq(0, 12.5, by = 2.5), limits = c(0,12.5)) +
  xlab(expression(paste("CD40L"^"+", "4-1BB"^"−", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  ylab(expression(paste("log2 of difference in antibody titer"))) +
  theme(aspect.ratio = 1,
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none")

correlation_plot_2 <-
  ggscatter(data = dplyr::filter(manuscript_data, Time_Point == 60), 
            x = "CD137pos_CD154neg_memory_CD4_out_memory_CD4", 
            y = "Antibody_titre_FC_365_0",
            fill = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"), 
            shape = 21, size = 4, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "#E7B800", fill = "lightgray", size = 1), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "spearman",label.sep = "\n", size = 5)) +
  scale_y_continuous(breaks=seq(0, 12.5, by = 2.5), limits = c(0,12.5)) +
  xlab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  ylab(expression(paste("log2 of difference in antibody titer"))) +
  theme(aspect.ratio = 1,
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none")

plot_legend <-
  get_legend(p = ggscatter(data = dplyr::filter(manuscript_data, Time_Point == 60), 
                           x = "CD137pos_CD154neg_memory_CD4_out_memory_CD4", 
                           y = "Antibody_titre_FC_365_0",
                           fill = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"), 
                           shape = 21, size = 4, # Points color, shape and size
                           add = "reg.line",  # Add regressin line
                           add.params = list(color = "#E7B800", fill = "lightgray", size = 1), # Customize reg. line
                           conf.int = TRUE, # Add confidence interval
                           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                           cor.coeff.args = list(method = "spearman",label.sep = "\n", size = 5)) +
               geom_label_repel(data = dplyr::filter(manuscript_data, Time_Point == 60), aes(label = Vaccinee), 
                                label.size = 0.05,
                                segment.size = 0.1,
                                box.padding   = 0.2, 
                                point.padding = 0.2,
                                color = 'red',
                                segment.color = 'red') + 
               xlab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ", "Memory CD4 T cells"))) +
               ylab(expression(paste("log2 of difference in antibody titer"))) +
               theme(aspect.ratio = 1,
                     axis.text  = element_text(size = 10),
                     axis.title = element_text(size = 10, face = "bold"),
                     legend.title = element_blank()))

Fig_S6 <- gridExtra::grid.arrange(as_ggplot(plot_legend), correlation_plot_1, correlation_plot_2, layout_matrix = matrix(data = c(1,1,1,1,2,2,3,3,2,2,3,3,2,2,3,3), nrow = 4, ncol = 4, byrow = TRUE),
                                 top = textGrob("Correlation of vaccine-specific memory CD4 T cells at day 60 and \n antibody titer difference between day 365 and day 0",gp=gpar(fontsize=18,fontface="bold")))
ggsave(plot = Fig_S6, file=paste(c(supp_fig, "Fig_S6_Antigen_specific_memory_CD4_out_memory_CD4_vs_Antibody_titer_FC_with_geom_text", ".png"), collapse = ""), height = 16, width = 23, dpi = 600, units = "cm")

########################################################################################################
# Fig_S7a: CD137neg_CD154pos_memory_CD4_out_memory_CD4 and CD137pos_CD154neg_memory_CD4_out_memory_CD4
########################################################################################################

my_comparisons <- list(c("0", "60"),
                       c("60", "180"),
                       c("180", "365"),
                       c("0", "180"),
                       c("0", "365"))

# CD137neg_CD154pos_memory_CD4_out_memory_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(CD137neg_CD154pos_memory_CD4_out_memory_CD4 ~ Time_Point, 
                                paired = TRUE, comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Time_Point")

plot_1 <-
  ggboxplot(data = dplyr::filter(manuscript_data), x = "Time_Point", y = "CD137neg_CD154pos_memory_CD4_out_memory_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.0001)  + # Add statistical test p-values
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"+", "4-1BB"^"−", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        aspect.ratio = 2)

# CD137pos_CD154neg_memory_CD4_out_memory_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(CD137pos_CD154neg_memory_CD4_out_memory_CD4 ~ Time_Point, 
                                paired = TRUE, comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Time_Point")

plot_2 <-
  ggboxplot(data = dplyr::filter(manuscript_data), x = "Time_Point", y = "CD137pos_CD154neg_memory_CD4_out_memory_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.0001)  + # Add statistical test p-values
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        aspect.ratio = 2)

Fig_S7a <- gridExtra::grid.arrange(plot_1, plot_2, ncol = 2, top = textGrob("Vaccine-specific memory CD4 T cells over time \n for early, late and non-converters",gp=gpar(fontsize=18,fontface="bold")))
ggsave(plot = Fig_S7a, file=paste(c(supp_fig, "Fig_S7a_Antigen_specific_memory_CD4_out_memory_CD4_per_group", ".png"), collapse = ""), height = 18, width = 36, dpi = 600, units = "cm")

########################################################################################################
# Fig_S7b: CD137neg_CD154pos_CD4_out_CD4 and CD137pos_CD154neg_CD4_out_CD4
########################################################################################################

my_comparisons <- list(c("0", "60"),
                       c("60", "180"),
                       c("180", "365"),
                       c("0", "180"),
                       c("0", "365"))

# CD137neg_CD154pos_CD4_out_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(CD137neg_CD154pos_CD4_out_CD4 ~ Time_Point, 
                                paired = TRUE, comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Time_Point")

plot_1 <-
  ggboxplot(data = filter(manuscript_data), x = "Time_Point", y = "CD137neg_CD154pos_CD4_out_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.0001)  + # Add statistical test p-values
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"+", "4-1BB"^"−", "CD4 T cells / ", "10"^"6", " Total CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        aspect.ratio = 2)

# CD137pos_CD154neg_CD4_out_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(CD137pos_CD154neg_CD4_out_CD4 ~ Time_Point, 
                                paired = TRUE, comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Time_Point")

plot_2 <-
  ggboxplot(data = filter(manuscript_data), x = "Time_Point", y = "CD137pos_CD154neg_CD4_out_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.0001)  + # Add statistical test p-values
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"−", "4-1BB"^"+", "CD4 T cells / ", "10"^"6", " Total CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        aspect.ratio = 2)

Fig_S7b <- gridExtra::grid.arrange(plot_1, plot_2, ncol = 2, top = textGrob("Vaccine-specific CD4 T cells over time \n for early, late and non-converters",gp=gpar(fontsize=18,fontface="bold")))
ggsave(plot = Fig_S7b, file=paste(c(supp_fig, "Fig_S7b_Antigen_specific_CD4_out_CD4_per_group", ".png"), collapse = ""), height = 18, width = 36, dpi = 600, units = "cm")

########################################################################################################
# Fig_S8
########################################################################################################

# CD137pos_CD45RAneg_Treg_out_CD4_freq
# AND
# CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq

my_comparisons <- list(c("Early-converter", "Late-converter"),
                       c("Early-converter", "Non-converter"),
                       c("Late-converter", "Non-converter"))

manuscript_data$CD137pos_CD45RAneg_Treg_out_CD4_freq_perc <- manuscript_data$CD137pos_CD45RAneg_Treg_out_CD4_freq/10000

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Time_Point) %>%
  rstatix::wilcox_test(CD137pos_CD45RAneg_Treg_out_CD4_freq_perc ~ Status_2, 
                                comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

plot_1 <- 
  ggboxplot(data = filter(manuscript_data), x = "Status_2", y = "CD137pos_CD45RAneg_Treg_out_CD4_freq_perc",
            color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
            add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", step.increase = 0.02, tip.length = 0.02)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("% 4-1BB"^"+", " CD45RA"^"−", "Treg / ", "CD4 T cells"))) +
  facet_grid(.~Time_Point)+
  theme(legend.title=element_blank(), 
        axis.text.x  = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom", aspect.ratio = 2) 

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  group_by(Time_Point) %>%
  rstatix::wilcox_test(CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq ~ Status_2, 
                       comparisons = my_comparisons, p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

plot_2 <- 
  ggboxplot(data = filter(manuscript_data), x = "Status_2", y = "CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq",
            color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
            add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", step.increase = 0.02, tip.length = 0.02)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("% 4-1BB"^"+", " CD45RA"^"−", "Treg / ", "CD45RA"^"−", " Treg"))) +
  ylim(0,65) +
  facet_grid(.~Time_Point)+
  theme(legend.title=element_blank(), 
        axis.text.x  = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom", aspect.ratio = 2)

Fig_S8 <-
  gridExtra::grid.arrange(plot_1, plot_2, nrow = 2,
                          top = textGrob("",gp=gpar(fontface="bold", fontsize=12)))

ggsave(plot = Fig_S8, file=paste(c(supp_fig, "Fig_S8_CD137pos_CD45RAneg_Treg_out_CD4_AND_out_CD45RAneg_Treg_freq_all_Time_Point", ".png"), collapse = ""), height = 20, width = 20, dpi = 600, units = "cm")


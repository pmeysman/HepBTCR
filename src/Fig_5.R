

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

main_fig <- 'figures/main figures'
dir.create(main_fig)

########################################################################################################
# load data
########################################################################################################

manuscript_data <- read_tsv(file.path("results/fcm/df_1.tsv"))
# order Status_2
manuscript_data$Status_2 <- factor(manuscript_data$Status_2, levels = c("Early-converter", "Late-converter", "Non-converter"))

########################################################################################################
# Fig. 5a
########################################################################################################

# CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Time_Point == 0) %>%
  rstatix::wilcox_test(CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq ~ Status_2 ,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

# push y.position upwards
stat.test$y.position <- stat.test$y.position + 5
  

plot_1 <- 
ggboxplot(data = dplyr::filter(manuscript_data, Time_Point == 0), x = "Status_2", y = "CD137pos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq",
          color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
          add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.1)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("%  4-1BB"^"+", " CD45RA"^"−", "Treg / ", "CD45RA"^"−", " Treg"))) +
  ylim(0,55) +
  theme(legend.title=element_blank(), 
        axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        legend.position = "none", aspect.ratio = 2)

Fig_5a <-
  gridExtra::grid.arrange(plot_1, ncol = 1,
                          top = textGrob("Frequency before vaccination",gp=gpar(fontface="bold", fontsize=12)))

ggsave(plot = Fig_5a, file=paste(c(main_fig, "Fig_5a_4-1BBpos_CD45RAneg_Treg_out_CD45RAneg_Treg_freq", ".png"), collapse = ""), height = 12, width = 8, dpi = 600, units = "cm")

########################################################################################################
# Fig. 5b
########################################################################################################

# keep columns with MFI
master_MFI <- dplyr::select(manuscript_data, c(Set:Status_365_0, contains("MFI")))

# fitler for Time_Point 0
master_MFI_NC_0 <- filter(master_MFI, Time_Point == 0)

# Remove Sample, Time_Point columns
master_MFI_NC_0 <- dplyr::select(master_MFI_NC_0, -c(Sample, Time_Point))

# first MFI plot
# subset
for_first_MFI_plot <- subset(master_MFI_NC_0, select = c(Set:Status_365_0, TH_CD137_MFI, TFH_CD137_MFI, Treg_CD137_MFI, TFR_CD137_MFI))
# convert to long formate
for_first_MFI_plot <- tidyr::gather(for_first_MFI_plot, key = "Variable", value = "MFI", -c(Set:Status_365_0))
# separate Variable by _
for_first_MFI_plot <- for_first_MFI_plot %>% separate(Variable, into = c("Subset1", "Subset2", "Marker", "Extra"), sep = "_", extra = "merge", fill = "left") %>%
  # unite Subset1 and Subset2
  unite("Subset", c("Subset1","Subset2"), sep = "_") %>% 
  # remove NA_ from Subset
  dplyr::mutate(Subset = stringr::str_replace_all(Subset, 'NA_?', '')) %>%
  # remove Extra column
  dplyr::select(-Extra)
# order Subset
for_first_MFI_plot$Subset <- factor(for_first_MFI_plot$Subset, levels = c("TH","TFH","Treg","TFR"))

my_comparisons <- list( c("TH", "Treg"),
                        c("TFH", "Treg"),
                        c("Treg", "TFR"))

# wilcox_test
stat.test <- 
  filter(for_first_MFI_plot, grepl("^T",Subset)) %>%
  rstatix::pairwise_wilcox_test(
    MFI ~ Subset, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "Subset", y.trans = log10)

plot_1 <-
  ggboxplot(data = filter(for_first_MFI_plot, grepl("^T",Subset)), 
            x = "Subset", y = "MFI",
               color = "Subset", add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", step.increase = 0.02, tip.length = 0.02)  + # Add statistical test p-values
  scale_x_discrete(labels=c("TH" = expression('T'[H]), 
                            "TFH" = expression('T'[FH]),
                            "Treg" = expression('T'[REG]),
                            "TFR" = expression('T'[FR]))) +
  scale_y_continuous(trans = "log10") +
  ylab(expression(paste("4-1BB MFI"))) +
  #stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme(axis.title.x =  element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        aspect.ratio = 2)

Fig_5b <-
  gridExtra::grid.arrange(plot_1, ncol = 1,
                          top = textGrob(expression(bold('4-1BB expression on CD4 T cell subsets')),gp=gpar(fontface="bold", fontsize=14)))
ggsave(plot = Fig_5b, file=paste(c(main_fig, "Fig_5b_MFI_4-1BB_CD4_subsets", ".png"), collapse = ""), height = 12, width = 10, dpi = 600, units = "cm")

########################################################################################################
# Fig. 5c
########################################################################################################

# second MFI plot
for_second_MFI_plot <- subset(master_MFI_NC_0, select = c(Set:Status_365_0, CD45RAneg_Treg_CD137_MFI, CD45RApos_Treg_CD137_MFI, CD45RAneg_Treg_CD25_MFI, CD45RApos_Treg_CD25_MFI))

# convert to long formate
for_second_MFI_plot <- tidyr::gather(for_second_MFI_plot, key = "Variable", value = "MFI", -c(Set:Status_365_0))

# separate Variable by _
for_second_MFI_plot <- for_second_MFI_plot %>% separate(Variable, into = c("Subset1", "Subset2", "Marker", "Extra"), sep = "_", extra = "merge", fill = "left") %>%
  # unite Subset1 and Subset2
  unite("Subset", c("Subset1","Subset2"), sep = "_") %>% 
  # remove NA_ from Subset
  dplyr::mutate(Subset = stringr::str_replace_all(Subset, 'NA_?', '')) %>%
  # remove Extra column
  dplyr::select(-Extra)

# order Subset
for_second_MFI_plot$Subset <- factor(for_second_MFI_plot$Subset, levels = c("CD45RAneg_Treg","CD45RApos_Treg"))

my_comparisons <- list( c("CD45RAneg_Treg", "CD45RApos_Treg"))

plot_1 <-
  ggboxplot(data = filter(for_second_MFI_plot, grepl("CD137",Marker))
               , x = "Subset", y = "MFI",
               color = "Subset", add = "jitter", add.params = list(alpha = 0.5),
               palette = c("#00e2e8","#047c80"))+
  scale_x_discrete(labels=c("CD45RAneg_Treg" = expression(paste("CD45RA"^"−", " T"[REG])), 
                            "CD45RApos_Treg" = expression(paste("CD45RA"^"+", " T"[REG])))) +
  scale_y_continuous(trans = "log10") +
  ylab(expression(paste("4-1BB MFI"))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "none", aspect.ratio = 2)

plot_2 <-
  ggboxplot(data = filter(for_second_MFI_plot, grepl("CD25",Marker))
               , x = "Subset", y = "MFI",
               color = "Subset",add = "jitter", add.params = list(alpha = 0.5),
               palette = c("#00e2e8","#047c80"))+
  scale_x_discrete(labels=c("CD45RAneg_Treg" = expression(paste("CD45RA"^"−", " T"[REG])), 
                            "CD45RApos_Treg" = expression(paste("CD45RA"^"+", " T"[REG])))) +
  scale_y_continuous(trans = "log10") +
  ylab(expression(paste("CD25 MFI"))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "none", aspect.ratio = 2)

Fig_5c <-
  gridExtra::grid.arrange(plot_1, plot_2, ncol = 2,
                          top = textGrob(expression(bold('CD25 and 4-1BB expression in T'[REG]*' cell subsets')),gp=gpar(fontface="bold", fontsize=12)))

ggsave(plot = Fig_5c, file=paste(c(main_fig, "Fig_5c_MFI_4-1BB_CD25_Treg_subsets", ".png"), collapse = ""), height = 12, width = 12, dpi = 600, units = "cm")

########################################################################################################
# Fig. 5d
########################################################################################################

# CD45RAneg_Treg_out_CD4_freq
my_comparisons <- list(c("Early-converter", "Late-converter"),
                       c("Late-converter", "Non-converter"),
                       c("Early-converter", "Non-converter"))
manuscript_data$Treg_out_CD4_perc <- manuscript_data$Treg_out_CD4_freq/10000

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Time_Point == 0) %>%
  rstatix::wilcox_test(Treg_out_CD4_perc ~ Status_2, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

p1 <- 
ggboxplot(data = filter(manuscript_data, Time_Point == 0), x = "Status_2", y = "Treg_out_CD4_perc",
          color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
          add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.1)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("% Treg / ", "Total CD4 T cells"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.y  = element_text(size = 8),
        legend.position = "none", aspect.ratio = 1.8)

manuscript_data$CD45RAneg_Treg_out_CD4_perc <- manuscript_data$CD45RAneg_Treg_out_CD4_freq/10000

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Time_Point == 0) %>%
  rstatix::wilcox_test(CD45RAneg_Treg_out_CD4_perc ~ Status_2, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

p2 <- 
ggboxplot(data = filter(manuscript_data, Time_Point == 0), x = "Status_2", y = "CD45RAneg_Treg_out_CD4_perc",
          color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
          add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.1)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("% CD45RA"^"−", "Treg / ", "Total CD4 T cells"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.y  = element_text(size = 8),
        legend.position = "none", aspect.ratio = 1.8)

manuscript_data$CD45RApos_Treg_out_CD4_perc <- manuscript_data$CD45RApos_Treg_out_CD4_freq/10000

# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Time_Point == 0) %>%
  rstatix::wilcox_test(CD45RApos_Treg_out_CD4_perc ~ Status_2, comparisons = my_comparisons,
                       p.adjust.method = "bonferroni") %>% 
  rstatix::add_xy_position(x = "Status_2")

p3 <- 
ggboxplot(data = filter(manuscript_data, Time_Point == 0), x = "Status_2", y = "CD45RApos_Treg_out_CD4_perc",
          color = "Status_2", palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"),
          add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.1)  + # Add statistical test p-values
  xlab("Status") +
  ylab(expression(paste("% CD45RA"^"+", "Treg / ", "Total CD4 T cells"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.y  = element_text(size = 8),
        legend.position = "none", aspect.ratio = 1.8)

Fig_5d <-
gridExtra::grid.arrange(p1, p2, p3, ncol = 3,
                        top = textGrob(expression(bold('Frequency of T'[REG]*' cells in total CD4 T cells')),gp=gpar(fontface="bold", fontsize=12)))

ggsave(plot = Fig_5d, file=paste(c(main_fig, "Fig_5d_Treg_out_CD4", ".png"), collapse = ""), height = 9, width = 15, dpi = 600, units = "cm")

########################################################################################################
# Fig. 5e
########################################################################################################

# Subsets of Treg cells by CD45RA and 4-1BB
# Plot Status_2
Treg_profile <- select(manuscript_data, Vaccinee, Status_2, Time_Point,
                       CD137neg_CD45RAneg_Treg_out_Treg_freq,
                       CD137pos_CD45RAneg_Treg_out_Treg_freq,
                       CD137neg_CD45RApos_Treg_out_Treg_freq,
                       CD137pos_CD45RApos_Treg_out_Treg_freq) %>%
  dplyr::filter(Time_Point == 0) %>%
  group_by(Status_2) %>%
  summarise(CD137neg_CD45RAneg_Treg_out_Treg_freq_median = median(CD137neg_CD45RAneg_Treg_out_Treg_freq),
            CD137pos_CD45RAneg_Treg_out_Treg_freq_median = median(CD137pos_CD45RAneg_Treg_out_Treg_freq),
            CD137neg_CD45RApos_Treg_out_Treg_freq_median = median(CD137neg_CD45RApos_Treg_out_Treg_freq),
            CD137pos_CD45RApos_Treg_out_Treg_freq_median = median(CD137pos_CD45RApos_Treg_out_Treg_freq)) %>%
  pivot_longer(cols = c(CD137neg_CD45RAneg_Treg_out_Treg_freq_median, 
                        CD137pos_CD45RAneg_Treg_out_Treg_freq_median,
                        CD137neg_CD45RApos_Treg_out_Treg_freq_median,
                        CD137pos_CD45RApos_Treg_out_Treg_freq_median),
               names_to = "subset",
               values_to = "freq_median")

# round freq_median
Treg_profile$freq_median <- round(Treg_profile$freq_median, digits = 1)
# convert subset to factor and reorder
Treg_profile$subset <- as.factor(Treg_profile$subset)
Treg_profile$subset <- factor(Treg_profile$subset, 
                              levels = c("CD137neg_CD45RApos_Treg_out_Treg_freq_median",
                                         "CD137pos_CD45RApos_Treg_out_Treg_freq_median",
                                         "CD137neg_CD45RAneg_Treg_out_Treg_freq_median",
                                         "CD137pos_CD45RAneg_Treg_out_Treg_freq_median"))


plot <- 
  ggplot(Treg_profile, aes(x = Status_2, y = freq_median)) +
  geom_bar(aes(fill = subset),
           stat = "identity", position = "fill", width = 0.6) +
  #geom_text(aes(label=freq_median), vjust=-1, size=6) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  scale_fill_manual(values = c("#047c80","#e87a00","#00e2e8","#e8006e"),
                    labels = c(expression(paste("4-1BB"^"−", " CD45RA"^"+", " T"[REG])),
                               expression(paste("4-1BB"^"+", " CD45RA"^"+", " T"[REG])),
                               expression(paste("4-1BB"^"−", " CD45RA"^"−", " T"[REG])),
                               expression(paste("4-1BB"^"+", " CD45RA"^"−", " T"[REG])))) +
  ylab(expression(paste("% in", " T"[REG])))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 8, face = "bold"),
        legend.title = element_blank(), 
        legend.text  = element_text(size = 8), 
        legend.position = "left", aspect.ratio = 2)

Fig_5e <-
  gridExtra::grid.arrange(plot, ncol = 1,
                          top = textGrob(expression(bold('Subsets of T'[REG]*' cells by CD45RA and 4-1BB')),gp=gpar(fontface="bold", fontsize=12)))

ggsave(plot = Fig_5e, file=paste(c(main_fig, "Fig_5e_Treg_composition", ".png"), collapse = ""), height = 10, width = 12, dpi = 600, units = "cm")


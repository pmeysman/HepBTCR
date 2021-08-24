
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
# visualize figure
########################################################################################################

########################################################################################################
# Fig_4a
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
  rstatix::pairwise_wilcox_test(
    CD137neg_CD154pos_memory_CD4_out_memory_CD4 ~ Time_Point, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "Time_Point")

plot_1 <- 
  ggboxplot(data = filter(manuscript_data), x = "Time_Point", y = "CD137neg_CD154pos_memory_CD4_out_memory_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.01)  + # Add statistical test p-values
  ylim(-100, 5000) + 
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"+", "4-1BB"^"−", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text  = element_text(size = 10),
        axis.title  = element_text(size = 10), 
        legend.position = "none",
        aspect.ratio = 2)

# CD137pos_CD154neg_CD4_out_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data) %>%
  rstatix::pairwise_wilcox_test(
    CD137pos_CD154neg_CD4_out_CD4 ~ Time_Point, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "Time_Point")

plot_2 <- 
  ggboxplot(data = filter(manuscript_data), x = "Time_Point", y = "CD137pos_CD154neg_CD4_out_CD4",
            color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
            add = "jitter", add.params = list(alpha = 0.5))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.01)  + # Add statistical test p-values
  ylim(-100, 5000) + 
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text  = element_text(size = 10),
        axis.title  = element_text(size = 10),
        legend.position = "none",
        aspect.ratio = 2)
Fig_4a <- gridExtra::grid.arrange(plot_1, plot_2, ncol = 2, top = textGrob("Vaccine-specific memory CD4 T cells over time",gp=gpar(fontsize=18,fontface="bold")))

ggsave(plot = Fig_4a, file=paste(c(main_fig, "Fig_4a_CD137neg_CD154pos_memory_CD4_and_CD137pos_CD154neg_memory_CD4_out_memory_CD4", ".png"), collapse = ""), height = 15, width = 15, dpi = 600, units = "cm")

########################################################################################################
# Fig_4b
########################################################################################################

correlation_plot_1 <-
  ggscatter(data = dplyr::filter(manuscript_data, Time_Point == 60), 
            x = "CD137neg_CD154pos_memory_CD4_out_memory_CD4", 
            y = "Antibody_titre_FC_365_0",
            fill = "#00AFBB", shape = 21, size = 4, # Points color, shape and size
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
        axis.title = element_text(size = 10, face = "bold"))

correlation_plot_2 <-
  ggscatter(data = dplyr::filter(manuscript_data, Time_Point == 60), 
            x = "CD137pos_CD154neg_memory_CD4_out_memory_CD4", 
            y = "Antibody_titre_FC_365_0",
            fill = "#00AFBB", shape = 21, size = 4, # Points color, shape and size
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
        axis.title = element_text(size = 10, face = "bold"))

Fig_4b <- gridExtra::grid.arrange(correlation_plot_1, correlation_plot_2, ncol = 2,
                                  top = textGrob("Correlation of vaccine-specific memory CD4 T cells at day 60 and \n antibody titer difference between day 365 and day 0",gp=gpar(fontsize=18,fontface="bold")))

ggsave(plot = Fig_4b, file=paste(c(main_fig, "Fig_4b_Antigen_specific_memory_CD4_out_memory_CD4_vs_Antibody_titer_FC", ".png"), collapse = ""), height = 16, width = 23, dpi = 600, units = "cm")

########################################################################################################
# Fig_4c
########################################################################################################

my_comparisons <- list(c("0", "60"),
                       c("60", "180"),
                       c("180", "365"),
                       c("0", "180"),
                       c("0", "365"))

# CD137neg_CD154pos_memory_CD4_out_memory_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Status_2 != "Non-converter") %>%
  dplyr::group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(
    CD137neg_CD154pos_memory_CD4_out_memory_CD4 ~ Time_Point, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "Time_Point")

plot_1 <-
ggboxplot(data = filter(manuscript_data, Status_2 != "Non-converter"), x = "Time_Point", y = "CD137neg_CD154pos_memory_CD4_out_memory_CD4",
          color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
          add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.00001, tip.length = 0.02)  + # Add statistical test p-values
  scale_y_continuous(limits = c(-100, 5000), expand = expansion(mult = c(0, 0.1))) + 
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"+", "4-1BB"^"−", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 10),
        axis.title  = element_text(size = 10),
        legend.position = "none",
        aspect.ratio = 2)

# CD137pos_CD154neg_memory_CD4_out_memory_CD4
# wilcox_test
stat.test <- 
  dplyr::filter(manuscript_data, Status_2 != "Non-converter") %>%
  dplyr::group_by(Status_2) %>%
  rstatix::pairwise_wilcox_test(
    CD137pos_CD154neg_memory_CD4_out_memory_CD4 ~ Time_Point, paired = TRUE, comparisons = my_comparisons,
    p.adjust.method = "bonferroni") %>%
  rstatix::add_xy_position(x = "Time_Point")

plot_2 <-
ggboxplot(data = filter(manuscript_data, Status_2 != "Non-converter"), x = "Time_Point", y = "CD137pos_CD154neg_memory_CD4_out_memory_CD4",
          color = "Time_Point", palette =c("#98b7dd", "#5d8ec9", "#3667a1", "#224167"),
          add = "jitter", add.params = list(alpha = 0.5), facet.by = "Status_2")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = FALSE, step.increase = 0.00001, tip.length = 0.02)  + # Add statistical test p-values
  scale_y_continuous(limits = c(-100, 5000), expand = expansion(mult = c(0, 0.1))) + 
  xlab("Time Point in Days") +
  ylab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 10),
        axis.title  = element_text(size = 10),
        legend.position = "none",
        aspect.ratio = 2)

Fig_4c <- gridExtra::grid.arrange(plot_1, plot_2, ncol = 2, top = textGrob("Vaccine-specific memory CD4 T cells over time \n for early and late-converters",gp=gpar(fontsize=18,fontface="bold")))

ggsave(plot = Fig_4c, file=paste(c(main_fig, "Fig_4c_Antigen_specific_memory_CD4_out_memory_CD4_per_group", ".png"), collapse = ""), height = 15, width = 24, dpi = 600, units = "cm")

########################################################################################################
# Fig_4d
########################################################################################################

my_comparisons <- list(c("Early-converter", "Late-converter"))

# CD137neg_CD154pos_memory_CD4_out_memory_CD4
plot_1 <-
ggboxplot(data = filter(manuscript_data, Status_2 != "Non-converter" & Time_Point == 0), x = "Status_2", y = "CD137neg_CD154pos_memory_CD4_out_memory_CD4",
          color = "Status_2", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", add.params = list(alpha = 0.5))+
  ylab(expression(paste("CD40L"^"+", "4-1BB"^"−", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  scale_y_continuous(limits = c(-100, 5000), expand = expansion(mult = c(0, 0.1))) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        aspect.ratio = 3) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = FALSE, label = "p.signif", label.y = 2750)

# CD137pos_CD154neg_memory_CD4_out_memory_CD4
plot_2 <- 
ggboxplot(data = filter(manuscript_data, Status_2 != "Non-converter" & Time_Point == 0), x = "Status_2", y = "CD137pos_CD154neg_memory_CD4_out_memory_CD4",
          color = "Status_2", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", add.params = list(alpha = 0.5))+
  ylab(expression(paste("CD40L"^"−", "4-1BB"^"+", "Memory CD4 T cells / ","10"^"6", " Memory CD4 T cells"))) +
  scale_y_continuous(limits = c(-100, 5000), expand = expansion(mult = c(0, 0.1))) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        aspect.ratio = 3) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", paired = FALSE, label = "p.signif", label.y = 2750)

Fig_4d <- gridExtra::grid.arrange(plot_1, plot_2, ncol = 2, top = textGrob("Vaccine-specific memory CD4 T cells \n prior to vaccination",gp=gpar(fontsize=10,fontface="bold")))

ggsave(plot = Fig_4d, file=paste(c(main_fig, "Fig_4d_Antigen_specific_memory_CD4_out_memory_CD4_day_0", ".png"), collapse = ""), height = 10, width = 8, dpi = 600, units = "cm")

########################################################################################################
# end
########################################################################################################


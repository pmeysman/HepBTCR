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

manuscript_data <- read_tsv(file="results/fcm/df_1.tsv")
# order Status_2
manuscript_data$Status_2 <- factor(manuscript_data$Status_2, levels = c("Early-converter", "Late-converter", "Non-converter"))

########################################################################################################
# Fig_1b
########################################################################################################

# Plot Status_2
converter_group <- manuscript_data %>%
  dplyr::filter(Time_Point == 0) %>%
  group_by(Status_2) %>%
  summarise(counts = n()) 

converter_group
# Status_2        counts
# 1 Early-converter     21
# 2 Late-converter       9
# 3 Non-converter        4

plot <- 
  ggplot(converter_group, aes(x = Status_2, y = counts)) +
  geom_bar(aes(fill = Status_2),
           stat = "identity", position = position_stack(), width = 0.6) +
  geom_text(aes(label=counts), vjust=-1, size=6) +
  scale_fill_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07")) +
  ylim(0,25) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "none", aspect.ratio = 2)

Fig_1b <-
  gridExtra::grid.arrange(plot, ncol = 1,
                          top = textGrob("Groups in the cohort determined by \n anti-HBs antibody titer serocoversion",gp=gpar(fontface="bold", fontsize=12)))


ggsave(plot = Fig_1b, file=paste(c(main_fig, "Fig_1b_converter_groups", ".png"), collapse = ""), height = 12, width = 8, dpi = 600, units = "cm")
########################################################################################################
########################################################################################################

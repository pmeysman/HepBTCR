#############################
# Fig. S3a
# Scatter plot of the DNA-based TCRβ reads for each vaccinee at each time point
#############################

p<-ggplot(tcr,aes(x=B0reads,y=B60reads, color=responder)) + geom_point(size = 3, alpha = 0.7)
p<-p + scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"), name = "Group")
p<-p + ggrepel::geom_label_repel(data = dplyr::filter(tcr, B60reads>2e+05 | B0reads > 125000 | B0reads < 50000), aes(label = Cohort_ID), show.legend = FALSE)
p<-p + xlab("Number of TCRβ reads at day 0") + ylab("Number of TCRβ reads at day 60")
p<-p + scale_x_continuous(labels = scales::scientific) + scale_y_continuous(labels = scales::scientific)
p<-p + theme_bw()  + theme(legend.position = "none", 
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 12, face = "bold"),
                           axis.title.y = element_text(size = 12, face = "bold"),
                           aspect.ratio = 1) 
p
ggsave(file.path(supp_fig, 'Fig_S3a_TCRreads.png'),p, height = 8, width = 9, units = "cm", dpi = 300)

legend <- p + theme(legend.position = "top")
# save legend separately 
my_legend <- get_legend(p + theme(legend.position = "top"))
legend <- ggpubr::as_ggplot(my_legend)
ggsave(file.path(supp_fig, 'Fig_S3a_and_b_legend.png'),legend, height = 2, width = 20, units = "cm", dpi = 300)

#############################
# Fig. S3b
# Scatter plot of number of unique TCRβ amino acid sequences for each vaccinee at the two time points, colored by the vaccinee group based on antibody response
#############################

p<-ggplot(tcr,aes(x=B0total,y=B60total,color=responder)) + geom_point(size=3, alpha = 0.7)
p<-p + xlab("Unique TCRβ at day 0") + ylab("Unique TCRβ at day 60")
p<-p + scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"),name = "Group")
p<-p + scale_x_continuous(labels = scales::scientific) + scale_y_continuous(labels = scales::scientific)
p<-p + theme_bw()  + theme(legend.position = "none", 
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 12, face = "bold"),
                           axis.title.y = element_text(size = 12, face = "bold"),
                           aspect.ratio = 1) 
p
ggsave(file.path(supp_fig, 'Fig_S3b_uniqueTCRs.png'),p, height = 8, width = 9, units = "cm", dpi = 300)

#############################
# Fig. S3c
# Overview of unique TCRβ amino acid sequences in the memory CD4 T cell repertoire of each vaccinee. 
# The bottom blue bar denotes TCRβ sequences that were found at both time points. 
# The green and red bars denote the number of unique TCRβ sequences at each time point. 
# The total bar height thus represents the total number of unique memory CD4 T cell clonotypes sequences for a specific vaccinee.
#############################

supp<-data.frame(vaccinee = rep(tcr$X,3))
supp$TCRcounts <- c(tcr$B0B60overlap, (tcr$B0total - tcr$B0B60overlap), ( tcr$B60total - tcr$B0B60overlap))
supp$label <- c(rep("Overlap",33),rep("Day0 Unique",33),rep("Day60 Unique",33))

p<-ggplot(supp,aes(x=vaccinee,y=TCRcounts,fill=label)) + geom_col(position = "fill") 
p<-p+scale_y_continuous(labels = scales::percent)
p<-p+ ylab("Number of unique TCRβ sequences") + xlab("Vaccinees")
p<-p + scale_fill_manual(values = c("#4DBBD5FF","#00A087FF","#E64B35FF"),
                         labels = c("unique to day 0", "unique to day 60", "shared between day 0 and day 60"))
p<-p + theme_bw()  + theme(legend.position = "right", 
                           legend.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 15, face = "bold"),
                           axis.title.y = element_text(size = 15, face = "bold"),
                           aspect.ratio = 1) 
p

ggsave(file.path(supp_fig, 'Fig_S3c_TCR_sharing_day_0_day_60 (geom_col fill).png'),p, height = 10, width = 16, units = "cm", dpi = 300)

p<-ggplot(supp,aes(x=vaccinee,y=TCRcounts,fill=label)) + geom_col(position = "stack") 
p<-p+ ylab("Number of unique TCRβ sequences") + xlab("Vaccinees")
p<-p + scale_fill_manual(values = c("#4DBBD5FF","#00A087FF","#E64B35FF"),
                         labels = c("unique to day 0", "unique to day 60", "shared between day 0 and day 60"))
p<-p + theme_bw()  + theme(legend.position = "right", 
                           legend.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 15, face = "bold"),
                           axis.title.y = element_text(size = 15, face = "bold"),
                           aspect.ratio = 1) 
p
ggsave(file.path(supp_fig, 'Fig_S3c_TCR_sharing_day_0_day_60 (geom_col stack).png'),p, height = 10, width = 16, units = "cm", dpi = 300)

#############################
# Fig. S3d
# Change in frequency of HBsAg-specific CD4 T-cells before and after vaccination. The (ns) mark denotes a non-significant paired Wilcoxon signed-rank test (p-value = 0.7577)
#############################

p<-ggplot(time,aes(y=PPIntfreq,x=day, color=day)) + geom_point(size = 3, alpha = 0.7) + geom_boxplot()
p<-p + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
p<-p+geom_line(aes(group = volunteer),alpha = 0.1)
p<-p+scale_y_continuous(labels = scales::scientific)
p<-p + scale_x_discrete(labels=c("Day0" = "0", "Day60" = "60"))
p<-p+ylab("% HBsAg-specific TCRβ / Total TCRβ") + xlab("Time Point (days)")
p<-p + theme_bw() + theme(aspect.ratio = 1.5, 
                          legend.title=element_blank(), 
                          legend.position = "none",
                          plot.title = element_text(size = 20, face = "bold"),
                          axis.text.x  = element_text(size = 10),
                          axis.title.x  = element_text(size = 10, face = "bold"), 
                          axis.title.y = element_text(size = 10, face = "bold"))
p
ggsave(file.path(supp_fig, 'Fig_S3d_PPIntfreq.png'),p, height = 8, width = 6, units = "cm", dpi = 300)

#############################
#############################

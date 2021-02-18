
#############################
#  load packages
#############################

require('ggplot2')
require('ggpubr')
require('gridExtra')
require('scales')
require('ggrepel')
require('ggpubr')
require('pROC')

#############################
#  load data
#############################

# Full repertoire analyses with peptide-specific TCRs
rep <- read.csv(file="results/peptide/loocvrep_ab.txt")

# load TCR data
tcr <- read.table('results/timescatter/volunteercounts.txt',header=TRUE,sep=',',row.names='volunteer')
# add Cohort_ID column (tcr is a data.frame and not a tibble)
tcr$Cohort_ID <- rownames(tcr)

# derive time data.frame from tcr
time<-data.frame(volunteer = rep(tcr$Cohort_ID,2),
                 entropy   = c(tcr$B0entropy,tcr$B60entropy),
                 total     = c(tcr$B0total,tcr$B60total),
                 reads     = c(tcr$B0reads,tcr$B60reads),
                 PPIntfreq = c(tcr$PPIntB0freq,tcr$PPIntB60freq),
                 day       = c(rep('Day0',33),rep('Day60',33)),
                 responder = rep(tcr$responder,2))

#############################
# create directory to save figures
#############################

fig_dir <- 'figures'
dir.create(fig_dir)

main_fig <- 'figures/main figures'
dir.create(main_fig)
supp_fig <- 'figures/supplementary figures'
dir.create(supp_fig)

#############################
# Fig. 2a
# Repertoire diversity and entropy between the two time points
#############################

pte<-ggboxplot(time[time['responder']=='Early-converter',],y="total",x="day", color="day", add = "point")
pte<-pte + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pte<-pte + geom_line(aes(group = volunteer),alpha = 0.1)
pte<-pte + scale_x_discrete(labels=c("Day0" = "0", "Day60" = "60")) + scale_y_continuous(labels = scales::scientific)
pte<-pte + xlab("Time Point (days)") + ylab("Number of unique CDR3s") + ggtitle("Early-converter") + theme(plot.title = element_text(hjust = 0.5))
pte<-pte + theme(aspect.ratio = 1.5,legend.position = "none") 
pte

ptl<-ggboxplot(time[time['responder']=='Late-converter',],y="total",x="day", color="day", add = "point")
ptl<-ptl + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
ptl<-ptl + geom_line(aes(group = volunteer),alpha = 0.1)
ptl<-ptl + scale_x_discrete(labels=c("Day0" = "0", "Day60" = "60")) + scale_y_continuous(labels = scales::scientific)
ptl<-ptl + xlab("Time Point (days)") + ylab("Number of unique CDR3s")  + ggtitle("Late-converter") + theme(plot.title = element_text(hjust = 0.5))
ptl<-ptl + theme(aspect.ratio = 1.5,legend.position = "none") 
ptl

pee<-ggboxplot(time[time['responder']=='Early-converter',],y="entropy",x="day", color="day", add = "point")
pee<-pee + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pee<-pee + geom_line(aes(group = volunteer),alpha = 0.1)
pee<-pee + scale_x_discrete(labels=c("Day0" = "0", "Day60" = "60"))
pee<-pee + ylim(6,14)
pee<-pee + xlab("Time Point (days)") + ylab("Entropy")
pee<-pee + theme(aspect.ratio = 1.5,legend.position = "none") 
pee

pel<-ggboxplot(time[time['responder']=='Late-converter',],y="entropy",x="day", color="day", add = "point")
pel<-pel + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pel<-pel + geom_line(aes(group = volunteer),alpha = 0.1)
pel<-pel + scale_x_discrete(labels=c("Day0" = "0", "Day60" = "60"))
pel<-pel + ylim(6,14)
pel<-pel + xlab("Time Point (days)") + ylab("Entropy")
pel

etplot<-grid.arrange(pte, ptl, pee, pel, nrow = 2, ncol= 2)
ggsave(file.path(main_fig, 'Fig_2a_entropy-total.png'),etplot, height = 13, width = 12, units = "cm", dpi = 300)

#############################
# 
#############################

wilcox.test(tcr$B0entropy,tcr$B60entropy,paired=TRUE)
wilcox.test(tcr$PPIntB0freq,tcr$PPIntB60freq,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$B0entropy,tcr[tcr['responder']=='Late-converter',]$B60entropy,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$B0total,tcr[tcr['responder']=='Late-converter',]$B60total,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$PPIntB0freq,tcr[tcr['responder']=='Late-converter',]$PPIntB60freq,paired=TRUE)

wilcox.test(tcr[tcr['responder']=='Early-converter',]$B0entropy,tcr[tcr['responder']=='Early-converter',]$B60entropy,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Early-converter',]$B0total,tcr[tcr['responder']=='Early-converter',]$B60total,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Early-converter',]$PPIntB0freq,tcr[tcr['responder']=='Early-converter',]$PPIntB60freq,paired=TRUE)

#############################
# Fig. 2b
# Frequency of unique vaccine-specific TCRβ sequences out of 
# total sequenced TCRβ sequences between two time points for all vaccinees colored by group
#############################

# HBV-associated TCR counts are stored in iPPB0 and iPPB60 (for day 0 and day 60)
# Normalize for the size of the repertoire (B0 and B60)
df <- data.frame(PP = c(rep$iPPB0/rep$B0,rep$iPPB60/rep$B60), response = rep$responder, time = c(rep('T0',33),rep('T60',33)),sample=rep(rownames(rep),2))
wt<-wilcox.test(PP~time,data=df,paired=TRUE)
p<-ggplot(df,aes(y=PP,x=time)) + geom_boxplot() + stat_summary(aes(group=sample,color=response),geom="line",fun=sum,alpha = 0.3) + geom_point(stat='summary', fun=sum, aes(group=sample,color=response))
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p+ scale_x_discrete(labels=c("T0" = "0", "T60" = "60"))
p<-p + ylab('% Unique HBsAg-specific TCRβ sequence') + xlab('Time Point (days)')
p<-p + labs(title="",subtitle=paste0("Paired Wilcox P-value = ",format(wt$p.value,digits=2)))
p<-p + theme_bw() + theme(aspect.ratio = 1.5, 
                          legend.title=element_blank(), 
                          legend.position = "none",
                          plot.title = element_text(size = 20, face = "bold"),
                          axis.text.x  = element_text(size = 10),
                          axis.title.x  = element_text(size = 10, face = "bold"), 
                          axis.title.y = element_text(size = 10, face = "bold"))
p
ggsave(file.path(main_fig, 'Fig_2b_tcrincrease.png'),p, scale = 1, height = 12, width = 8, units = "cm", dpi = 300)

# now extract the legend
legend <- as_ggplot(get_legend(p))
ggsave(file.path(main_fig, 'Fig_2_legend.png'),legend, scale = 1, height = 4, width = 4, units = "cm", dpi = 300)

#############################
# Fig. ADDED
#############################
# "Early-converter"
df_subset <- dplyr::filter(df, response=="Early-converter")
wt<-wilcox.test(PP~time,data=df_subset,paired=TRUE)
p1<-ggplot(df_subset,aes(y=PP,x=time)) + geom_boxplot() + stat_summary(aes(group=sample,color=response),geom="line",fun=sum,alpha = 0.3) + geom_point(stat='summary', fun=sum, aes(group=sample,color=response))
p1<-p1+scale_color_manual(values=c("#00AFBB"))
p1<-p1+ ylim(c(0, 8e-04)) + scale_x_discrete(labels=c("T0" = "0", "T60" = "60"))
p1<-p1 + ylab('% Unique HBsAg-specific TCRβ sequence') + xlab('Time Point (days)')
p1<-p1 + labs(title="",subtitle=paste0("Paired Wilcox P-value = ",format(wt$p.value,digits=2)))
p1<-p1 + theme_bw() + theme(aspect.ratio = 1.5, 
                          legend.title=element_blank(), 
                          legend.position = "none",
                          plot.title = element_text(size = 20, face = "bold"),
                          axis.text.x  = element_text(size = 10),
                          axis.title.x  = element_text(size = 10, face = "bold"), 
                          axis.title.y = element_text(size = 10, face = "bold"))
p1

df_subset <- dplyr::filter(df, response=="Late-converter")
wt<-wilcox.test(PP~time,data=df_subset,paired=TRUE)
p2<-ggplot(df_subset,aes(y=PP,x=time)) + geom_boxplot() + stat_summary(aes(group=sample,color=response),geom="line",fun=sum,alpha = 0.3) + geom_point(stat='summary', fun=sum, aes(group=sample,color=response))
p2<-p2+scale_color_manual(values=c("#E7B800"))
p2<-p2+ ylim(c(0, 8e-04)) + scale_x_discrete(labels=c("T0" = "0", "T60" = "60"))
p2<-p2 + ylab('% Unique HBsAg-specific TCRβ sequence') + xlab('Time Point (days)') + theme_bw() + theme(aspect.ratio = 1.5, legend.position = "none", legend.title = element_blank()) 
p2<-p2 + labs(title="",subtitle=paste0("Paired Wilcox P-value = ",format(wt$p.value,digits=2)))
p2<-p2 + theme_bw() + theme(aspect.ratio = 1.5, 
                            legend.title=element_blank(), 
                            legend.position = "none",
                            plot.title = element_text(size = 20, face = "bold"),
                            axis.text.x  = element_text(size = 10),
                            axis.title.x  = element_text(size = 10, face = "bold"), 
                            axis.title.y = element_text(size = 10, face = "bold"))
p2


df_subset <- dplyr::filter(df, response=="Non-converter")
wt<-wilcox.test(PP~time,data=df_subset,paired=TRUE)
p3<-ggplot(df_subset,aes(y=PP,x=time)) + geom_boxplot() + stat_summary(aes(group=sample,color=response),geom="line",fun=sum,alpha = 0.3) + geom_point(stat='summary', fun=sum, aes(group=sample,color=response))
p3<-p3+scale_color_manual(values=c("#FC4E07"))
p3<-p3+ ylim(c(0, 8e-04)) + scale_x_discrete(labels=c("T0" = "0", "T60" = "60"))
p3<-p3 + ylab('% Unique HBsAg-specific TCRβ sequence') + xlab('Time Point (days)') + theme_bw() + theme(aspect.ratio = 1.5, legend.position = "none", legend.title = element_blank()) 
p3<-p3 + labs(title="",subtitle=paste0("Paired Wilcox P-value = ",format(wt$p.value,digits=2)))
p3<-p3 + theme_bw() + theme(aspect.ratio = 1.5, 
                            legend.title=element_blank(), 
                            legend.position = "none",
                            plot.title = element_text(size = 20, face = "bold"),
                            axis.text.x  = element_text(size = 10),
                            axis.title.x  = element_text(size = 10, face = "bold"), 
                            axis.title.y = element_text(size = 10, face = "bold"))
p3

p <- grid.arrange(p1, p2, p3, nrow = 1)
# ggsave(file.path(main_fig, 'Fig_ADDED_tcrincrease_per_group.png'),p, scale = 1, height = 10, width = 22, units = "cm", dpi = 300)

#############################
# Fig. 2d
# Frequency of vaccine-specific TCRβ sequences within memory CD4 T cell repertoire at time point 60 for each vaccinee, 
# normalized by amount of HBsAg-specific TCRβ found for each vaccinee.
#############################

# HBV-associated TCR vounts are stored iniPPB60 (for day 60)
# Normalize for number of HBV-associated TCRs found in peptide pool experiments

p<-ggplot(rep,aes(y=iPPB60/PP,x=responder, color=responder)) + geom_boxplot() + geom_point(size = 3, alpha = 0.7)
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")))
p<-p + ggtitle("Day 60") + xlab("") + ylab("% of Normalized HBsAg-specific TCRβ")
p<-p+ylim(0,0.7)
p<-p + theme_bw() + theme(legend.title=element_blank(), 
                          plot.title = element_text(size = 20, face = "bold"),
                          axis.text.x  = element_blank(), 
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 10, face = "bold"),
                          legend.position = "none", aspect.ratio = 1.5)
p
ggsave(file.path(main_fig, 'Fig_2d_ppnorm_day60.png'),p, height = 10, width = 7, units = "cm", dpi = 300)


# Same for Day 0
p<-ggplot(rep,aes(y=iPPB0/PP,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")))
p<-p + ggtitle("Day 0") + xlab("") + ylab("% of Normalized HBsAg-specific TCRβ")
p<-p+ylim(0,0.7)
p<-p + theme_bw() + theme(legend.title=element_blank(), 
                          plot.title = element_text(size = 20, face = "bold"),
                          axis.text.x  = element_blank(), 
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 10, face = "bold"),
                          legend.position = "none", aspect.ratio = 1.5)
p
#ggsave(file.path(main_fig, 'Fig_2d_ppnorm_day0.png') ,p, height = 10, width = 7, units = "cm", dpi = 300)

#############################
# Fig. 3b
# Scatter plot with the percentage predicted epitope-specific and bystander TCRβ sequences. Predictions done as a leave-one-out cross-validation. Each dot represents a vaccinee with the color denoting the responder status (blue: early-converter, yellow: late-converter, red: non-converter)
#############################

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB0/B0,x=PPnrB0/B0, color=responder)) + geom_point(size = 3, alpha = 0.7)
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + xlab("% Predicted bystander TCRβ") + ylab("% Predicted HBsAg-specific TCRβ") + ggtitle("TCRβ data from day 0")
p<-p +  xlim(3.5e-07, 1.25e-06) + ylim(3.8e-07, 6.6e-07)
p<-p + theme_bw()  + theme(legend.position = "none", 
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 12, face = "bold"),
                           axis.title.y = element_text(size = 12, face = "bold"),
                           aspect.ratio = 1) 
p
ggsave(file.path(main_fig, 'Fig_3b_tcrspecific_bystander_day0.png'),p, height = 8, width = 8, units = "cm", dpi = 300)

# same but for day 60
p<-ggplot(rep,aes(y=PSB60/B60,x=PPnrB60/B60, color=responder)) + geom_point(size = 3, alpha = 0.7)
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + xlab("% Predicted bystander TCRβ") + ylab("% Predicted HBsAg-specific TCRβ") + ggtitle("TCRβ data from day 60")
p<-p +  xlim(3.5e-07, 1.25e-06) + ylim(3.8e-07, 6.6e-07)
p<-p + theme_bw()  + theme(legend.position = "none", 
                           plot.title   = element_text(size = 15, face = "bold"),
                           axis.title.x = element_text(size = 12, face = "bold"),
                           axis.title.y = element_text(size = 12, face = "bold"),
                           aspect.ratio = 1) 

p
#ggsave(file.path(main_fig, 'Fig_3b_tcrspecific_bystander_day60.png'),p, height = 8, width = 8, units = "cm", dpi = 300)

#############################
# Fig. 3c
#############################

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB60/PPnrB60,x=responder, color=responder)) + geom_boxplot() + geom_point(size = 5, alpha = 0.7)
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
p<-p + xlab("") + ylab("HBsAg-predictive ratio R(hbs)") + ggtitle("Day 60")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            plot.title = element_text(size = 20, face = "bold"),
                            axis.text.x  = element_blank(), 
                            axis.title.x = element_blank(),
                            axis.text.y  = element_text(size = 14),
                            axis.title.y = element_text(size = 18, face = "bold"),
                            legend.position = "none", aspect.ratio = 1.5)
p
ggsave(file.path(main_fig, 'Fig_3c_tcrspecific_day60_bystandernorm.png'),p, height = 10, width = 8, units = "cm", dpi = 300)

#############################
# Fig. 3d
#############################

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB0/PPnrB0,x=responder, color=responder)) + geom_boxplot() + geom_point(size = 5, alpha = 0.7)
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
p<-p + xlab("") + ylab("HBsAg-predictive ratio R(hbs)") + ggtitle("Day 0")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            plot.title = element_text(size = 20, face = "bold"),
                            axis.text.x  = element_blank(), 
                            axis.title.x = element_blank(),
                            axis.text.y  = element_text(size = 14),
                            axis.title.y = element_text(size = 18, face = "bold"),
                            legend.position = "none", aspect.ratio = 1.5)
p
ggsave(file.path(main_fig, 'Fig_3c_tcrspecific_day0_bystandernorm.png'),p, height = 10, width = 8, units = "cm", dpi = 300)

#############################
# Fig. 3e
#############################

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires
# ROC plot early vs all

rep$PSB0divnr = rep$PSB0/rep$PPnrB0
r.psb0divnr<-roc(responder ~ PSB0divnr, rep[!rep$responder=='Non-converter',], plot=TRUE, grid=TRUE,print.auc=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE, conf.level=0.95, levels=c('Late-converter','Early-converter'))
p<-ggroc(r.psb0divnr, colour = 'steelblue', size = 2) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
p<-p+xlab("Specificity") + ylab("Sensitivity")
p<-p + labs(title="ROC for HBsAg-predictive ratio R(hbs) at day 0",subtitle=paste0("AUC = ",format(r.psb0divnr$ci,digits=2)[2], " (95% CI: ",format(r.psb0divnr$ci,digits=2)[1],"-",format(r.psb0divnr$ci,digits=2)[3],")"))
p<-p + theme_bw() +   theme(plot.title = element_text(size = 12, face = "bold"),
                            plot.subtitle = element_text(size = 10),
                            axis.title.x = element_text(size = 20, face = "bold"),
                            axis.title.y = element_text(size = 20, face = "bold"),
                            legend.position = "none", aspect.ratio = 1)
p
ggsave(file.path(main_fig, 'Fig_3e_tcrspecific_day0_bystandernorm_ROC.png'), p, height = 13, width = 12, units = "cm", dpi = 300)

#############################
# Fig. 3f
#############################
# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires
# ROC plot CD154 vs none

rep$PSB0divnr = rep$PSB0/rep$PPnrB0
r.cd154b0divnr<-roc(CD154>0 ~ PSB0divnr, rep, plot=TRUE, grid=TRUE,print.auc=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE, conf.level=0.95)
p<-ggroc(r.cd154b0divnr, colour = 'steelblue', size = 2) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
p<-p+xlab("Specificity")+ylab("Sensitivity")
p<-p + labs(title="CD40L ROC for HBsAg-predictive ratio R(hbs) at day 0",subtitle=paste0("AUC = ",format(r.cd154b0divnr$ci,digits=2)[2], " (95% CI: ",format(r.cd154b0divnr$ci,digits=2)[1],"-",format(r.cd154b0divnr$ci,digits=2)[3],")"))
p<-p + theme_bw() +   theme(plot.title = element_text(size = 12, face = "bold"),
                            plot.subtitle = element_text(size = 10),
                            axis.title.x = element_text(size = 20, face = "bold"),
                            axis.title.y = element_text(size = 20, face = "bold"),
                            legend.position = "none", aspect.ratio = 1)
p
ggsave(file.path(main_fig, 'Fig_3f_cd40l_day0_bystandernorm_ROC.png'), p, height = 15, width = 12, units = "cm", dpi = 300)
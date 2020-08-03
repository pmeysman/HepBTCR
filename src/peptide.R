require('ggplot2')
require('ggpubr')
require('pROC')

# Full repertoire analyses with peptide-specific TCRs

rep <- read.csv(file="../results/peptide/loocvrep_ab.txt")

# HBV-associated TCR vounts are stored in iPPB0 and iPPB60 (for day 0 and day 60)
# Normalize for the size of the repertoire (B0 and B60)

df <- data.frame(PP = c(rep$iPPB0/rep$B0,rep$iPPB60/rep$B60), response = rep$responder, time = c(rep('T0',33),rep('T60',33)),sample=rep(rownames(rep),2))
wt<-wilcox.test(PP~time,data=df,paired=TRUE)
p<-ggplot(df,aes(y=PP,x=time)) + geom_boxplot() + stat_summary(aes(group=sample,color=response),geom="line",fun.y=sum,alpha = 0.3) + geom_point(stat='summary', fun.y=sum, aes(group=sample,color=response))
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p+ scale_x_discrete(labels=c("T0" = "Day0", "T60" = "Day60"))
p<-p + ylab('HBV-associated TCR%') + xlab('Time point') + theme_bw() + theme(legend.title = element_blank()) 
p<-p + labs(title="HBV TCR increase after vaccination",subtitle=paste0("Paired Wilcox P-value = ",format(wt$p.value,digits=2)))
p
ggsave('../results/peptide/tcrincrease.pdf',p)

# HBV-associated TCR vounts are stored iniPPB60 (for day 60)
# Normalize for number of HBV-associated TCRs found in peptide pool experiments

p<-ggplot(rep,aes(y=iPPB60/PP,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")))
p<-p + xlab("") + ylab("Normalized HBV-associated TCR%") + ggtitle("HBV TCR frequency at day 60")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
                            axis.title.x = element_blank(),
                            legend.position = "none", aspect.ratio = 2)
p
ggsave('../results/peptide/ppnorm_day60.pdf',p)

# Same for Day 0

p<-ggplot(rep,aes(y=iPPB0/PP,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")))
p<-p + xlab("") + ylab("Normalized HBV-associated TCR%") + ggtitle("HBV TCR frequency at day 0")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
                            axis.title.x = element_blank(),
                            legend.position = "none", aspect.ratio = 2)
p
ggsave('../results/peptide/ppnorm_day0.pdf',p)

# HBV-specific prediction (epitope specific)
# Same as before, normalize by repertoire size

p<-ggplot(rep,aes(y=PSB0/B0,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + xlab("") + ylab("Predicted HBV-specific TCR%") + ggtitle("Cross-validated HBV-specific TCR frequency at day 0")
p<-p + theme_bw() + theme(legend.position = "none")
p
ggsave('../results/peptide/tcrspecific_day0.pdf',p)

# HBV-specific prediction (epitope specific)
# Same as before, normalize by repertoire size
# ROC plot early vs all

rep$PSB0divB0 = rep$PSB0/rep$B0
r.psb0divb0<-roc(responder ~ PSB0divB0, rep[!rep$responder=='Non-converter',], plot=TRUE, grid=TRUE,print.auc=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE, conf.level=0.95, levels=c('Late-converter','Early-converter'))
p<-ggroc(r.psb0divb0) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
p<-p+xlab("Specificity")+ylab("Sensitivity")+theme_bw()
p<-p + labs(title="ROC for cross-validated HBV-specific TCR frequency at day 0",subtitle=paste0("AUC = ",format(r.psb0divb0$ci,digits=2)[2], " (95% CI: ",format(r.psb0divb0$ci,digits=2)[1],"-",format(r.psb0divb0$ci,digits=2)[3],")"))
p
ggsave('../results/peptide/tcrspecific_day0_ROC.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB60/PPnrB60,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + xlab("") + ylab("HBsAg-predictive ratio R(hbs)") + ggtitle("TCR HBsAg predictions at day 60")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
                            axis.title.x = element_blank(),
                            legend.position = "none", aspect.ratio = 2)
p
ggsave('../results/peptide/tcrspecific_day60_bystandernorm.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB0/PPnrB0,x=responder, color=responder)) + geom_boxplot() + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + xlab("") + ylab("HBsAg-predictive ratio R(hbs)") + ggtitle("TCR HBsAg predictions at day 0")
p<-p + theme_bw() +   theme(legend.title=element_blank(), 
                            axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
                            axis.title.x = element_blank(),
                            legend.position = "none", aspect.ratio = 2)
p
ggsave('../results/peptide/tcrspecific_day0_bystandernorm.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB60/B60,x=PPnrB60/B60, color=responder)) + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + xlab("Predicted bystander TCR%") + ylab("Predicted HBsAg-specific TCR%") + ggtitle("Cross-validated HBsAg-specific / Bystander TCR frequency at day 60")
p<-p + theme_bw()  + theme(legend.title = element_blank()) 
p
ggsave('../results/peptide/tcrspecific_bystander_day60.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires

p<-ggplot(rep,aes(y=PSB0/B0,x=PPnrB0/B0, color=responder)) + geom_point()
p<-p+scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07"))
p<-p + xlab("Predicted bystander TCR%") + ylab("Predicted HBsAg-specific TCR%") + ggtitle("Cross-validated HBV-specific / Bystander TCR frequency at day 0")
p<-p + theme_bw()  + theme(legend.title = element_blank()) 
p
ggsave('../results/peptide/tcrspecific_bystander_day0.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires
# ROC plot early vs all

rep$PSB0divnr = rep$PSB0/rep$PPnrB0
r.psb0divnr<-roc(responder ~ PSB0divnr, rep[!rep$responder=='Non-converter',], plot=TRUE, grid=TRUE,print.auc=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE, conf.level=0.95, levels=c('Late-converter','Early-converter'))
p<-ggroc(r.psb0divnr) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
p<-p+xlab("Specificity")+ylab("Sensitivity")+theme_bw()
p<-p + labs(title="ROC for HBsAg-predictive ratio R(hbs) at day 0",subtitle=paste0("AUC = ",format(r.psb0divnr$ci,digits=2)[2], " (95% CI: ",format(r.psb0divnr$ci,digits=2)[1],"-",format(r.psb0divnr$ci,digits=2)[3],")"))
p
ggsave('../results/peptide/tcrspecific_day0_bystandernorm_ROC.pdf',p)

# CD154 comparison

rep[is.na(rep$CD154),'CD154'] <- 0

p<-ggboxplot(rep,y="CD154",x="responder", color="responder",palette =c("#00AFBB", "#E7B800", "#FC4E07", "#28fc07")) + geom_point(aes(color=responder))
p<-p + xlab("") + ylab("CD40L assay") + ggtitle("CD40L assay for Ab-defined categories")
p<-p + stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),comparisons=list(c("Early-converter", "Late-converter"),c("Early-converter", "Non-converter"),c("Late-converter", "Non-converter")),method.args = list(alternative = "greater"))
p<-p + theme_bw() + theme(legend.title=element_blank(), 
                           axis.text.x  = element_text(size = 12, angle = 45, hjust = 1), 
                           axis.title.x = element_blank(),
                           legend.position = "none", aspect.ratio = 2)
p
ggsave('../results/peptide/cd40l_ab_comparison.pdf',p)

# HBV-specific prediction (epitope specific)
# Normalize with bystander repertoires
# ROC plot CD154 vs none

rep$PSB0divnr = rep$PSB0/rep$PPnrB0
r.cd154b0divnr<-roc(CD154>0 ~ PSB0divnr, rep, plot=TRUE, grid=TRUE,print.auc=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE, conf.level=0.95)
p<-ggroc(r.cd154b0divnr) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
p<-p+xlab("Specificity")+ylab("Sensitivity")+theme_bw()
p<-p + labs(title="CD40L ROC for HBsAg-predictive ratio R(hbs) at day 0",subtitle=paste0("AUC = ",format(r.cd154b0divnr$ci,digits=2)[2], " (95% CI: ",format(r.cd154b0divnr$ci,digits=2)[1],"-",format(r.cd154b0divnr$ci,digits=2)[3],")"))
p
ggsave('../results/peptide/cd40l_day0_bystandernorm_ROC.pdf',p)



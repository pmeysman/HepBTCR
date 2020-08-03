
#Get TCR data
tcr <- read.table('../results/timescatter/volunteercounts.txt',header=TRUE,sep=',',row.names='volunteer')

#Make plots
require('ggplot2')
require('ggpubr')
require('gridExtra')

time<-data.frame(volunteer=rep(row.names(tcr),2),
                 entropy=c(tcr$B0entropy,tcr$B60entropy),
                 total=c(tcr$B0total,tcr$B60total),
                 reads=c(tcr$B0reads,tcr$B60reads),
                 PPIntfreq=c(tcr$PPIntB0freq,tcr$PPIntB60freq),
                 day=c(rep('Day0',33),rep('Day60',33)),
                 responder=rep(tcr$responder,2))

p<-ggplot(tcr,aes(x=B0total,y=B60total,shape=responder)) + geom_point() + theme_classic()
p<-p + xlab("Day 0 unique TCRs") + ylab("Day 60 unique TCRs")
p
ggsave('../results/uniqueTCRs.pdf',p)

p<-ggplot(tcr,aes(x=B0reads,y=B60reads,label=row.names(tcr))) + geom_text() + theme_classic()
p<-p + xlab("Day 0 TCRs reads") + ylab("Day 60 TCRs reads")
p
ggsave('../results/timescatter/TCRreads.pdf',p)

p<-ggboxplot(time,y="PPIntfreq",x="day", color="day")
p<-p + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
p<-p+geom_line(aes(group = volunteer),alpha = 0.1) + theme(legend.position = "none") 
p<-p+ylab("HBV-specific T-cell frequency")+xlab("")
p
ggsave('../results/timescatter/PPIntfreq.pdf',p)

pte<-ggboxplot(time[time['responder']=='Early-converter',],y="total",x="day", color="day")
pte<-pte + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pte<-pte +geom_line(aes(group = volunteer),alpha = 0.1)
pte<-pte + xlab("") + ylab("Unique CDR3s") + ggtitle("Early-converter") + theme(plot.title = element_text(hjust = 0.5))
pte<-pte + theme(legend.position = "none") 

ptl<-ggboxplot(time[time['responder']=='Late-converter',],y="total",x="day", color="day")
ptl<-ptl + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
ptl<-ptl +geom_line(aes(group = volunteer),alpha = 0.1)
ptl<-ptl + xlab("") + ylab("")  + ggtitle("Late-converter") + theme(plot.title = element_text(hjust = 0.5))
ptl<-ptl + theme(legend.position = "none") 

pee<-ggboxplot(time[time['responder']=='Early-converter',],y="entropy",x="day", color="day")
pee<-pee + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pee<-pee +geom_line(aes(group = volunteer),alpha = 0.1)
pee<-pee + xlab("") + ylab("Entropy") + ylim(6,8)
pee<-pee + theme(legend.position = "none") 

pel<-ggboxplot(time[time['responder']=='Late-converter',],y="entropy",x="day", color="day")
pel<-pel + stat_compare_means(method = "wilcox.test", paired = TRUE, aes(label = ..p.signif..),label.x = 1.5)
pel<-pel +geom_line(aes(group = volunteer),alpha = 0.1)
pel<-pel + xlab("") + ylab("") +ylim(6,8)
pel<-pel + theme(legend.position = "none") 

etplot<-grid.arrange(pte, ptl, pee, pel, nrow = 2, ncol= 2)
ggsave('../results/entropy-total.pdf',etplot)

wilcox.test(tcr$B0entropy,tcr$B60entropy,paired=TRUE)
wilcox.test(tcr$PPIntB0freq,tcr$PPIntB60freq,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$B0entropy,tcr[tcr['responder']=='Late-converter',]$B60entropy,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$B0total,tcr[tcr['responder']=='Late-converter',]$B60total,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Late-converter',]$PPIntB0freq,tcr[tcr['responder']=='Late-converter',]$PPIntB60freq,paired=TRUE)


wilcox.test(tcr[tcr['responder']=='Early-converter',]$B0entropy,tcr[tcr['responder']=='Early-converter',]$B60entropy,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Early-converter',]$B0total,tcr[tcr['responder']=='Early-converter',]$B60total,paired=TRUE)
wilcox.test(tcr[tcr['responder']=='Early-converter',]$PPIntB0freq,tcr[tcr['responder']=='Early-converter',]$PPIntB60freq,paired=TRUE)


# 


# Supplemental figure

supp<-data.frame(vaccinee = rep(tcr$X,3))
supp$TCRcounts <- c(tcr$B0B60overlap, (tcr$B0total - tcr$B0B60overlap), ( tcr$B60total - tcr$B0B60overlap))
supp$label <- c(rep("Overlap",33),rep("Day0 Unique",33),rep("Day60 Unique",33))

s.unique<-ggplot(supp,aes(x=vaccinee,y=TCRcounts,fill=label)) + geom_col() 
s.unique<-s.unique+ ylab("Unique TCRs") + theme_classic() + theme(legend.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Vaccinees")
s.unique
ggsave('../results/s_unique.pdf',s.unique)

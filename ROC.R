#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
wd=args[1]
print(wd)
setwd(wd)

library(ggplot2)
#library(dplyr)


# - - - - - - - - - - - - - - - ROC CURVES - - - - - - - - - - - - - - - - - - - - - - - 
prep4rocs<-function(DF) {
  p<-rep(1, nrow(DF))
  n<-rep(0, nrow(DF))
  lab<-c(p,n)
  sco<-c(DF$V1, DF$V2)
  newDF<-data.frame(lab, sco)
  return(newDF)
}

# from https://blog.revolutionanalytics.com/2016/08/roc-curves-in-two-lines-of-code.html
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

plotRocs<-function(PWM, AUCpwm){

df1<-read.table(PWM, header=FALSE)
d1<-prep4rocs(df1)
#sr1<-simple_roc(d1$lab, d1$sco) %>% mutate(Model=paste("PWM (", AUCpwm, ")"))
sr1<-simple_roc(d1$lab, d1$sco)
print(head(sr1))


SR<-rbind(sr1)
#print(head(SR))
p<-ggplot(data = SR, 
       aes(x = FPR, y = TPR, color = "red")) +
    #geom_point() + 
	geom_line(aes(color="red"), lwd=.6)
	#scale_linetype_manual(values=c("twodash", "dotted", "solid")) +
	#theme(axis.text=element_text(size=TS), axis.title=element_text(size=TS), legend.position = c(0.56, 0.27), legend.text=element_text(size=8),
	#		legend.title=element_text(size=10), plot.title = element_text(size=TS, hjust = 0.5))
return(p)
}

roc<-plotRocs("scores/tab_pfm.tsv", 0.05)


pdf("Figure.pdf", width=12, height=10)
roc
dev.off()


#1system("convert -density 144 -background white Figure4.pdf Figure4.png")

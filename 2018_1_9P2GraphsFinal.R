

setwd("/Users/cansav091/Desktop/Current Projects /Western Blots/P2E4AdultDieldrinWesterns")
std.e <- function(x) sd(x)/sqrt(length(x))
dat=read.csv("AllTargetsCleanData.csv")
dat=dat[-which(dat$Treatment=="Standard"),]
badsamples=c("P2-E1-027")
dat=dat[is.na(match(dat$Sample,badsamples)),]

target=c("TH","DAT","VMAT2")

for(ii in 1:length(target)){

### Make a boxplot using the regression estimated values
  png(paste0("P2BarPlot",target[ii],".png"))
  par(mar=c(5,7,5,5))
  #### Do a ANOVA; make a boxplot d
  xx=which(dat$Target==unique(dat$Target)[ii])### Only graph the data for a particular target
  target.data=dat[xx,]
  tx.groups=factor(as.character(target.data$Treatment),unique(as.character(target.data$Treatment)))
  avg=tapply(target.data$EstimatVals.RelRevert,tx.groups,mean)
  bars=barplot(avg,names=c("0","0.3","3"),cex.main=3,cex.names=1.5,cex.lab=1.5,cex.axis = 1.5,ylim=c(0,20),xlab="Dieldrin (mg/kg)",ylab=paste0("Striatal ",target[ii],"\n (Relative Units)"),col="black")
  axis(side=1,at=c(0,.70,1.9,3.1,4),lwd=3,labels=FALSE)
  axis(side=2,at=c(0,5,10,15,20),lwd=3,labels=FALSE)
  
  xx=which(tx.groups=="Control")
  yy=which(tx.groups=="Low")
  zz=which(tx.groups=="High")
  
  ses=c(std.e(target.data$EstimatVals.RelRevert[xx]),std.e(target.data$EstimatVals.RelRevert[yy]),std.e(target.data$EstimatVals.RelRevert[zz]))
  
  
  segments(bars, avg - ses * 2, bars,  avg + ses * 2, lwd = 3)
  arrows(bars, avg - ses * 2, bars,  avg + ses * 2, lwd = 3, angle = 90, code = 3, length = 0.2)
  dev.off()
  
  }

getwd()



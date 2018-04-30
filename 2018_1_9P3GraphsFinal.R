

setwd("/Users/cansav091/Desktop/Current Projects /Western Blots/P3E4DevDieldrinWesterns/")
std.e <- function(x) sd(x)/sqrt(length(x))
files=grep("CleanData",dir(),value=TRUE)

dat=rbind(read.csv(files[1]),read.csv(files[2]))
dat=cbind(dat,as.factor(c(rep("Female",54),rep("Male",54))))
colnames(dat)[16]="Sex"
dat=dat[-which(dat$Treatment=="Standard"),]
#badsamples=c("P2-E1-027")
#dat=dat[is.na(match(dat$Sample,badsamples)),]

target=c("TH","DAT","VMAT2")

aovs=c()
rev.aovs=c()
t.tests=matrix(ncol=5)
colnames(t.tests)=c("df","t stat","pval","ControlAvg","DieldrinAvg")
sexes=c("Female","Male")
rnames=c()

for(ii in 1:length(target)){
  for(jj in 1:length(sexes)){
    rnames=c(rnames,paste0(sexes[jj],target[ii]))
    ### Make a boxplot using the regression estimated values
    png(paste0("P3BarPlot",target[ii],sexes[jj],".png"))
    par(mar=c(7,7,5,5))
    #### Do a ANOVA; make a boxplot d
    xx=which(dat$Target==unique(dat$Target)[ii])### Only graph the data for a particular target
    target.data=dat[xx,]
    sex=which(target.data$Sex==sexes[jj])
    target.data=target.data[sex,]

    tx.groups=factor(as.character(target.data$Treatment),unique(as.character(target.data$Treatment)))
    xx=which(tx.groups=="Control")
    ttest=t.test(target.data$EstimatVals.RelRevert[xx],target.data$EstimatVals.RelRevert[-xx]) 
    t.tests=rbind(t.tests,c(ttest$parameter,ttest$statistic, ttest$p.value,ttest$estimate))

    avg=tapply(target.data$EstimatVals.RelRevert,tx.groups,mean)
    
    bars=barplot(avg,names=c("0","0.3"),mgp=c(2.5,1,0),cex.main=3,border=NA,cex.names=1.5,cex.lab=1.25,cex.sub=1,cex.axis = 1.5,cex.lab=1.2,ylim=c(0,20),xlab="Dieldrin (mg/kg)",ylab=paste0("Striatal ",target[ii],"\n (Relative Units)"),col="black")
    axis(side=1,at=c(0,.70,1.9,3.1,4.3,5),lwd=3,labels=FALSE)
    axis(side=2,at=c(0,5,10,15,20),lwd=3,labels=FALSE)
    

    ses=c(std.e(target.data$EstimatVals.RelRevert[xx]),std.e(target.data$EstimatVals.RelRevert[-xx]))

    segments(bars, avg - ses * 2, bars,  avg + ses * 2, lwd = 3)
    arrows(bars, avg - ses * 2, bars,  avg + ses * 2, lwd = 3, angle = 90, code = 3, length = 0.2)
    dev.off()
    
}
}
t.tests=t.tests[-1,]
rownames(t.tests)=rnames
write.csv(t.tests,"t-test results P3.csv")








rad.date="7/19/17"
half.life=4537


rad.date=as.Date(rad.date,"%m/%d/%y")
rad.date=abs(as.numeric(rad.date-Sys.Date()))
percent=.5^(rad.date/half.life)

micromolar=(percent*1/.55)*1000

######## Analyze Data 
home="/Users/cansav091/Desktop/Current Projects /UptakeAssays"
setwd(home)

files=select.list(dir(),multiple=TRUE, title="Which output file(s) do you want to quantify? Choose 0 if you want the most recently added output file.")
files=dir()[-7]
alldata=matrix(ncol=10,nrow=0)
rows=c()
cols=c()

for(ii in 1:length(files)){
read=suppressWarnings(grep("COLUMNS",readLines(files[ii])))
end=suppressWarnings(grep("PLATES:",readLines(files[ii])))
data=readLines(files[ii],n=end)
xx=which(data=="========")+1
data=data[xx:(length(data)-2)]

data=unlist(strsplit(data,"\t"))
data=data[-which(data==" ")]
data=data[!is.na(data)]

matrix=matrix(ncol=3,nrow=(length(data)/4)-1)
matrix[,1]=data[seq(from=6,to=length(data),by=4)]
matrix[,2]=data[seq(from=7,to=length(data),by=4)]
matrix[,3]=data[seq(from=8,to=length(data),by=4)]
matrix=cbind(matrix,c(rep(paste0("Plate",ii),nrow(matrix))))

rows=c(rows,data[seq(from=5,to=length(data),by=4)])
cols=data[2:10]

rowsplate=substr(matrix[,1],1,1)
colsplate=substr(matrix[,1],2,4)

row.avg=tapply(as.numeric(matrix[,2]),rowsplate,mean)
col.avg=tapply(as.numeric(matrix[-which(rowsplate=="D"),2]),colsplate[-which(rowsplate=="D")],mean)

row.sd=tapply(as.numeric(matrix[,2]),rowsplate,sd)
col.sd=tapply(as.numeric(matrix[-which(rowsplate=="D"),2]),colsplate[-which(rowsplate=="D")],sd)


matrix=cbind(matrix,rep(row.avg,summary(as.factor(rowsplate))))
matrix=cbind(matrix,rep(col.avg,summary(as.factor(colsplate))))

matrix=cbind(matrix,rep(row.sd,summary(as.factor(rowsplate))))
matrix=cbind(matrix,rep(col.sd,summary(as.factor(colsplate))))

matrix.stat=apply(matrix[,c(2,5:8)],2,as.numeric)

rowzscore=(matrix.stat[,1]-matrix.stat[,2])/matrix.stat[,4]
colzscore=(matrix.stat[,1]-matrix.stat[,3])/matrix.stat[,5]

matrix=cbind(matrix,rowzscore,colzscore)
alldata=rbind(alldata,matrix)
}

alldata=cbind(substr(alldata[,1],1,1),substr(alldata[,1],2,4),alldata[,2:10])
colnames(alldata)=c("Row","Col",cols[2:10],"Plate")
rownames(alldata)=rows

alldata[(nrow(alldata)-23):nrow(alldata),5]=rep("EndPlate",24)
TBZ=which(alldata[,1]=="D")

CCPMs=as.numeric(alldata[,3])

xx=which(alldata[,5]=="Plate1")
CCPMs.avg=tapply(CCPMs[-xx],as.factor(rownames(alldata))[-xx],mean)
xx=which(alldata[,5]=="Plate2")
CCPMs.avg=c(CCPMs.avg,tapply(CCPMs[-xx],as.factor(rownames(alldata))[-xx],mean))

xx=which(alldata[,5]=="EndPlate")
std.curve=rep(rep(c(30,10,3,1,.3,.1),1,each=4),2)
plot(CCPMs[-xx],log(std.curve))

plot(CCPMs.avg,log(c(30,10,3,1,.3,.1)))
#convert CCPMs to fmol/mg protein/min
log(std.curve)

DA=CCPMs[which(rownames(alldata)=="Unk_21")]


write.csv(alldata,file="CleanedUptakeData.csv")




wells=paste0(alldata[,1],alldata[,2])
alldata[order(wells),]
corr()



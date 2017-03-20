
# Plot frequency of most common base in CpG sites

library(data.table)
mydata<-fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_A01.1.txt")
mydata[, total := nA + nC + nG + nT]

hist(100*mydata[ref=="C"]$nC/mydata[ref=="C"]$total, main="%C in ref=C", xlab="%C")
hist(100*mydata[ref=="C"]$nC/mydata[ref=="C"]$total, main="%C in ref=C", xlab="%C", ylim=c(0,10000), breaks=100)
hist(100*mydata[ref=="G"]$nG/mydata[ref=="G"]$total, main="%G in ref=G", xlab="%G")
hist(100*mydata[ref=="G"]$nG/mydata[ref=="G"]$total, main="%G in ref=G", xlab="%G", ylim=c(0,10000), breaks=100)



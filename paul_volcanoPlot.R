res=read.table("atac.diffreps.annotated", sep="\t",header=T)

#down:23466 up:15360 

pdf("pie_Open_in_NT.pdf")
tdown=table(res$Feature[res$Event=="Down"])
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
pie(tdown, main="Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)")
dev.off()

pdf("pie_Open_in_SH.pdf")
tup=table(res$Feature[res$Event=="Up"])
names(tup) = paste(names(tup)," ",round(tup/sum(tup)*100,digits=2),"%",sep="")
pie(tup, main="Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)")
dev.off()

#########
#annotatePeaks.pl Down_NT.bed hg19 -annStats Down_NT.annStats > Down_NT.anno
#annotatePeaks.pl Up_SH.bed hg19 -annStats Up_SH.annStats > Up_SH.anno

pdf("pie_homer_Open_in_NT.pdf")
par(mar=c(5.1,4.1,4.1,7))
res=read.table(pipe("more Down_NT.annStats |grep -v '0.0'|cut -f1,2,4|tail -n +10"), sep="\t",header=F)
tdown = res[,2]
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>16]
pie(tdown, main="Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)")
dev.off()

pdf("pie_homer_Open_in_SH.pdf")
res=read.table(pipe("more Up_SH.annStats |grep -v '0.0'|cut -f1,2,4|tail -n +10"), sep="\t",header=F)
tup = res[,2]
names(tup) = res[,1]
names(tup) = paste(names(tup)," ",round(tup/sum(tup)*100,digits=2),"%",sep="")
tup = tup[tup>16]
pie(tup, main="Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)")
dev.off()

#####
pdf("enrichment_homer_in_NT.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Down_NT.annStats |grep -v '0.0'|cut -f1,2,4|tail -n +10"), sep="\t",header=F)
tdown = res[,3]
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[res[,2]>16]
barplot(sort(tdown),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
dev.off()

pdf("enrichment_homer_in_SH.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Up_SH.annStats |grep -v '0.0'|cut -f1,2,4|tail -n +10"), sep="\t",header=F)
tup = res[,3]
names(tup) = res[,1]
names(tup) = paste(names(tup)," ",round(tup/sum(tup)*100,digits=2),"%",sep="")
tup = tup[res[,2]>16]
barplot(sort(tup),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
dev.off()

###########################################################################################
################################################################################################################
pdf("pie_homer_Open_in_NT.pdf")
res=read.table(pipe("more Down_NT.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8)
title("Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)", cex.main=.9)
dev.off()

pdf("pie_homer_Open_in_SH.pdf")
res=read.table(pipe("more Up_SH.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8)
title("Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)", cex.main=.9)
dev.off()

####################################################################################
#####
pdf("enrichment_homer_in_NT.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Down_NT.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")

enrich = as.numeric(as.character(res[tdown>200,3]))
tdown = tdown[tdown>200]
names(enrich) = names(tdown)
barplot(sort(enrich),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
abline(h=0)
dev.off()

pdf("enrichment_homer_in_SH.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Up_SH.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")

enrich = as.numeric(as.character(res[tdown>200,3]))
tdown = tdown[tdown>200]
names(enrich) = names(tdown)
barplot(sort(enrich),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
abline(h=0)
dev.off()

###########################################################################################
################################################################################################################
pdf("pie_homer_Open_in_NT_noThreshold.pdf")
res=read.table(pipe("more Down_NT.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
pie(sort(tdown), main=,cex=.8)
title("Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)", cex.main=.9)
dev.off()

pdf("pie_homer_Open_in_SH_noThreshold.pdf")
res=read.table(pipe("more Up_SH.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
pie(sort(tdown), main=,cex=.8)
title("Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)", cex.main=.9)
dev.off()

####################################################################################
#####
pdf("enrichment_homer_in_NT_noThreshold.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Down_NT.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")

enrich = as.numeric(as.character(res[,3]))
tdown = tdown
names(enrich) = names(tdown)
barplot(sort(enrich),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
abline(h=0)
dev.off()

pdf("enrichment_homer_in_SH_noThreshold.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more Up_SH.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")

enrich = as.numeric(as.character(res[,3]))
tdown = tdown
names(enrich) = names(tdown)
barplot(sort(enrich),las=2,ylim=c(-2.5,5),ylab="Log2 Enrichment",col="lightblue3")
abline(h=0)
dev.off()


###########################################################################################
################################################################################################################
pdf("pie_homer_Open_in_NT_white.pdf")
res=read.table(pipe("more Down_NT.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8,col=NA)
title("Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)", cex.main=.9)
dev.off()

pdf("pie_homer_Open_in_SH_white.pdf")
res=read.table(pipe("more Up_SH.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8,col=NA)
title("Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)", cex.main=.9)
dev.off()

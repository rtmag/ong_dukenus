##!R

system("annotatePeaks.pl /root/ong_dukenus/h2az/shNT_IP_broad_-log10_5.bed hg19 -annStats shNT_IP_broad_-log10_5.annStats > shNT_IP_broad_-log10_5.anno")

pdf("H2AZ_-log10_5_pie.pdf")
res=read.table(pipe("more shNT_IP_broad_-log10_5.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("ATAC-down_log2.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more ATAC-down.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,3]))
names(tdown) = res[,1]
tdown = tdown[as.numeric(as.character(res[,2]))>50]
barplot(sort(tdown),las=2,ylim=c(-4,6),ylab="Log2 Enrichment over random genomic background",col="lightblue3")
abline(h=0)
dev.off()

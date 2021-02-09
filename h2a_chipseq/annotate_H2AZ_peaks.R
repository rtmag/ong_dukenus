##!R

system("annotatePeaks.pl /root/ong_dukenus/h2az/shNT_IP_broad_-log10_5.bed hg19 -annStats shNT_IP_broad_-log10_5.annStats > shNT_IP_broad_-log10_5.anno")

pdf("H2AZ_-log10_5_pie.pdf",width=9)
res=read.table(pipe("more shNT_IP_broad_-log10_5.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
ix = tdown>500
combined = sum(tdown[!ix])
names(combined) = paste("Others ",round(combined/sum(tdown)*100,digits=2),"%",sep="")
merged = c(tdown[ix],combined)
pie(sort(merged),cex=1.2, main=paste0("Genomic distribution of ",sum(tdown)," H2AZ peaks"))
dev.off()

pdf("H2AZ_-log10_5_enrichmentBars.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more shNT_IP_broad_-log10_5.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,3]))
names(tdown) = res[,1]
tdown = tdown[as.numeric(as.character(res[,2]))>500]
barplot(sort(tdown),las=2,ylim=c(-4,6),ylab="Log2 Enrichment over a random genomic background",col="lightblue3",yaxt='n')
axis(2, at=c(-4, -2, 0, 2,4,6), labels=c("-4", "-2", "0","2","4","6"), cex.axis=1, srt=45, col.ticks = "black")
abline(h=0)
dev.off()

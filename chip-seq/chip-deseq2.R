
#bedtools intersect -v -a mergedpeaks.bed -b ~/resources/hg19_consensusBlacklist.bed > mergedpeaks_noBL.bed



library(Rsubread)

x=read.table('/root/ong_dukenus/chip-DIFF/macs/mergedpeaks_noBL.bed',sep="\t",stringsAsFactors=F)

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c('/root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam',
'/root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam',
'/root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam',
'/root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam',
'/root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam',
'/root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam')

fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)
countData=fc_SE$counts
colnames(countData)=c("shNT_1","shNT_2","sh143_1","sh143_2","sh400_1","sh400_2")
saveRDS(countData,'chip_deseq2_counts.rds')
##
#

countData=readRDS('chip_deseq2_counts.rds')

require(DESeq2)
colData <- data.frame(group=c("shNT","shNT","shH2","shH2","shH2","shH2"))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dLRT_res <- results(dLRT)
dLRT_res$padj[is.na(dLRT_res$padj)]=1

saveRDS(dLRT_res,"deseq2_res.rds")

dLRT_res$padj<0.05
#########
grep "Up" diffChip_batch2_m.bed|annotatePeaks.pl - hg19 -annStats chip_up_forGREAT.annStats > chip_up_forGREAT.bed.anno

pdf("chip_up_forGREAT.pdf")
res=read.table(pipe("more chip_up_forGREAT.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
title("ChIP peaks gained after knock-down\n(9,431 regions)", cex.main=.9)
dev.off()


pdf("enrichment_chip_up_forGREAT.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more chip_up_forGREAT.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,3]))
names(tdown) = res[,1]
tdown = tdown[as.numeric(as.character(res[,2]))>50]
barplot(sort(tdown),las=2,ylim=c(-4,6),ylab="Log2 Enrichment over random genomic background",col="lightblue3")
abline(h=0)
dev.off()

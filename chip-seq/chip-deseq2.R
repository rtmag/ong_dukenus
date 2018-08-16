
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

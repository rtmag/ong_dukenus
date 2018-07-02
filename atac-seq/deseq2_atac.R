
library(Rsubread)

x=read.table('/root/ong_dukenus/ATAC-SEQ/macs2/merged_peaks.bed',sep="\t",stringsAsFactors=F)

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c("/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam",
"/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam",
"/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam",
"/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam",
"/root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam",
"/root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam" )




fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)

countData=fc_SE$counts

colnames(countData)=c ( "shH2_I_1","shH2_I_2","shH2_II_1","shH2_II_2","shNT_1","shNT_2" )

saveRDS(countData,'atac_countdata.rds')
##
#

countData=readRDS('atac_countdata.rds')


#

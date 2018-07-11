
library(Rsubread)

x=read.table('/root/ong_dukenus/ATAC-SEQ/macs2/merged_peaks_blacklisted.bed',sep="\t",stringsAsFactors=F)

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
######################################################################################


countData=readRDS('atac_countdata.rds')
require(DESeq2)
colData <- data.frame(group=  c ( "shH2","shH2","shH2","shH2","shNT","shNT" ) )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
design = ~ group)
#####################################################################################

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dLRT_res <- results(dLRT)
dLRT_res$padj[is.na(dLRT_res$padj)]=1

pdf("Diagnostic_design_pca2.pdf")
plotPCA(dLRT_vsd,ntop=90000,intgroup=c('group'))
dev.off()
#####################################################################################

design<-data.frame(cells = c ( "shH2","shH2","shH2","shH2","shNT","shNT" ) )
dds <- DESeqDataSetFromMatrix(countData = countData, colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","shH2","shNT"))
write.csv(res,"shH2_vs_shNT.csv")

#####################################################################################
res = read.csv("shH2_vs_shNT.csv",row.names=1)
res$padj[is.na(res$padj)]=1

#####################################################################################
regionFinder = function(res,chr,start,end){
       coord = row.names(res)
       coord = t(matrix(unlist(strsplit(coord,"\\_\\!\\_")),nrow=3))
       coord = data.frame( chr=coord[,1],start=as.numeric(coord[,2]),end=as.numeric(coord[,3]) )
       ix = coord$chr==chr & (coord$start<=start & coord$end>=start) & 

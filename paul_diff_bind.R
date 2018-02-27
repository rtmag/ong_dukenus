
library(Rsubread)
#
x=read.table('/root/ong_dukenus/paul_peakcalls/atac_merged_broadPeak.bed',sep="\t",stringsAsFactors=F)
ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c("/root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam",
"/root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam",
"/root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam",
"/root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam",
"/root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam",
"/root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam")

fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)
countData=fc_SE$counts
colnames(countData)=c("1_NT","2_143","3_400","4_NT","5_143","6_400")
saveRDS(countData,'atac_broadPeak_countdata.rds')
########################################################################################################
########################################################################################################
x=read.table('/root/ong_dukenus/paul_peakcalls/atac_merged_narrowPeak.bed',sep="\t",stringsAsFactors=F)
ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')
fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)
countData=fc_SE$counts
colnames(countData)=c("1_NT","2_143","3_400","4_NT","5_143","6_400")
saveRDS(countData,'atac_narrowPeak_countdata.rds')

########################################################################################################
########################################################################################################

countData=readRDS('atac_broadPeak_countdata.rds')

require(DESeq2)

colData <- data.frame(group=c("NT","143","400","NT","143","400") )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)

pdf("Diagnostic_design_pca_broadPeak.pdf")
plotPCA(dLRT_vsd,ntop=136800,intgroup=c('group'))
dev.off()

countData=readRDS('atac_narrowPeak_countdata.rds')

require(DESeq2)

colData <- data.frame(group=c("NT","143","400","NT","143","400") )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)

pdf("Diagnostic_design_pca_broadPeak.pdf")
plotPCA(dLRT_vsd,ntop=136800,intgroup=c('group'))
dev.off()




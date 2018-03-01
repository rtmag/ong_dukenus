library(Rsubread)
#
x=read.table('/root/ong_dukenus/peakcalls/star_merged_broad_noBlackList.bed',sep="\t",stringsAsFactors=F)
ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c("/root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.rmdup.out.bam")

fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)
countData=fc_SE$counts
colnames(countData)=c("1_NT","2_143","3_400","4_NT","5_143","6_400")
saveRDS(countData,'atac_broadPeak_countdata.rds')

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

########################################################################################################
########################################################################################################

require(DESeq2)
library(ggplot2)

design<-data.frame(cells = c("NT","143","NT","143") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,4,5)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","NT","143"))
bed_NT = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
bed_143 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(bed_NT,"NT_over_143.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(bed_143,"143_over_NT.bed",sep="\t",quote=F,row.names=F,col.names=F)


pdf("Volcano_NT_vs_143.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( NT / 143 ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("NT :",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("143 :",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n") 
dev.off()
#
design<-data.frame(cells = c("NT","400","NT","400") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,3,4,6)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","NT","400"))
bed_NT = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
bed_400 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(bed_NT,"NT_over_400.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(bed_400,"400_over_NT.bed",sep="\t",quote=F,row.names=F,col.names=F)


pdf("Volcano_NT_vs_400.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( NT / 400 ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("NT :",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("400 :",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n") 
dev.off()
#

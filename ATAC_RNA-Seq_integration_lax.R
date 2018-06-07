####################################################################
library(graphics)

deseq=read.csv("deseq2_results.csv")
atac = read.table("integration.txt",sep="\t")

atac = atac[atac[,4]<0.05 & abs(atac[,3])>.6 & (atac[,1]>60 | atac[,2]>60),]

ix = match(atac[,5],deseq[,8])
atac = atac[!is.na(ix),]
deseq=deseq[ix[!is.na(ix)],]

table = cbind(atac[,3],deseq$log2FoldChange)

table = table[complete.cases(table),]

smoothScatter(table,ylab="RNA log2FC",xlab="ATAC log2FC")

###


atac = atac[atac[,4]<0.05 & abs(atac[,3])>.6 & (atac[,1]>60 | atac[,2]>60),]

ix = match(atac[,5],deseq[,8])
atac = atac[!is.na(ix),]
deseq=deseq[ix[!is.na(ix)],]

table = cbind(atac[,3],deseq$log2FoldChange)

table = table[complete.cases(table),]
pdf("integration_atac_rna.pdf")
smoothScatter(table,ylab="RNA log2FC",xlab="ATAC log2FC")
dev.off()

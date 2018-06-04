options(scipen=999)

remap=read.table("ong_atac_down.txt",sep="\t",header=T)
deseq=read.csv("deseq2_results.csv")


deseq = deseq[match(toupper(remap[,1]), deseq$gene_symbol),]
deseq = deseq[!is.na(deseq$gene_symbol),]
remap = remap[!is.na(match(toupper(remap[,1]), deseq$gene_symbol)),]

x = cbind(deseq[,8],deseq[,c(3,7)],remap[,c(4,2,3)] )
colnames(x) = c("GeneName","log2FoldChange","padj","log10Evalue","ObservedOverlap","ExpectedOverlap")

x = x[order(x[,4],decreasing=T),]

write.table(x,"remap_strict_atacseq_shNT.txt",sep="\t",quote=F,row.names=F)

write.table(x[which(x[,3]<0.05 & x[,4]>110),],"remap_strict_atacseq_shNT_threshold.txt",sep="\t",quote=F,row.names=F)

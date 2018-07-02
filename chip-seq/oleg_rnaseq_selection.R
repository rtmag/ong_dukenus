
hg19=read.table('~/resources/hg19_geneBody.bed',sep="\t",stringsAsFactors=F)

genes_nt=read.table('oleg_downregulated.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"oleg_downregulated_tss.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
            
genes_h2afv=read.table('oleg_upregulated.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"oleg_upregulated_tss.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
##

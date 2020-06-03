rna = readRDS("/Users/wone/CSI/ong/atacseq_redo/MASTER_RNASEQ_TABLE_OLEG.rds")
x<-read.table("sp1.reg",stringsAsFactors=F)

ix <- match(x[,2],rna[,1])

 mat=data.frame(x,log2RNA=-log2(rna[ix, 'PostFC']))
 
 mat<-na.omit(mat)
 
 write.table(mat,"spi1_network.txt",sep="\t",quote=F,row.names=F,col.names=F)

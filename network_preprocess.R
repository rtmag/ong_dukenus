rna = readRDS("/Users/wone/CSI/ong/atacseq_redo/MASTER_RNASEQ_TABLE_OLEG.rds")
log2RNA<-(-log2(rna$PostFC))
names(log2RNA) <- rna$X

rep<-read.table("5_rep_TFs_network_edges.txt",header=T)
act<-read.table("8_tfs_prom_ass_edges.txt",header=T)

rep <- data.frame(rep,log2RNA = log2RNA[match(rep[,2],names(log2RNA))])
act <- data.frame(act,log2RNA = log2RNA[match(act[,2],names(log2RNA))])

write.table(rep,"rep.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(act,"act.txt",sep="\t",quote=F,row.names=F,col.names=F)

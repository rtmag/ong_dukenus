#######################################################################################################

atac = read.table("atac_diffreps_blacklist.tsv",header=T)
atac = atac[atac[,7]>100 | atac[,8]>100,]
write.table(atac,'atac_diffreps_blacklist_100reads.tsv',sep="\t",quote=F,row.names=F)

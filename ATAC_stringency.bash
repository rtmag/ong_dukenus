x = read.table(pipe('more atac.diffreps|grep -v "#"|cut -f1,2,3,7,8,11,12,14'),header=T,sep="\t")

x1=x[x$Control.avg>40 & x$Treatment.avg>40 & abs(x$log2FC)>1,]

write.table( x1[x1$Event=="Up",1:3], "atac_diffreps_Up_log2fc1_avg40.bed", sep="\t",quote=F,col.names=F,row.names=F)
write.table( x1[x1$Event=="Down",1:3], "atac_diffreps_Down_log2fc1_avg40.bed", sep="\t",quote=F,col.names=F,row.names=F)

x1=x[x$Control.avg>100 & x$Treatment.avg>100 & abs(x$log2FC)>1,]

write.table( x1[x1$Event=="Up",1:3], "atac_diffreps_Up_log2fc1_avg100.bed", sep="\t",quote=F,col.names=F,row.names=F)
write.table( x1[x1$Event=="Down",1:3], "atac_diffreps_Down_log2fc1_avg100.bed", sep="\t",quote=F,col.names=F,row.names=F)


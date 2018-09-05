more atac_diffreps_w100|grep -v "#"|grep -v "Treatment.avg"|awk -F"\t" '{if($14<0.05){print$0}}'|cat diffreps_head.txt  - \
> atac_diffreps_w100_FDR.bed
##!R
x = read.table("atac_diffreps_w100_FDR.bed",sep="\t",header=T)
x = x[(x$Control.avg>50 | x$Treatment.avg>50) & abs(x$log2FC)>.5,] 
write.table( x[,1:3], "atac_diffreps_w100_FDR5_50R_l2FC5.bed", sep="\t",quote=F,col.names=F,row.names=F)
##!R


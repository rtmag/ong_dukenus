###!SHELL
more atac_diffreps_w100|grep -v "#"|grep -v "Treatment.avg"|awk -F"\t" '{if($14<0.05){print$0}}'|cat diffreps_head.txt  - \
> atac_diffreps_w100_FDR.bed
###!SHELL
##!R
x = read.table("atac_diffreps_w100_FDR.bed",sep="\t",header=T)
x = x[(x$Control.avg>50 | x$Treatment.avg>50) & abs(x$log2FC)>.5,] 
write.table( x[x[,11]=="Up",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Up.bed", sep="\t",quote=F,col.names=F,row.names=F)
write.table( x[x[,11]=="Down",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Down.bed", sep="\t",quote=F,col.names=F,row.names=F)

write.table( x[x[,11]=="Up",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Up_FullStats.txt", sep="\t",quote=F,col.names=T,row.names=F)
write.table( x[x[,11]=="Down",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Down_FullStats.txt", sep="\t",quote=F,col.names=T,row.names=F)
##!R
##!SHELL
cat atac_diffreps_w100_FDR5_50R_l2FC5_Up.bed > /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
echo "#shH2AFV specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
cat atac_diffreps_w100_FDR5_50R_l2FC5_Down.bed |cut -f1-3 >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
echo "#shNT specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
###!SHELL
annotatePeaks.pl atac_diffreps_w100_FDR5_50R_l2FC5_Down.bed hg19 -annStats ATAC-down.annStats > ATAC-down.anno &
annotatePeaks.pl atac_diffreps_w100_FDR5_50R_l2FC5_Up.bed hg19 -annStats ATAC-up.annStats > ATAC-up.bed.anno &
more ATAC-down.anno|grep -i "promoter" > ATAC-down_promoter.anno
more ATAC-down.anno|grep -i -v "promoter" > ATAC-down_NOTpromoter.anno


##!R
pdf("ATAC-down_pie.pdf")
res=read.table(pipe("more ATAC-down.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("ATAC-down_log2.pdf")

##!R

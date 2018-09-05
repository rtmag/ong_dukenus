more atac_diffreps_w100|grep -v "#"|grep -v "Treatment.avg"|awk -F"\t" '{if($14<0.05){print$0}}'|cat diffreps_head.txt  - \
> atac_diffreps_w100_FDR.bed
##!R
x = read.table("atac_diffreps_w100_FDR.bed",sep="\t",header=T)
x = x[(x$Control.avg>50 | x$Treatment.avg>50) & abs(x$log2FC)>.5,] 
write.table( x[x[,11]=="Up",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Up.bed", sep="\t",quote=F,col.names=F,row.names=F)
write.table( x[x[,11]=="Down",1:3], "atac_diffreps_w100_FDR5_50R_l2FC5_Down.bed", sep="\t",quote=F,col.names=F,row.names=F)
##!R
cat atac_diffreps_w100_FDR5_50R_l2FC5_Up.bed > /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
echo "#shH2AFV specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
cat atac_diffreps_w100_FDR5_50R_l2FC5_Down.bed |cut -f1-3 >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed
echo "#shNT specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_w100_FDR5_50R_l2FC5.bed


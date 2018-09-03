more atac.diffreps|grep -v "#"|grep -v "Treatment.cnt"|bedtools intersect -a - -b ../macs2/merged_peaks_blacklisted.bed > atac_diffreps_macs2.bed
more atac.diffreps|grep "Treatment.cnt" > diffreps_head.txt

cat diffreps_head.txt atac_diffreps_macs2.bed > atac_diffreps_macs2.tmp
rm atac_diffreps_macs2.bed 
mv atac_diffreps_macs2.tmp atac_diffreps_macs2.bed 

#########!R
x=read.table("atac_diffreps_macs2.bed",sep="\t",header=T,stringsAsFactors=F)
ix = (x[,7]>100|x[,8]>100) & (abs(x[,12])>0.66) & (x[,14]<0.05)
filtered= x[ix,]


table(filtered[,11])

#########!R

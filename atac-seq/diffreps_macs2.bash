more atac.diffreps|grep -v "#"|grep -v "Treatment.cnt"|bedtools intersect -a - -b ../macs2/merged_peaks_blacklisted.bed > atac_diffreps_macs2.bed
more atac.diffreps|grep "Treatment.cnt" > diffreps_head.txt

cat diffreps_head.txt atac_diffreps_macs2.bed > atac_diffreps_macs2.tmp
rm atac_diffreps_macs2.bed 
mv atac_diffreps_macs2.tmp atac_diffreps_macs2.bed 

#########!R
x=read.table("atac_diffreps_macs2.bed",sep="\t",header=T,stringsAsFactors=F)
ix = (x[,7]>100|x[,8]>100) & (abs(x[,12])>.5) & (x[,14]<0.05)
filtered= x[ix,]


table(filtered[,11])

#########!R

more atac.diffreps_w100_s100_cell|grep -v "#"|grep -v "Treatment.cnt"|bedtools intersect -a - -b ../macs2/merged_peaks_blacklisted.bed > atac.diffreps_w100_s100_macs2.bed

cat diffreps_head.txt atac_diffreps_macs2.bed > atac_diffreps_macs2.tmp
rm atac_diffreps_macs2.bed 
mv atac_diffreps_macs2.tmp atac_diffreps_macs2.bed 

#########!R
x=read.table("atac_diffreps_macs2.bed",sep="\t",header=T,stringsAsFactors=F)
ix =  (x[,14]<0.05)
filtered= x[ix,]


table(filtered[,11])

#########!R

more atac.diffreps_w100_s100_cell| grep -v "#"|grep -v "Treatment.cnt"|awk -F"\t" '{if($14<0.05){print $1"\t"$2"\t"$3"\t"$12}}'

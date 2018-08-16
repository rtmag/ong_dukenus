macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam &

#
macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh143_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/sh143_input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh143_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_rmdup.bam &

#
macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh400_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/sh400-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh400_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_rmdup.bam &
#
# DIFFREPS
bamToBed -i /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/shNT_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/shNT_IP_2_rmdup.bed &
bamToBed -i /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh143_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh143_IP_2_rmdup.bed &
bamToBed -i /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh400_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh400_IP_2_rmdup.bed &


diffReps.pl --treatment ./bed/sh143_IP_1_rmdup.bed ./bed/sh143_IP_2_rmdup.bed \
./bed/sh400_IP_1_rmdup.bed ./bed/sh400_IP_2_rmdup.bed \
--control ./bed/shNT_IP_1_rmdup.bed ./bed/shNT_IP_2_rmdup.bed \
--meth nb --gname hg19 --report /root/ong_dukenus/chip-DIFF/diffreps/diffChip_both --frag 0 --nproc 30 &

diffReps.pl --treatment ./bed/sh143_IP_2_rmdup.bed \
./bed/sh400_IP_2_rmdup.bed \
--control ./bed/shNT_IP_2_rmdup.bed \
--meth gt --gname hg19 --report /root/ong_dukenus/chip-DIFF/diffreps/diffChip_batch2 --frag 0 --nproc 30 &
##################v
more ../diffreps/diffChip_both|grep -v "#"|grep -w Up |cut -f1,2,3 > diffboth.bed
echo "#shH2AFV-up" >> diffboth.bed
more ../diffreps/diffChip_both|grep -v "#"|grep -w Down |cut -f1,2,3 >> diffboth.bed
echo "#shNT-up" >> diffboth.bed
 
more ../diffreps/diffChip_batch2|grep -v "#"|grep -w Up |cut -f1,2,3 > diffbatch2.bed
echo "#shH2AFV-up" >> diffbatch2.bed
more ../diffreps/diffChip_batch2|grep -v "#"|grep -w Down |cut -f1,2,3 >> diffbatch2.bed
echo "#shNT-up" >> diffbatch2.bed

computeMatrix reference-point -S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R diffboth.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 -bs 20 -a 2000 -b 2000 -p max -out diffboth.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m diffboth.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out diffboth.pdf

computeMatrix reference-point -S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R diffbatch2.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 -bs 20 -a 2000 -b 2000 -p max -out diffbatch2.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m diffbatch2.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out diffbatch2.pdf
###############################################################
more diffChip_batch2 |grep -v "#"|grep -v "Treatment.avg"| \
bedtools intersect -v -a - -b ~/resources/hg19_consensusBlacklist.bed -wa|cut -f1,2,3,11,12 > diffChip_batch2.bed


more diffChip_both |grep -v "#"|grep -v "Treatment.avg"| \
bedtools intersect -v -a - -b ~/resources/hg19_consensusBlacklist.bed -wa|cut -f1,2,3,11,12 > diffChip_both.bed

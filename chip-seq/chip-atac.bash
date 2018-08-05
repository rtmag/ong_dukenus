
computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 -bs 20 -a 2000 -b 2000 -p 40 -out chip_atac.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m chip_atac.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out chip_atac.pdf

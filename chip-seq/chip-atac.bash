
computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 -bs 20 -a 2000 -b 2000 -p 40 -out chip_atac.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m chip_atac.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out chip_atac.pdf



computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/ong_dukenus/mnase_analysis/ATAC_shNT_only.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 -bs 20 -a 2000 -b 2000 -p 40 -out chip_atac_shNT.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m chip_atac_shNT.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out chip_atac_shNT.pdf
########################################################################################################################
computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/ong_dukenus/chip-seq/heatmaps/oleg_rna_chip_integration.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 -bs 20 -a 2000 -b 2000 -p 40 -out /root/ong_dukenus/mnase_analysis/chip_rna.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m /root/ong_dukenus/mnase_analysis/chip_rna.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out /root/ong_dukenus/mnase_analysis/chip_rna.pdf
##############################
computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_batch2/danpos2/shH2.Fnor.smooth.bw \
/root/ong_dukenus/mnase_batch2/danpos2/shNT.Fnor.smooth.bw \
-R /root/ong_dukenus/chip-seq/heatmaps/oleg_rna_chip_integration.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p max -out /root/ong_dukenus/mnase_analysis/mnasedanpose_rna.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC" --colorMap Blues \
-m /root/ong_dukenus/mnase_analysis/mnasedanpose_rna.mat  \
 --samplesLabel "shH2" "shNT"  \
-out /root/ong_dukenus/mnase_analysis/mnasedanpose_rna.pdf
##############################
computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
-R /root/ong_dukenus/chip-seq/heatmaps/oleg_rna_chip_integration.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 5 6 -bs 20 -a 2000 -b 2000 -p max -out /root/ong_dukenus/mnase_analysis/atac_rna.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/mnase_analysis/atac_rna.mat \
 --samplesLabel "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" "shNT" "shNT" \
-out /root/ong_dukenus/mnase_analysis/atac_rna.pdf
###

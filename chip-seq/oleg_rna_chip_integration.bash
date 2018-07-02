cat oleg_upregulated_tss.bed > oleg_rna_chip_integration.bed
echo "#UpRegulated" >> oleg_rna_chip_integration.bed
cat oleg_downregulated_tss.bed >> oleg_rna_chip_integration.bed
echo "#DownRegulated" >> oleg_rna_chip_integration.bed
###############################################################################################################################

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
-R oleg_rna_chip_integration.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 -bs 20 -a 2000 -b 2000 -p 40 -out oleg_rna_chip_integration.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m oleg_rna_chip_integration.mat \
 --samplesLabel "shNT" "shH2AFV I" "shH2AFV II" \
-out oleg_rna_chip_integration.pdf

##############################################################################################################################


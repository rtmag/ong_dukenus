
computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
-R /root/ong_dukenus/paul_diffreps/atac_diffreps_log2fc1_avg40.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/chip-seq/heatmaps/chipseq_on_atac_diffreps_log2fc1_avg40.mat

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
-R /root/ong_dukenus/paul_diffreps/atac_diffreps_log2fc1_avg100.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/chip-seq/heatmaps/chipseq_on_atac_diffreps_log2fc1_avg100.mat

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
-R /root/ong_dukenus/chip-seq/macs2/highconfidence.broadpeak --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/chip-seq/heatmaps/chipseq_comparison.mat

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p 40 -out /root/ong_dukenus/chip-seq/heatmaps/chipseq_meta_tss_h19_allKnownCanonicalGenes.mat
###############################
plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/chip-seq/heatmaps/chipseq_on_atac_diffreps_log2fc1_avg40.mat \
 --samplesLabel "shH2AFV#1-rep1" "shH2AFV#2-rep1" "shNT-rep1" \
-out /root/ong_dukenus/paul_heatmap/chipseq_on_atac_diffreps_log2fc1_avg40.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/chip-seq/heatmaps/chipseq_on_atac_diffreps_log2fc1_avg100.mat \
 --samplesLabel "shH2AFV#1-rep1" "shH2AFV#2-rep1" "shNT-rep1" \
-out /root/ong_dukenus/paul_heatmap/chipseq_on_atac_diffreps_log2fc1_avg100.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFV Peak" --colorMap Blues \
-m /root/ong_dukenus/chip-seq/heatmaps/chipseq_comparison.mat \
 --samplesLabel "shH2AFV#1-rep1" "shH2AFV#2-rep1" "shNT-rep1" \
-out /root/ong_dukenus/paul_heatmap/chipseq_comparison.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m /root/ong_dukenus/chip-seq/heatmaps/chipseq_meta_tss_h19_allKnownCanonicalGenes.mat \
 --samplesLabel "shH2AFV#1-rep1" "shH2AFV#2-rep1" "shNT-rep1" \
-out /root/ong_dukenus/paul_heatmap/chipseq_meta_tss_h19_allKnownCanonicalGenes.pdf
######

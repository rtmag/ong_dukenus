cut -f 1,2,3 atac_diffreps_Up_log2fc1_avg40.bed > atac_diffreps_log2fc1_avg40.bed
echo "# shH2AFV specific" >> atac_diffreps_log2fc1_avg40.bed
cut -f 1,2,3 atac_diffreps_Down_log2fc1_avg40.bed >> atac_diffreps_log2fc1_avg40.bed
echo "# shNT specific" >> atac_diffreps_log2fc1_avg40.bed

cut -f 1,2,3 atac_diffreps_Up_log2fc1_avg100.bed > atac_diffreps_log2fc1_avg100.bed
echo "# shH2AFV specific" >> atac_diffreps_log2fc1_avg100.bed
cut -f 1,2,3 atac_diffreps_Down_log2fc1_avg100.bed >> atac_diffreps_log2fc1_avg100.bed
echo "# shNT specific" >> atac_diffreps_log2fc1_avg100.bed
##
computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R /root/ong_dukenus/paul_diffreps/atac_diffreps_log2fc1_avg40.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg40.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg40.mat \
 --samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg40.pdf
##
computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R /root/ong_dukenus/paul_diffreps/atac_diffreps_log2fc1_avg100.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg100.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg100.mat \
 --samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out /root/ong_dukenus/paul_heatmap/atac_diffreps_log2fc1_avg100.pdf
##

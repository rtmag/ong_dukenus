cut -f 1,2,3 /root/ong_dukenus/paul_deseq2/143_over_NT.bed > /root/ong_dukenus/paul_heatmap/143_vs_NT.bed
echo "# 143 specific" >> /root/ong_dukenus/paul_heatmap/143_vs_NT.bed
cut -f 1,2,3 /root/ong_dukenus/paul_deseq2/NT_over_143.bed >> /root/ong_dukenus/paul_heatmap/143_vs_NT.bed
echo "# NT specific" >> /root/ong_dukenus/paul_heatmap/143_vs_NT.bed


computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R /root/ong_dukenus/paul_heatmap/143_vs_NT.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/paul_heatmap/143_vs_NT.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap RdBu_r \
-m /root/ong_dukenus/paul_heatmap/143_vs_NT.mat -out /root/ong_dukenus/paul_heatmap/143_vs_NT.pdf

#


perl -pe 's/ /\t/g' tss_genes_names_high_expression_shH2AFV.bed |cut -f1-3 > tss_atac_rnaseq_heatmap.bed
echo "# shH2AFV-specific" >> tss_atac_rnaseq_heatmap.bed
perl -pe 's/ /\t/g' tss_genes_names_high_expression_shNT.bed |cut -f1-3 >> tss_atac_rnaseq_heatmap.bed
echo "# shNT-specific" >> tss_atac_rnaseq_heatmap.bed


computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R tss_atac_rnaseq_heatmap.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p 40 -out tss_atac_rnaseq_heatmap.mat

plotHeatmap --xAxisLabel "" --refPointLabel "TSS" --colorMap Blues -m tss_atac_rnaseq_heatmap.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out tss_atac_rnaseq_heatmap.pdf

plotProfile -m tss_atac_rnaseq_heatmap.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colors "#ffb3ba" "#ff6961" "#bae1ff" "#aec6cf" "#77dd77" "#baffc9" \
              -out tss_atac_rnaseq_onlyProfile.pdf --perGroup 
              
              

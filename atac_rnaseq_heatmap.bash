
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
              
              
###############
# body

hg19=read.table('~/resources/hg19_geneBody.bed',sep="\t",stringsAsFactors=F)

genes_nt=read.table('genes_high_expression_shNT_log2fc.6.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"body_genes_names_high_expression_shNT.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
            
genes_h2afv=read.table('genes_high_expression_shH2AFV_log2fc.6.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"body_genes_names_high_expression_shH2AFV.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
##

perl -pe 's/ /\t/g' body_genes_names_high_expression_shH2AFV.bed > body_atac_rnaseq_heatmap.bed
echo "# Upregulated in shH2AFV" >> body_atac_rnaseq_heatmap.bed
perl -pe 's/ /\t/g' body_genes_names_high_expression_shNT.bed >> body_atac_rnaseq_heatmap.bed
echo "# Upregulated in shNT" >> body_atac_rnaseq_heatmap.bed


computeMatrix scale-regions \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R body_atac_rnaseq_heatmap.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--sortRegions descend -bs 100 -p max \
-out atac_genesbody.mat


plotHeatmap -m atac_genesbody.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colorMap Blues -out atac_genesbody.pdf


plotProfile -m  atac_genesbody.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colors "#ffb3ba" "#ff6961" "#bae1ff" "#aec6cf" "#77dd77" "#baffc9" \
              -out  atac_genesbody_onlyProfile.pdf --perGroup 

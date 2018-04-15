awk -F"\t" '{if($6>0){print $2 } }' KOBAYASHI_EGFR_SIGNALING_24HR_DN.xls | grep -v "PROBE" > EGFR_up.txt
awk -F"\t" '{if($6<0){print $2 } }' KOBAYASHI_EGFR_SIGNALING_24HR_DN.xls | grep -v "PROBE" > EGFR_dw.txt

awk -F"\t" '{if($6>0){print $2 } }' SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP.xls | grep -v "PROBE" > EMT_up.txt
awk -F"\t" '{if($6<0){print $2 } }' SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP.xls | grep -v "PROBE" > EMT_dw.txt

awk -F"\t" '{if($6>0){print $2 } }' VERHAAK_GLIOBLASTOMA_PRONEURAL.xls | grep -v "PROBE" > glioblastoma_up.txt
awk -F"\t" '{if($6<0){print $2 } }' VERHAAK_GLIOBLASTOMA_PRONEURAL.xls | grep -v "PROBE" > glioblastoma_dw.txt
#

##################################################################################################################
hg19=read.table('~/resources/hg19_tss_knownCanonical_noUnasembled.bed',sep="\t",stringsAsFactors=F)
#
genes_nt=read.table('EGFR_dw.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"EGFR_shNT.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
            
genes_h2afv=read.table('EGFR_up.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"EGFR_shH2AFV.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
#
genes_nt=read.table('EMT_dw.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"EMT_shNT.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
            
genes_h2afv=read.table('EMT_up.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"EMT_shH2AFV.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
#glioblastoma
genes_nt=read.table('glioblastoma_dw.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"glioblastoma_shNT.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
            
genes_h2afv=read.table('glioblastoma_up.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"glioblastoma_shH2AFV.bed",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
#
##################################################################################################################
cat EGFR_shH2AFV.bed > EGFR.bed
echo "# shH2AFV open" >> EGFR.bed
cat EGFR_shNT.bed >> EGFR.bed
echo "# shNT Open" >> EGFR.bed

cat EMT_shH2AFV.bed > EMT.bed
echo "# shH2AFV open" >> EMT.bed
cat EMT_shNT.bed >> EMT.bed
echo "# shNT Open" >> EMT.bed

cat glioblastoma_shH2AFV.bed > glioblastoma.bed
echo "# shH2AFV open" >> glioblastoma.bed
cat glioblastoma_shNT.bed >> glioblastoma.bed
echo "# shNT Open" >> glioblastoma.bed
##################################################################################################################

# EGFR
computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R EGFR.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out EGFR.mat

plotHeatmap --xAxisLabel "" --refPointLabel "TSS" --colorMap Reds -m EGFR.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out EGFR.pdf


plotProfile -m EGFR.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colors "#ffb3ba" "#ff6961" "#bae1ff" "#aec6cf" "#77dd77" "#baffc9" \
              -out EGFR_onlyProfile.pdf --perGroup 

# EMT

computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R EMT.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out EMT.mat

plotHeatmap --xAxisLabel "" --refPointLabel "TSS" --colorMap Reds -m EMT.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out EMT.pdf


plotProfile -m EMT.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colors "#ffb3ba" "#ff6961" "#bae1ff" "#aec6cf" "#77dd77" "#baffc9" \
              -out EMT_onlyProfile.pdf --perGroup 

# GLIOBLASTOMA

computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R glioblastoma.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out glioblastoma.mat

plotHeatmap --xAxisLabel "" --refPointLabel "TSS" --colorMap Reds -m glioblastoma.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
-out glioblastoma.pdf


plotProfile -m glioblastoma.mat \
--samplesLabel "shNT-rep1" "shNT-rep2" "shH2AFV#1-rep1" "shH2AFV#1-rep2" "shH2AFV#2-rep1" "shH2AFV#2-rep2" \
--colors "#ffb3ba" "#ff6961" "#bae1ff" "#aec6cf" "#77dd77" "#baffc9" \
              -out glioblastoma_onlyProfile.pdf --perGroup 

#####


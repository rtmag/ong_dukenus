more hg19_tss_knownCanonical.bed|awk -F"\t" '{
    if($5=="+"){print $1"\t"$2-1"\t"$2"\t"$4"\t"999"\t"$5} \
    if($5=="-"){print $1"\t"$3"\t"$3+1"\t"$4"\t"999"\t"$5} \
    }'|grep -E -   "_random|17_ctg5|chrUn_gl|6_ssto_hap|6_mann_hap|ctg9_hap1|qbl_hap|mcf_hap|dbb_hap|cox_hap|_hap" \
    > hg19_tss_knownCanonical_noUnasembled.bed

computeMatrix reference-point \
-S \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p 40 -out meta_tss_h19_allKnownCanonicalGenes.mat

plotHeatmap --xAxisLabel "" --refPointLabel "TSS" --colorMap Blues -m meta_tss_h19_allKnownCanonicalGenes.mat \
--samplesLabel "NT 1" "NT 4" "143 2" "143 5" "400 3" "400 6" \
-out meta_tss_h19_allKnownCanonicalGenes.pdf

python /home/rtm/myprograms/danpos-2.2.2/danpos.py dpos --paired 1 -o res shH2:shNT &> danpos.log

python /home/rtm/myprograms/danpos-2.2.2/danpos.py dpos --paired 1 -o res shH2:shNT &> danpos.log

~/myPrograms/kentUtils/bin/wigToBigWig shH2.Fnor.smooth.wig ../hg19.chrom.sizes shH2.Fnor.smooth.bw &
~/myPrograms/kentUtils/bin/wigToBigWig shNT.Fnor.smooth.wig ../hg19.chrom.sizes shNT.Fnor.smooth.bw &
#########################################################################################################
#########################################################################################################
computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_analysis/combined/shH2.Fnor.smooth.bw \
/root/ong_dukenus/mnase_analysis/combined/shNT.Fnor.smooth.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p max -out /root/ong_dukenus/mnase_analysis/combined/mnasedanpose_tss.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m /root/ong_dukenus/mnase_analysis/combined/mnasedanpose_tss.mat --regionsLabel "genes" \
 --samplesLabel "shH2" "shNT"  \
-out /root/ong_dukenus/mnase_analysis/combined/mnasedanpose_tss_combined.pdf

computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_analysis/batch1/shH2.Fnor.smooth.bw \
/root/ong_dukenus/mnase_analysis/batch1/shNT.Fnor.smooth.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p max -out /root/ong_dukenus/mnase_analysis/batch1/mnasedanpose_tss.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m /root/ong_dukenus/mnase_analysis/batch1/mnasedanpose_tss.mat --regionsLabel "genes" \
 --samplesLabel "shH2" "shNT"  \
-out /root/ong_dukenus/mnase_analysis/batch1/mnasedanpose_tss_batch1.pdf

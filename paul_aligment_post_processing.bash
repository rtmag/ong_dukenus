java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9.bam \
O=/root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam \
M=/root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10.bam \
O=/root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam \
M=/root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11.bam \
O=/root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam \
M=/root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12.bam \
O=/root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam \
M=/root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13.bam \
O=/root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam \
M=/root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14.bam \
O=/root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam \
M=/root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14.mfile

#

samtools index /root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam &
samtools index /root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam &
samtools index /root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam &
samtools index /root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam &
samtools index /root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam &
samtools index /root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam &


macs2 callpeak -t /root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 1_3502DukeNus_TS543-NT-031117_hs_i9_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 1_3502DukeNus_TS543-NT-031117_hs_i9_narrow &

macs2 callpeak -t /root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 2_3502DukeNus_TS543-143-031117_hg19_i10_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 2_3502DukeNus_TS543-143-031117_hg19_i10_narrow &

macs2 callpeak -t /root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 3_3502DukeNus_TS543-400-031117_hg19_i11_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 3_3502DukeNus_TS543-400-031117_hg19_i11_narrow &

macs2 callpeak -t /root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 4_3502DukeNus_TS543-NT-241117_hg19_i12_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 4_3502DukeNus_TS543-NT-241117_hg19_i12_narrow &

macs2 callpeak -t /root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 5_3502DukeNus_TS543-143-241117_hg19_i13_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 5_3502DukeNus_TS543-143-241117_hg19_i13_narrow &

macs2 callpeak -t /root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 6_3502DukeNus_TS543-400-241117_hg19_i14_broad &

macs2 callpeak -t /root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/paul_peakcalls -n 6_3502DukeNus_TS543-400-241117_hg19_i14_narrow &
#
cat 1_3502DukeNus_TS543-NT-031117_hs_i9_broad_peaks.broadPeak \
2_3502DukeNus_TS543-143-031117_hg19_i10_broad_peaks.broadPeak \
3_3502DukeNus_TS543-400-031117_hg19_i11_broad_peaks.broadPeak \
4_3502DukeNus_TS543-NT-241117_hg19_i12_broad_peaks.broadPeak \
5_3502DukeNus_TS543-143-241117_hg19_i13_broad_peaks.broadPeak \
6_3502DukeNus_TS543-400-241117_hg19_i14_broad_peaks.broadPeak |
sort -k1,1 -k2,2n |bedtools merge -i - > atac_merged_broadPeak.bed

cat 1_3502DukeNus_TS543-NT-031117_hs_i9_narrow_peaks.narrowPeak \
2_3502DukeNus_TS543-143-031117_hg19_i10_narrow_peaks.narrowPeak \
3_3502DukeNus_TS543-400-031117_hg19_i11_narrow_peaks.narrowPeak \
4_3502DukeNus_TS543-NT-241117_hg19_i12_narrow_peaks.narrowPeak \
5_3502DukeNus_TS543-143-241117_hg19_i13_narrow_peaks.narrowPeak \
6_3502DukeNus_TS543-400-241117_hg19_i14_narrow_peaks.narrowPeak |
sort -k1,1 -k2,2n |bedtools merge -i - > atac_merged_narrowPeak.bed

# filter blacklisted h19
intersectBed -v -a atac_merged_broadPeak.bed -b ~/resources/hg19consensusBlacklist.bed \
> atac_merged_broadPeak_noBlackList.bed 

intersectBed -v -a atac_merged_narrowPeak.bed -b ~/resources/hg19consensusBlacklist.bed \
> atac_merged_narrowPeak_noBlackList.bed 

#
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam -o /root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam -o /root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam -o /root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam -o /root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam -o /root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam -o /root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw
##
#

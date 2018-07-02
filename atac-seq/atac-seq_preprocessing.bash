#########################################################################################
# Trimmed

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/1_3502DukeNus_TS543-NT-031117_hs_i9_r1.fastq.gz \
/root/ong_dukenus/Dereck/1_3502DukeNus_TS543-NT-031117_hs_i9_r2.fastq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/2_3502DukeNus_TS543-143-031117_hs_i10_r1.fastq.gz \
/root/ong_dukenus/Dereck/2_3502DukeNus_TS543-143-031117_hs_i10_r2.fastq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/3_3502DukeNus_TS543-400-031117_hs_i11_r1.fastq.gz \
/root/ong_dukenus/Dereck/3_3502DukeNus_TS543-400-031117_hs_i11_r2.fastq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/4_3502DukeNus_TS543-NT-241117_hs_i12_r1.fastq.gz \
/root/ong_dukenus/Dereck/4_3502DukeNus_TS543-NT-241117_hs_i12_r2.fastq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/5_3502DukeNus_TS543-143-241117_hs_i13_r1.fastq.gz \
/root/ong_dukenus/Dereck/5_3502DukeNus_TS543-143-241117_hs_i13_r2.fastq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/ATAC-SEQ/trimmed/ \
/root/ong_dukenus/Dereck/6_3502DukeNus_TS543-400-241117_hs_i14_r1.fastq.gz \
/root/ong_dukenus/Dereck/6_3502DukeNus_TS543-400-241117_hs_i14_r2.fastq.gz &

#########################################################################################
# Aligment

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/1_3502DukeNus_TS543-NT-031117_hs_i9_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/1_3502DukeNus_TS543-NT-031117_hs_i9_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shNT_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/2_3502DukeNus_TS543-143-031117_hs_i10_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/2_3502DukeNus_TS543-143-031117_hs_i10_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/3_3502DukeNus_TS543-400-031117_hs_i11_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/3_3502DukeNus_TS543-400-031117_hs_i11_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/4_3502DukeNus_TS543-NT-241117_hs_i12_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/4_3502DukeNus_TS543-NT-241117_hs_i12_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shNT_2_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/5_3502DukeNus_TS543-143-241117_hs_i13_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/5_3502DukeNus_TS543-143-241117_hs_i13_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/ATAC-SEQ/trimmed/6_3502DukeNus_TS543-400-241117_hs_i14_r1_val_1.fq.gz \
/root/ong_dukenus/ATAC-SEQ/trimmed/6_3502DukeNus_TS543-400-241117_hs_i14_r2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_
#########################################################################################
# remove dup
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shNT_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shNT_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shNT_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shNT_2.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam \
M=/root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2.mfile

#######################################################################
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam &
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam &
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam &
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam &
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam &
samtools index /root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam &

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM \
-b /root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam \
-o /root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw
#######################################################################

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_I_1_narrow &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_I_2_narrow &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_II_1_narrow &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_II_2_narrow &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shNT_1_narrow &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shNT_2_narrow &
#

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_1_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_I_1_broad &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_I_2_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_I_2_broad &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_1_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_II_1_broad &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shH2_II_2_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shH2_II_2_broad &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shNT_1_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shNT_1_broad &

macs2 callpeak -t /root/ong_dukenus/ATAC-SEQ/bam/shNT_2_rmdup.bam --broad \
-f BAMPE --keep-dup all --nomodel -g hs -q 0.05 --outdir /root/ong_dukenus/ATAC-SEQ/macs2/ -n shNT_2_broad &
#
cat *narrowPeak|sort -k1,1 -k2,2n|mergeBed -i - > merged_peaks.bed
#
bedtools intersect -a merged_peaks.bed -b ~/resources/hg19_consensusBlacklist.bed -v > merged_peaks_blacklisted.bed

######################################################################################################################

bamToBed -i shH2_I_1_rmdup.bam > ../bed/shH2_I_1.bed &
bamToBed -i shH2_I_2_rmdup.bam > ../bed/shH2_I_2.bed &
bamToBed -i shH2_II_1_rmdup.bam > ../bed/shH2_II_1.bed &
bamToBed -i shH2_II_2_rmdup.bam > ../bed/shH2_II_2.bed &
bamToBed -i shNT_1_rmdup.bam > ../bed/shNT_1.bed &
bamToBed -i shNT_2_rmdup.bam > ../bed/shNT_2.bed &


diffReps.pl --treatment shH2_I_1.bed shH2_I_2.bed shH2_II_1.bed shH2_II_2.bed \
--control shNT_1.bed shNT_2.bed \
--meth nb --gname hg19 --report atac.diffreps_w600_nsd20 --frag 0 --nproc 50 --window 600 --nsd 20

#####################################################################################################################
grep "Up" /root/ong_dukenus/ATAC-SEQ/bed/atac_diffreps_blacklist.bed|cut -f1-3 > /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.bed
echo "#shH2AFV specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.bed
grep "Down" /root/ong_dukenus/ATAC-SEQ/bed/atac_diffreps_blacklist.bed|cut -f1-3 >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.bed
echo "#shNT specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.bed


computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.mat \
 --samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt.pdf
#####################################################################################################################
grep "Up" atac_diffreps_blacklist_50reads.tsv|cut -f1-3 > /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.bed
echo "#shH2AFV specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.bed
grep "Down" atac_diffreps_blacklist_50reads.tsv|cut -f1-3 >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.bed
echo "#shNT specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.bed

computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.mat \
 --samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_50reads.pdf
#####################################################################################################################
#####################################################################################################################
# 100 filter
grep "Up" atac_diffreps_blacklist_100reads.tsv|cut -f1-3 > /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed
echo "#shH2AFV specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed
grep "Down" atac_diffreps_blacklist_100reads.tsv|cut -f1-3 >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed
echo "#shNT specific" >> /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed

computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues \
-m /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.mat \
 --samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out /root/ong_dukenus/ATAC-SEQ/heatmap/h2_vs_nt_100reads.pdf
#####################################################################################################################

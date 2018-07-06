trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-input/FCHTL3TBBXX_L6_CHKPEI85218050142_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-input/FCHTL3TBBXX_L6_CHKPEI85218050142_2.fq.gz &

trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-IP/FCHTL3TBBXX_L6_CHKPEI85218050139_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-IP/FCHTL3TBBXX_L6_CHKPEI85218050139_2.fq.gz &

trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-input/FCHTL3TBBXX_L6_CHKPEI85218050143_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-input/FCHTL3TBBXX_L6_CHKPEI85218050143_2.fq.gz &

trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-IP/FCHTL3TBBXX_L6_CHKPEI85218050140_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-IP/FCHTL3TBBXX_L6_CHKPEI85218050140_2.fq.gz &

trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/shNT-input/FCHTL3TBBXX_L6_CHKPEI85218050141_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/shNT-input/FCHTL3TBBXX_L6_CHKPEI85218050141_2.fq.gz &

trim_galore --illumina -q 20 --paired --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/shNT-IP/FCHTL3TBBXX_L6_CHKPEI85218050138_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/shNT-IP/FCHTL3TBBXX_L6_CHKPEI85218050138_2.fq.gz &

wait

##########################################################################################      
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050142_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050142_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_
##############
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050139_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050139_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_
##############
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050143_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050143_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_
##############
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050140_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050140_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_
##############
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050141_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050141_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_
##############
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050138_1_val_1.fq.gz \
/root/ong_dukenus/mnase_batch2/trimmed/FCHTL3TBBXX_L6_CHKPEI85218050138_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_

##########################################################################################
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2.mfile &

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2.mfile &

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2.mfile &
################################
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/sh143_IP_2.mfile &

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/sh400_IP_2.mfile &

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam \
M=/root/ong_dukenus/mnase_batch2/bam/shNT_IP_2.mfile &
##########################################################################################

wait

samtools index  /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam &
samtools index  /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam &
samtools index  /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam &
samtools index  /root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_rmdup.bam &
samtools index  /root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_rmdup.bam &
samtools index  /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam &

wait
##########################################################################################

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/sh143_mnase_2.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/sh400_mnase_2.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam \
-o /root/ong_dukenus/mnase_batch2/bw/shNT_mnase_2.bw
##########################################################################################

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_IP_2_broad_.001 --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam &
######################################

computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p max -out /root/ong_dukenus/mnase_batch2/bw/chipseq_tss2.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m /root/ong_dukenus/mnase_batch2/bw/chipseq_tss2.mat --regionsLabel "genes" \
 --samplesLabel "sh143_2" "sh400_2" "shNT_2" \
-out /root/ong_dukenus/mnase_batch2/bw/chipseq_tss2.pdf


computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_batch2/bw/sh143_IP_2.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_IP_2.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_IP_2.bw \
-R /root/ong_dukenus/chip-seq/macs2/shNT_IP_2_broad_.001_peaks.broadPeak --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p max -out chipseq_peak2.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFZ peak" --colorMap Blues \
-m chipseq_peak2.mat --regionsLabel "peaks" \
 --samplesLabel "sh143_2" "sh400_2" "shNT_2" \
-out /root/ong_dukenus/mnase_batch2/bw/chipseq_peak2.pdf

#

computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_batch2/bw/sh143_mnase_2.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_mnase_2.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_mnase_2.bw \
-R /root/ong_dukenus/chip-seq/macs2/shNT_IP_2_broad_.001_peaks.broadPeak --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p max -out mnase_peak2.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFZ peak" --colorMap Blues \
-m mnase_peak2.mat --regionsLabel "peaks" \
 --samplesLabel "sh143_2" "sh400_2" "shNT_2" \
-out /root/ong_dukenus/mnase_batch2/bw/mnase_peak2.pdf



computeMatrix reference-point \
-S \
/root/ong_dukenus/mnase_batch2/bw/sh143_mnase_2.bw \
/root/ong_dukenus/mnase_batch2/bw/sh400_mnase_2.bw \
/root/ong_dukenus/mnase_batch2/bw/shNT_mnase_2.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p max -out mnase_tss2.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFZ peak" --colorMap Blues \
-m mnase_tss2.mat --regionsLabel "peaks" \
 --samplesLabel "sh143_2" "sh400_2" "shNT_2" \
-out /root/ong_dukenus/mnase_batch2/bw/mnase_tss2.pdf

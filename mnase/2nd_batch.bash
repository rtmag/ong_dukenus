trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-input/FCHTL3TBBXX_L6_CHKPEI85218050142_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-input/FCHTL3TBBXX_L6_CHKPEI85218050142_2.fq.gz &

trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-IP/FCHTL3TBBXX_L6_CHKPEI85218050139_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh143-IP/FCHTL3TBBXX_L6_CHKPEI85218050139_2.fq.gz &

trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-input/FCHTL3TBBXX_L6_CHKPEI85218050143_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-input/FCHTL3TBBXX_L6_CHKPEI85218050143_2.fq.gz &

trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-IP/FCHTL3TBBXX_L6_CHKPEI85218050140_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/Sh400-IP/FCHTL3TBBXX_L6_CHKPEI85218050140_2.fq.gz &

trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
/root/ong_dukenus/mnase_batch2/Clean/shNT-input/FCHTL3TBBXX_L6_CHKPEI85218050141_1.fq.gz \
/root/ong_dukenus/mnase_batch2/Clean/shNT-input/FCHTL3TBBXX_L6_CHKPEI85218050141_2.fq.gz &

trim_galore --illumina -q 20 --fastqc -o /root/ong_dukenus/mnase_batch2/trimmed/ \
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

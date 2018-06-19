trim_galore --illumina --paired -o /root/ong_dukenus/chip-seq/trimmed/ \
c_1_1.fq.gz \
c_1_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/SH143_gDNA/FCHTM5YBBXX_L7_CHKPEI85218050078_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/SH143_gDNA/FCHTM5YBBXX_L7_CHKPEI85218050078_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh143_IP/FCHTM5YBBXX_L7_CHKPEI85218050075_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh143_IP/FCHTM5YBBXX_L7_CHKPEI85218050075_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh143_input/FCHTM5YBBXX_L7_CHKPEI85218050072_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh143_input/FCHTM5YBBXX_L7_CHKPEI85218050072_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-IP/FCHTM5YBBXX_L7_CHKPEI85218050076_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-IP/FCHTM5YBBXX_L7_CHKPEI85218050076_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-gDNA/FCHTM5YBBXX_L7_CHKPEI85218050079_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-gDNA/FCHTM5YBBXX_L7_CHKPEI85218050079_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-input/FCHTM5YBBXX_L7_CHKPEI85218050073_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/sh400-input/FCHTM5YBBXX_L7_CHKPEI85218050073_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-IP/FCHTM5YBBXX_L7_CHKPEI85218050074_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-IP/FCHTM5YBBXX_L7_CHKPEI85218050074_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-gDNA/FCHTM5YBBXX_L7_CHKPEI85218050077_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-gDNA/FCHTM5YBBXX_L7_CHKPEI85218050077_2.fq.gz

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chip-seq/trimmed/ \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-input/FCHTM5YBBXX_L7_CHKPEI85218050071_1.fq.gz \
/root/ong_dukenus/chip-seq/MNase-seqDATA_1st_batch/shNT-input/FCHTM5YBBXX_L7_CHKPEI85218050071_2.fq.gz

##############################################################################################################
      
      
      
STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050071_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050071_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/shNT-input_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050072_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050072_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/sh143_input_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050073_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050073_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/sh400-input_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050074_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050074_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/shNT-IP_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050075_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050075_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/sh143_IP_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050076_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050076_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/sh400-IP_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050077_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050077_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/shNT-gDNA_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050078_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050078_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/SH143_gDNA_1_

STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterMismatchNoverLmax 0.09 \
--alignMatesGapMax 2000 \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--readFilesIn \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050079_1_val_1.fq.gz \
/root/ong_dukenus/chip-seq/trimmed/FCHTM5YBBXX_L7_CHKPEI85218050079_2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chip-seq/bam/sh400-gDNA_1_

#########Aligned.sortedByCoord.out.bam
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/shNT-input_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/shNT-input_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/sh143_input_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/sh143_input_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/sh143_input_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/sh400-input_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/sh400-input_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/sh400-input_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/shNT-IP_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/shNT-IP_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/sh143_IP_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/sh143_IP_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/sh400-IP_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/sh400-IP_1.mfile


java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/shNT-gDNA_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/shNT-gDNA_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/shNT-gDNA_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/SH143_gDNA_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/SH143_gDNA_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/SH143_gDNA_1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/sh400-gDNA_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/sh400-gDNA_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/sh400-gDNA_1.mfile
########################
samtools index /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/shNT-gDNA_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/SH143_gDNA_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/sh400-gDNA_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/sh143_input_1_rmdup.bam &
samtools index /root/ong_dukenus/chip-seq/bam/sh400-input_1_rmdup.bam &

######################

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/shNT-gDNA_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/shNT-gDNA_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/SH143_gDNA_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/SH143_gDNA_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/sh400-gDNA_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/sh400-gDNA_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam  \
-o /root/ong_dukenus/chip-seq/bw/shNT-input_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/sh143_input_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/sh143_input_1.bw

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/sh400-input_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/sh400-input_1.bw
###
##
#

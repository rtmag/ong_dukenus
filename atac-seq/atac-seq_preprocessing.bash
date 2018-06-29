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


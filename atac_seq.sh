#STAR --runThreadN 35 --runMode genomeGenerate --genomeDir /root/resources/hg38_noanno/ \
#--genomeFastaFiles /root/resources/hg38_allchr.fasta

##
/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/1_3502DukeNus_TS543-NT-031117_hs_i9_r1.fastq.gz \
/root/ong_dukenus/Dereck/1_3502DukeNus_TS543-NT-031117_hs_i9_r2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/1_NT_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/2_3502DukeNus_TS543-143-031117_hs_i10_r1.fastq.gz \
/root/ong_dukenus/Dereck/2_3502DukeNus_TS543-143-031117_hs_i10_r2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/2_143_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/3_3502DukeNus_TS543-400-031117_hs_i11_r1.fastq.gz \
/root/ong_dukenus/Dereck/3_3502DukeNus_TS543-400-031117_hs_i11_r2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/3_400_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/4_3502DukeNus_TS543-NT-241117_hs_i12_r1.fastq.gz \
/root/ong_dukenus/Dereck/4_3502DukeNus_TS543-NT-241117_hs_i12_r1.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/4_NT_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/5_3502DukeNus_TS543-143-241117_hs_i13_r1.fastq.gz \
/root/ong_dukenus/Dereck/5_3502DukeNus_TS543-143-241117_hs_i13_r2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/5_143_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/Dereck/6_3502DukeNus_TS543-400-241117_hs_i14_r1.fastq.gz \
/root/ong_dukenus/Dereck/6_3502DukeNus_TS543-400-241117_hs_i14_r2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/bam/6_400_

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.out.mfile

samtools index /root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.rmdup.out.bam 
samtools index /root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.rmdup.out.bam 
samtools index /root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.rmdup.out.bam 

bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/1_NT.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/2_143.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/3_400.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/4_NT.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/5_143.bw
bamCoverage -p max -bs 1 --normalizeUsingRPKM -b /root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.rmdup.out.bam -o /root/ong_dukenus/bam/6_400.bw
##

samtools merge -f -h /root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam /root/ong_dukenus/bam/NT.bam \
/root/ong_dukenus/bam/1_NT_Aligned.sortedByCoord.rmdup.out.bam\
/root/ong_dukenus/bam/4_NT_Aligned.sortedByCoord.rmdup.out.bam &

samtools merge -f -h /root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam /root/ong_dukenus/bam/143.bam \
/root/ong_dukenus/bam/2_143_Aligned.sortedByCoord.rmdup.out.bam \
/root/ong_dukenus/bam/5_143_Aligned.sortedByCoord.rmdup.out.bam &

samtools merge -f -h /root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam /root/ong_dukenus/bam/400.bam \
/root/ong_dukenus/bam/3_400_Aligned.sortedByCoord.rmdup.out.bam \
/root/ong_dukenus/bam/6_400_Aligned.sortedByCoord.rmdup.out.bam &

wait

##
samtools sort /root/ong_dukenus/bam/NT.bam /root/ong_dukenus/bam/NT_sort &
samtools sort /root/ong_dukenus/bam/143.bam /root/ong_dukenus/bam/143_sort &
samtools sort /root/ong_dukenus/bam/400.bam /root/ong_dukenus/bam/400_sort &

wait

macs2 callpeak -t /root/ong_dukenus/bam/NT_sort.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/peakcalls -n NT &

macs2 callpeak -t /root/ong_dukenus/bam/143_sort.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/peakcalls -n 143 &

macs2 callpeak -t /root/ong_dukenus/bam/400_sort.bam \
-f BAMPE --keep-dup all --nomodel --broad -g hs -q 0.05 --outdir /root/ong_dukenus/peakcalls -n 400 &

###
##
#

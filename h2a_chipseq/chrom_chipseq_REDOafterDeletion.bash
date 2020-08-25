trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq/ \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_2_FT.fq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq \
/root/ong_dukenus/chrom_chipseq/fastq/TS543_H3K27ac_CHKSEI85218090013_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/TS543-H3K27ac_CHKSEI85218090013_2_FT.fq.gz &
#######################################
/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/chrom_chipseq/fastq/TS543_H3K27ac_CHKSEI85218090013_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/TS543-H3K27ac_CHKSEI85218090013_2_FT.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_2_FT.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chrom_chipseq/bam/input_
#######################################
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chrom_chipseq/bam/input_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
M=/root/ong_dukenus/chrom_chipseq/bam/input.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
M=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac.mfile

samtools index input_rmdup.bam &
samtools index H3K27ac_rmdup.bam
#######################################

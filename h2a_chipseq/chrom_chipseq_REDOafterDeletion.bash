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
#######################################
scp -P 60057 root@172.18.149.78:/root/ong_dukenus/chrom_chipseq/bam/*bam .

java -jar /home/rtm/myprograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=input_Aligned.sortedByCoord.out.bam \
O=input_rmdup.bam \
M=input.mfile

java -jar /home/rtm/myprograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=H3K27ac_Aligned.sortedByCoord.out.bam \
O=H3K27ac_rmdup.bam \
M=H3K27ac.mfile

samtools index input_rmdup.bam &
samtools index H3K27ac_rmdup.bam

scp -P 60057 *rmdup.bam* root@172.18.149.78:/root/ong_dukenus/chrom_chipseq/bam/

macs2 callpeak -f BAMPE -g hs -q 0.01 --keep-dup auto -n H3K27ac_gbm --outdir ./ -t H3K27ac_rmdup.bam -c input_rmdup.bam

scp -P 60057 *narrowPeak root@172.18.149.78:/root/ong_dukenus/chrom_chipseq/bam/

macs2 callpeak -f BAMPE --broad -g hs -q 0.01 --keep-dup auto -n H3K27ac_gbm_broad --outdir ./ -t H3K27ac_rmdup.bam -c input_rmdup.bam

scp -P 60057 *broadPeak root@172.18.149.78:/root/ong_dukenus/chrom_chipseq/bam/

#######################################
cut -f1-3 /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks.narrowPeak > /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks.bed
cut -f1-6 /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks.narrowPeak > /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks_6Columns.bed
cut -f1-6 /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_broad_peaks.broadPeak > /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_broad_6Columns.bed


python2.7 /root/myPrograms/rose/ROSE_main_hg38.py \
-r /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
-c /root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
-i /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks_6Columns.bed \
-g HG19 -t 2500 \
-o /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_SE_hg38 &> H3K27ac_gbm_SE_hg38.log 


python2.7 /root/myPrograms/rose/ROSE_main.py \
-r /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
-c /root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
-i /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks_6Columns.bed \
-g HG19 -t 2500 \
-o /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_SE &> H3K27ac_gbm_SE.log 


python2.7 /root/myPrograms/rose/ROSE_main.py \
-r /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
-c /root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
-i /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_gbm_peaks_6Columns.bed \
-g HG19 -t 2500 \
-o /root/ong_dukenus/chrom_chipseq/rose_test &> H3K27ac_gbm_SE_test.log 


##
# GOOD ONE! TRICK IS TO USE GFF INSTEAD OF BED, bug in conversion...
python2.7 /root/myPrograms/rose/ROSE_main.py -g HG19 -i /root/ong_dukenus/chrom_chipseq/rose_test/gff/H3K27ac_gbm_peaks_6Columns.gff \
-r /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
-c /root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
-t 2500 \
-o /root/ong_dukenus/chrom_chipseq/rose_test
##


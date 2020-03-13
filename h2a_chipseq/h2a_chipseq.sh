~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR9102954 > H2A.Z1_rep1_U2OS_SRR9102954.fastq.gz 
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR9102955 > H2A.Z2_rep1_U2OS_SRR9102955.fastq.gz
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR9102956 > H2A.Z1_rep2_U2OS_SRR9102956.fastq.gz
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR9102957 > H2A.Z2_rep2_U2OS_SRR9102957.fastq.gz 
############

trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z1_rep1_U2OS_SRR9102954.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z2_rep1_U2OS_SRR9102955.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z1_rep2_U2OS_SRR9102956.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z2_rep2_U2OS_SRR9102957.fastq.gz &

############

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep1_U2OS_SRR9102954_trimmed.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z1_rep1_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z2_rep1_U2OS_SRR9102955_trimmed.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z2_rep1_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep2_U2OS_SRR9102956_trimmed.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z1_rep2_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z2_rep2_U2OS_SRR9102957_trimmed.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z2_rep2_

############

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=H2A.Z1_rep1_Aligned.sortedByCoord.out.bam \
O=H2A.Z1_rep1_rmdup.bam \
M=H2A.Z1_rep1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=H2A.Z2_rep1_Aligned.sortedByCoord.out.bam \
O=H2A.Z2_rep1_rmdup.bam \
M=H2A.Z2_rep1.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=H2A.Z1_rep2_Aligned.sortedByCoord.out.bam \
O=H2A.Z1_rep2_rmdup.bam \
M=H2A.Z1_rep2.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=H2A.Z2_rep2_Aligned.sortedByCoord.out.bam \
O=H2A.Z2_rep2_rmdup.bam \
M=H2A.Z2_rep2.mfile

############
samtools index H2A.Z1_rep1_rmdup.bam &
samtools index H2A.Z2_rep1_rmdup.bam &
samtools index H2A.Z1_rep2_rmdup.bam
samtools index H2A.Z2_rep2_rmdup.bam
############

bamCoverage -p max -bs 1 --normalizeUsing CPM -b H2A.Z1_rep1_rmdup.bam -o H2AZ1_rep1_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsing CPM -b H2A.Z2_rep1_rmdup.bam -o H2AZ2_rep1_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsing CPM -b H2A.Z1_rep2_rmdup.bam -o H2AZ1_rep2_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsing CPM -b H2A.Z2_rep2_rmdup.bam -o H2AZ2_rep2_rmdup.bw

#############

macs2 callpeak -f BAM -g hs -q 0.01 --broad --keep-dup auto -n H2AZ1_broad --outdir ./ -t H2A.Z1_rep1_rmdup.bam H2A.Z1_rep2_rmdup.bam &
macs2 callpeak -f BAM -g hs -q 0.01 --broad --keep-dup auto -n H2AZ2_broad --outdir ./ -t H2A.Z2_rep1_rmdup.bam H2A.Z2_rep2_rmdup.bam &

macs2 callpeak -f BAM -g hs -q 0.01 --call-summits --keep-dup auto -n H2AZ1_narrow --outdir ./ -t H2A.Z1_rep1_rmdup.bam H2A.Z1_rep2_rmdup.bam &
macs2 callpeak -f BAM -g hs -q 0.01 --call-summits --keep-dup auto -n H2AZ2_narrow --outdir ./ -t H2A.Z2_rep1_rmdup.bam H2A.Z2_rep2_rmdup.bam &
##################
more H2AZ2_broad_peaks.broadPeak|awk -F"\t" '{if($5>200 && $7>2){print $0}}' > H2AZ2_filtered.broadPeak
more H2AZ1_broad_peaks.broadPeak|awk -F"\t" '{if($5>200 && $7>2){print $0}}' > H2AZ1_filtered.broadPeak
##################
bedtools intersect -a H2AZ1_hg19.bed -b H2AZ2_hg19.bed -v > H2AZ1_hg19_only.bed
bedtools intersect -b H2AZ1_hg19.bed -a H2AZ2_hg19.bed -v > H2AZ2_hg19_only.bed
##################
computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R H2AZ1_hg19.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out H2AZ1_peak_atac.mat \
--outFileNameMatrix H2AZ1_peak_atac.rmat

computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R H2AZ2_hg19.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out H2AZ2_peak_atac.mat \
--outFileNameMatrix H2AZ2_peak_atac.rmat
##################

computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R H2AZ1_hg19_only.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out H2AZ1_only_atac.mat \
--outFileNameMatrix H2AZ1_only_atac.rmat

computeMatrix reference-point \
-S \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shNT_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_I_2_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_1_rmdup.bw \
/root/ong_dukenus/ATAC-SEQ/bw/shH2_II_2_rmdup.bw \
-R H2AZ2_hg19_only.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 1000 -b 1000 -p 40 -out H2AZ2_only_atac.mat \
--outFileNameMatrix H2AZ2_only_atac.rmat
##################

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AZ1 Peaks" --colorMap Blues \
-m H2AZ1_peak_atac.mat \
--samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out H2AZ1_peak_atac.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AZ2 Peaks" --colorMap Blues \
-m H2AZ2_peak_atac.mat \
--samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out H2AZ2_peak_atac.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AZ1 Only" --colorMap Blues \
-m H2AZ1_only_atac.mat \
--samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out H2AZ1_only_atac.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AZ1 Only" --colorMap Blues \
-m H2AZ2_only_atac.mat \
--samplesLabel "shNT" "shNT" "shH2AFV I" "shH2AFV I" "shH2AFV II" "shH2AFV II" \
-out H2AZ2_only_atac.pdf

############################################
grep "Down" /root/ong_dukenus/ATAC-SEQ/heatmap/atac_diffreps_blacklist_100reads.tsv|cut -f1-3 > Down_atacseq.bed
# liftover and filezzilla

computeMatrix reference-point \
-S \
/root/ong_dukenus/h2a_chipseq/fastq/H2AZ1_rep1_rmdup.bw \
/root/ong_dukenus/h2a_chipseq/fastq/H2AZ1_rep2_rmdup.bw \
/root/ong_dukenus/h2a_chipseq/fastq/H2AZ2_rep1_rmdup.bw \
/root/ong_dukenus/h2a_chipseq/fastq/H2AZ2_rep2_rmdup.bw \
-R Down_atacseq_hg38.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 3 4 -bs 20 -a 2000 -b 2000 -p 40 -out H2AZ_ON_ATACDOWN.mat \
--outFileNameMatrix H2AZ_ON_ATACDOWN.rmat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC-Lost" --colorMap Reds \
-m H2AZ_ON_ATACDOWN.mat \
--samplesLabel "H2AZ1-I" "H2AZ1-II" "H2AZ2-I" "H2AZ2-II" \
-out H2AZ_ON_ATACDOWN.pdf

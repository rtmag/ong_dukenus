############
trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq/ \
/root/ong_dukenus/chrom_chipseq/TS543-H3K27ac/TS543_H3K27ac_CHKSEI85218090013_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/TS543-H3K27ac/TS543-H3K27ac_CHKSEI85218090013_2_FT.fq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq/ \
/root/ong_dukenus/chrom_chipseq/TS543-H3K4me3/H3K4me3_CHKSEI85218090014_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/TS543-H3K4me3/H3K4me3_CHKSEI85218090014_2_FT.fq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq/ \
/root/ong_dukenus/chrom_chipseq/TS543_input/T543_input_CHKSEI85218090010_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/TS543_input/T543_input_CHKSEI85218090010_2_FT.fq.gz &

############
/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/chrom_chipseq/fastq/H3K4me3_CHKSEI85218090014_1_FT_val_1.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/H3K4me3_CHKSEI85218090014_2_FT_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chrom_chipseq/bam/H3K4me3_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_1_FT_val_1.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_2_FT_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chrom_chipseq/bam/input_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg19_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn /root/ong_dukenus/chrom_chipseq/fastq/TS543_H3K27ac_CHKSEI85218090013_1_FT_val_1.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/TS543-H3K27ac_CHKSEI85218090013_2_FT_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ong_dukenus/chrom_chipseq/bam/H3K27ac_
############
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chrom_chipseq/bam/H3K4me3_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chrom_chipseq/bam/H3K4me3_rmdup.bam \
M=/root/ong_dukenus/chrom_chipseq/bam/H3K4me3.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chrom_chipseq/bam/input_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bam \
M=/root/ong_dukenus/chrom_chipseq/bam/input.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bam \
M=/root/ong_dukenus/chrom_chipseq/bam/H3K27ac.mfile
############
samtools index H3K4me3_rmdup.bam &
samtools index input_rmdup.bam
samtools index H3K27ac_rmdup.bam
#############
bamCoverage -p max -bs 1 --normalizeUsing CPM -b H3K4me3_rmdup.bam -o H3K4me3_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsing CPM -b input_rmdup.bam -o input_rmdup.bw
bamCoverage -p max -bs 1 --normalizeUsing CPM -b H3K27ac_rmdup.bam -o H3K27ac_rmdup.bw
#############
cut -f 1,2,3 /root/ong_dukenus/h2a_chipseq/heatmap/Up_atacseq.bed > /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff.DPbed
echo "#Up ATAC" >> /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff.DPbed
cut -f 1,2,3 /root/ong_dukenus/h2a_chipseq/heatmap/Down_atacseq.bed >> /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff.DPbed
echo "#Down ATAC" >> /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff.DPbed


computeMatrix reference-point \
-S \
/root/ong_dukenus/chrom_chipseq/bam/H3K4me3_rmdup.bw \
/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bw \
/root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bw \
-R /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff.DPbed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 4000 -b 4000 -p 40 -out /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom.mat \
--outFileNameMatrix /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom.rmat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --colorMap "Greens" "Reds" "Blues" \
-m /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom.mat \
--samplesLabel "H3K4me3" "H3K27ac" "Input" \
-out /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom.pdf

####

computeMatrix reference-point \
-S \
/root/ong_dukenus/chrom_chipseq/bam/H3K4me3_rmdup.bw \
/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bw \
/root/ong_dukenus/chrom_chipseq/bam/input_rmdup.bw \
-R /root/ong_dukenus/h2a_chipseq/heatmap/Down_atacseq.bed --referencePoint center \
--sortRegions descend --sortUsingSamples 1 2 -bs 20 -a 4000 -b 4000 -p 40 -out /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom_down.mat \
--outFileNameMatrix /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom_down.rmat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --colorMap "Greens" "Reds" "Blues" \
-m /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom_down.mat \
--samplesLabel "H3K4me3" "H3K27ac" "Input" \
-out /root/ong_dukenus/h2a_chipseq/heatmap/atac_diff_chrom_down.pdf


################
# multipanel chips

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chrom_chipseq/bam/H3K4me3_rmdup.bw \
/root/ong_dukenus/chrom_chipseq/bam/H3K27ac_rmdup.bw \
-R /root/ong_dukenus/h2a_chipseq/heatmap/Down_atacseq.bed --referencePoint center \
--sortRegions descend -bs 20 -a 4000 -b 4000 -p 40 -out /root/ong_dukenus/h2a_chipseq/heatmap/multipanel_chip_chrom.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --colorMap "Purples" "Reds" "Greens" \
-m /root/ong_dukenus/h2a_chipseq/heatmap/multipanel_chip_chrom.mat \
--samplesLabel "H2AZ" "H3K4me3" "H3K27ac" --refPointLabel "ATAC-Down" --regionsLabel "ATAC-Down" --zMax .3 \
-out /root/ong_dukenus/h2a_chipseq/heatmap/multipanel_chip_chrom.pdf
################

macs2 callpeak -f BAM -g hs -q 0.01 --broad --keep-dup auto -n H3K4me3 --outdir ./ -t H3K4me3_rmdup.bam -c input_rmdup.bam &
macs2 callpeak -f BAM -g hs -q 0.01 --broad --keep-dup auto -n H3K27ac --outdir ./ -t H3K27ac_rmdup.bam -c input_rmdup.bam &
#####

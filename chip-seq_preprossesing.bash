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
macs2 callpeak -f BAMPE -g hs -q 0.05 --broad --keep-dup auto -n shNT_IP_broad --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.05 --call-summits --keep-dup auto -n shNT_IP_narrow --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &
############################
macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_IP_broad_.001 --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &
######################################

sudo python -m pip install rpy2==2.8.6
#################################


multiBigwigSummary bins -b \
/root/ong_dukenus/paul_bw/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bw \
/root/ong_dukenus/paul_bw/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bw \
/root/ong_dukenus/paul_bw/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bw \
/root/ong_dukenus/paul_bw/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bw \
/root/ong_dukenus/paul_bw/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bw \
/root/ong_dukenus/paul_bw/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bw \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/SH143_gDNA_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-gDNA_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-gDNA_1.bw \
/root/ong_dukenus/chip-seq/bw/sh143_input_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-input_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-input_1.bw \
 -p max -bs 200 -o multiBigwig.npz

plotPCA --corData multiBigwig.npz -o multiBigwig_pca.pdf
plotCorrelation --whatToPlot heatmap --corData multiBigwig.npz -c spearman -o multiBigwig_plotcorrelation.pdf
plotCorrelation --whatToPlot scatterplot --corData multiBigwig.npz -c spearman -o multiBigwig_plotcorrelation_scatter.pdf
#######################################################################################################

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
-R /root/ong_dukenus/chip-seq/macs2/highconfidence.broadpeak --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out geochipseq_comparison.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFV Peak" --colorMap Blues \
-m geochipseq_comparison.mat \
 --samplesLabel "SKmel147-H2AZ" "SKMel147-EGFP_Z1" "SKMel147-EGFP_Z2" "Melanocytes-H2AZ" \
-out geochipseq_comparison.pdf

multiBigwigSummary bins -b \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
 -p max -bs 200 -o multiBigwig200.npz

multiBigwigSummary bins -b \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
 -p max -bs 5000 -o multiBigwig5000.npz
 
 multiBigwigSummary bins -b \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
 -p max -bs 10000 -o multiBigwig10000.npz

plotCorrelation --whatToPlot heatmap --corData multiBigwig200.npz -c spearman -o multiBigwig_plotcorrelation200.pdf
plotCorrelation --whatToPlot scatterplot --corData multiBigwig200.npz -c spearman -o multiBigwig_plotcorrelation_scatter200.pdf
plotCorrelation --whatToPlot heatmap --corData multiBigwig5000.npz -c spearman -o multiBigwig_plotcorrelation5000.pdf
plotCorrelation --whatToPlot scatterplot --corData multiBigwig5000.npz -c spearman -o multiBigwig_plotcorrelation_scatter5000.pdf
plotCorrelation --whatToPlot heatmap --corData multiBigwig10000.npz -c spearman -o multiBigwig_plotcorrelation10000.pdf
plotCorrelation --whatToPlot scatterplot --corData multiBigwig10000.npz -c spearman -o multiBigwig_plotcorrelation_scatter10000.pdf


computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 2000 -b 2000 -p max -out mnaseseq.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m mnaseseq.mat --regionsLabel "genes" \
 --samplesLabel "SKmel147-H2AZ" "SKMel147-EGFP_Z1" "SKMel147-EGFP_Z2" "Melanocytes-H2AZ" \
-out mnaseseq2.pdf
###############################
##########
##########
##########
multiBigwigSummary BED-file --BED /root/ong_dukenus/chip-seq/macs2/highconfidence.broadpeak -b \
/root/ong_dukenus/chip-seq/bw/sh143_IP_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-IP_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665991_SKmel147-H2AZ-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665993_SKMel147-EGFP_Z1-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665994_SKMel147-EGFP_Z2-FE_hg19.bigWig \
/root/ong_dukenus/chip-seq/ong_geo/GSM1665995_Melanocytes-H2AZ-FE_hg19.bigWig \
 -p max -bs 5000 -o multiBig_h2peaks.npz
 
plotCorrelation --whatToPlot heatmap --corData multiBig_h2peaks.npz -c spearman -o multiBig_h2peaks_heatmap.pdf
plotCorrelation --whatToPlot scatterplot --corData multiBig_h2peaks.npz -c spearman -o multiBig_h2peaks_scatter.pdf
#################
computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_input_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-input_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-input_1.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p max -out mnaseseq.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m mnaseseq.mat --regionsLabel "genes" \
 --samplesLabel "sh143" "sh400" "shNT" \
-out mnaseseq20.pdf

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_input_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-input_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-input_1.bw \
-R /root/resources/hg19_tss_knownCanonical_noUnasembled.bed --referencePoint center \
--sortRegions descend -bs 1 -a 1000 -b 1000 -p max -out mnaseseq.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "TSS" --colorMap Blues \
-m mnaseseq.mat --regionsLabel "genes" \
 --samplesLabel "sh143" "sh400" "shNT" \
-out mnaseseq1.pdf
###########
###
###


computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_input_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-input_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-input_1.bw \
-R /root/ong_dukenus/chip-seq/macs2/summit_highconfidence.bed --referencePoint center \
--sortRegions descend -bs 20 -a 500 -b 500 -p max -out mnaseseqpeak20.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFV peak" --colorMap Blues \
-m mnaseseqpeak20.mat --regionsLabel "peaks" \
 --samplesLabel "sh143" "sh400" "shNT" \
-out mnaseseqpeak20.pdf

computeMatrix reference-point \
-S \
/root/ong_dukenus/chip-seq/bw/sh143_input_1.bw \
/root/ong_dukenus/chip-seq/bw/sh400-input_1.bw \
/root/ong_dukenus/chip-seq/bw/shNT-input_1.bw \
-R /root/ong_dukenus/chip-seq/macs2/summit_highconfidence.bed --referencePoint center \
--sortRegions descend -bs 1 -a 500 -b 500 -p max -out mnaseseqpeak1.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "H2AFV peak" --colorMap Blues \
-m mnaseseqpeak1.mat --regionsLabel "peaks" \
 --samplesLabel "sh143" "sh400" "shNT" \
-out mnaseseqpeak1.pdf
###########

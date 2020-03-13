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
more H2AZ2_broad_peaks.broadPeak|awk -F"\t" '{if($5>100 && $7>2){print $0}}' > H2AZ2_filtered.broadPeak
more H2AZ1_broad_peaks.broadPeak|awk -F"\t" '{if($5>100 && $7>2){print $0}}' > H2AZ1_filtered.broadPeak
##################


~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files -Z SRR9102954 > H2A.Z1_rep1_U2OS_SRR9102954.fastq &
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files -Z SRR9102955 > H2A.Z2_rep1_U2OS_SRR9102955.fastq &
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files -Z SRR9102956 > H2A.Z1_rep2_U2OS_SRR9102956.fastq &
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files -Z SRR9102957 > H2A.Z2_rep2_U2OS_SRR9102957.fastq &
############

trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z1_rep1_U2OS_SRR9102954.fastq &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z2_rep1_U2OS_SRR9102955.fastq &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z1_rep2_U2OS_SRR9102956.fastq &
trim_galore --illumina -q 20 --fastqc -o ./ H2A.Z2_rep2_U2OS_SRR9102957.fastq &

############

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep1_U2OS_SRR9102954.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z1_rep1_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep1_U2OS_SRR9102954.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z2_rep1_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep1_U2OS_SRR9102954.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z1_rep2_

/root/myPrograms/STAR/bin/STAR --genomeDir /root/resources/hg38_noanno/ \
--readFilesCommand zcat \
--runThreadN 35 \
--alignIntronMax 1 \
--outFilterScoreMinOverLread 0.40 \
--outFilterMatchNminOverLread 0.40 \
--alignEndsType EndToEnd \
--readFilesIn H2A.Z1_rep1_U2OS_SRR9102954.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix H2A.Z2_rep2_

############

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/chip-seq/bam/shNT-input_1_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam \
M=/root/ong_dukenus/chip-seq/bam/shNT-input_1.mfile

############

bamCoverage -p max -bs 1 --normalizeUsing CPM -b /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam \
-o /root/ong_dukenus/chip-seq/bw/shNT-IP_1.bw

#############


macs2 callpeak -f BAMPE -g hs -q 0.05 --broad --keep-dup auto -n shNT_IP_broad --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.05 --call-summits --keep-dup auto -n shNT_IP_narrow --outdir /root/ong_dukenus/chip-seq/macs2/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &

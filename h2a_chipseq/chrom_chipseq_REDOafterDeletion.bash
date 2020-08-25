trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq/ \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/T543_input_CHKSEI85218090010_2_FT.fq.gz &

trim_galore --illumina --paired -q 20 --fastqc -o /root/ong_dukenus/chrom_chipseq/fastq \
/root/ong_dukenus/chrom_chipseq/fastq/TS543_H3K27ac_CHKSEI85218090013_1_FT.fq.gz \
/root/ong_dukenus/chrom_chipseq/fastq/TS543-H3K27ac_CHKSEI85218090013_2_FT.fq.gz &
#######################################

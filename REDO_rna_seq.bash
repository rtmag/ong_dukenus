ls /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/*/*L001_R1*| perl -pe 's/(.+)L001_R1_001.fastq.gz/cat $1L001_R1_001.fastq.gz $1L002_R1_001.fastq.gz $1L003_R1_001.fastq.gz $1L004_R1_001.fastq.gz \> $1.fastq.gz\n/g'

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-1_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-2_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-3_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_1-109818737/sh202-1_S8_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_1-109818737/sh202-1_S8_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_1-109818737/sh202-1_S8_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_1-109818737/sh202-1_S8_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/shNT-1_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_2-109819732/sh202-2_S7_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_2-109819732/sh202-2_S7_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_2-109819732/sh202-2_S7_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_2-109819732/sh202-2_S7_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/shNT-2_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_3-109827732/sh202-3_S5_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_3-109827732/sh202-3_S5_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_3-109827732/sh202-3_S5_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh202_3-109827732/sh202-3_S5_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/shNT-3_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh400-1_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh400-2_R1.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L001_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L002_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L003_R1_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L004_R1_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh400-3_R1.fastq.gz &


ls /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/*/*L001_R1*| perl -pe 's/(.+)L001_R1_001.fastq.gz/cat $1L001_R2_001.fastq.gz $1L002_R2_001.fastq.gz $1L003_R2_001.fastq.gz $1L004_R2_001.fastq.gz \n/g'

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_1-109826732/sh143-1_S6_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-1_R2.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_2-109814733/sh143-2_S3_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-2_R2.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh143_3-109813729/sh143-3_S9_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/fastq/sh143-3_R2.fastq.gz &

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_1-109823740/sh400-1_S2_.fastq.gz

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_2-109815738/sh400-2_S4_.fastq.gz

cat /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L001_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L002_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L003_R2_001.fastq.gz /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_L004_R2_001.fastq.gz > /home/oleg/DATA/H2AFV_PROJECT/gsc_543_RNA-seq_raw/sh400_3-109823741/sh400-3_S1_.fastq.gz


#######
ls *L001_R1*| perl -pe 's/(.+)L001_R1_001.fastq.gz/cat $1L001_R1_001.fastq.gz $1L002_R1_001.fastq.gz $1L003_R1_001.fastq.gz $1L004_R1_001.fastq.gz \> $1R1.fastq.gz/g'

cat sh143-1_S6_L001_R1_001.fastq.gz sh143-1_S6_L002_R1_001.fastq.gz sh143-1_S6_L003_R1_001.fastq.gz sh143-1_S6_L004_R1_001.fastq.gz > sh143-1_S6_R1.fastq.gz
cat sh143-2_S3_L001_R1_001.fastq.gz sh143-2_S3_L002_R1_001.fastq.gz sh143-2_S3_L003_R1_001.fastq.gz sh143-2_S3_L004_R1_001.fastq.gz > sh143-2_S3_R1.fastq.gz
cat sh143-3_S9_L001_R1_001.fastq.gz sh143-3_S9_L002_R1_001.fastq.gz sh143-3_S9_L003_R1_001.fastq.gz sh143-3_S9_L004_R1_001.fastq.gz > sh143-3_S9_R1.fastq.gz
cat sh202-1_S8_L001_R1_001.fastq.gz sh202-1_S8_L002_R1_001.fastq.gz sh202-1_S8_L003_R1_001.fastq.gz sh202-1_S8_L004_R1_001.fastq.gz > sh202-1_S8_R1.fastq.gz
cat sh202-2_S7_L001_R1_001.fastq.gz sh202-2_S7_L002_R1_001.fastq.gz sh202-2_S7_L003_R1_001.fastq.gz sh202-2_S7_L004_R1_001.fastq.gz > sh202-2_S7_R1.fastq.gz
cat sh202-3_S5_L001_R1_001.fastq.gz sh202-3_S5_L002_R1_001.fastq.gz sh202-3_S5_L003_R1_001.fastq.gz sh202-3_S5_L004_R1_001.fastq.gz > sh202-3_S5_R1.fastq.gz
cat sh400-1_S2_L001_R1_001.fastq.gz sh400-1_S2_L002_R1_001.fastq.gz sh400-1_S2_L003_R1_001.fastq.gz sh400-1_S2_L004_R1_001.fastq.gz > sh400-1_S2_R1.fastq.gz
cat sh400-2_S4_L001_R1_001.fastq.gz sh400-2_S4_L002_R1_001.fastq.gz sh400-2_S4_L003_R1_001.fastq.gz sh400-2_S4_L004_R1_001.fastq.gz > sh400-2_S4_R1.fastq.gz
cat sh400-3_S1_L001_R1_001.fastq.gz sh400-3_S1_L002_R1_001.fastq.gz sh400-3_S1_L003_R1_001.fastq.gz sh400-3_S1_L004_R1_001.fastq.gz > sh400-3_S1_R1.fastq.gz


ls *L001_R1*| perl -pe 's/(.+)L001_R1_001.fastq.gz/cat $1L001_R2_001.fastq.gz $1L002_R2_001.fastq.gz $1L003_R2_001.fastq.gz $1L004_R2_001.fastq.gz \> $1R2.fastq.gz/g'



cat sh143-1_S6_L001_R2_001.fastq.gz sh143-1_S6_L002_R2_001.fastq.gz sh143-1_S6_L003_R2_001.fastq.gz sh143-1_S6_L004_R2_001.fastq.gz > sh143-1_S6_R2.fastq.gz
cat sh143-2_S3_L001_R2_001.fastq.gz sh143-2_S3_L002_R2_001.fastq.gz sh143-2_S3_L003_R2_001.fastq.gz sh143-2_S3_L004_R2_001.fastq.gz > sh143-2_S3_R2.fastq.gz
cat sh143-3_S9_L001_R2_001.fastq.gz sh143-3_S9_L002_R2_001.fastq.gz sh143-3_S9_L003_R2_001.fastq.gz sh143-3_S9_L004_R2_001.fastq.gz > sh143-3_S9_R2.fastq.gz
cat sh202-1_S8_L001_R2_001.fastq.gz sh202-1_S8_L002_R2_001.fastq.gz sh202-1_S8_L003_R2_001.fastq.gz sh202-1_S8_L004_R2_001.fastq.gz > sh202-1_S8_R2.fastq.gz
cat sh202-2_S7_L001_R2_001.fastq.gz sh202-2_S7_L002_R2_001.fastq.gz sh202-2_S7_L003_R2_001.fastq.gz sh202-2_S7_L004_R2_001.fastq.gz > sh202-2_S7_R2.fastq.gz
cat sh202-3_S5_L001_R2_001.fastq.gz sh202-3_S5_L002_R2_001.fastq.gz sh202-3_S5_L003_R2_001.fastq.gz sh202-3_S5_L004_R2_001.fastq.gz > sh202-3_S5_R2.fastq.gz
cat sh400-1_S2_L001_R2_001.fastq.gz sh400-1_S2_L002_R2_001.fastq.gz sh400-1_S2_L003_R2_001.fastq.gz sh400-1_S2_L004_R2_001.fastq.gz > sh400-1_S2_R2.fastq.gz
cat sh400-2_S4_L001_R2_001.fastq.gz sh400-2_S4_L002_R2_001.fastq.gz sh400-2_S4_L003_R2_001.fastq.gz sh400-2_S4_L004_R2_001.fastq.gz > sh400-2_S4_R2.fastq.gz
cat sh400-3_S1_L001_R2_001.fastq.gz sh400-3_S1_L002_R2_001.fastq.gz sh400-3_S1_L003_R2_001.fastq.gz sh400-3_S1_L004_R2_001.fastq.gz > sh400-3_S1_R2.fastq.gz



###################################################
mv sh143-1_S6_R1.fastq.gz RNA_shH2_I_rep1_R1.fastq.gz
mv sh143-1_S6_R2.fastq.gz RNA_shH2_I_rep1_R2.fastq.gz
mv sh143-2_S3_R1.fastq.gz RNA_shH2_I_rep2_R1.fastq.gz
mv sh143-2_S3_R2.fastq.gz RNA_shH2_I_rep2_R2.fastq.gz
mv sh143-3_S9_R1.fastq.gz RNA_shH2_I_rep3_R1.fastq.gz
mv sh143-3_S9_R2.fastq.gz RNA_shH2_I_rep3_R2.fastq.gz
mv sh202-1_S8_R1.fastq.gz RNA_shNT_rep1_R1.fastq.gz
mv sh202-1_S8_R2.fastq.gz RNA_shNT_rep1_R2.fastq.gz
mv sh202-2_S7_R1.fastq.gz RNA_shNT_rep2_R1.fastq.gz
mv sh202-2_S7_R2.fastq.gz RNA_shNT_rep2_R2.fastq.gz
mv sh202-3_S5_R1.fastq.gz RNA_shNT_rep3_R1.fastq.gz
mv sh202-3_S5_R2.fastq.gz RNA_shNT_rep3_R2.fastq.gz
mv sh400-1_S2_R1.fastq.gz RNA_shH2_II_rep1_R1.fastq.gz
mv sh400-1_S2_R2.fastq.gz RNA_shH2_II_rep1_R2.fastq.gz
mv sh400-2_S4_R1.fastq.gz RNA_shH2_II_rep2_R1.fastq.gz
mv sh400-2_S4_R2.fastq.gz RNA_shH2_II_rep2_R2.fastq.gz
mv sh400-3_S1_R1.fastq.gz RNA_shH2_II_rep3_R1.fastq.gz
mv sh400-3_S1_R2.fastq.gz RNA_shH2_II_rep3_R2.fastq.gz
################################################

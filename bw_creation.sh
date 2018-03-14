samtools merge -f -h sh400-1.bam \
sh400.bam \
sh400-1.bam \
sh400-2.bam \
sh400-3.bam &

samtools merge -f -h sh_202_1.bam \
sh_202.bam \
sh_202_1.bam \
sh_202_2.bam \
sh_202_3.bam &

######

samtools sort sh400.bam > sh400_sort.bam &
samtools sort sh400-1.bam > sh400-1_sort.bam &
samtools sort sh400-2.bam > sh400-2_sort.bam &
samtools sort sh400-3.bam > sh400-3_sort.bam &

samtools sort sh_202.bam > sh_202_sort.bam &
samtools sort sh_202_1.bam > sh_202_1_sort.bam &
samtools sort sh_202_2.bam > sh_202_2_sort.bam &
samtools sort sh_202_3.bam > sh_202_3_sort.bam &

######

samtools index sh400-1.bam &
samtools index sh400-2.bam &
samtools index sh400-3.bam &
samtools index sh400.bam &
samtools index sh202-1.bam &
samtools index sh202-2.bam &
samtools index sh202-3.bam &
samtools index sh202.bam &


#####

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh400-1.bam \
-o sh400-1.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh400-2.bam \
-o sh400-2.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh400-3.bam \
-o sh400-3.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh202-1.bam \
-o sh202-1.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh202-1.bam \
-o sh202-1.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh202-1.bam \
-o sh202-1.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh400.bam \
-o sh400.bw

bamCoverage -p 10 -bs 1 --normalizeUsingRPKM -b sh202.bam \
-o sh202.bw

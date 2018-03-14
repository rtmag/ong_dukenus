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

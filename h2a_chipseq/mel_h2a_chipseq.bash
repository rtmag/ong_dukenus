~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR1993915 > Skmel147_H2AZ_ChIP.fastq.gz 
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR1993916 > Skmel147_H2A-GFP.fastq.gz 
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR1993917 > Skmel147_H2AZ.1-GFP.fastq.gz 
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR1993918 > Skmel147_H2AZ.2-GFP.fastq.gz 
~/myPrograms/sra-tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip -Z SRR1993919 > Melanocytes_H2AZ_ChIP.fastq.gz 
############
trim_galore --illumina -q 20 --fastqc -o ./ Skmel147_H2AZ_ChIP.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ Skmel147_H2A-GFP.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ Skmel147_H2AZ.1-GFP.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ Skmel147_H2AZ.2-GFP.fastq.gz &
trim_galore --illumina -q 20 --fastqc -o ./ Melanocytes_H2AZ_ChIP.fastq.gz &
############

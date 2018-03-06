library(BiocInstaller)
biocLite(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
           "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene",
           "phastCons100way.UCSC.hg19"))
           
library(ATACseqQC)
## input the bamFile from the ATACseqQC package 
bamfile <- system.file("extdata", "GL1.bam", 
                        package="ATACseqQC", mustWork=TRUE)
bamfile.labels <- gsub(".bam", "", basename(bamfile))

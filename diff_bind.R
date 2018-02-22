
# cat 1_3502DukeNus_TS543-NT-031117_hs_i9_peaks.broadPeak \
# 3_3502DukeNus_TS543-400-031117_hs_i11_peaks.broadPeak \
# 4_3502DukeNus_TS543-NT-241117_hs_i12_peaks.broadPeak \
# 5_3502DukeNus_TS543-143-241117_hs_i13_peaks.broadPeak \
# 6_3502DukeNus_TS543-400-241117_hs_i14_peaks.broadPeak | \
# sort -k1,1 -k2,2n |bedtools merge -i - > dukenus_merged_peaks.bed



bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

#
require(csaw)

param <- readParam(minq=10, pe='both')

regions=bed_to_granges("~/ayako/ayako_dejavu/peakcall/CD41_merged_peaks.bed")


counts <- regionCounts(bam.files, regions, param=param)
##
#
library(Rsubread)

x=read.table('~/ayako/ayako_dejavu/peakcall/CD41_merged_peaks.bed',sep="\t",stringsAsFactors=F)

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c('/root/ayako/ayako_dejavu/bam/CD41+_untr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_3_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/cd41-_untreated/bam/CD41-_untr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/cd41-_untreated/bam/CD41-_untr_2_Aligned_rmdup.sortedByCoord.out.bam')




fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)

#
##

countData=fc_SE$counts

colnames(countData)=gsub('X.root.ayako.ayako_dejavu.bam.',"",colnames(countData))

colnames(countData)=gsub('_Aligned_rmdup.sortedByCoord.out.bam',"",colnames(countData))

saveRDS(countData,'atac_countdata.rds')

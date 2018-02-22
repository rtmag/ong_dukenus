
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

regions=bed_to_granges("/root/ong_dukenus/peakcalls/dukenus_merged_peaks.bed")


counts <- regionCounts(bam.files, regions, param=param)
##
#
library(Rsubread)

x=read.table('/root/ong_dukenus/peakcalls/dukenus_merged_peaks.bed',sep="\t",stringsAsFactors=F)

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c("/root/ong_dukenus/bam/1_3502DukeNus_TS543-NT-031117_hs_i9_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/3_3502DukeNus_TS543-400-031117_hs_i11_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/5_3502DukeNus_TS543-143-241117_hs_i13_Aligned.sortedByCoord.rmdup.out.bam",
"/root/ong_dukenus/bam/6_3502DukeNus_TS543-400-241117_hs_i14_Aligned.sortedByCoord.rmdup.out.bam")




fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)

#
##

countData=fc_SE$counts

colnames(countData)=c("1_NT","3_400","4_NT","5_143","6_400")

saveRDS(countData,'atac_countdata.rds')

countData=readRDS('atac_countdata.rds')

require(DESeq2)

colData <- data.frame(group=colnames(countData) )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)


pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=136500,intgroup=c('group'))
dev.off()


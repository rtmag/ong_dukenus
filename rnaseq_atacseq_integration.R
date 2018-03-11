#
hg19=read.table('~/resources/hg19_tss_knownCanonical_noUnasembled.bed',sep="\t",stringsAsFactors=F)

genes_nt=read.table('genes_names_high_expression_shNT.txt',sep="\t",stringsAsFactors=F)
tss_nt = hg19[(which(hg19[,4] %in% genes_nt[,1])),]
tss_nt = tss_nt[!duplicated(tss_nt[,4]),]
write.table(tss_nt,"tss_genes_names_high_expression_shNT.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
            
genes_h2afv=read.table('genes_names_high_expression_shH2AFV.txt',sep="\t",stringsAsFactors=F)
tss_h2afv = hg19[(which(hg19[,4] %in% genes_h2afv[,1])),]
tss_h2afv = tss_h2afv[!duplicated(tss_h2afv[,4]),]
write.table(tss_h2afv,"tss_genes_names_high_expression_shH2AFV.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

##################################################################################################################
hg19=read.table('~/resources/hg19_tss_knownCanonical_noUnasembled.bed',sep="\t",stringsAsFactors=F)

gene_symbol=read.table('gene_symbol.txt',sep="\t",stringsAsFactors=F)
gene_tss = hg19[(which(hg19[,4] %in% gene_symbol[,1])),]
gene_tss = gene_tss[!duplicated(gene_tss[,4]),]
gene_tss_1kb=data.frame(gene_tss[,1],gene_tss[,2]-500,gene_tss[,3]+500,gene_tss[,4:6])

write.table(gene_tss_1kb,"gene_tss_1kb.bed",sep="\t",quote=F,col.names=F,row.names=F)
#######################################################################################################


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
require(DESeq2)


blacklist=bed_to_granges("../hg19_blacklist.bed")
param <- readParam(minq=10,discard=blacklist,pe="both")
bam.files <- c(
"1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam",
"2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam",
"3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam",
"4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam",
"5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam",
"6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam"
)

binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)

normfacs <- normOffsets(binned)
saveRDS(normfacs,"normfacs.rds")
##
regions=bed_to_granges("../gene_tss_1kb.bed")
counts <- regionCounts(bam.files, regions, ext=0, param=param)
countData=assay(counts)
colnames(countData)=c("shNT_1","shH2AFV#1_1","shH2AFV#2_1","shNT_2","shH2AFV#1_2","shH2AFV#2_2")

saveRDS(countData,"atac_counts_gene_tss_1kb.rds")

colData <- data.frame(group=colnames(countData))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- normfacs

saveRDS(dds,"dds.rds")

dds_vsd = varianceStabilizingTransformation(dds,fitType='local')
vsd = assay(dds_vsd)
atac_log2fc = rowMeans(vsd[,c(2,3,5,6)]) - rowMeans(vsd[,c(1,4)])

names(atac_log2fc) = regions$id

saveRDS(atac_log2fc,"atac_log2fc_shH2AFV_vs_shNT.rds")

#######################################################################################################
rna = read.csv("../4_out_tables/deseq2_results.csv")
atac = readRDS("atac_log2fc_shH2AFV_vs_shNT.rds")

 rna[rna$gene_symbol %in% names(atac),]

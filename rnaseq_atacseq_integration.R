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


blacklist=bed_to_granges("~/resources/hg19consensusBlacklist.bed")
param <- readParam(minq=10,discard=blacklist,pe="both")
bam.files <- c(
"/root/ong_dukenus/paul_bam/1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam",
"/root/ong_dukenus/paul_bam/2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam",
"/root/ong_dukenus/paul_bam/3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam",
"/root/ong_dukenus/paul_bam/4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam",
"/root/ong_dukenus/paul_bam/5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam",
"/root/ong_dukenus/paul_bam/6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam"
)

binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)

normfacs <- normOffsets(binned)
saveRDS(normfacs,"normfacs.rds")
##
regions=bed_to_granges("/root/ong_dukenus/rnaseq_atacseq/gene_tss_1kb.bed")
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
atac = readRDS("atac_log2fc_shH2AFV_vs_shNT.rds")
atac = atac[order(names(atac))]

rna = read.csv("../2_log2fc_1/2_3_out_tables/deseq2_results.csv")
rna = rna[rna$gene_symbol %in% names(atac),]
rna = rna[!duplicated(rna[,8]),]
rna = rna[order(rna[,8]),]

library(graphics)
pdf("rna_atac_tss_scatterplot_thresholds.pdf")
smoothScatter(rna$log2FoldChange, atac,ylim=c(-3,3),xlim=c(-3,3),
             xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),
             ylab=expression('ATAC-Seq Log'[2]*' Fold Change TSS ( shH2AFV / shNT )'))
abline(v=0)
abline(h=0)

length(which(rna$log2FoldChange<(-0.6) & atac<(-0.6)))
lines(c(-3,-.6),c(-.6,-.6),col="red",lty=2)
lines(c(-3,-.6),c(-3,-3),col="red",lty=2)
lines(c(-.6,-.6),c(-3,-.6),col="red",lty=2)
lines(c(-3,-3),c(-3,-.6),col="red",lty=2)

legend(-2.7,-2.3, paste("Genes:",length(which(rna$log2FoldChange<(-0.6) & atac<(-0.6))) ), bty="n") 

dev.off()

selected_genes = (cbind(rna[which(rna$log2FoldChange<(-0.6) & atac<(-0.6)),],atac[rna$log2FoldChange<(-0.6) & atac<(-0.6)]))

selected_genes = selected_genes[,c(1,8,3,9,7)]
colnames(selected_genes) = c("esemblid","gene_symbol","RNA_log2FC","ATAC_log2FC","PADJ")

write.table(selected_genes,"selected_genes_3rdquadrant_log2FC.6.txt",quote=F,row.names=F,col.names=T,sep="\t")

#######################################################################################################
# MA PLOTS
require(DESeq2)
options(scipen=999)

# MA plot ATAC
countData = readRDS("atac_counts_gene_tss_1kb.rds")
normfacs = readRDS("normfacs.rds")
atac = readRDS("atac_log2fc_shH2AFV_vs_shNT.rds")

rownames(countData) = names(atac)

colData <- data.frame(group=c("shNT","shH2AFV","shH2AFV","shNT","shH2AFV","shH2AFV") )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- normfacs
dds <- DESeq(dds)
#dds_vsd <- varianceStabilizingTransformation(dds)
dds_res <- results(dds,contrast=c("group","shH2AFV","shNT"))
pdf("maplot_atac_tss.pdf")
plotMA(dds_res)
dev.off()

# MA plot RNA
rna = read.csv("deseq2_results.csv",row.names=1)
#pdf("maplot_rna.pdf")
#plotMA(data.frame(rna$baseMean,rna$log2FoldChange,rna$padj<0.1))
#abline(h=0,col="red")
#dev.off()
#
# MA plot ATAC integration RNA-seq
#dim(rna) #58243
#dim(dds_res) #13,172

ix = match(rownames(dds_res), rna[,7])

 library(RColorBrewer)
colors <- colorRampPalette(c("blue","red"))(13172)

col_rna_log2fc <- colors[as.numeric(cut(rna[ix,2],breaks = 257))]

rbPal <- colorRampPalette(c('blue','grey','red'))

col_rna_log2fc <- rbPal(200)[as.numeric(cut(rna[ix,2],breaks = 200))]

pdf("maplot_integration.pdf")
plot(dds_res$baseMean, dds_res$log2FoldChange, ylim=c(-2.5,2.5), pch = 20, log="x", col=col_rna_log2fc,
    xlab="BaseMeans",ylab="Log2FC (shH2AFV / NT)")
dev.off()

pdf("maplot_integration_up_.6.pdf")
plot(dds_res$baseMean[rna[ix,2]>(.6)], dds_res$log2FoldChange[rna[ix,2]>(.6)], ylim=c(-2.5,2.5), pch = 20, log="x", col=col_rna_log2fc[rna[ix,2]>(.6)],
    xlab="BaseMeans",ylab="Log2FC (shH2AFV / NT)")
dev.off()

pdf("maplot_integration_down.6.pdf")
plot(dds_res$baseMean[rna[ix,2]<(-.6)], dds_res$log2FoldChange[rna[ix,2]<(-.6)], ylim=c(-2.5,2.5), pch = 20, log="x", col=col_rna_log2fc[rna[ix,2]<(-.6)],
    xlab="BaseMeans",ylab="Log2FC (shH2AFV / NT)")
dev.off()

library(ggplot2)
dat <- data.frame(Normalized_Mean_Counts = log10(dds_res$baseMean),ATAC_log2FC_shH2AFV_vs_NT = dds_res$log2FoldChange )
RNASeq_log2FC=rna[ix,2]
pdf("maplot_ggplot_white.pdf")
qplot(Normalized_Mean_Counts, ATAC_log2FC_shH2AFV_vs_NT, data=dat, colour=RNASeq_log2FC) + scale_colour_gradientn(colours=c("blue","white","red"))

dev.off()


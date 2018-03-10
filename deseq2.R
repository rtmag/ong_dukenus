countData = read.table(pipe("grep -v '__' 202_143_400_gene_counts.gff"),sep=" ",row.names=1 )
colnames(countData) = c("202_1","202_2","202_3","143_1","143_2","143_3","400_1","400_2","400_3")

options(scipen=999)
library(DESeq2)

#############################################################################
colData <- data.frame(condition = c("shNT","shNT","shNT","shH2AFV#1","shH2AFV#1","shH2AFV#1","shH2AFV#2","shH2AFV#2","shH2AFV#2") )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ condition)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)

pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=60000,intgroup=c('condition'))
dev.off()
#############################################################################
# shH2AFV VS shNT
design<-data.frame(experiment=colnames(countData), 
                   sh = c("shNT","shNT","shNT","shH2AFV","shH2AFV","shH2AFV","shH2AFV","shH2AFV","shH2AFV") )

dds <- DESeqDataSetFromMatrix(countData = countData, colData = design, design = ~ sh  )
dds <- DESeq(dds)
dds_vsd <- varianceStabilizingTransformation(dds)
dds_res <- results(dds,contrast=c("sh","shH2AFV","shNT"))

sig_vsd=dds_vsd[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
sig_vsd = assay(sig_vsd)
colnames(sig_vsd) = c("shNT","shNT","shNT","shH2AFV#1","shH2AFV#1","shH2AFV#1","shH2AFV#2","shH2AFV#2","shH2AFV#2")
#############################################################################
library(gplots)
library(factoextra)
      
 library(RColorBrewer)
colors <- colorRampPalette(c("blue","white","red"))(45)

pdf("heatmap_differentially_expressed_genes.pdf")
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="spearman"),srtCol=25,
          labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression")
dev.off()

pdf("volcano_plot_differentially_expressed_genes.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( shH2AFV / shNT )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red")
legend("topright", paste("shH2AFV:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("shNT:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 
dev.off()

write.csv(dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),],"differentially_expressed_genes.csv")

up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
#############################################################################
# Annotation

library(ensembldb)
library("EnsDb.Hsapiens.v86")   
edb <- EnsDb.Hsapiens.v86
tx <- genes(edb, columns=c("gene_id", "gene_name"))
tx=data.frame(gene_id=tx$gene_id, gene_name=tx$gene_name)
#############################################################################
rankedlist[["ensid"]] <- tx[ match(rankedlist[['ensid']], tx[['gene_id']] ) , 'gene_name']

write.table(rankedlist,"genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)

sig_res = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),c(2,6)]
sig_res = sig_res[order(sig_res$log2FoldChange),]
sig_res = data.frame(sig_res,gene_symbol=rownames(sig_res))
sig_res[["gene_symbol"]] <- tx[ match(sig_res[['gene_symbol']], tx[['gene_id']] ) , 'gene_name']
sig_res = data.frame(sig_res,ensembl_id=rownames(sig_res))
write.table(sig_res,"genes_res_FC_FDR.txt", sep="\t", quote=F,col.names=F,row.names=F)

dds_res = data.frame(dds_res,gene_symbol=rownames(dds_res))
dds_res[["gene_symbol"]] <- tx[ match(dds_res[['gene_symbol']], tx[['gene_id']] ) , 'gene_name']

write.csv(dds_res[order(dds_res$log2FoldChange),],"deseq2_results.csv")
#############################################################################
## qpcr
pdf("QC_genes.pdf", width=15, height=6)

vsd = assay(dds_vsd)
rownames(vsd) = tx[ match(rownames(vsd), tx[['gene_id']] ) , 'gene_name']

notch=dds_res[grep("NOTCH",rownames(vsd)),2]
names(notch) = as.character(rownames(vsd[grep("NOTCH",rownames(vsd)),]))
notch = notch[!is.na(notch)]

par(mar=c(5.1,4.1,4.1,2.1))
barplot(notch,ylim=c(-1.5,1.5),col="#bae1ff",ylab="Log2 Fold Change")
abline(h=0)

barplot(t(vsd[rownames(vsd) %in% names(notch),]),beside=T,ylim=c(0,14),ylab="Log2 Normalized Counts",
                                                                       col=c("#ffb3ba","#ffb3ba","#ffb3ba", 
                                                                            "#baffc9","#baffc9","#baffc9",
                                                                            "#bae1ff","#bae1ff","#bae1ff"))
legend("topright",c("202","143","400"),fill=c("#ffb3ba","#baffc9","#bae1ff") )
abline(h=0)

###

h2a=dds_res[grep("h2a",rownames(vsd),ignore.case=T),2]
names(h2a) = as.character( rownames( vsd[grep("h2a",rownames(vsd),ignore.case=T),]) )
h2a = h2a[!is.na(h2a)]
h2a = h2a[abs(h2a)>.2]

# bottom, left, top and right margins respectively (5.1,4.1,4.1,2.1
par(mar=c(5.1,4.1,4.1,2.1))
barplot(h2a,ylim=c(-1.5,1.5),col="#966FD6",ylab="Log2 Fold Change")
abline(h=0)

barplot(t(vsd[rownames(vsd) %in% names(h2a),]),beside=T,ylim=c(0,15),ylab="Log2 Normalized Counts",col=c("#ffb3ba","#ffb3ba","#ffb3ba",
                                                                            "#baffc9","#baffc9","#baffc9",
                                                                            "#bae1ff","#bae1ff","#bae1ff"))
legend("topright",c("202","143","400"),fill=c("#ffb3ba","#baffc9","#bae1ff") )
abline(h=0)
dev.off()

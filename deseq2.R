countData = read.table(pipe("grep -v '__' 202_143_400_gene_counts.gff"),sep=" ",row.names=1 )
colnames(countData) = c("202_1","202_2","202_3","143_1","143_2","143_3","400_1","400_2","400_3")

options(scipen=999)
library(DESeq2)


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

library(gplots)
library(factoextra)
      
 library(RColorBrewer)
colors <- colorRampPalette(c("blue","white","red"))(45)

pdf("heatmap_differentially_expressed_genes.pdf")
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="spearman"),srtCol=25,
          labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression")
dev.off()

pdf(volcano_plot_differentially_expressed_genes.pdf)
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


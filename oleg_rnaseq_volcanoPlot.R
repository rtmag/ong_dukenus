data=x[,c('PPEE','PostFC','X')]
data = data[data$X!="N/A",]
data = data[complete.cases(data),]
data[is.na(data[,1]),1] = 1

data[data[,1]==0,1] = 1e-16
data[,2] = log2(data[,2])

dds_res = data.frame(log2FoldChange=data[,2],padj=data[,1])

library(graphics)
library(gplots)
library(ggplot2)
pdf("oleg_rnaseq_volcanoPlot.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( shNT / shH2AFV )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5))
abline(v=-.5,lty = 2,col="grey")
abline(v=.5,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>.6 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>.6 & dds_res$padj<0.05],
      col="red")
legend("topright", paste("shNT:",length(which(dds_res$log2FoldChange>.6 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("shH2AFV:",length(which(dds_res$log2FoldChange<(-.6) & dds_res$padj<0.05))), bty="n") 
dev.off()

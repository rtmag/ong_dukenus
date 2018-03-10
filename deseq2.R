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

#pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=60000,intgroup=c('condition'))
#dev.off()



design<-data.frame(experiment=colnames(countData), sh = c("r1","r1","r1","r2","r2","r2","r3","r3","r3"),
                                            condition = c("NT","NT","NT","SH","SH","SH","SH","SH","SH") )

dLRT <- DESeqDataSetFromMatrix(countData = countData, colData = design, design = ~ sh + condition )
dLRT <- DESeq(dLRT, test="LRT",full= ~ sh + condition , reduced=~ sh )
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dDif_res <- results(dLRT,contrast=c("condition","SH","NT"))


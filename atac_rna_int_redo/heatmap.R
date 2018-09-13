options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

x = read.csv("Master_tab_o_202_143_400_R.csv")
vsd=x[,16:24]
colnames(vsd) = c("shNT 1","shNT 2","shNT 3",
"shH2AFV #I 1","shH2AFV #I 2","shH2AFV #I 3",
"shH2AFV #II 1","shH2AFV #II 2","shH2AFV #II 3")
rownames(vsd) = make.unique(as.character(x[,4]), sep = "_") 

pdf("dreamComplex.pdf")
sig_vsd = vsd[rownames(vsd) %in% c("LIN9","LIN54","LIN52","LIN37","FOXM1","MYB","MYBL2","KCNIP3","RB1","KIF2A","E2F2",
                                    "DYRK1A","TP53","CDKN1A","E2F4","RBBP4","TFDP1","RBL2","RBL1","E2F5","TFDP2","MYBL1",
                                    "NOLC1","FLNA","BIRC5","PLK1","CDC25C","CLTC","LATS2","SPTA1","MEIS2","SPTBN1","CLTA",
                                    "FAM111A","NEIL3","E2F1","E2F3","E2F6","E2F7","E2F8"),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="DREAM complex Genes",key.title="Gene expression",cexCol=.65,cexRow=.6)
dev.off()

pdf("dreamComplex_rowSCALED.pdf")
sig_vsd = vsd[rownames(vsd) %in% c("LIN9","LIN54","LIN52","LIN37","FOXM1","MYB","MYBL2","KCNIP3","RB1","KIF2A","E2F2",
                                    "DYRK1A","TP53","CDKN1A","E2F4","RBBP4","TFDP1","RBL2","RBL1","E2F5","TFDP2","MYBL1",
                                    "NOLC1","FLNA","BIRC5","PLK1","CDC25C","CLTC","LATS2","SPTA1","MEIS2","SPTBN1","CLTA",
                                    "FAM111A","NEIL3","E2F1","E2F3","E2F6","E2F7","E2F8"),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="DREAM complex Genes",key.title="Gene expression",cexCol=.65,cexRow=.6)
dev.off()

gmt = read.table("HallMark_e2f_targets.txt",sep = "\t", head=T,stringsAsFactors=F)

pdf("e2f_targets.pdf")
sig_vsd = vsd[(rownames(vsd) %in% gmt[,1]),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="E2F Target Genes",key.title="Gene expression",cexCol=.65,cexRow=.2)
dev.off()

pdf("e2f_targets_rowSCALED.pdf")
sig_vsd = vsd[(rownames(vsd) %in% gmt[,1]),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="E2F Target Genes",key.title="Gene expression",cexCol=.65,cexRow=.2)
dev.off()

pdf("TF.pdf")
sig_vsd = vsd[rownames(vsd) %in% c("E2F4", "CEBPD", "FOXM1", "REST", "RBPJ"),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="TF Genes",key.title="Gene expression",cexCol=.65,cexRow=.6)
dev.off()

pdf("TF_rowSCALED.pdf")
sig_vsd = vsd[rownames(vsd) %in% c("E2F4", "CEBPD", "FOXM1", "REST", "RBPJ"),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="TF Genes",key.title="Gene expression",cexCol=.65,cexRow=.6)
dev.off()

##############################################################################################################################
atac_down = read.table(pipe("cut -f 16 /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno"),sep="\t",header=T)
rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")
gmt = read.table("HallMark_e2f_targets.txt",sep = "\t", head=T,stringsAsFactors=F)


sig_vsd = vsd[which( (rna[,1] %in% atac_down[,1]) & rna$PPEE<0.05 & abs(log2(rna$PostFC))>.37 & (rna[,1] %in% gmt[,1])),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=1, key=FALSE)
dev.off()
##############################################################################################################################
atac_down = read.table(pipe("cut -f 16 /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno"),sep="\t",header=T)
rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")
gmt = read.table("~/Downloads/CEBPD_targets.human.tsv",sep = "\t", head=T,stringsAsFactors=F)

sig_vsd = vsd[which( rna$PPEE<0.05 & (rna[,1] %in% gmt[,2])),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("CEBPD_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=1, key=FALSE)
dev.off()
##############################################################################################################################
atac_down = read.table(pipe("cut -f 16 /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno"),sep="\t",header=T)
rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")
gmt = read.table("~/Downloads/REST_targets.human.tsv",sep = "\t", head=T,stringsAsFactors=F)

sig_vsd = vsd[which( rna$PPEE<0.05 & (rna[,1] %in% gmt[,2])),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("REST_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=1, key=FALSE)
dev.off()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# prediction
more /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno |cut -f2,3,4,16|grep -v "Gene Name"| \
bedtools intersect -a - -b /Users/wone/Downloads/CEBPD.bed|cut -f4 > putative_cebpd_targets.txt

more /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno |cut -f2,3,4,16|grep -v "Gene Name"| \
bedtools intersect -a - -b /Users/wone/Downloads/REST.bed|cut -f4 > putative_rest_targets.txt
# CEBPD
atac_down = read.table(pipe("cut -f 16 /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno"),sep="\t",header=T)
rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")
gmt = read.table("putative_cebpd_targets.txt",sep = "\t", head=T,stringsAsFactors=F)


sig_vsd = vsd[which( (rna[,1] %in% atac_down[,1]) & rna$PPEE<0.05 & (-log2(rna$PostFC))<(-.37) & (rna[,1] %in% gmt[,1])),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("CEBPD_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=.4, key=FALSE,Colv="none")
dev.off()
# REST
atac_down = read.table(pipe("cut -f 16 /Users/wone/CSI/ong/atacseq_redo/variantWindows/anno/ATAC-down_promoter.anno"),sep="\t",header=T)
rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")
gmt = read.table("putative_rest_targets.txt",sep = "\t", head=T,stringsAsFactors=F)


sig_vsd = vsd[which( (rna[,1] %in% atac_down[,1]) & rna$PPEE<0.05 & (-log2(rna$PostFC))>(.37) & (rna[,1] %in% gmt[,1])),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("REST_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=.4, key=FALSE,Colv="none")
dev.off()

##########################
# CEBPD

gmt = rbind(
"FASN",
"HMGB2",
"VDAC1",
"PTTG1",
"PYGL",
"SOCS2",
"SOD2",
"ACTN1",
"CTNNAL1",
"JMJD6")


sig_vsd = vsd[which( rownames(vsd) %in% gmt[,1]),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("CEBPD_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=.6, key=FALSE,Colv="none")
dev.off()
# REST

gmt = rbind('MAP2',
'ARPC3',
'MXD4',
'MTSS1',
'HOMER1',
'PBX3',
'DUSP16',
'TULP3')

sig_vsd = vsd[which( rownames(vsd) %in% gmt[,1]),]
sig_vsd = sig_vsd[complete.cases(sig_vsd),]
sig_vsd = sig_vsd/rowMeans(sig_vsd)
sig_vsd = sig_vsd[complete.cases(sig_vsd),]

pdf("REST_TARGET_GENES.pdf")
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(as.matrix(sig_vsd),dendrogram='none',col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="",key.title="Gene expression",cexCol=.6,cexRow=.4, key=FALSE,Colv="none")
dev.off()

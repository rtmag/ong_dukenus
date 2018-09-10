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

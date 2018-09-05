# MATCH ID
x = read.csv("/Users/wone/CSI/ong/RNA_SEQ_analysis/2_log2fc_1/2_3_out_tables/deseq2_results.csv")
vsd = readRDS("ong_vsd.rds")
vsd = assay(vsd)
table(rownames(vsd[match(x[,1],rownames(vsd)),])==x[,1])
vsd = vsd[match(x[,1],rownames(vsd)),]
rownames(vsd) = x$gene_symbol
colnames(vsd) = c("shNT","shNT","shNT","shH2AFV#1","shH2AFV#1","shH2AFV#1","shH2AFV#2","shH2AFV#2","shH2AFV#2")
saveRDS(vsd,"ong_vsd.rds")
# END

options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

vsd = readRDS("ong_vsd.rds")
sox2 = read.table("~/Downloads/SOX2_targets.human.tsv",sep = "\t", head=F,stringsAsFactors=F)

pdf("sox2_targets.pdf")
sig_vsd = vsd[rownames(vsd) %in% sox2[,2],]
ix = apply(sig_vsd,1,sd)!=0
sig_vsd = sig_vsd[ix,]
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  
track = sox2[match(rownames(sig_vsd),sox2[,2]),3]
track[track=="Activation"]=1
track[track=="Repression"]=2
track[track=="Unknown"]=3
track=as.numeric(track)
                               
colores=c("#30d24b","#d3303f","#606360")
  heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="SOX2 Target Genes",key.title="Gene expression",cexCol=.65,cexRow=.65,RowSideColors=colores[track])
  legend("topright",c("Activation","Repression","Unkown"),fill=c("#30d24b","#d3303f","#606360"), border=T, bty="n")

dev.off()

sox2 = read.table("~/Downloads/SOX10_targets.human.tsv",sep = "\t", head=F,stringsAsFactors=F)

pdf("sox10_targets.pdf")
sig_vsd = vsd[rownames(vsd) %in% sox2[,2],]
ix = apply(sig_vsd,1,sd)!=0
sig_vsd = sig_vsd[ix,]
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  
track = sox2[match(rownames(sig_vsd),sox2[,2]),3]
track[track=="Activation"]=1
track[track=="Repression"]=2
track[track=="Unknown"]=3
track=as.numeric(track)
                               
colores=c("#30d24b","#d3303f","#606360")
  heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
  xlab="", ylab="SOX10 Target Genes",key.title="Gene expression",cexCol=.65,cexRow=.65,RowSideColors=colores[track])
  legend("topright",c("Activation","Repression","Unkown"),fill=c("#30d24b","#d3303f","#606360"), border=T, bty="n")

dev.off()

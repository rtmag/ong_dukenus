options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

signature <- read.csv("/home/rtm/CSI/ong/sulconazole/3_pathway_signature.csv",stringsAsFactors=FALSE)
signature[signature$Geneset_for_CMAP_analysis=="KOBAYASHI_EGFR_SIGNALING_24HR_DN",2] <- "EGFR"
signature[signature$Geneset_for_CMAP_analysis=="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",2] <- "Mesenchymal"
signature[signature$Geneset_for_CMAP_analysis=="generation of neurons",2] <- "GenerationOfNeurons"


GBMLGG <- read.table("/home/rtm/CSI/ong/sulconazole/public/2019-08-21_TCGA_GBMLGG_expression.txt", 
                     header=TRUE, row.names=1)
GBMLGG <- t(GBMLGG)
GBMLGG_sig <- GBMLGG[rownames(GBMLGG) %in% signature$Gene.symbol,]

track= signature[match(rownames(GBMLGG_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

GBMLGG.pheno <- read.table("/home/rtm/CSI/ong/sulconazole/public/2019-08-21_TCGA_GBMLGG_pheno.txt", 
                     header=TRUE, row.names=1)

ctrack = as.character(GBMLGG.pheno$Histology)
ctrack[ctrack=="Astrocytoma"]=1
ctrack[ctrack=="GBM"]=2
ctrack[ctrack=="Oligoastrocytoma"]=3
ctrack[ctrack=="Oligodendroglioma"]=4
colores=c("grey","black","#ffb347","#966fd6")
clab=as.character(colores[as.numeric(ctrack)])

png("heatmap_diff_TUMOR_VS_NORMAL_FDR5p.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)

GBMLGG_sig_centered = GBMLGG_sig - rowMeans(GBMLGG_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered > 7] = 7
GBMLGG_sig_centered[GBMLGG_sig_centered < (-7)] = -7

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(29))

heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="GBM-LGG Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)

legend("topright",legend=c("EGFR","Mesenchymal","GenerationOfNeurons",""),
       fill=c("#ffb3ba","#baffc9","#bae1ff","white"), border=T, bty="n" )

dev.off()

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

png("heatmap_TCGA_GBMLG_Signature.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)

GBMLGG_sig_centered = GBMLGG_sig - rowMeans(GBMLGG_sig)

GBMLGG_sig_centered[GBMLGG_sig_centered > 6] = 6
GBMLGG_sig_centered[GBMLGG_sig_centered < (-6)] = -6

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(11))

x=heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="GBM-LGG Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()

pdf("heatmap_TCGA_GBMLG_labels.pdf")
plot.new()
legend("center",legend=c("EGFR","Mesenchymal","GenerationOfNeurons",
                           "Astrocytoma","GBM","Oligoastrocytoma","Oligodendroglioma"),
       fill=c("#ffb3ba","#baffc9","#bae1ff","grey","black","#ffb347","#966fd6"), border=T, bty="n" )

dev.off()

hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=4 )

track=as.numeric(groups)
colores=c("red","blue","green","orange")
clab=(colores[track])

png("heatmap_TCGA_GBMLG_Signature_K4Cut.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)

heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="GBM-LGG Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()

library(RTCGA.clinical)

clinical <- data.frame(times = GBMLGG.pheno$survival,
                       bcr_patient_barcode = rownames(GBMLGG.pheno),
                       patient.vital_status = as.numeric(GBMLGG.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1
clinical$bcr_patient_barcode == colnames(GBMLGG_sig_centered)

pdf("survival_GBMLG_Signature_K4Cut_pval.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=FALSE,palette = c("red","blue","green","orange"))
dev.off()

pdf("survival_GBMLG_Signature_K4Cut.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=FALSE,palette = c("red","blue","green","orange"))
dev.off()

pdf("survival_GBMLG_Signature_K4Cut_pval_table.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=TRUE,palette = c("red","blue","green","orange"))
dev.off()

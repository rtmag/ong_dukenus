options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

# READ AND PARSE SIGNATURES
signature <- read.csv("3_pathway_signature.csv",stringsAsFactors=FALSE)
signature[signature$Geneset_for_CMAP_analysis=="KOBAYASHI_EGFR_SIGNALING_24HR_DN",2] <- "EGFR"
signature[signature$Geneset_for_CMAP_analysis=="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",2] <- "Mesenchymal"
signature[signature$Geneset_for_CMAP_analysis=="generation of neurons",2] <- "GenerationOfNeurons"

#READ EXPRESSION TABLE
GBMLGG <- read.table("2019-08-21_TCGA_GBMLGG_expression.txt", header=TRUE, row.names=1)
GBMLGG <- t(GBMLGG)

#READ PHENO TABLE
GBMLGG.pheno <- read.table("2019-08-21_TCGA_GBMLGG_pheno.txt", header=TRUE, row.names=1)

#CHECK GRADE OF ASTRO AND GBM
ctrack = data.frame(hist=as.character(GBMLGG.pheno$Histology),grade=as.character(GBMLGG.pheno$Grade))
table(ctrack[ctrack$hist=="Astrocytoma",2])
table(ctrack[ctrack$hist=="GBM",2])

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(GBMLGG.pheno$Histology)!="Oligoastrocytoma" & 
     as.character(GBMLGG.pheno$Histology)!="Oligodendroglioma" & 
     !is.na(GBMLGG.pheno$Grade) & !is.na(GBMLGG.pheno$Histology)
GBMLGG_astro_gbm <- GBMLGG[,ix]
GBMLGG_astro_gbm.pheno <- GBMLGG.pheno[ix,]

#GET SIGNATURE GENES ONLY
GBMLGG_astro_gbm_sig <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% signature$Gene.symbol,]

track= signature[match(rownames(GBMLGG_astro_gbm_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

#create label track with grade
ctrack = as.character(GBMLGG_astro_gbm.pheno$Grade)
ctrack[ctrack=="II"]=1
ctrack[ctrack=="III"]=2
ctrack[ctrack=="IV"]=3
colores=c("black","grey","red")
clab=as.character(colores[as.numeric(ctrack)])

#heatmap
png("heatmap_TCGA_GBM_ASTRO_Signature.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
GBMLGG_sig_centered = GBMLGG_astro_gbm_sig - rowMeans(GBMLGG_astro_gbm_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 6] = 6
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-6)] = -6

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(11))


x=heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="TCGA GBM-LGG Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()

hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=3 )

track=as.numeric(groups)
colores=c("purple","orange","blue")
clab=(colores[track])

png("heatmap_TCGA_GBM_ASTRO_Signature_K3Cut.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="TCGA GBM-LGG Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()

pdf("heatmap_TCGA_GBMLG_GBM_ASTRO_labels.pdf")
plot.new()
legend("center",legend=c("EGFR","Mesenchymal","GenerationOfNeurons",
                           "Grade II","Grade III","Grade IV"),
       fill=c("#ffb3ba","#baffc9","#bae1ff","black","grey","red"), border=T, bty="n" )
dev.off()


library(RTCGA.clinical)

clinical <- data.frame(times = GBMLGG_astro_gbm.pheno$survival,
                       bcr_patient_barcode = rownames(GBMLGG_astro_gbm.pheno),
                       patient.vital_status = as.numeric(GBMLGG_astro_gbm.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

pdf("survival_GBM_ASTRO_GBMLG_Signature_K3Cut_pval.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue"))
dev.off()

pdf("survival_GBM_ASTRO_GBMLG_Signature_K3Cut.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
## ONLY GBM

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
GBMLGG_sig_centered = GBMLGG - rowMeans(GBMLGG)

ix = as.character(GBMLGG.pheno$Histology)!="Oligoastrocytoma" & 
     as.character(GBMLGG.pheno$Histology)!="Oligodendroglioma" & 
     as.character(GBMLGG.pheno$Histology)!="Astrocytoma" & 
     !is.na(GBMLGG.pheno$Grade) & !is.na(GBMLGG.pheno$Histology)
GBMLGG_astro_gbm <- GBMLGG_sig_centered[,ix]
GBMLGG_astro_gbm.pheno <- GBMLGG.pheno[ix,]

#GET SIGNATURE GENES ONLY
GBMLGG_astro_gbm_sig <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% signature$Gene.symbol,]

track= signature[match(rownames(GBMLGG_astro_gbm_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

#create label track with grade
ctrack = as.character(GBMLGG_astro_gbm.pheno$Grade)
ctrack[ctrack=="II"]=1
ctrack[ctrack=="III"]=2
ctrack[ctrack=="IV"]=3
colores=c("black","grey","red")
clab=as.character(colores[as.numeric(ctrack)])

#heatmap
png("heatmap_TCGA_GBMonly_Signature.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
GBMLGG_astro_gbm_sig <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% signature$Gene.symbol,]
GBMLGG_astro_gbm_sig[GBMLGG_astro_gbm_sig >= 9] = 9
GBMLGG_astro_gbm_sig[GBMLGG_astro_gbm_sig <= (-9)] = -9

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(11))


x=heatmap.2(as.matrix(GBMLGG_astro_gbm_sig),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="TCGA GBM Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()


hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=3 )

track=as.numeric(groups)
colores=c("purple","orange","blue")
clab=(colores[track])


png("heatmap_TCGA_GBMonly_Signature_K3Cut.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(GBMLGG_astro_gbm_sig),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="TCGA GBM Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()

library(RTCGA.clinical)

clinical <- data.frame(times = GBMLGG_astro_gbm.pheno$survival,
                       bcr_patient_barcode = rownames(GBMLGG_astro_gbm.pheno),
                       patient.vital_status = as.numeric(GBMLGG_astro_gbm.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

pdf("survival_GBMonly_GBMLG_Signature_K3Cut_pval.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue"))
dev.off()

pdf("survival_GBMonly_GBMLG_Signature_K3Cut.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue"))
dev.off()


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
## REMB


remb <- read.table("2019-08-21_Rembrandt_expression.txt", 
                     header=TRUE, row.names=1)
remb <- t(remb)

remb.pheno <- read.table("2019-08-21_Rembrandt_pheno.txt", 
                     header=TRUE, row.names=1)

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(remb.pheno$Histology)!="Mixed glioma" & 
     as.character(remb.pheno$Histology)!="Oligodendroglioma" & 
     as.character(remb.pheno$Histology)!="Non-tumor" & 
     as.character(remb.pheno$Histology)!="Unknown" & 
     as.character(remb.pheno$Grade)!="I" & 
     !is.na(remb.pheno$Grade) & !is.na(remb.pheno$Histology)
remb_imp <- remb[,ix]
remb_imp.pheno <- remb.pheno[ix,]

#GET SIGNATURE GENES ONLY
remb_imp_sig <- remb_imp[rownames(remb_imp) %in% signature$Gene.symbol,]

track= signature[match(rownames(remb_imp_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

#create label track with grade
ctrack = as.character(remb_imp.pheno$Grade)
ctrack[ctrack=="II"]=1
ctrack[ctrack=="III"]=2
ctrack[ctrack=="IV"]=3
colores=c("black","grey","red")
clab=as.character(colores[as.numeric(ctrack)])

#heatmap
png("heatmap_Rembrandt_GBM_ASTRO_Signature.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
GBMLGG_sig_centered = remb_imp_sig - rowMeans(remb_imp_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 4] = 4
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-4)] = -4

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(11))
x=heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="Rembrandt Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()


hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=4 )

track=as.numeric(groups)
colores=c("purple","orange","blue","darkgreen")
clab=(colores[track])

png("heatmap_Rembrandt_GBM_ASTRO_Signature_K4Cut.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="Rembrandt Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()


library(RTCGA.clinical)

clinical <- data.frame(times = remb_imp.pheno$survival,
                       bcr_patient_barcode = rownames(remb_imp.pheno),
                       patient.vital_status = as.numeric(remb_imp.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

pdf("survival_GBM_ASTRO_Rembrandt_Signature_K4Cut_pval.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue","darkgreen"))
dev.off()

pdf("survival_GBM_ASTRO_Rembrandt_Signature_K4Cut.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue","darkgreen"))
dev.off()

###
## REMB


remb <- read.table("2019-08-21_Rembrandt_expression.txt", 
                     header=TRUE, row.names=1)
remb <- t(remb)

remb.pheno <- read.table("2019-08-21_Rembrandt_pheno.txt", 
                     header=TRUE, row.names=1)

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(remb.pheno$Histology)!="Mixed glioma" & 
     as.character(remb.pheno$Histology)!="Oligodendroglioma" & 
     as.character(remb.pheno$Histology)!="Non-tumor" & 
     as.character(remb.pheno$Histology)!="Unknown" & 
     as.character(remb.pheno$Grade)!="I" & 
     !is.na(remb.pheno$Grade) & !is.na(remb.pheno$Histology)
remb_imp <- remb[,ix]
remb_imp.pheno <- remb.pheno[ix,]

#GET SIGNATURE GENES ONLY
remb_imp_sig <- remb_imp[rownames(remb_imp) %in% signature$Gene.symbol,]

track= signature[match(rownames(remb_imp_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

#create label track with grade
ctrack = as.character(remb_imp.pheno$Grade)
ctrack[ctrack=="II"]=1
ctrack[ctrack=="III"]=2
ctrack[ctrack=="IV"]=3
colores=c("black","grey","red")
clab=as.character(colores[as.numeric(ctrack)])

#heatmap
png("heatmap_Rembrandt_GBM_ASTRO_Signature.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
GBMLGG_sig_centered = remb_imp_sig - rowMeans(remb_imp_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 4] = 4
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-4)] = -4

colors <- rev(colorRampPalette( (brewer.pal(11, "RdBu")) )(11))
x=heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="Rembrandt Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()


hc <- as.hclust( x$colDendrogram )
groups=cutree( hc, k=4 )

track=as.numeric(groups)
colores=c("purple","orange","blue","darkgreen")
clab=(colores[track])

png("heatmap_Rembrandt_GBM_ASTRO_Signature_K4Cut.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(GBMLGG_sig_centered),col=colors,scale="none", trace="none",
              distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = "",labCol = "",xlab="Rembrandt Patient Sample", ylab="Signature Genes",key.title="",
         RowSideColors=rlab,ColSideColors=clab)
dev.off()


library(RTCGA.clinical)

clinical <- data.frame(times = remb_imp.pheno$survival,
                       bcr_patient_barcode = rownames(remb_imp.pheno),
                       patient.vital_status = as.numeric(remb_imp.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

pdf("survival_GBM_ASTRO_Rembrandt_Signature_K4Cut_pval.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue","darkgreen"))
dev.off()

pdf("survival_GBM_ASTRO_Rembrandt_Signature_K4Cut.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=FALSE,palette = c("purple","orange","blue","darkgreen"))
dev.off()

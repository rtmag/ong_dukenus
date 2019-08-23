options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)

signature <- read.csv("3_pathway_signature.csv",stringsAsFactors=FALSE)
signature[signature$Geneset_for_CMAP_analysis=="KOBAYASHI_EGFR_SIGNALING_24HR_DN",2] <- "EGFR"
signature[signature$Geneset_for_CMAP_analysis=="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",2] <- "Mesenchymal"
signature[signature$Geneset_for_CMAP_analysis=="generation of neurons",2] <- "GenerationOfNeurons"

GBMLGG <- read.table("2019-08-21_TCGA_GBMLGG_expression.txt", header=TRUE, row.names=1)
GBMLGG <- t(GBMLGG)
GBMLGG_sig <- GBMLGG[rownames(GBMLGG) %in% signature$Gene.symbol,]

track= signature[match(rownames(GBMLGG_sig),signature[,1]),2]
track[track=="EGFR"]=1
track[track=="Mesenchymal"]=2
track[track=="GenerationOfNeurons"]=3
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
rlab=as.character(colores[track])

GBMLGG.pheno <- read.table("2019-08-21_TCGA_GBMLGG_pheno.txt", header=TRUE, row.names=1)

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("org.Hs.eg.db")
options(scipen=999)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats) ## for reordering the factor


genes<-read.table("H2AFV_correlatedGenes.txt")
genes<-as.character(genes[,1])

gene.df <- bitr(genes, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Hs.eg.db)

geneEntrez <- list(genes = gene.df$ENTREZID)

x=compareCluster(geneEntrez, fun='enrichGO',
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP")
pdf("dotplot_H2AFV_correlated_GO_biologicalProcess.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun="enrichPathway", organism = "human")

pdf("dotplot_H2AFV_correlated_enrichPathway.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun="enrichKEGG", organism = "human")

pdf("dotplot_H2AFV_correlated_enrichKEGG.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

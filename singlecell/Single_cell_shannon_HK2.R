library(gplots)
library(ggplot2)
library(graphics)
library(scales)
library(vegan)
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
library("xlsx")
library(beeswarm)

# Shannon Index
shannonIndex = function ( box ){
    samples = as.character(unique(box[,1]))
    shannonIndex = rep(0, length(samples))
    names(shannonIndex) = samples
    for( i in 1:length(samples) ){ 
        ix = box[,1] == samples[i]
        shannonIndex[i] = diversity(box[ix,2], index = "shannon", MARGIN = 1, base = exp(1))
        }
    return(shannonIndex)
    }


# MGH's  #OLIG2 - cancer stem cell marker
data = read.table("GSE57872_GBM_data_matrix.txt",sep="\t",header=T,row.names=1)
data = data[,grep("_",colnames(data))]
data = data[,grep("MGH",colnames(data))]
#
H2AFV = data[rownames(data)=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=2^(as.numeric(H2AFV)) )

#H2AFZ = data[rownames(data)=="H2AFZ",]
#H2AFZ = data.frame(cell = gsub("\\_.+","",colnames(H2AFZ),perl=TRUE), gene=2^(as.numeric(H2AFZ)) )

H2AFY = data[rownames(data)=="H2AFY",]
H2AFY = data.frame(cell = gsub("\\_.+","",colnames(H2AFY),perl=TRUE), gene=2^(as.numeric(H2AFY)) )
#
pdf("MGH_jitter.pdf")
par(mfrow = c(2,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

#stripchart(gene ~ cell, vertical = TRUE, data = H2AFZ, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFZ'),
#    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = H2AFY, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFY'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()
#
pdf("MGH_beeSwarm.pdf")
par(mfrow = c(2,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)
#beeswarm(gene ~ cell, vertical = TRUE, data = H2AFZ,method = "swarm",pch = 16,xlab="",
#         ylab = expression('Single Cell RNA-Seq H2AFZ'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFY,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFY'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

####################################################################################################
# Shannon Index
shannonIndex = function ( box ){
    samples = as.character(unique(box[,1]))
    shannonIndex = rep(0, length(samples))
    names(shannonIndex) = samples
    for( i in 1:length(samples) ){ 
        ix = box[,1] == samples[i]
        shannonIndex[i] = diversity(box[ix,2], index = "shannon", MARGIN = 1, base = exp(1))
        }
    return(shannonIndex)
    }

#summary(c(H2AFV[,2],H2AFZ[,2],OLIG1[,2]))

#H2AFV[,2] = H2AFV[,2]+7
#H2AFZ[,2] = H2AFZ[,2]+7
#H2AFY[,2] = H2AFY[,2]+7

shannon_H2AFV = shannonIndex(H2AFV)
#shannon_H2AFZ = shannonIndex(H2AFZ)
shannon_H2AFY = shannonIndex(H2AFY)

sig = rbind(shannon_H2AFV,shannon_H2AFY)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_MGH.pdf")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig, dendrogram = "none",
          cellnote=sig,Rowv = FALSE,
          scale="column",  trace="none", 
          notecex=0.9, distfun = function(x) get_dist(x,method="pearson"),
          notecol="black",col=colors,
          na.color=par("bg"),cexCol=.8, cexRow=.9,key=FALSE,
          sepwidth=c(0.01,0.005),
           sepcolor="white",
           colsep=1:ncol(sig),
           rowsep=1:nrow(sig) )

dev.off()

#################################################################################################################################

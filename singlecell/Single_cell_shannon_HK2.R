library(gplots)
library(ggplot2)
library(graphics)
library(scales)
library(vegan)
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
library(openxlsx)
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
#
# RPL5, RUVBL1, RUVBL2, ACTL6A, YEATS4 .
####################################################################################################

# MGH's  #OLIG2 - cancer stem cell marker
data = read.table("GSE57872_GBM_data_matrix.txt",sep="\t",header=T,row.names=1)
data = data[,grep("_",colnames(data))]
data = data[,grep("MGH",colnames(data))]
#
H2AFV = data[rownames(data)=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=2^(as.numeric(H2AFV)) )

RUVBL1 = data[rownames(data)=="RUVBL1",]
RUVBL1 = data.frame(cell = gsub("\\_.+","",colnames(RUVBL1),perl=TRUE), gene=2^(as.numeric(RUVBL1)) )

RUVBL2 = data[rownames(data)=="RUVBL2",]
RUVBL2 = data.frame(cell = gsub("\\_.+","",colnames(RUVBL2),perl=TRUE), gene=2^(as.numeric(RUVBL2)) )

ACTL6A = data[rownames(data)=="ACTL6A",]
ACTL6A = data.frame(cell = gsub("\\_.+","",colnames(ACTL6A),perl=TRUE), gene=2^(as.numeric(ACTL6A)) )
#
pdf("MGH_jitter.pdf")
par(mfrow = c(4,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL1, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL1'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL2, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL2'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = ACTL6A, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTL6A'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()
#
pdf("MGH_beeSwarm.pdf")
par(mfrow = c(4,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL1,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL1'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL2,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL2'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = ACTL6A,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTL6A'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_RUVBL1 = shannonIndex(RUVBL1)
shannon_RUVBL2 = shannonIndex(RUVBL2)
shannon_ACTL6A = shannonIndex(ACTL6A)

sig = rbind(shannon_H2AFV,shannon_RUVBL1,shannon_RUVBL2,shannon_ACTL6A)
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



# BREAST #ITG86 - cancer stem cell marker

data = read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",sep="\t",header=T)

# BC01-BC02, estrogen receptor positive (ER+); 
# BC03, double positive (ER+ and HER2+); 
# BC04-BC06, human epidermal growth factor receptor 2 positive (HER2+); 
# BC07-BC11, triple-negative breast cancer (TNBC); 
# BC03LN, lymph node metastasis of BC03; 
# BC07LN, lymph node metastasis of BC07)

tpm = data[,18:dim(data)[2]]
rownames(tpm) = make.names(data[,2],unique=T)
colnames(tpm) = gsub("\\_.+","",colnames(tpm),perl=T)

#remove LN
tpm = tpm[,grep("LN",colnames(tpm),invert=T)]
colnames(tpm) = gsub("\\..+","",colnames(tpm),perl=T)


colnames(tpm)[colnames(tpm)=="BC01"] = "ER+ #1"
colnames(tpm)[colnames(tpm)=="BC02"] = "ER+ #2"
colnames(tpm)[colnames(tpm)=="BC03"] = "ER+_HER2+"
colnames(tpm)[colnames(tpm)=="BC03LN"] = "ER+_HER2+_LN"
colnames(tpm)[colnames(tpm)=="BC04"] = "HER2+ #1"
colnames(tpm)[colnames(tpm)=="BC05"] = "HER2+ #2"
colnames(tpm)[colnames(tpm)=="BC06"] = "HER2+ #3"
colnames(tpm)[colnames(tpm)=="BC07"] = "TNBC #1"
colnames(tpm)[colnames(tpm)=="BC07LN"] = "TNBC_LN"
colnames(tpm)[colnames(tpm)=="BC08"] = "TNBC #2"
colnames(tpm)[colnames(tpm)=="BC09"] = "TNBC #3"
colnames(tpm)[colnames(tpm)=="BC10"] = "TNBC #4"
colnames(tpm)[colnames(tpm)=="BC11"] = "TNBC #5"

# RPL5, RUVBL1, RUVBL2, ACTL6A, YEATS4 .

H2AFV = tpm[rownames(tpm)=="H2AFV",]
H2AFV = data.frame(cell = colnames(tpm), gene=(as.numeric(H2AFV)) )

RPL5 = tpm[rownames(tpm)=="RPL5",]
RPL5 = data.frame(cell = colnames(tpm), gene=(as.numeric(RPL5)) )

RUVBL1 = tpm[rownames(tpm)=="RUVBL1",]
RUVBL1 = data.frame(cell = colnames(tpm), gene=(as.numeric(RUVBL1)) )

RUVBL2 = tpm[rownames(tpm)=="RUVBL2",]
RUVBL2 = data.frame(cell = colnames(tpm), gene=(as.numeric(RUVBL2)) )

ACTL6A = tpm[rownames(tpm)=="ACTL6A",]
ACTL6A = data.frame(cell = colnames(tpm), gene=(as.numeric(ACTL6A)) )

YEATS4 = tpm[rownames(tpm)=="YEATS4",]
YEATS4 = data.frame(cell = colnames(tpm), gene=(as.numeric(YEATS4)) )

pdf("breast_jitter.pdf",width=13)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL5, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL5'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL1, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL1'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL2, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL2'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = ACTL6A, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTL6A'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = YEATS4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq YEATS4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("breast_beeSwarm.pdf",width=13)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL5'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL1,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL1'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL2,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL2'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = ACTL6A,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTL6A'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = YEATS4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq YEATS4'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_RPL5 = shannonIndex(RPL5)
shannon_RUVBL1 = shannonIndex(RUVBL1)
shannon_RUVBL2 = shannonIndex(RUVBL2)
shannon_ACTL6A = shannonIndex(ACTL6A)
shannon_YEATS4 = shannonIndex(YEATS4)

sig = rbind(shannon_H2AFV,shannon_RPL5,shannon_RUVBL1,shannon_RUVBL2,shannon_ACTL6A,shannon_YEATS4)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_BREAST.pdf")
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
# MELANOMA
data = read.table("GSE72056_melanoma_single_cell_revised_v2.txt",sep="\t",header=T)
genenames = as.character( data[ 4:dim(data)[1], 1 ] )
cells = as.character(data[1,2:dim(data)[2]])
data = data[4:dim(data)[1], 2:dim(data)[2]]
####
# RPL5, RUVBL1, RUVBL2, ACTL6A, YEATS4 .
H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = cells, gene=(as.numeric(H2AFV)) )

RPL5 = data[genenames=="RPL5",]
RPL5 = data.frame(cell = cells, gene=(as.numeric(RPL5)) )

RUVBL1 = data[genenames=="RUVBL1",]
RUVBL1 = data.frame(cell = cells, gene=(as.numeric(RUVBL1)) )

RUVBL2 = data[genenames=="RUVBL2",]
RUVBL2 = data.frame(cell = cells, gene=(as.numeric(RUVBL2)) )

ACTL6A = data[genenames=="ACTL6A",]
ACTL6A = data.frame(cell = cells, gene=(as.numeric(ACTL6A)) )

YEATS4 = data[genenames=="YEATS4",]
YEATS4 = data.frame(cell = cells, gene=(as.numeric(YEATS4)) )
####
pdf("melanoma_jitter.pdf",width=23)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL5, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL5'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL1, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL1'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL2, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL2'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = ACTL6A, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTL6A'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = YEATS4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq YEATS4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("melanoma_beeSwarm.pdf",width=33)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL5'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL1,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL1'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL2,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL2'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = ACTL6A,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTL6A'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = YEATS4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq YEATS4'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_RPL5 = shannonIndex(RPL5)
shannon_RUVBL1 = shannonIndex(RUVBL1)
shannon_RUVBL2 = shannonIndex(RUVBL2)
shannon_ACTL6A = shannonIndex(ACTL6A)
shannon_YEATS4 = shannonIndex(YEATS4)

sig = rbind(shannon_H2AFV,shannon_RPL5,shannon_RUVBL1,shannon_RUVBL2,shannon_ACTL6A,shannon_YEATS4)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_MELANOMA.pdf")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig, dendrogram = "none",
          cellnote=sig,Rowv = FALSE,
          scale="column",  trace="none", 
          notecex=0.7, distfun = function(x) get_dist(x,method="pearson"),
          notecol="black",col=colors,
          na.color=par("bg"),cexCol=.8, cexRow=.9,key=FALSE,
          sepwidth=c(0.01,0.005),
           sepcolor="white",
           colsep=1:ncol(sig),
           rowsep=1:nrow(sig) )

dev.off()
#################################################################################################################################
# LUNG
data = read.table("GSE69405_PROCESSED_GENE_TPM_ALL.txt",sep="\t",header=T)
genenames = as.character(data[,2])
data = data[, grep("_SC",colnames(data))]
######
# RPL5, RUVBL1, RUVBL2, ACTL6A, YEATS4 .

H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=(as.numeric(H2AFV)) )

RPL5 = data[genenames=="RPL5",]
RPL5 = data.frame(cell = gsub("\\_.+","",colnames(RPL5),perl=TRUE), gene=(as.numeric(RPL5)) )

RUVBL1 = data[genenames=="RUVBL1",]
RUVBL1 = data.frame(cell = gsub("\\_.+","",colnames(RUVBL1),perl=TRUE), gene=(as.numeric(RUVBL1)) )

RUVBL2 = data[genenames=="RUVBL2",]
RUVBL2 = data.frame(cell = gsub("\\_.+","",colnames(RUVBL2),perl=TRUE), gene=(as.numeric(RUVBL2)) )

ACTL6A = data[genenames=="ACTL6A",]
ACTL6A = data.frame(cell = gsub("\\_.+","",colnames(ACTL6A),perl=TRUE), gene=(as.numeric(ACTL6A)) )

YEATS4 = data[genenames=="YEATS4",]
YEATS4 = data.frame(cell = gsub("\\_.+","",colnames(YEATS4),perl=TRUE), gene=(as.numeric(YEATS4)) )
######
#
pdf("LUNG_jitter.pdf")
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL5, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL5'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL1, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL1'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RUVBL2, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RUVBL2'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = ACTL6A, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTL6A'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = YEATS4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq YEATS4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()
#
pdf("LUNG_beeSwarm.pdf")
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL5'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL1,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL1'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RUVBL2,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RUVBL2'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = ACTL6A,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTL6A'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = YEATS4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq YEATS4'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()


shannon_H2AFV = shannonIndex(H2AFV)
shannon_RPL5 = shannonIndex(RPL5)
shannon_RUVBL1 = shannonIndex(RUVBL1)
shannon_RUVBL2 = shannonIndex(RUVBL2)
shannon_ACTL6A = shannonIndex(ACTL6A)
shannon_YEATS4 = shannonIndex(YEATS4)

sig = rbind(shannon_H2AFV,shannon_RPL5,shannon_RUVBL1,shannon_RUVBL2,shannon_ACTL6A,shannon_YEATS4)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_LUNG.pdf")
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

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
# UBC, RPL4, HIST1H2AC, HIST2H2BE, H2AFV
####################################################################################################

# MGH's  #OLIG2 - cancer stem cell marker
data = read.table("GSE57872_GBM_data_matrix.txt",sep="\t",header=T,row.names=1)
data = data[,grep("_",colnames(data))]
data = data[,grep("MGH",colnames(data))]
#
H2AFV = data[rownames(data)=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=2^(as.numeric(H2AFV)) )

UBC = data[rownames(data)=="UBC",]
UBC = data.frame(cell = gsub("\\_.+","",colnames(UBC),perl=TRUE), gene=2^(as.numeric(UBC)) )

RPL4 = data[rownames(data)=="RPL4",]
RPL4 = data.frame(cell = gsub("\\_.+","",colnames(RPL4),perl=TRUE), gene=2^(as.numeric(RPL4)) )

HIST1H2AC = data[rownames(data)=="HIST1H2AC",]
HIST1H2AC = data.frame(cell = gsub("\\_.+","",colnames(HIST1H2AC),perl=TRUE), gene=2^(as.numeric(HIST1H2AC)) )

HIST2H2BE = data[rownames(data)=="HIST2H2BE",]
HIST2H2BE = data.frame(cell = gsub("\\_.+","",colnames(HIST2H2BE),perl=TRUE), gene=2^(as.numeric(HIST2H2BE)) )
## UBC, RPL4, HIST1H2AC, HIST2H2BE, H2AFV

pdf("MGH_jitter.pdf")
par(mfrow = c(4,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = UBC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq UBC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST1H2AC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST1H2AC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
    
stripchart(gene ~ cell, vertical = TRUE, data = HIST2H2BE, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST2H2BE'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)


dev.off()
## UBC, RPL4, HIST1H2AC, HIST2H2BE, H2AFV
pdf("MGH_beeSwarm.pdf")
par(mfrow = c(4,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = UBC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq UBC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL4'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST1H2AC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST1H2AC'),col = alpha(colour='red',alpha=.8),cex = .8)
         
beeswarm(gene ~ cell, vertical = TRUE, data = HIST2H2BE,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST2H2BE'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()
## UBC, RPL4, HIST1H2AC, HIST2H2BE, H2AFV

shannon_H2AFV = shannonIndex(H2AFV)
shannon_UBC = shannonIndex(UBC)
shannon_RPL4 = shannonIndex(RPL4)
shannon_HIST1H2AC = shannonIndex(HIST1H2AC)
shannon_HIST2H2BE = shannonIndex(HIST2H2BE)

sig = rbind(shannon_H2AFV,shannon_UBC,shannon_RPL4,shannon_HIST1H2AC,shannon_HIST2H2BE)
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

# UBC, RPL4, HIST1H2AC, HIST2H2BE

H2AFV = tpm[rownames(tpm)=="H2AFV",]
H2AFV = data.frame(cell = colnames(tpm), gene=(as.numeric(H2AFV)) )

UBC = tpm[rownames(tpm)=="UBC",]
UBC = data.frame(cell = colnames(tpm), gene=(as.numeric(UBC)) )

RPL4 = tpm[rownames(tpm)=="RPL4",]
RPL4 = data.frame(cell = colnames(tpm), gene=(as.numeric(RPL4)) )

HIST1H2AC = tpm[rownames(tpm)=="HIST1H2AC",]
HIST1H2AC = data.frame(cell = colnames(tpm), gene=(as.numeric(HIST1H2AC)) )

HIST2H2BE = tpm[rownames(tpm)=="HIST2H2BE",]
HIST2H2BE = data.frame(cell = colnames(tpm), gene=(as.numeric(HIST2H2BE)) )

# UBC, RPL4, HIST1H2AC, HIST2H2BE

pdf("breast_jitter.pdf",width=13)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = UBC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq UBC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST1H2AC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST1H2AC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST2H2BE, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST2H2BE'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()

# UBC, RPL4, HIST1H2AC, HIST2H2BE

pdf("breast_beeSwarm.pdf",width=13)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = UBC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq UBC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL4'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST1H2AC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST1H2AC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST2H2BE,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST2H2BE'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_UBC = shannonIndex(UBC)
shannon_RPL4 = shannonIndex(RPL4)
shannon_HIST1H2AC = shannonIndex(HIST1H2AC)
shannon_HIST2H2BE = shannonIndex(HIST2H2BE)

sig = rbind(shannon_H2AFV,shannon_UBC,shannon_RPL4,shannon_HIST1H2AC,shannon_HIST2H2BE)
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
# UBC, RPL4, HIST1H2AC, HIST2H2BE
H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = cells, gene=(as.numeric(H2AFV)) )

UBC = data[genenames=="UBC",]
UBC = data.frame(cell = cells, gene=(as.numeric(UBC)) )

RPL4 = data[genenames=="RPL4",]
RPL4 = data.frame(cell = cells, gene=(as.numeric(RPL4)) )

HIST1H2AC = data[genenames=="HIST1H2AC",]
HIST1H2AC = data.frame(cell = cells, gene=(as.numeric(HIST1H2AC)) )

HIST2H2BE = data[genenames=="HIST2H2BE",]
HIST2H2BE = data.frame(cell = cells, gene=(as.numeric(HIST2H2BE)) )

####
# UBC, RPL4, HIST1H2AC, HIST2H2BE
pdf("melanoma_jitter.pdf",width=23)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = UBC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq UBC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST1H2AC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST1H2AC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST2H2BE, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST2H2BE'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()

# UBC, RPL4, HIST1H2AC, HIST2H2BE
pdf("melanoma_beeSwarm.pdf",width=33)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = UBC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq UBC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL4'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST1H2AC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST1H2AC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST2H2BE,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST2H2BE'),col = alpha(colour='red',alpha=.8),cex = .8)

dev.off()


shannon_H2AFV = shannonIndex(H2AFV)
shannon_UBC = shannonIndex(UBC)
shannon_RPL4 = shannonIndex(RPL4)
shannon_HIST1H2AC = shannonIndex(HIST1H2AC)
shannon_HIST2H2BE = shannonIndex(HIST2H2BE)

sig = rbind(shannon_H2AFV,shannon_UBC,shannon_RPL4,shannon_HIST1H2AC,shannon_HIST2H2BE)
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
# UBC, RPL4, HIST1H2AC, HIST2H2BE

H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=(as.numeric(H2AFV)) )

RPL5 = data[genenames=="UBC",]
RPL5 = data.frame(cell = gsub("\\_.+","",colnames(UBC),perl=TRUE), gene=(as.numeric(UBC)) )

RUVBL1 = data[genenames=="RPL4",]
RUVBL1 = data.frame(cell = gsub("\\_.+","",colnames(RPL4),perl=TRUE), gene=(as.numeric(RPL4)) )

RUVBL2 = data[genenames=="HIST1H2AC",]
RUVBL2 = data.frame(cell = gsub("\\_.+","",colnames(HIST1H2AC),perl=TRUE), gene=(as.numeric(HIST1H2AC)) )

ACTL6A = data[genenames=="HIST2H2BE",]
ACTL6A = data.frame(cell = gsub("\\_.+","",colnames(HIST2H2BE),perl=TRUE), gene=(as.numeric(HIST2H2BE)) )

######
## UBC, RPL4, HIST1H2AC, HIST2H2BE
pdf("LUNG_jitter.pdf")
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = UBC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq UBC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = RPL4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq RPL4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST1H2AC, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST1H2AC'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HIST2H2BE, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HIST2H2BE'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()
### UBC, RPL4, HIST1H2AC, HIST2H2BE
pdf("LUNG_beeSwarm.pdf")
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = UBC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq UBC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = RPL4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq RPL4'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST1H2AC,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST1H2AC'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HIST2H2BE,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HIST2H2BE'),col = alpha(colour='red',alpha=.8),cex = .8)

dev.off()



shannon_H2AFV = shannonIndex(H2AFV)
shannon_UBC = shannonIndex(UBC)
shannon_RPL4 = shannonIndex(RPL4)
shannon_HIST1H2AC = shannonIndex(HIST1H2AC)
shannon_HIST2H2BE = shannonIndex(HIST2H2BE)

sig = rbind(shannon_H2AFV,shannon_UBC,shannon_RPL4,shannon_HIST1H2AC,shannon_HIST2H2BE)
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

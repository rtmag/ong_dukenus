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
#################################################################################################################################
# GBM 1 - 5 #OLIG2 - cancer stem cell marker
data = read.xlsx("Count_RNA_scope_tissue.xlsx", sheetIndex = 1,stringsAsFactors=F)
colnames(data) = data[1,]
data = data[2:dim(data)[1],2:dim(data)[2]]
h2afv = data.frame(patient = c( rep(colnames(data)[1],dim(data)[1]),
                                rep(colnames(data)[2],dim(data)[1]),
                                rep(colnames(data)[3],dim(data)[1]),
                                rep(colnames(data)[4],dim(data)[1]),
                                rep(colnames(data)[5],dim(data)[1]),
                                rep(colnames(data)[6],dim(data)[1]) ),
           signal = c( as.numeric(data[,1]), as.numeric(data[,2]), as.numeric(data[,3]),
                       as.numeric(data[,4]), as.numeric(data[,5]), as.numeric(data[,6]) ) 
                  )

pdf("RNA_scope_jitter.pdf")
stripchart(signal ~ patient, vertical = TRUE, data = h2afv, jitter = 0.3, ylab = "H2AFV RNA Scope Singal",
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("RNA_scope_beeSwarm.pdf")
beeswarm(signal ~ patient, vertical = TRUE, data = h2afv,method = "hex",pch = 16,xlab="",
         ylab = "H2AFV RNA Scope Singal",col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()


#################################################################################################################################
# MGH's  #OLIG2 - cancer stem cell marker
data = read.table("GSE57872_GBM_data_matrix.txt",sep="\t",header=T,row.names=1)
data = data[,grep("_",colnames(data))]
#
H2AFV = data[rownames(data)=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=(as.numeric(H2AFV)) )

H2AFZ = data[rownames(data)=="H2AFZ",]
H2AFZ = data.frame(cell = gsub("\\_.+","",colnames(H2AFZ),perl=TRUE), gene=(as.numeric(H2AFZ)) )

OLIG1 = data[rownames(data)=="OLIG1",]
OLIG1 = data.frame(cell = gsub("\\_.+","",colnames(OLIG1),perl=TRUE), gene=(as.numeric(OLIG1)) )
#
pdf("MGH_jitter.pdf")
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = H2AFZ, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = OLIG1, jitter = 0.3, ylab = expression('Single Cell RNA-Seq OLIG1'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()
#
pdf("MGH_beeSwarm.pdf")
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFZ,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFZ'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = OLIG1,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq OLIG1'),col = alpha(colour='red',alpha=.8),cex = .8)
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

summary(c(H2AFV[,2],H2AFZ[,2],OLIG1[,2]))

H2AFV[,2] = H2AFV[,2]+7
H2AFZ[,2] = H2AFZ[,2]+7
OLIG1[,2] = OLIG1[,2]+7

shannon_H2AFV = shannonIndex(H2AFV)
shannon_H2AFZ = shannonIndex(H2AFZ)
shannon_OLIG1 = shannonIndex(OLIG1)

sig = rbind(shannon_H2AFV,shannon_H2AFZ,shannon_OLIG1)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_MGH.pdf")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig, dendrogram = "none",
          cellnote=sig,
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

H2AFV = tpm[rownames(tpm)=="H2AFV",]
H2AFV = data.frame(cell = colnames(tpm), gene=(as.numeric(H2AFV)) )


H2AFZ = tpm[rownames(tpm)=="H2AFZ",]
H2AFZ = data.frame(cell = colnames(tpm), gene=(as.numeric(H2AFZ)) )

ITGA6 = tpm[rownames(tpm)=="ITGA6",]
ITGA6 = data.frame(cell = colnames(tpm), gene=(as.numeric(ITGA6)) )

pdf("breast_jitter.pdf",width=13)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell H2AFV TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = H2AFZ, jitter = 0.3, ylab = expression('Single Cell H2AFV TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = ITGA6, jitter = 0.3, ylab = expression('Single Cell ITGA6 TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("breast_beeSwarm.pdf",width=13)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFZ,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFZ TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = ITGA6,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ITGA6'),col = alpha(colour='red',alpha=.8),cex = .8)
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

shannon_H2AFV = shannonIndex(H2AFV)
shannon_H2AFZ = shannonIndex(H2AFZ)
shannon_ITGA6 = shannonIndex(ITGA6)

sig = rbind(shannon_H2AFV,shannon_H2AFZ,shannon_ITGA6)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_BREAST.pdf")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig, dendrogram = "none",
          cellnote=sig,
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



# LUNG

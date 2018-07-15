

#################################################################################################################################
# MGH's  #OLIG2 - cancer stem cell marker
data = read.table("GSE57872_GBM_data_matrix.txt",sep="\t",header=T,row.names=1)
data = data[,grep("_",colnames(data))]
data = data[,grep("MGH",colnames(data))]
#
H2AFV = data[rownames(data)=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=2^(as.numeric(H2AFV)) )

TBP = data[rownames(data)=="TBP",]
TBP = data.frame(cell = gsub("\\_.+","",colnames(TBP),perl=TRUE), gene=2^(as.numeric(TBP)) )

ACTB = data[rownames(data)=="ACTB",]
ACTB = data.frame(cell = gsub("\\_.+","",colnames(ACTB),perl=TRUE), gene=2^(as.numeric(ACTB)) )

HSPA4 = data[rownames(data)=="HSPA4",]
HSPA4 = data.frame(cell = gsub("\\_.+","",colnames(HSPA4),perl=TRUE), gene=2^(as.numeric(HSPA4)) )

HSPA5 = data[rownames(data)=="HSPA5",]
HSPA5 = data.frame(cell = gsub("\\_.+","",colnames(HSPA5),perl=TRUE), gene=2^(as.numeric(HSPA5)) )

GAPDH = data[rownames(data)=="GAPDH",]
GAPDH = data.frame(cell = gsub("\\_.+","",colnames(GAPDH),perl=TRUE), gene=2^(as.numeric(GAPDH)) )
#
pdf("MGH_jitter_houseGenes.pdf")
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = TBP, jitter = 0.3, ylab = expression('Single Cell RNA-Seq TBP'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = ACTB, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTB'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HSPA4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HSPA4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = HSPA5, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HSPA5'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

#stripchart(gene ~ cell, vertical = TRUE, data = GAPDH, jitter = 0.3, ylab = expression('Single Cell RNA-Seq GAPDH'),
#    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()
#
pdf("MGH_beeSwarm_houseGenes.pdf")
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = TBP,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq TBP'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = ACTB,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTB'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HSPA4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA4'),col = alpha(colour='red',alpha=.8),cex = .8)

beeswarm(gene ~ cell, vertical = TRUE, data = HSPA5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA5'),col = alpha(colour='red',alpha=.8),cex = .8)

#beeswarm(gene ~ cell, vertical = TRUE, data = GAPDH,method = "swarm",pch = 16,xlab="",
#         ylab = expression('Single Cell RNA-Seq GAPDH'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_TBP = shannonIndex(TBP)
shannon_ACTB = shannonIndex(ACTB)
shannon_HSPA4 = shannonIndex(HSPA4)
shannon_HSPA5 = shannonIndex(HSPA5)
#shannon_GAPDH = shannonIndex(GAPDH)

sig = rbind(shannon_H2AFV,shannon_TBP,shannon_ACTB,shannon_HSPA4,shannon_HSPA5)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_MGH_houseGenes.pdf")
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

############
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

TBP = tpm[rownames(tpm)=="TBP",]
TBP = data.frame(cell = colnames(tpm), gene=(as.numeric(TBP)) )

ACTB = tpm[rownames(tpm)=="ACTB",]
ACTB = data.frame(cell = colnames(tpm), gene=(as.numeric(ACTB)) )

HSPA4 = tpm[rownames(tpm)=="HSPA4",]
HSPA4 = data.frame(cell = colnames(tpm), gene=(as.numeric(HSPA4)) )

HSPA5 = tpm[rownames(tpm)=="HSPA5",]
HSPA5 = data.frame(cell = colnames(tpm), gene=(as.numeric(HSPA5)) )

GAPDH = tpm[rownames(tpm)=="GAPDH",]
GAPDH = data.frame(cell = colnames(tpm), gene=(as.numeric(GAPDH)) )
#

pdf("breast_jitter_houseGenes.pdf",width=13)
par(mfrow = c(3,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell H2AFV TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = TBP, jitter = 0.3, ylab = expression('Single Cell TBP TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = ACTB, jitter = 0.3, ylab = expression('Single Cell ACTB TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA4, jitter = 0.3, ylab = expression('Single Cell HSPA4 TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA5, jitter = 0.3, ylab = expression('Single Cell HSPA5 TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = GAPDH, jitter = 0.3, ylab = expression('Single Cell GAPDH TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("breast_beeSwarm_houseGenes.pdf",width=13)
par(mfrow = c(3,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = TBP,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq TBP TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = ACTB,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTB'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA4 TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA5 TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = GAPDH,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq GAPDH'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()

shannon_H2AFV = shannonIndex(H2AFV)
shannon_TBP = shannonIndex(TBP)
shannon_ACTB = shannonIndex(ACTB)
shannon_HSPA4 = shannonIndex(HSPA4)
shannon_HSPA5 = shannonIndex(HSPA5)
shannon_GAPDH = shannonIndex(GAPDH)

sig = rbind(shannon_H2AFV,shannon_TBP,shannon_ACTB,shannon_HSPA4,shannon_HSPA5,shannon_GAPDH)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_BREAST_houseGenes.pdf")
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

H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = cells, gene=(as.numeric(H2AFV)) )

TBP = data[genenames=="TBP",]
TBP = data.frame(cell = cells, gene=(as.numeric(TBP)) )

ACTB = data[genenames=="ACTB",]
ACTB = data.frame(cell = cells, gene=(as.numeric(ACTB)) )

HSPA4 = data[genenames=="HSPA4",]
HSPA4 = data.frame(cell = cells, gene=(as.numeric(HSPA4)) )

HSPA5 = data[genenames=="HSPA5",]
HSPA5 = data.frame(cell = cells, gene=(as.numeric(HSPA5)) )

GAPDH = data[genenames=="GAPDH",]
GAPDH = data.frame(cell = cells, gene=(as.numeric(GAPDH)) )
####
pdf("melanoma_jitter_houseGenes.pdf",width=23)
par(mfrow = c(2,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell H2AFV TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = TBP, jitter = 0.3, ylab = expression('Single Cell TBP TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = ACTB, jitter = 0.3, ylab = expression('Single Cell ACTB TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA4, jitter = 0.3, ylab = expression('Single Cell HSPA4 TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA5, jitter = 0.3, ylab = expression('Single Cell HSPA5 TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = GAPDH, jitter = 0.3, ylab = expression('Single Cell GAPDH TPM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()

pdf("melanoma_beeSwarm_houseGenes.pdf",width=33)
par(mfrow = c(2,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = TBP,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq TBP TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = ACTB,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTB TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA4 TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA5 TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = GAPDH,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq GAPDH TPM'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()


shannon_H2AFV = shannonIndex(H2AFV)
shannon_TBP = shannonIndex(TBP)
shannon_ACTB = shannonIndex(ACTB)
shannon_HSPA4 = shannonIndex(HSPA4)
shannon_HSPA5 = shannonIndex(HSPA5)
shannon_GAPDH = shannonIndex(GAPDH)

sig = rbind(shannon_H2AFV,shannon_TBP,shannon_ACTB,shannon_HSPA4,shannon_HSPA5,shannon_GAPDH)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_MELANOMA_houseGenes.pdf")
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
H2AFV = data[genenames=="H2AFV",]
H2AFV = data.frame(cell = gsub("\\_.+","",colnames(H2AFV),perl=TRUE), gene=(as.numeric(H2AFV)) )

TBP = data[genenames=="TBP",]
TBP = data.frame(cell = gsub("\\_.+","",colnames(TBP),perl=TRUE), gene=(as.numeric(TBP)) )

ACTB = data[genenames=="ACTB",]
ACTB = data.frame(cell = gsub("\\_.+","",colnames(ACTB),perl=TRUE), gene=(as.numeric(ACTB)) )

HSPA4 = data[genenames=="HSPA4",]
HSPA4 = data.frame(cell = gsub("\\_.+","",colnames(HSPA4),perl=TRUE), gene=(as.numeric(HSPA4)) )

HSPA5 = data[genenames=="HSPA5",]
HSPA5 = data.frame(cell = gsub("\\_.+","",colnames(HSPA5),perl=TRUE), gene=(as.numeric(HSPA5)) )

GAPDH = data[genenames=="GAPDH",]
GAPDH = data.frame(cell = gsub("\\_.+","",colnames(GAPDH),perl=TRUE), gene=(as.numeric(GAPDH)) )
#
######
#
pdf("LUNG_jitter_houseGenes.pdf")
par(mfrow = c(2,1) )
stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = expression('Single Cell RNA-Seq H2AFV'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = TBP, jitter = 0.3, ylab = expression('Single Cell RNA-Seq TBP'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = ACTB, jitter = 0.3, ylab = expression('Single Cell RNA-Seq ACTB'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA4, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HSPA4'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = HSPA5, jitter = 0.3, ylab = expression('Single Cell RNA-Seq HSPA5'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
stripchart(gene ~ cell, vertical = TRUE, data = GAPDH, jitter = 0.3, ylab = expression('Single Cell RNA-Seq GAPDH'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)
dev.off()
#
pdf("LUNG_beeSwarm_houseGenes.pdf")
par(mfrow = c(2,1) )
beeswarm(gene ~ cell, vertical = TRUE, data = H2AFV,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq H2AFV'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = TBP,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq TBP'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = ACTB,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq ACTB'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA4,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA4'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = HSPA5,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq HSPA5'),col = alpha(colour='red',alpha=.8),cex = .8)
beeswarm(gene ~ cell, vertical = TRUE, data = GAPDH,method = "swarm",pch = 16,xlab="",
         ylab = expression('Single Cell RNA-Seq GAPDH'),col = alpha(colour='red',alpha=.8),cex = .8)
dev.off()


shannon_H2AFV = shannonIndex(H2AFV)
shannon_TBP = shannonIndex(TBP)
shannon_ACTB = shannonIndex(ACTB)
shannon_HSPA4 = shannonIndex(HSPA4)
shannon_HSPA5 = shannonIndex(HSPA5)
shannon_GAPDH = shannonIndex(GAPDH)

sig = rbind(shannon_H2AFV,shannon_TBP,shannon_ACTB,shannon_HSPA4,shannon_HSPA5,shannon_GAPDH)
rownames(sig) = gsub("shannon\\_","",rownames(sig))
sig = round(sig,digits=2)

pdf("shannon_diversity_index_LUNG_houseGenes.pdf")
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



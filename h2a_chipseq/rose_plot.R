#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================
options(bitmapType="cairo")

# [4] "../samtools1.8/H3K27ac_gbm_peaks_6Columns_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt" # this is enhancerFile
# [5] "H3K27ac_gbm_peaks_6Columns" # enhancerName
# [6] "input_rmdup.bam" # wceName

wceName = "present"

#Read enhancer regions with closestGene columns
stitched_regions <- read.delim(file= "H3K27ac_gbm_peaks_6Columns_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt",sep="\t")
annot<- read.table(pipe("cut -f1,2,3,4,16 H3K27ac_gbm_SuperEnhancers.anno"),sep="\t",header = T,stringsAsFactors=F)

stitched_regions_x<- cbind( stitched_regions, annot[match(stitched_regions[,1], annot[,1]),5])

#perform WCE subtraction. Using pipeline table to match samples to proper background.
rankBy_factor = colnames(stitched_regions)[7]
prefix = unlist(strsplit(rankBy_factor,'_'))[1]

if(wceName == 'NONE'){
	rankBy_vector = as.numeric(stitched_regions[,7])
}else{
	wceName = colnames(stitched_regions)[8]
	print('HERE IS THE WCE NAME')
	print(wceName)
	rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])
}	

#SETTING NEGATIVE VALUES IN THE rankBy_vector to 0
rankBy_vector[rankBy_vector < 0] <- 0

#FIGURING OUT THE CUTOFF
cutoff_options <- list()
cutoff_options$absolute <- 8223.0995

#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)

#MAKING HOCKEY STICK PLOT
pdf("superEnhancer_H3K27ac_names_Plot_points.pdf")
signalOrder = order(rankBy_vector,decreasing=TRUE)

plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab="Super enhancers ranked by H3K27ac signal",ylab="H3K27ac ChIP-Seq signal - Input ",pch=19,cex=2)
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')

se_labels<-as.character(stitched_regions_x[signalOrder,9][1:20])
library(wordcloud)
nc=wordlayout( (length(rankBy_vector):1)[1:20], rankBy_vector[signalOrder][1:20],words=se_labels)
nc[,1] <- nc[,1]-500
nc[5:dim(nc)[1],2] <- nc[5:dim(nc)[1],2]+50000
text(nc[,1],nc[,2],label=se_labels,cex=1)

segments(nc[,1]+(nc[,3]/2),nc[,2],(length(rankBy_vector):1)[1:20], rankBy_vector[signalOrder][1:20])

dev.off()



####
# BROAD
options(bitmapType="cairo")
wceName = "present"

stitched_regions <- read.delim(file= "H3K27ac_gbm_broad_6Columns_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt",sep="\t")
annot<- read.table(pipe("cut -f1,2,3,4,16 H3K27ac_gbm_SuperEnhancers_broad.anno"),sep="\t",header = T,stringsAsFactors=F)

stitched_regions_x<- cbind( stitched_regions, annot[match(stitched_regions[,1], annot[,1]),5])

#perform WCE subtraction. Using pipeline table to match samples to proper background.
rankBy_factor = colnames(stitched_regions)[7]
prefix = unlist(strsplit(rankBy_factor,'_'))[1]

wceName = colnames(stitched_regions)[8]
rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])

#SETTING NEGATIVE VALUES IN THE rankBy_vector to 0
rankBy_vector[rankBy_vector < 0] <- 0

#FIGURING OUT THE CUTOFF
cutoff_options <- list()
cutoff_options$absolute <- 9202.9535

#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)

#MAKING HOCKEY STICK PLOT
pdf("superEnhancer_H3K27ac_names_Plot_points_broad.pdf")
signalOrder = order(rankBy_vector,decreasing=TRUE)

plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab="Super enhancers ranked by H3K27ac signal",ylab="H3K27ac ChIP-Seq signal - Input ",pch=19,cex=2)
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')

se_labels<-as.character(stitched_regions_x[signalOrder,9][1:20])
library(wordcloud)
nc=wordlayout( (length(rankBy_vector):1)[1:20], rankBy_vector[signalOrder][1:20],words=se_labels)
nc[,1] <- nc[,1]-500
nc[5:dim(nc)[1],2] <- nc[5:dim(nc)[1],2]+50000
text(nc[,1],nc[,2],label=se_labels,cex=1)

segments(nc[,1]+(nc[,3]/2),nc[,2],(length(rankBy_vector):1)[1:20], rankBy_vector[signalOrder][1:20])

dev.off()



##############################
dReg_genes<-read.table("downr_genes_superEnhancers_expressionSulconazol.txt")
dReg_genes<-as.character(dReg_genes[,1])
dReg_genes<-dReg_genes[dReg_genes!="CCDC26"]

#Read enhancer regions with closestGene columns
stitched_regions <- read.delim(file= "H3K27ac_gbm_peaks_6Columns_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt",sep="\t")
annot<- read.table(pipe("cut -f1,2,3,4,16 H3K27ac_gbm_SuperEnhancers.anno"),sep="\t",header = T,stringsAsFactors=F)

stitched_regions_x<- cbind( stitched_regions, annot[match(stitched_regions[,1], annot[,1]),5])

#perform WCE subtraction. Using pipeline table to match samples to proper background.
rankBy_factor = colnames(stitched_regions)[7]
prefix = unlist(strsplit(rankBy_factor,'_'))[1]

wceName = colnames(stitched_regions)[8]
rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])

#SETTING NEGATIVE VALUES IN THE rankBy_vector to 0
rankBy_vector[rankBy_vector < 0] <- 0

#FIGURING OUT THE CUTOFF
cutoff_options <- list()
cutoff_options$absolute <- 8223.0995

#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)

signalOrder = order(rankBy_vector,decreasing=TRUE)

#MAKING HOCKEY STICK PLOT
pdf("superEnhancer_H3K27ac_names_Plot_points_withDownRegGenes.pdf")

pdf("ROSE_allDownreg_genes.pdf")
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab="Super enhancers ranked by H3K27ac signal",ylab="H3K27ac ChIP-Seq signal - Input ",pch=19,cex=2)
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')

se_labels<-as.character(stitched_regions_x[signalOrder,9][stitched_regions_x[signalOrder,9] %in% dReg_genes][1:10])
library(wordcloud)
nc=wordlayout( (length(rankBy_vector):1)[stitched_regions_x[signalOrder,9] %in% dReg_genes], rankBy_vector[signalOrder][stitched_regions_x[signalOrder,9] %in% dReg_genes],words=se_labels)
nc[,1] <- nc[,1]-500
nc[5:dim(nc)[1],2] <- nc[5:dim(nc)[1],2]+50000
text(nc[,1],nc[,2],label=se_labels,cex=1)

segments(nc[,1]+(nc[,3]/2),nc[,2],(length(rankBy_vector):1)[stitched_regions_x[signalOrder,9] %in% dReg_genes], rankBy_vector[signalOrder][stitched_regions_x[signalOrder,9] %in% dReg_genes])

dev.off()


# top 10 genes downreg
#MAKING HOCKEY STICK PLOT
pdf("superEnhancer_H3K27ac_names_Plot_points_withDownRegGenes.pdf")
signalOrder = order(rankBy_vector,decreasing=TRUE)

pdf("ROSE_allDownreg_genes_TOP10.pdf")
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab="Super enhancers ranked by H3K27ac signal",ylab="H3K27ac ChIP-Seq signal - Input ",pch=19,cex=3,
    cex.lab=1.5)
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')

se_labels<-as.character(stitched_regions_x[signalOrder,9][stitched_regions_x[signalOrder,9] %in% dReg_genes][c(1,2,3,5,6,8,9,10,11,12,14)])
library(wordcloud)
nc=wordlayout( (length(rankBy_vector):1)[stitched_regions_x[signalOrder,9] %in% dReg_genes][c(1,2,3,5,6,8,9,10,11,12,14)], rankBy_vector[signalOrder][stitched_regions_x[signalOrder,9] %in% dReg_genes][c(1,2,3,5,6,8,9,10,11,12,14)],words=se_labels,cex=2)
nc[,1] <- nc[,1]-500
nc[2:dim(nc)[1],2] <- nc[2:dim(nc)[1],2]+80000
text(nc[,1],nc[,2],label=se_labels,cex=2)

segments(nc[,1]+(nc[,3]/2),nc[,2],(length(rankBy_vector):1)[stitched_regions_x[signalOrder,9] %in% dReg_genes][c(1,2,3,5,6,8,9,10,11,12,14)], rankBy_vector[signalOrder][stitched_regions_x[signalOrder,9] %in% dReg_genes][c(1,2,3,5,6,8,9,10,11,12,14)])

dev.off()

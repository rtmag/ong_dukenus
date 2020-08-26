#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================
options(bitmapType="cairo")

# [4] "../samtools1.8/H3K27ac_gbm_peaks_6Columns_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt" # this is enhancerFile
# [5] "H3K27ac_gbm_peaks_6Columns" # enhancerName
# [6] "input_rmdup.bam" # wceName

wceName = "present"

#Read enhancer regions with closestGene columns
stitched_regions <- read.delim(file= enhancerFile,sep="\t")

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
cutoff_options <- ""
cutoff_options$absolute <- 8223.0995

#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)

#MAKING HOCKEY STICK PLOT
plotFileName = paste('superEnhancer_H3K27ac_names_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
if(wceName == 'NONE'){
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal'),pch=19,cex=2)	
	
}else{
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal','- ',wceName),pch=19,cex=2)
}
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')
#text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)

dev.off()

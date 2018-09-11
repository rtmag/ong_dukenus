rna = readRDS("MASTER_RNASEQ_TABLE_OLEG.rds")

atac = read.table(pipe('more ATAC_TSS_500bp_20bp.rmat|grep -v "#"|grep -v "shNT_1_rmdup"'),sep="\t",header=F)
bed = read.table("hg19_tss_knownCanonical_noUnasembled.bed",sep="\t",stringsAsFactors=F)
rownames(atac) = make.unique(as.character(bed[,4]), sep = "_")

nt1 = rowSums(atac[,1:50])
nt2 = rowSums(atac[,51:100])
h11 = rowSums(atac[,101:150])
h12 = rowSums(atac[,151:200])
h21 = rowSums(atac[,201:250])
h22 = rowSums(atac[,251:300])

atacfc = log2(rowMeans(cbind(h11,h12,h21,h22))+0.001) - log2(rowMeans(cbind(nt1,nt2))+0.001)


ix = match(names(atacfc),rna[,1])
atacfc = atacfc[!is.na(ix)]
rna=rna[ix[!is.na(ix)],]

jx = rowSums(rna[,6:14]>5)>0
rna = rna[jx,]
atacfc = atacfc[jx]


options(scipen=999)
library(graphics)
library(RColorBrewer)
 
buylrd <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026",
            "#A50026", "#A50026", "#A50026", "#A50026","#A50026", "#A50026", "#A50026", "#A50026","#A50026", "#A50026", "#A50026",
           "#A50026", "#A50026", "#A50026", "#A50026","#A50026", "#A50026", "#A50026", "#A50026","#A50026", "#A50026", "#A50026")

plot(1, type="n",xlab="RNA-Seq Log2 Fold Change",ylab="ATAC-Seq Log2 Fold Change",main = "",xlim=c(-4,4), ylim=c(-3,3))
rect(xleft=-5, ybottom=-5, xright=5, ytop=5,col="#313695",border=NA)
par(new=T)
smoothScatter(-log2(rna$PostFC),atacfc, nbin=1000, colramp = colorRampPalette(c(buylrd)), 
              nrpoints=Inf,pch="", cex=.7, transformation = function(x) x^.6, col="black",axes=F, ann=F,xlim=c(6,13), ylim=c(6,13))
abline(h=0)
abline(v=0)
par(new=NULL)
#
 
smoothScatter(-log2(rna$PostFC),atacfc,nrpoints=0)

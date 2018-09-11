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

atacfc = log2(rowMeans(cbind(h11,h12,h21,h22))) - log2(rowMeans(cbind(nt1,nt2)))


ix = match(names(atacfc),rna[,1])
atacfc = atacfc[!is.na(ix)]
rna=rna[ix[!is.na(ix)],]

jx = rowSums(rna[,6:14]>5)>0
rna = rna[jx,]
atacfc = atacfc[jx]

library(graphics)

smoothScatter(-log2(rna$PostFC),atacfc)

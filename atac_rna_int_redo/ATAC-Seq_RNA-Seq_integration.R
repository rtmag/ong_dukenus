rna = readRDS("/Users/wone/CSI/ong/atacseq_redo/MASTER_RNASEQ_TABLE_OLEG.rds")

atac = read.table(pipe('more /Users/wone/CSI/ong/atacseq_redo/ATAC_TSS_500bp_20bp.rmat|grep -v "#"|grep -v "shNT_1_rmdup"'),sep="\t",header=F)
bed = read.table("/Users/wone/CSI/ong/atacseq_redo/hg19_tss_knownCanonical_noUnasembled.bed",sep="\t",stringsAsFactors=F)
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

pdf("ATAC-RNA_integration_TSS.pdf")
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),
             )
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()
 
targets=read.table("~/CSI/sjlab/deseq2/HallMark_e2f_targets.txt")
targets<-as.character(unique(targets[,1]))
 
pdf("ATAC-RNA_integration_TSS_E2F_targets.pdf")
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),
             )
abline(h=0,lty=3)
abline(v=0,lty=3)
points(-log2(rna$PostFC)[rna[,1] %in% targets],
       atacfc[names(atacfc) %in% targets],
      col="red",pch=20)
dev.off()
 
scp -P 60035 remap2018_REST_nr_macs2_hg38_v1_2.bed root@172.18.149.78:./root/ong_dukenus/TF_targets
annotatePeaks.pl remap2018_REST_nr_macs2_hg38_v1_2.bed hg38 -annStats REST_remap.annStats > REST_remap.anno 

more remap2018_REST_all_macs2_hg38_v1_2.bed |grep -i "Neural" > rest_neural_remap.bed
scp -P 60035 rest_neural_remap.bed root@172.18.149.78:/root/ong_dukenus/TF_targets
annotatePeaks.pl rest_neural_remap.bed hg38 -annStats rest_neural_remap.annStats > rest_neural_remap.anno 

library(TFregulomeR)

rest <- loadPeaks(id = "GTRD-EXP010296_HSA_U87_REST", includeMotifOnly = TRUE)
write.table(rest,"rest_gbm_quy.bed",sep="\t",quote=F,row.names=F,col.names=F)


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
points(-log2(rna$PostFC)[rna[,1] %in% targets],
       atacfc[names(atacfc) %in% targets],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()
 
#scp -P 60035 remap2018_REST_nr_macs2_hg38_v1_2.bed root@172.18.149.78:./root/ong_dukenus/TF_targets
#annotatePeaks.pl remap2018_REST_nr_macs2_hg38_v1_2.bed hg38 -annStats REST_remap.annStats > REST_remap.anno 

#more remap2018_REST_all_macs2_hg38_v1_2.bed |grep -i "Neural" > rest_neural_remap.bed
#scp -P 60035 rest_neural_remap.bed root@172.18.149.78:/root/ong_dukenus/TF_targets
#annotatePeaks.pl rest_neural_remap.bed hg38 -annStats rest_neural_remap.annStats > rest_neural_remap.anno 
################################################################################################################
################################################################################################################
library(TFregulomeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

TFBS_brain <- dataBrowser(organ = "brain",species="HUMAN")
TFBS <- dataBrowser(disease_state="tumor",species="HUMAN",tf="E2F1")
########################################################################
foxm1 <- loadPeaks(id = "MM1_HSA_K562_FOXM1", includeMotifOnly = TRUE)
foxm1_anno <- genomeAnnotate(peaks = foxm1,return_annotation = TRUE)
foxm1_target <- sort(unique(foxm1_anno[foxm1_anno$annotation=="promoter-TSS",'geneName']))
foxm1_target <- foxm1_target[foxm1_target != ""]
########################################################################
cebpb <- loadPeaks(id = "MM1_HSA_K562_CEBPB", includeMotifOnly = TRUE)
cebpb_anno <- genomeAnnotate(peaks = cebpb,return_annotation = TRUE)
cebpb_target <- sort(unique(cebpb_anno[cebpb_anno$annotation=="promoter-TSS",'geneName']))
cebpb_target <- cebpb_target[cebpb_target != ""]
########################################################################
hif1a <- loadPeaks(id = "GTRD-EXP035004_HSA_U2OS_HIF1A", includeMotifOnly = FALSE)
hif1a_anno <- genomeAnnotate(peaks = hif1a,return_annotation = TRUE)
hif1a_target <- sort(unique(hif1a_anno[hif1a_anno$annotation=="promoter-TSS",'geneName']))
hif1a_target <- hif1a_target[hif1a_target != ""]
########################################################################
rest <- loadPeaks(id = "GTRD-EXP010296_HSA_U87_REST", includeMotifOnly = TRUE)
rest_anno <- genomeAnnotate(peaks = rest,return_annotation = TRUE)
rest_target <- sort(unique(rest_anno[rest_anno$annotation=="promoter-TSS",'geneName']))
rest_target <- rest_target[rest_target != ""]
########################################################################
twist1 <- loadPeaks(id = "GTRD-EXP036335_HSA_SK-N-BE2-C_TWIST1", includeMotifOnly = TRUE)
twist1_anno <- genomeAnnotate(peaks = twist1,return_annotation = TRUE)
twist1_target <- sort(unique(twist1_anno[twist1_anno$annotation=="promoter-TSS",'geneName']))
twist1_target <- twist1_target[twist1_target != ""]
########################################################################
e2f4 <- loadPeaks(id = "MM1_HSA_K562_E2F4", includeMotifOnly = TRUE)
e2f4_anno <- genomeAnnotate(peaks = e2f4,return_annotation = TRUE)
e2f4_target <- sort(unique(e2f4_anno[e2f4_anno$annotation=="promoter-TSS",'geneName']))
e2f4_target <- e2f4_target[e2f4_target != ""]
########################################################################
cebpd <- loadPeaks(id = "MM1_HSA_HepG2_CEBPD", includeMotifOnly = TRUE)
cebpd_anno <- genomeAnnotate(peaks = cebpd,return_annotation = TRUE)
cebpd_target <- sort(unique(cebpd_anno[cebpd_anno$annotation=="promoter-TSS",'geneName']))
cebpd_target <- cebpd_target[cebpd_target != ""]
########################################################################
nfkb1 <- loadPeaks(id = "GTRD-EXP034345_HSA_L1236_NFKB1", includeMotifOnly = TRUE)
nfkb1_anno <- genomeAnnotate(peaks = nfkb1,return_annotation = TRUE)
nfkb1_target <- sort(unique(nfkb1_anno[nfkb1_anno$annotation=="promoter-TSS",'geneName']))
nfkb1_target <- nfkb1_target[nfkb1_target != ""]
########################################################################
e2f1 <- loadPeaks(id = "MM1_HSA_K562_E2F1", includeMotifOnly = TRUE)
e2f1_anno <- genomeAnnotate(peaks = e2f1,return_annotation = TRUE)
e2f1_target <- sort(unique(e2f1_anno[e2f1_anno$annotation=="promoter-TSS",'geneName']))
e2f1_target <- e2f1_target[e2f1_target != ""]
########################################################################
pdf("RNA_ATAC_chipseq_targetID.pdf")
par(mfrow=c(3,3))
#1
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab="",
              ylab="",
             main = "FOXM1 targets")
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_target],
       atacfc[names(atacfc) %in% foxm1_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)

#2
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPB targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_target],
       atacfc[names(atacfc) %in% cebpb_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#3
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "HIF1A targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_target],
       atacfc[names(atacfc) %in% hif1a_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#4
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "REST targets",xlab="",
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'))
points(-log2(rna$PostFC)[rna[,1] %in% rest_target],
       atacfc[names(atacfc) %in% rest_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#5
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "TWIST1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% twist1_target],
       atacfc[names(atacfc) %in% twist1_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#6
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F4 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_target],
       atacfc[names(atacfc) %in% e2f4_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#7
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPD targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_target],
       atacfc[names(atacfc) %in% cebpd_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#8
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "NFKB1 targets",
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_target],
       atacfc[names(atacfc) %in% nfkb1_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#9
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_target],
       atacfc[names(atacfc) %in% e2f1_target],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()

##########################
# TRRUST targets
trrust <- read.table("trrust_rawdata.human.tsv",header=T,stringsAsFactors=FALSE)

foxm1_a <- trrust[trrust$Source=="FOXM1" & trrust$Effect=="Activation",2]
foxm1_r <- trrust[trrust$Source=="FOXM1" & trrust$Effect=="Repression",2]
foxm1_u <- trrust[trrust$Source=="FOXM1" & trrust$Effect=="Unknown",2]

cebpb_a <- trrust[trrust$Source=="CEBPB" & trrust$Effect=="Activation",2]
cebpb_r <- trrust[trrust$Source=="CEBPB" & trrust$Effect=="Repression",2]
cebpb_u <- trrust[trrust$Source=="CEBPB" & trrust$Effect=="Unknown",2]

hif1a_a <- trrust[trrust$Source=="HIF1A" & trrust$Effect=="Activation",2]
hif1a_r <- trrust[trrust$Source=="HIF1A" & trrust$Effect=="Repression",2]
hif1a_u <- trrust[trrust$Source=="HIF1A" & trrust$Effect=="Unknown",2]

rest_a <- trrust[trrust$Source=="REST" & trrust$Effect=="Activation",2]
rest_r <- trrust[trrust$Source=="REST" & trrust$Effect=="Repression",2]
rest_u <- trrust[trrust$Source=="REST" & trrust$Effect=="Unknown",2]

twist1_a <- trrust[trrust$Source=="TWIST1" & trrust$Effect=="Activation",2]
twist1_r <- trrust[trrust$Source=="TWIST1" & trrust$Effect=="Repression",2]
twist1_u <- trrust[trrust$Source=="TWIST1" & trrust$Effect=="Unknown",2]

e2f4_a <- trrust[trrust$Source=="E2F4" & trrust$Effect=="Activation",2]
e2f4_r <- trrust[trrust$Source=="E2F4" & trrust$Effect=="Repression",2]
e2f4_u <- trrust[trrust$Source=="E2F4" & trrust$Effect=="Unknown",2]

cebpd_a <- trrust[trrust$Source=="CEBPD" & trrust$Effect=="Activation",2]
cebpd_r <- trrust[trrust$Source=="CEBPD" & trrust$Effect=="Repression",2]
cebpd_u <- trrust[trrust$Source=="CEBPD" & trrust$Effect=="Unknown",2]

nfkb1_a <- trrust[trrust$Source=="NFKB1" & trrust$Effect=="Activation",2]
nfkb1_r <- trrust[trrust$Source=="NFKB1" & trrust$Effect=="Repression",2]
nfkb1_u <- trrust[trrust$Source=="NFKB1" & trrust$Effect=="Unknown",2]

e2f1_a <- trrust[trrust$Source=="E2F1" & trrust$Effect=="Activation",2]
e2f1_r <- trrust[trrust$Source=="E2F1" & trrust$Effect=="Repression",2]
e2f1_u <- trrust[trrust$Source=="E2F1" & trrust$Effect=="Unknown",2]

pdf("RNA_ATAC_trrust_targetID.pdf")
par(mfrow=c(3,3))
#1
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab="",
              ylab="",
             main = "FOXM1 targets")
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_a],
       atacfc[names(atacfc) %in% foxm1_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_r],
       atacfc[names(atacfc) %in% foxm1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_u],
       atacfc[names(atacfc) %in% foxm1_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)

#2
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPB targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_a],
       atacfc[names(atacfc) %in% cebpb_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_r],
       atacfc[names(atacfc) %in% cebpb_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_u],
       atacfc[names(atacfc) %in% cebpb_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#3
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "HIF1A targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_a],
       atacfc[names(atacfc) %in% hif1a_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_r],
       atacfc[names(atacfc) %in% hif1a_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_u],
       atacfc[names(atacfc) %in% hif1a_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#4
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "REST targets",xlab="",
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'))
points(-log2(rna$PostFC)[rna[,1] %in% rest_a],
       atacfc[names(atacfc) %in% rest_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% rest_r],
       atacfc[names(atacfc) %in% rest_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% rest_u],
       atacfc[names(atacfc) %in% rest_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#5
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "TWIST1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% twist1_a],
       atacfc[names(atacfc) %in% twist1_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% twist1_r],
       atacfc[names(atacfc) %in% twist1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% twist1_u],
       atacfc[names(atacfc) %in% twist1_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#6
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F4 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_a],
       atacfc[names(atacfc) %in% e2f4_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_r],
       atacfc[names(atacfc) %in% e2f4_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_u],
       atacfc[names(atacfc) %in% e2f4_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#7
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPD targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_a],
       atacfc[names(atacfc) %in% cebpd_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_r],
       atacfc[names(atacfc) %in% cebpd_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_u],
       atacfc[names(atacfc) %in% cebpd_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#8
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "NFKB1 targets",
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_a],
       atacfc[names(atacfc) %in% nfkb1_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_r],
       atacfc[names(atacfc) %in% nfkb1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_u],
       atacfc[names(atacfc) %in% nfkb1_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#9
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_a],
       atacfc[names(atacfc) %in% e2f1_a],
      col="green",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_r],
       atacfc[names(atacfc) %in% e2f1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_u],
       atacfc[names(atacfc) %in% e2f1_u],
      col="grey",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()


##############################################################################################################################
##############################################################################################################################

pdf("RNA_ATAC_trrust_targetID_allRED.pdf")
par(mfrow=c(3,3))
#1
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab="",
              ylab="",
             main = "FOXM1 targets")
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_a],
       atacfc[names(atacfc) %in% foxm1_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_r],
       atacfc[names(atacfc) %in% foxm1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_u],
       atacfc[names(atacfc) %in% foxm1_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)

#2
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPB targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_a],
       atacfc[names(atacfc) %in% cebpb_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_r],
       atacfc[names(atacfc) %in% cebpb_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_u],
       atacfc[names(atacfc) %in% cebpb_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#3
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "HIF1A targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_a],
       atacfc[names(atacfc) %in% hif1a_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_r],
       atacfc[names(atacfc) %in% hif1a_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_u],
       atacfc[names(atacfc) %in% hif1a_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#4
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "REST targets",xlab="",
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'))
points(-log2(rna$PostFC)[rna[,1] %in% rest_a],
       atacfc[names(atacfc) %in% rest_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% rest_r],
       atacfc[names(atacfc) %in% rest_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% rest_u],
       atacfc[names(atacfc) %in% rest_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#5
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "TWIST1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% twist1_a],
       atacfc[names(atacfc) %in% twist1_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% twist1_r],
       atacfc[names(atacfc) %in% twist1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% twist1_u],
       atacfc[names(atacfc) %in% twist1_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#6
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F4 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_a],
       atacfc[names(atacfc) %in% e2f4_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_r],
       atacfc[names(atacfc) %in% e2f4_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_u],
       atacfc[names(atacfc) %in% e2f4_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#7
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPD targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_a],
       atacfc[names(atacfc) %in% cebpd_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_r],
       atacfc[names(atacfc) %in% cebpd_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_u],
       atacfc[names(atacfc) %in% cebpd_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#8
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "NFKB1 targets",
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_a],
       atacfc[names(atacfc) %in% nfkb1_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_r],
       atacfc[names(atacfc) %in% nfkb1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_u],
       atacfc[names(atacfc) %in% nfkb1_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
#9
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_a],
       atacfc[names(atacfc) %in% e2f1_a],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_r],
       atacfc[names(atacfc) %in% e2f1_r],
      col="red",pch=20)
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_u],
       atacfc[names(atacfc) %in% e2f1_u],
      col="red",pch=20)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()

#######################

pdf("RNA_ATAC_chipseq_targetID_PPEElow5.pdf")
par(mfrow=c(3,3))
#1
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab="",
              ylab="",
             main = "FOXM1 targets")
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% foxm1_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)

#2
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPB targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% cebpb_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#3
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "HIF1A targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% hif1a_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#4
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "REST targets",xlab="",
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'))
points(-log2(rna$PostFC)[rna[,1] %in% rest_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% rest_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#5
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "TWIST1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% twist1_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% twist1_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#6
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F4 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% e2f4_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#7
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPD targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% cebpd_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#8
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "NFKB1 targets",
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% nfkb1_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#9
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_target & rna$PPEE<0.05],
       atacfc[names(atacfc) %in% e2f1_target & rna$PPEE<0.05],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()
#############################

pdf("RNA_ATAC_chipseq_targetID_PPEElow5_l2FC5.pdf")
par(mfrow=c(3,3))
#1
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,
              xlab="",
              ylab="",
             main = "FOXM1 targets")
points(-log2(rna$PostFC)[rna[,1] %in% foxm1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% foxm1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)

#2
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPB targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpb_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% cebpb_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#3
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "HIF1A targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% hif1a_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% hif1a_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#4
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "REST targets",xlab="",
              ylab=expression('ATAC-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'))
points(-log2(rna$PostFC)[rna[,1] %in% rest_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% rest_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#5
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "TWIST1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% twist1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% twist1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#6
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F4 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f4_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% e2f4_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#7
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "CEBPD targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% cebpd_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% cebpd_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#8
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "NFKB1 targets",
              xlab=expression('RNA-Seq Log'[2]*' Fold Change ( shH2AFV / shNT )'),ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% nfkb1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% nfkb1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
#9
smoothScatter(-log2(rna$PostFC),atacfc,xlim=c(-3,3), ylim=c(-2,2),nrpoints=0,main = "E2F1 targets",xlab="",ylab="")
points(-log2(rna$PostFC)[rna[,1] %in% e2f1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
       atacfc[names(atacfc) %in% e2f1_target & rna$PPEE<0.05 & abs(-log2(rna$PostFC))>.5],
      col="red",pch=20,cex=.4)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()

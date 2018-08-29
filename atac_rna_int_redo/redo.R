x=read.csv("Master_tab_o_202_143_400_R.csv")
x=x[,c(2,10,12)]

#
#more atac.diffreps |grep -v "Treatment.avg"|grep -v "#"|bedtools intersect -a - -b hg19_tss_knownCanonical_noUnasembled.bed -wa -wb|cut -f12,23 > ATAC_TSS_lo2FC.txt#
#

atac=read.table("/Users/wone/CSI/ong/atacseq_redo/atac_rna_int/ATAC_TSS_lo2FC.txt",sep="\t")

match(x[,1],atac[,2])

ix = match(atac[,2],x[,1])
atac = atac[!is.na(ix),]
deseq=x[ix[!is.na(ix)],]

table = cbind(deseq,atac[,1])
table[,3] = log2(table[,3])
colnames(table) = c("GeneName","RNA_FDR","RNA_log2FC","ATAC_TSS_log2FC")

saveRDS(table,"RNA_atac_int.rds")
write.csv(table,"RNA_atac_int.csv")


stringent = table[table[,2]<0.05 & table[,3]<(-.6),]
stringent = stringent[complete.cases(stringent),]
stringent= stringent[order(stringent[,4]),]

write.csv(stringent,"stringent_RNA_atac_integration.csv")

relaxed = table[table[,2]<0.05 & table[,3]<(-.3),]
relaxed = relaxed[complete.cases(relaxed),]
relaxed= relaxed[order(relaxed[,4]),]

write.csv(relaxed,"relaxed_RNA_atac_integration.csv")

######
pdf("RNA_atac_int.pdf")
library(graphics)
smoothScatter(table[,3],table[,4],ylab="ATAC-Seq TSS log2FC",xlab="RNA-Seq log2FC",nrpoints=0,main="ATAC-RNA Seq Integration")
abline(h=0);abline(v=0)
dev.off()

pdf("stringent_RNA_atac_integration.pdf")
smoothScatter(stringent[,3],stringent[,4],ylab="ATAC-Seq TSS log2FC",xlab="RNA-Seq log2FC",nrpoints=0,main="ATAC-RNA Seq Stringent")
abline(h=0);abline(v=0)
dev.off()

pdf("relaxed_RNA_atac_integration.pdf")
smoothScatter(relaxed[,3],relaxed[,4],ylab="ATAC-Seq TSS log2FC",xlab="RNA-Seq log2FC",nrpoints=0,main="ATAC-RNA Seq Relaxed")
abline(h=0);abline(v=0)
dev.off()


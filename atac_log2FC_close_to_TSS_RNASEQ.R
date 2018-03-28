
more ~/resources/hg19_tss_knownCanonical_noUnasembled.bed| \
awk -F"\t" '{ if($6=="+"){ print $1"\t"$2-5000"\t"$3+1000"\t"$4"\t"$5"\t"$6 }  
if($6=="-"){ print $1"\t"$2-1000"\t"$3+5000"\t"$4"\t"$5"\t"$6 }
}' | perl -pe 's/\-\d+\t/1\t/g' > hg19_tss_knownCanonical_Upstream5kb_Downstream1kb.bed

grep -P -v "#|Treatment.cnt" atac.diffreps| \
intersectBed -a hg19_tss_knownCanonical_Upstream5kb_Downstream1kb.bed -b - -wa -wb|cut -f4,7,8,9,18,20 \
> atac_diffreps_integration.txt
########################################################################################################
options(scipen=999)
rna = read.csv("deseq2_results.csv")
atac = read.table("atac_diffreps_integration.txt",sep="\t",stringsAsFactors=F)

inte = data.frame(atac,rna[match(atac[,1],rna$gene_symbol),c(3,7)])
colnames(inte) = c("gene_symbol","atac_chr","atac_start","atac_end","atac_log2FC","atac_FDR","rnaseq_log2FC","rnaseq_FDR")

inte$rnaseq_FDR = round(inte$rnaseq_FDR,digits=14)
inte = inte[order(inte$rnaseq_log2FC),]

write.csv(inte,"integration_rnaseq_atac_Upstream5kb_Downstream1kb.csv")

inte = inte[!is.na(inte$rnaseq_log2FC),]
write.csv(inte,"integration_rnaseq_atac_Upstream5kb_Downstream1kb_ZeroCountRemoved.csv")

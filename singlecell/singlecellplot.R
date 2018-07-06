data = read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",sep="\t",header=T)

# BC01-BC02, estrogen receptor positive (ER+); 
# BC03, double positive (ER+ and HER2+); 
# BC04-BC06, human epidermal growth factor receptor 2 positive (HER2+); 
# BC07-BC11, triple-negative breast cancer (TNBC); 
# BC03LN, lymph node metastasis of BC03; 
# BC07LN, lymph node metastasis of BC07)

library(ggplot2)

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



############################################################################################################
H2AFV = tpm[rownames(tpm)=="H2AFV",]
H2AFV = data.frame(cell = colnames(tpm), gene=as.numeric(H2AFV))


H2AFZ = tpm[rownames(tpm)=="H2AFZ",]
H2AFZ = data.frame(cell = colnames(tpm), gene=as.numeric(H2AFZ))

############################################################################################################

pdf("breast_primary_tumor.pdf",width=14)

par(mfrow = c(2,1) )

stripchart(gene ~ cell, vertical = TRUE, data = H2AFV, jitter = 0.3, ylab = "Single Cell H2AFV TPM",
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

stripchart(gene ~ cell, vertical = TRUE, data = H2AFZ, jitter = 0.3, ylab = "Single Cell H2AFZ TPM",
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2)

dev.off()


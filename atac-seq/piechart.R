#annotatePeaks.pl highconfidence.broadpeak hg19 -annStats highconfidence.annStats > highconfidence.bed.anno
#############################################################################################
pdf("shNT_h2afv_highconfidence.pdf")
res=read.table(pipe("more highconfidence.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8)
title("Distribution of H2AFV peaks\n(41,294 regions)", cex.main=.9)
dev.off()

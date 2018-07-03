#annotatePeaks.pl h2_vs_nt_100reads_down.bed hg19 -annStats h2_vs_nt_100reads_down.annStats > h2_vs_nt_100reads_down.bed.anno
#############################################################################################
pdf("h2_vs_nt_100reads_down.pdf")
res=read.table(pipe("more h2_vs_nt_100reads_down.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
title("Distribution of open chromatin in shNT\n(5,584 regions)", cex.main=.9)
dev.off()

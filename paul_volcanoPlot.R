res=read.table("atac.diffreps.annotated", sep="\t",header=T)

#down:23466 up:15360 

pdf("pie_Open_in_NT.pdf")
tdown=table(res$Feature[res$Event=="Down"])
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
pie(tdown, main="Distribution of chromatin regions open in NT compared to SH\n(23,466 regions)")
dev.off()

pdf("pie_Open_in_SH.pdf")
tup=table(res$Feature[res$Event=="Up"])
names(tup) = paste(names(tup)," ",round(tup/sum(tup)*100,digits=2),"%",sep="")
pie(tup, main="Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)")
dev.off()

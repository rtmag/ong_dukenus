h2az = read.table(pipe("grep -v '#' h2_vs_nt_100reads_tss.rmat|grep -v 'genes'"),sep="\t")

#########################################################################################################################

# 95% confidence interval
nt1=colMeans(h2az[,1:200],na.rm=T)
nt2=colMeans(h2az[,201:400],na.rm=T)
sh11=colMeans(h2az[,401:600],na.rm=T)
sh12=colMeans(h2az[,601:800],na.rm=T)
sh21=colMeans(h2az[,801:1000],na.rm=T)
sh22=colMeans(h2az[,1001:1200],na.rm=T)

nt1_conf=apply(h2az[1:200],2,function(x) t.test(x)$conf.int)
nt2_conf=apply(h2az[201:400],2,function(x) t.test(x)$conf.int)
sh11_conf=apply(h2az[401:600],2,function(x) t.test(x)$conf.int)
sh12_conf=apply(h2az[601:800],2,function(x) t.test(x)$conf.int)
sh21_conf=apply(h2az[801:1000],2,function(x) t.test(x)$conf.int)
sh22_conf=apply(h2az[1001:1200],2,function(x) t.test(x)$conf.int)

#

library(gplots)
library(ggplot2)
library(Cairo)

c.min=min(nt1_conf,nt2_conf,sh11_conf,sh12_conf,sh21_conf,sh22_conf)
c.max=max(nt1_conf,nt2_conf,sh11_conf,sh12_conf,sh21_conf,sh22_conf)
pdf("ATAC_seq_TSS.pdf")
plot(nt1,type="l",ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="CPM Normalized ATAC-Seq Tags",col='#ef4040',lwd=3)
Axis(side=1, labels=c("-2 KB","TSS","2 KB"),at=c(1,99,200))
#polygon( c(1:200,rev(c(1:200)) ),c(nt1_conf[1,],rev(nt1_conf[2,])), col = alpha('#ef4040',.3), border = NA)
                
lines(nt2,col="#f26565",lwd=3)
#polygon( c(1:200,rev(c(1:200)) ),c(nt2_conf[1,],rev(nt2_conf[2,])), col = alpha('#f26565',.3), border = NA)
                
lines(sh11,col="#77dd77",lwd=3)
#polygon( c(1:200,rev(c(1:200)) ),c(sh11_conf[1,],rev(sh11_conf[2,])), col = alpha('#77dd77',.3), border = NA)
                
lines(sh12,col="#03c03c",lwd=3)
#polygon( c(1:200,rev(c(1:200)) ),c(sh12_conf[1,],rev(sh12_conf[2,])), col = alpha('#03c03c',.3), border = NA)
                
lines(sh21,col="#779ecb",lwd=3)
#polygon( c(1:200,rev(c(1:200)) ),c(sh21_conf[1,],rev(sh21_conf[2,])), col = alpha('#779ecb',.3), border = NA)
                
lines(sh22,col="#1195c6",lwd=3)
#polygon( c(1:200,rev(c(1:200)) ),c(sh22_conf[1,],rev(sh22_conf[2,])), col = alpha('#1195c6',.3), border = NA)
                
legend("topright", legend=c("shNT","shH2AFV #1","shH2AFV #2"), fill=c('#ef4040','#77dd77','#779ecb'), bty = "n",cex = 1.3)
abline(v=99,lty=2)
dev.off()

require(xlsx)

rna = read.xlsx("TFs_downreg_prom.xlsx", sheetName = "Sheet1")
rnae =  rna[,c(1,5)]
rnae[,1] = toupper(rnae[,1])
rnae[,1] = gsub("\\_.+","",rnae[,1],perl=T)

atac = read.table("tab_1287439080", sep="\t",header=T)
atac[,1] = toupper(atac[,1])


rnae100 = head(rnae[order(rnae[,2]),],100)
atac100 = head(atac[order(atac[,4],decreasing=T),],100)

table(rnae100 %in% atac100)

unique(rnae100[rnae100[,1] %in% atac100[,1],1])

library(VennDiagram)

pdf("vennDiagram.pdf")
grid.newpage()
draw.pairwise.venn(area1 = 100, area2 =100, cross.area = 16, 
category = c("Top 100 enriched TF\nin RNA-Seq", "Top 100 enriched TF\nin ATAC-Seq"), 
lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.05, 2),cat.cex = c(1.8,1.8),cex= 2.7)
    dev.off()
    

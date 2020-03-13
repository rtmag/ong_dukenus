library(ComplexHeatmap)
options(scipen=999)
library(graphics)
library(gplots)
library(factoextra)
library(RColorBrewer)
library(survival)
library(survminer)
library(gridExtra)
library(grid)


#READ EXPRESSION TABLE
GBMLGG <- read.table("2019-08-21_TCGA_GBMLGG_expression.txt", header=TRUE, row.names=1)
GBMLGG <- t(GBMLGG)

#READ PHENO TABLE
GBMLGG.pheno <- read.table("2019-08-21_TCGA_GBMLGG_pheno.txt", header=TRUE, row.names=1)

#CHECK GRADE OF ASTRO AND GBM
ctrack = data.frame(hist=as.character(GBMLGG.pheno$Histology),grade=as.character(GBMLGG.pheno$Grade))
table(ctrack[ctrack$hist=="Astrocytoma",2])
table(ctrack[ctrack$hist=="GBM",2])

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(GBMLGG.pheno$Histology)!="Oligoastrocytoma" & 
     as.character(GBMLGG.pheno$Histology)!="Oligodendroglioma" & 
     !is.na(GBMLGG.pheno$Grade) & !is.na(GBMLGG.pheno$Histology)
GBMLGG_astro_gbm <- GBMLGG[,ix]
GBMLGG_astro_gbm.pheno <- GBMLGG.pheno[ix,]


#GET SIGNATURE GENES ONLY
GBMLGG_astro_gbm_sig <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% signature$Gene.symbol,]

track= signature[match(rownames(GBMLGG_astro_gbm_sig),signature[,1]),2]

#create label track with grade
ctrack = as.character(GBMLGG_astro_gbm.pheno$Grade)
ctrack2 = as.character(GBMLGG_astro_gbm.pheno$Subtype.original)

GBMLGG_sig_centered = GBMLGG_astro_gbm_sig - rowMeans(GBMLGG_astro_gbm_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 6] = 6
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-6)] = -6

column_ha = HeatmapAnnotation(Grade = ctrack,Subtype = ctrack2,
                              col = list(Grade = c("II" = "black", "III" = "grey", "IV" = "red"),
                                      Subtype= c("Classical"="darkgreen","G-CIMP"="gray","IDHmut-codel"="brown",
                                                 "IDHmut-non-codel"="darkblue",
                                                  "IDHwt"="darkred","Mesenchymal"="purple","Neural"="orange","Proneural"="black")))
row_ha = rowAnnotation(Signature = track,show_annotation_name = FALSE,
              col = list(Signature = c("GenerationOfNeurons" = "#bae1ff","EGFR" = "#ffb3ba", "Mesenchymal" = "#baffc9")))

ht=Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",row_split =track)


hc <- as.hclust( column_dend(ht) )
groups=cutree( hc, k=3 )

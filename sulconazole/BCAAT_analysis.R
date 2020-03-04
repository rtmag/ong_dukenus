
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



# READ AND PARSE SIGNATURES
signature <- read.csv("3_pathway_signature.csv",stringsAsFactors=FALSE)
signature[signature$Geneset_for_CMAP_analysis=="KOBAYASHI_EGFR_SIGNALING_24HR_DN",2] <- "EGFR"
signature[signature$Geneset_for_CMAP_analysis=="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",2] <- "Mesenchymal"
signature[signature$Geneset_for_CMAP_analysis=="generation of neurons",2] <- "GenerationOfNeurons"

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


codel =  as.character(GBMLGG_astro_gbm.pheno$Chr.1p_19q.codeletion)
cogain =  as.character(GBMLGG_astro_gbm.pheno$Chr.19_20.co_gain)
IDH.status =  as.character(GBMLGG_astro_gbm.pheno$IDH.status)

codel[is.na(codel)] = "NA"
cogain[is.na(cogain)] = "NA"
IDH.status[is.na(IDH.status)] = "NA"

cogain[cogain=="Gain chr 19/20"] = "gain"
cogain[cogain=="No chr 19/20 gain"] = "no-gain"


column_ha = HeatmapAnnotation(Cluster = as.character(groups),
                              coDeletion.1p_19q=codel,
                              IDH.status=IDH.status,
                              Grade = ctrack,
                              col = list( Cluster = c("1"="purple","2"="orange","3"="blue"),
                                         coDeletion.1p_19q = c("codel"="red","non-codel"="grey","NA"="black"),
                                         IDH.status = c("WT"="grey","Mutant"="red","NA"="black"),
                                         Grade = c("II" = "grey", "III" = "orange", "IV" = "red")
                                        )
                             )

ht1 = Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_split = 3,
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)
################################################################################################################################
# BCAAT TCGA
bcca<-read.table("~/Desktop/ong/BCCA_catabolic_enzyme_genes.txt")
bcca<-as.character(bcca[,1])

bcca_mat <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% bcca,]

bcca_mat_centered = bcca_mat - rowMeans(bcca_mat)
bcca_mat_centered[bcca_mat_centered >= 6] = 6
bcca_mat_centered[bcca_mat_centered <= (-6)] = -6


column_ha_bcca = HeatmapAnnotation(Cluster = as.character(groups),
                              coDeletion.1p_19q=codel,
                              IDH.status=IDH.status,
                              Grade = ctrack,
                              col = list( Cluster = c("1"="purple","2"="orange","3"="blue"),
                                         coDeletion.1p_19q = c("codel"="red","non-codel"="grey","NA"="black"),
                                         IDH.status = c("WT"="grey","Mutant"="red","NA"="black"),
                                         Grade = c("II" = "grey", "III" = "orange", "IV" = "red")
                                        )
                             )


pdf("~/Desktop/ong/TCGA_BCCA.pdf",width=6,height=8)
ht1 = Heatmap(GBMLGG_sig_centered,height =22,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_split = 3,
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)

ht2 = Heatmap(bcca_mat_centered,height = 36,
show_row_names = TRUE,show_column_names = FALSE,name = "BCCA genes",row_dend_reorder = T,
column_title_side = "bottom", row_title="", row_title_side = "right",row_names_gp = gpar(fontsize = 8),
        clustering_distance_rows = "pearson")

ht_list = ht1 %v% ht2 
draw(ht_list)
dev.off()

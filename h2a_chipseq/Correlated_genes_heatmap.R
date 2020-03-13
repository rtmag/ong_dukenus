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
GBMLGG <- read.table("~/CSI/ong/public/2019-08-21_TCGA_GBMLGG_expression.txt", header=TRUE, row.names=1)
GBMLGG <- t(GBMLGG)

#READ PHENO TABLE
GBMLGG.pheno <- read.table("~/CSI/ong/public/2019-08-21_TCGA_GBMLGG_pheno.txt", header=TRUE, row.names=1)

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
GBMLGG_astro_gbm_sig <- GBMLGG_astro_gbm[rownames(GBMLGG_astro_gbm) %in% c("H2AFV","TRIM24","ATAD2","EZH2"),]

#create label track with grade
ctrack = as.character(GBMLGG_astro_gbm.pheno$Grade)
ctrack2 = as.character(GBMLGG_astro_gbm.pheno$Subtype.original)

GBMLGG_sig_centered = GBMLGG_astro_gbm_sig - rowMeans(GBMLGG_astro_gbm_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 6] = 6
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-6)] = -6


ht=Heatmap(t(scale(t(GBMLGG_astro_gbm_sig),center=F)), 
           row_dend_reorder = T, 
           column_dend_reorder = T,
           clustering_distance_columns = "pearson", 
           clustering_distance_rows = "pearson")


hc <- as.hclust( column_dend(ht) )
groups=cutree( hc, k=2 )


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
                              col = list( Cluster = c("1"="#3399cc","2"="#CC5500"),
                                         coDeletion.1p_19q = c("codel"="red","non-codel"="grey","NA"="black"),
                                         IDH.status = c("WT"="grey","Mutant"="red","NA"="black"),
                                         Grade = c("II" = "grey", "III" = "orange", "IV" = "red")
                                        )
                             )

pdf("H2AFV_correlated_epigeneticRegulators.pdf",width=6,height=3)
 Heatmap(t(scale(t(GBMLGG_astro_gbm_sig),center=F)),
show_row_names = TRUE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = T,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="", row_title_side = "right",
bottom_annotation = column_ha, clustering_distance_columns = "pearson",column_split = 2,
        clustering_distance_rows = "pearson",show_row_dend = FALSE)
dev.off()


svg("H2AFV_correlated_epigeneticRegulators.svg",width=6,height=3)
 Heatmap(t(scale(t(GBMLGG_astro_gbm_sig),center=F)),
show_row_names = TRUE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = T,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="", row_title_side = "right",
bottom_annotation = column_ha, clustering_distance_columns = "pearson",column_split = 2,
        clustering_distance_rows = "pearson",show_row_dend = FALSE)
dev.off()


########################################################################################################################
########################################################################################################################
# Survival curves

library(RTCGA.clinical)

clinical <- data.frame(times = GBMLGG_astro_gbm.pheno$survival,
                       bcr_patient_barcode = rownames(GBMLGG_astro_gbm.pheno),
                       patient.vital_status = as.numeric(GBMLGG_astro_gbm.pheno$status),
                       signature = as.factor(groups)) # from the tree cut above.
# alive=0 and dead=1
pdf("survival_TCGA-GBM-LGG_H2AFV_correlated_epigeneticRegulators.pdf")
kmTCGA(clinical, explanatory.names="signature",  pval = TRUE,conf.int = FALSE, risk.table=TRUE ,palette = c("#3399cc","#CC5500"))
dev.off()

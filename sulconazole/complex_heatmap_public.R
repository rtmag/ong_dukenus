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

grid.ftable <- function(d, padding = unit(4, "mm"), ...) {

  nc <- ncol(d)
  nr <- nrow(d)

  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))

  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")

  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)

  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)

  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding

  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))

  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")

  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)

  grid.draw(g)
  invisible(g)
}


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

#column_ha = HeatmapAnnotation(Grade = ctrack,Subtype = ctrack2,Cluster = as.character(groups),
#                              col = list(Grade = c("II" = "black", "III" = "grey", "IV" = "red"),
#                                      Subtype= c("Classical"="darkgreen","G-CIMP"="gray","IDHmut-codel"="brown",
#                                                 "IDHmut-non-codel"="darkblue",
#                                                  "IDHwt"="darkred","Mesenchymal"="purple","Neural"="orange","Proneural"="black"),
#                                         Cluster = c("1"="purple","2"="orange","3"="blue")
#                                        )
#                             )

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

#pdf("TCGA_complexHeatmap.pdf")
png("TCGA_complexHeatmap.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="TCGA GBM-LGG Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_split = 3,
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)
dev.off()

clinical <- data.frame(times = GBMLGG_astro_gbm.pheno$survival,
                       bcr_patient_barcode = rownames(GBMLGG_astro_gbm.pheno),
                       patient.vital_status = as.numeric(GBMLGG_astro_gbm.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1


res <- pairwise_survdiff(Surv(times, patient.vital_status) ~ signature, data = clinical)


#pdf("Survival_curves_TCGA.pdf")
png("Survival_curves_TCGA.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=TRUE,return.survfit = F,
       palette = c("purple","orange","blue"))
dev.off()

pdf("Survival_pvalues_TCGA.pdf")
 grid.ftable(formatC(res$p.value, format = "e", digits = 2),x=.14,y=.14,
           gp = gpar(fill = rep(c("grey90", "grey95"), each = 6)))
dev.off()
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

remb <- read.table("2019-08-21_Rembrandt_expression.txt", 
                     header=TRUE, row.names=1)
remb <- t(remb)

remb.pheno <- read.table("2019-08-21_Rembrandt_pheno.txt", 
                     header=TRUE, row.names=1)

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(remb.pheno$Histology)!="Mixed glioma" & 
     as.character(remb.pheno$Histology)!="Oligodendroglioma" & 
     as.character(remb.pheno$Histology)!="Non-tumor" & 
     as.character(remb.pheno$Histology)!="Unknown" & 
     as.character(remb.pheno$Grade)!="I" & 
     !is.na(remb.pheno$Grade) & !is.na(remb.pheno$Histology)
remb_imp <- remb[,ix]
remb_imp.pheno <- remb.pheno[ix,]

#GET SIGNATURE GENES ONLY
remb_imp_sig <- remb_imp[rownames(remb_imp) %in% signature$Gene.symbol,]

track= signature[match(rownames(remb_imp_sig),signature[,1]),2]

ctrack = as.character(remb_imp.pheno$Grade)
ctrack2 = as.character(remb_imp.pheno$Subtype_Verhaak_2010)
ctrack2[is.na(ctrack2)]="NA"

column_ha = HeatmapAnnotation(Grade = ctrack,Subtype = ctrack2,
                              col = list(Grade = c("II" = "black", "III" = "orange", "IV" = "red"),
                                      Subtype= c("Classical"="darkgreen",
                                                "Mesenchymal"="purple","Neural"="darkred","Proneural"="darkblue","NA"="grey")))
row_ha = rowAnnotation(Signature = track,show_annotation_name = FALSE,
              col = list(Signature = c("GenerationOfNeurons" = "#bae1ff","EGFR" = "#ffb3ba", "Mesenchymal" = "#baffc9")))

GBMLGG_sig_centered = remb_imp_sig - rowMeans(remb_imp_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 4] = 4
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-4)] = -4


ht=Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Rembrandt patient samples", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson")

hc <- as.hclust( column_dend(ht) )
groups=cutree( hc, k=4 )

column_ha = HeatmapAnnotation(Cluster = as.character(groups),Grade = ctrack,
                              col = list( Cluster = c("1"="purple","2"="orange","3"="blue","4"="darkgreen"),
                                         Grade = c("II" = "grey", "III" = "orange", "IV" = "red")
                                        )
                             )
#pdf("Rembrandt_complexHeatmap.pdf")
png("Rembrandt_complexHeatmap.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Rembrandt patient samples", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_split = 4,
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)
dev.off()


library(RTCGA.clinical)

clinical <- data.frame(times = remb_imp.pheno$survival,
                       bcr_patient_barcode = rownames(remb_imp.pheno),
                       patient.vital_status = as.numeric(remb_imp.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=TRUE,
       palette = c("purple","orange","blue","darkgreen"))


res <- pairwise_survdiff(Surv(times, patient.vital_status) ~ signature, data = clinical)


#pdf("Survival_curves_Rembrandt.pdf")
png("Survival_curves_Rembrandt.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=TRUE,return.survfit = F,
       palette = c("purple","orange","blue","darkgreen"))
dev.off()

pdf("Survival_pvalues_Rembrandt.pdf")
 grid.ftable(formatC(res$p.value, format = "e", digits = 2),x=.14,y=.14,
           gp = gpar(fill = rep(c("grey90", "grey95"), each = 6)))
dev.off()
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

#READ EXPRESSION TABLE
grav <- read.table("2019-08-27_Gravendeel_expression.txt", header=TRUE, row.names=1)
grav <- t(grav)

#READ PHENO TABLE
grav.pheno <- read.table("2019-08-27_Gravendeel_pheno.txt", header=TRUE, row.names=1)

#CHECK GRADE OF ASTRO AND GBM
ctrack = data.frame(hist=as.character(grav.pheno$Histology),grade=as.character(grav.pheno$Grade))
table(ctrack[ctrack$hist=="Astrocytoma",2])
table(ctrack[ctrack$hist=="GBM",2])

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(grav.pheno$Histology)!="Mixed glioma" & 
     as.character(grav.pheno$Histology)!="Oligodendroglioma" & 
     as.character(grav.pheno$Histology)!="Non-tumor" & 
     as.character(grav.pheno$Histology)!="Pilocytic Astrocytoma" & 
     as.character(grav.pheno$Grade)!="I" & 
     !is.na(grav.pheno$Grade) & !is.na(grav.pheno$Histology)
grav.imp <- grav[,ix]
grav.imp.pheno <- grav.pheno[ix,]

#GET SIGNATURE GENES ONLY
grav_sig <- grav.imp[rownames(grav.imp) %in% signature$Gene.symbol,]

track= signature[match(rownames(grav_sig),signature[,1]),2]

ctrack = as.character(grav.imp.pheno$Grade)
ctrack2 = as.character(grav.imp.pheno$Subtype_Verhaak_2010)
ctrack2[is.na(ctrack2)]="NA"

LOH.1p =  as.character(grav.imp.pheno$LOH_1p)
LOH.19q =  as.character(grav.imp.pheno$LOH_19q)
IDH1.status =  as.character(grav.imp.pheno$IDH1_status)

LOH.1p[is.na(LOH.1p)] = "NA"
LOH.19q[is.na(LOH.19q)] = "NA"
IDH.status[is.na(IDH.status)] = "NA"

LOH.1p[LOH.1p=="PARTIAL"] = "YES"
LOH.19q[LOH.19q=="PARTIAL"] = "YES"

IDH1.status[IDH1.status=="Wild_type"] = "WT"
IDH1.status[IDH1.status=="Mut"] = "Mutant"

column_ha = HeatmapAnnotation(Grade = ctrack,Subtype = ctrack2,
                              col = list(Grade = c("II" = "black", "III" = "grey", "IV" = "red"),
                                      Subtype= c("Classical"="darkgreen",
                                                "Mesenchymal"="purple","Neural"="darkred","Proneural"="darkblue","NA"="grey")))
row_ha = rowAnnotation(Signature = track,show_annotation_name = FALSE,
              col = list(Signature = c("EGFR" = "#ffb3ba", "Mesenchymal" = "#baffc9", "GenerationOfNeurons" = "#bae1ff")))

GBMLGG_sig_centered = grav_sig - rowMeans(grav_sig)
GBMLGG_sig_centered[GBMLGG_sig_centered >= 4] = 4
GBMLGG_sig_centered[GBMLGG_sig_centered <= (-4)] = -4


ht=Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Gravendeel patient samples", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)

hc <- as.hclust( column_dend(ht) )
groups=cutree( hc, k=3 )

column_ha = HeatmapAnnotation(Cluster = as.character(groups),LOH.1p=LOH.1p,LOH.19q=LOH.19q,IDH1.status=IDH1.status,Grade = ctrack,
                              col = list( Cluster = c("1"="purple","2"="orange","3"="blue","4"="darkgreen"),
                                         LOH.1p = c("YES"="red","NO"="grey","NA"="black"),
                                         LOH.19q = c("YES"="red","NO"="grey","NA"="black"),
                                         IDH1.status = c("WT"="grey","Mutant"="red","NA"="black"),
                                         Grade = c("II" = "grey", "III" = "orange", "IV" = "red")
                                        )
                             )
pdf("complexHeatmap_Gravendeel.pdf")
Heatmap(GBMLGG_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Rembrandt patient samples", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_split = 3,
        clustering_distance_rows = "pearson",row_split =track,show_row_dend = FALSE)
dev.off()


library(RTCGA.clinical)

clinical <- data.frame(times = remb_imp.pheno$survival,
                       bcr_patient_barcode = rownames(remb_imp.pheno),
                       patient.vital_status = as.numeric(remb_imp.pheno$status),
                       signature = as.factor(groups))
# alive=0 and dead=1

kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=TRUE,
       palette = c("purple","orange","blue","darkgreen"))


res <- pairwise_survdiff(Surv(times, patient.vital_status) ~ signature, data = clinical)


pdf("Survival_curves_Rembrandt.pdf")

kmTCGA(clinical, explanatory.names="signature",  pval = FALSE,conf.int = FALSE, risk.table=TRUE,return.survfit = F,
       palette = c("purple","orange","blue","darkgreen"))
dev.off()

pdf("Survival_pvalues_Rembrandt.pdf")
 grid.ftable(formatC(res$p.value, format = "e", digits = 2),x=.14,y=.14,
           gp = gpar(fill = rep(c("grey90", "grey95"), each = 6)))
dev.off()
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##########################

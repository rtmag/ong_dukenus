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
###

#READ EXPRESSION TABLE
dataset <- read.table("2019-09-04_Freije_expression.txt", header=TRUE, row.names=1)
dataset <- t(dataset)

#READ PHENO TABLE
dataset.pheno <- read.table("2019-09-04_Freije_pheno.txt", header=TRUE, row.names=1)

#CHECK GRADE OF ASTRO AND GBM
ctrack = data.frame(hist=as.character(dataset.pheno$Histology),grade=as.character(dataset.pheno$Grade))
table(ctrack[ctrack$hist=="Astrocytoma",2])
table(ctrack[ctrack$hist=="GBM",2])

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = as.character(dataset.pheno$Histology)!="Mixed glioma" & 
     as.character(dataset.pheno$Histology)!="Oligodendroglioma" & 
     !is.na(dataset.pheno$Grade) 
dataset_astro_gbm <- dataset[,ix]
dataset_astro_gbm.pheno <- dataset.pheno[ix,]

#GET SIGNATURE GENES ONLY
dataset_astro_gbm_sig <- dataset_astro_gbm[rownames(dataset_astro_gbm) %in% signature$Gene.symbol,]

track= signature[match(rownames(dataset_astro_gbm_sig),signature[,1]),2]

#create label track with grade
ctrack = as.character(dataset_astro_gbm.pheno$Grade)
ctrack2 = as.character(dataset_astro_gbm.pheno$Subtype.original)

dataset_sig_centered = dataset_astro_gbm_sig - rowMeans(dataset_astro_gbm_sig)
dataset_sig_centered[dataset_sig_centered >= 6] = 6
dataset_sig_centered[dataset_sig_centered <= (-6)] = -6

column_ha = HeatmapAnnotation(Grade = ctrack,
                              col = list(Grade = c("II" = "black", "III" = "grey", "IV" = "red")
                                      ))
row_ha = rowAnnotation(Signature = track,show_annotation_name = FALSE,
              col = list(Signature = c("GenerationOfNeurons" = "#bae1ff","EGFR" = "#ffb3ba", "Mesenchymal" = "#baffc9")))

png("complexHeatmap_Freije_Forced3Groups_Kmeans3.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(dataset_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Freije Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_km = 3, column_km_repeats = 100,
        clustering_distance_rows = "pearson",row_split =track)
dev.off()

png("complexHeatmap_Freije_noKmeans.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(dataset_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Freije Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",row_split =track)
dev.off()
################################################################################################
################################################################################################
################################################################################################
# Phillips

#READ EXPRESSION TABLE
dataset <- read.table("2019-09-04_Phillips_expression.txt", header=TRUE, row.names=1)
dataset <- t(dataset)

#READ PHENO TABLE
dataset.pheno <- read.table("2019-09-04_Phillips_pheno.txt", header=TRUE, row.names=1)

#CHECK GRADE OF ASTRO AND GBM
ctrack = data.frame(hist=as.character(dataset.pheno$Histology),grade=as.character(dataset.pheno$Grade))
table(ctrack[ctrack$hist=="Astrocytoma",2])
table(ctrack[ctrack$hist=="GBM",2])

#REMOVE other than ASTRO AND GBM AND keep non NA from expression and pheno tables
ix = !is.na(dataset.pheno$Grade) 
dataset_astro_gbm <- dataset[,ix]
dataset_astro_gbm.pheno <- dataset.pheno[ix,]

#GET SIGNATURE GENES ONLY
dataset_astro_gbm_sig <- dataset_astro_gbm[rownames(dataset_astro_gbm) %in% signature$Gene.symbol,]

track= signature[match(rownames(dataset_astro_gbm_sig),signature[,1]),2]

#create label track with grade
ctrack = as.character(dataset_astro_gbm.pheno$Grade)
ctrack2 = as.character(dataset_astro_gbm.pheno$Subtype.original)

dataset_sig_centered = dataset_astro_gbm_sig - rowMeans(dataset_astro_gbm_sig)
dataset_sig_centered[dataset_sig_centered >= 6] = 6
dataset_sig_centered[dataset_sig_centered <= (-6)] = -6

column_ha = HeatmapAnnotation(Grade = ctrack,
                              col = list(Grade = c("II" = "black", "III" = "grey", "IV" = "red")
                                      ))
row_ha = rowAnnotation(Signature = track,show_annotation_name = FALSE,
              col = list(Signature = c("GenerationOfNeurons" = "#bae1ff","EGFR" = "#ffb3ba", "Mesenchymal" = "#baffc9")))

png("complexHeatmap_Phillips_Forced3Groups_Kmeans3.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(dataset_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Phillips Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",column_km = 3, column_km_repeats = 100,
        clustering_distance_rows = "pearson",row_split =track)
dev.off()

png("complexHeatmap_Phillips_noKmeans.png",width= 7.25,
  height= 7.25,units="in",
  res=1200,pointsize=4)
Heatmap(dataset_sig_centered,
show_row_names = FALSE,show_column_names = FALSE,name = "Expression",row_dend_reorder = T, column_dend_reorder = F,
column_title="Phillips Patients", column_title_side = "bottom", row_title="Gene Signature", row_title_side = "right",
bottom_annotation = column_ha, right_annotation = row_ha,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",row_split =track)
dev.off()

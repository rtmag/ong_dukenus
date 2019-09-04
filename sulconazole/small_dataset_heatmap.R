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

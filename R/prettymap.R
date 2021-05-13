library(colorspace)
## ===========pretty t-SNE===========
prettyTSNE <- function(tsne_matrix, object, fcol = "markers",
                       main = "", orgOrder = getMarkerClasses(object),
                       mainCol = getStockcol(), 
                       outlineCol = darken(getStockcol()), 
                       pchUn = 21, lwd = 1,
                       ...) {
  fData(object)[, "tmp"] <- fData(object)[, fcol]
  setUnknowncol(NULL)
  setUnknowncol(paste0(getUnknowncol(), 90))
  # setStockcol(NULL)
  plot2D(tsne_matrix, method = "none",  methargs = list(object),
         fcol = NULL, cex = 1, pch = pchUn, grid = FALSE, main = main, 
         cex.axis = 1.5,
         cex.lab = 1.5, lwd = lwd,
         bg = getUnknowncol(), col = darken(getUnknowncol()), ...)
  cl <- getMarkerClasses(object, fcol = "tmp")
  names(cl) <- getMarkerClasses(object, fcol = "tmp")
  cl <- cl[orgOrder]
  col1 <- outlineCol
  col2 <- mainCol
  names(col1) <- cl
  names(col2) <- cl
  for (i in seq(cl)) {
    ind <- which(fData(object)[, "tmp"] == cl[i])
    points(tsne_matrix[ind, ], col = col1[orgOrder][i], pch = 21, 
           bg = col2[orgOrder][i], cex = 1.5, lwd = lwd)
  }
}

prettyTSNE_overlay <- function(tsne_matrix, object, fcol = "markers",
                               main = "", orgOrder = getMarkerClasses(object),
                               mainCol = paste0(getStockcol(), 70), 
                               outlineCol = darken(getStockcol()), ...) {
  fData(object)[, "tmp"] <- fData(object)[, fcol]
  setUnknowncol(NULL)
  setUnknowncol(paste0(getUnknowncol(), 70))
  plot2D(tsne_matrix, method = "none",  methargs = list(object),
         fcol = NULL, pch = 19, grid = FALSE, main = main, 
         cex.axis = 1.5,
         cex.lab = 1.5, 
         col = getUnknowncol(), ...)
  cl <- getMarkerClasses(object, fcol = "tmp")
  names(cl) <- getMarkerClasses(object, fcol = "tmp")
  cl <- cl[orgOrder]
  col1 <- outlineCol
  col2 <- mainCol
  names(col1) <- cl
  names(col2) <- cl  
  for (i in seq(cl)) {
    ind <- which(fData(object)[, "tmp"] == cl[i])
    points(tsne_matrix[ind, ], col = col1[orgOrder][i], pch = 21, 
           bg = col2[orgOrder][i], cex = 1.5)
  }
}
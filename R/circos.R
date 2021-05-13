library("circlize")

## Create circos plot of translocations
customChord <- function(msnset1, 
                        msnset2, 
                        fcol1 = "localisation.pred", 
                        fcol2 = "localisation.pred", 
                        mrkCol1 = "markers",
                        mrkCol2 = "markers",
                        onlyMovers = TRUE,
                        labels = TRUE,
                        cols,
                        ...) {
                   
  ## make data.frame of translocations
  df <- makedf(msnset1, msnset2, 
               fcol1, fcol2, 
               mrkCol1,
               mrkCol2,
               onlyMovers)
  ## add colour scheme if not provided
  if (missing(cols)) {
    ll <- unique(c(levels(df[,1]), levels(df[,2])))
    grid.col <- segcols <- setNames(rainbow(length(ll)), ll)
  } else {
    grid.col <- cols
  }
  ## create circos
  chordDiagram(df, annotationTrack = "grid", 
               preAllocateTracks = 1, 
               grid.col = grid.col,
               directional = 1, 
               direction.type = c("diffHeight", "arrows"), 
               link.arr.type = "big.arrow", ...)
  ## annotate tracking regions and customise
  if (labels) {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), 
                  ylim[1] + .1, 
                  sector.name, 
                  facing = "clockwise",
                  niceFacing = TRUE, 
                  adj = c(-0.5, 0.5),
                  cex = 1.2,
                  col=grid.col[sector.name],
                  font = 2)
      circos.axis(h = "top",
                  labels.cex = .8,
                  major.tick.length = 1,
                  sector.index = sector.name,
                  track.index = 2)
    }, bg.border = NA)
  } else {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.axis(h = "top",
                  labels.cex = .8,
                  major.tick.length = 1,
                  sector.index = sector.name,
                  track.index = 2)
    }, bg.border = NA)
  }
  
  circos.clear()
}

## get count data (for organelle transfers between conditions) 
## from msnset
makedf <- function(msnset1, msnset2, 
                   fcol1, fcol2, 
                   mrkCol1 = "markers",
                   mrkCol2 = "markers",
                   # includeMarkers = FALSE, ## include markers
                   onlyMovers = TRUE) {  
  
  ## get levels to convert localisation info to factors
  ## needed for dplyr and .drop = FALSE
  .fct.lev1 <- getMarkerClasses(msnset1, mrkCol1)
  .fct.lev2 <- getMarkerClasses(msnset2, mrkCol2)
  fct.lev <- union(.fct.lev1, .fct.lev2)
  .locnam1 <- names(table(getMarkers(msnset1, fcol1, verbose = FALSE)))
  .locnam2 <- names(table(getMarkers(msnset2, fcol2, verbose = FALSE)))
  locnam <- union(.locnam1, .locnam2)
  locnam <- union(locnam, fct.lev)
  msnset1 <- unknownMSnSet(msnset1, fcol = mrkCol1)
  msnset2 <- unknownMSnSet(msnset2, fcol = mrkCol2)
  dat1 <- factor(fData(msnset1)[, fcol1], locnam)
  dat2 <- factor(fData(msnset2)[, fcol2], locnam)
  stopifnot(length(dat1) == length(dat2))
  dat <- data.frame(x = dat1, y = dat2)
  dat$z <- 1
  datdf <- dat %>% group_by(x, y, .drop = FALSE) %>% 
    dplyr:::summarise(count=sum(z), .groups = NULL)
  if (onlyMovers) {
    torm <- which(datdf$x == datdf$y)
    datdf <- datdf[-torm, ]
  }
  datdf <- as.data.frame(datdf)
}
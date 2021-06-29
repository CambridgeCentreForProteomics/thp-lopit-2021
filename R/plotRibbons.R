tcol <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

ribbonPlot <- function(profs, col = "red", q = c(0, 1), size = 1.6, ...) {
  # par(mar = c(5, 4, 1, 1), oma = c(0, 0, 0, 0))
  ylim <- range(profs)
  n <- nrow(profs)
  m <- ncol(profs)
  fracInds <- as.numeric(colnames(profs))
  # fracs <- colnames(profs)
  
  ## get quantiles for each fraction 
  quants <- apply(profs, MARGIN = 2, function(x) quantile(x, q))  # max and min for unknowns
  meanProf <- apply(profs, MARGIN = 2, function(x) mean(x))  # max and min for unknowns
  
  bound_low <- quants[1, ]
  bound_high <- quants[2, ]
  
  ## make polygon plots
  plot(0, ylim = ylim, xlim = c(0, 24),
       type = "n", xaxt = "n", yaxt = "n",  
       ylab = "Normalised intensities", xlab = "Time (hours)", 
       cex.axis = size, cex.lab = size,
       cex.lab = size, ...)
  v_x <- axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2), cex.axis = size)
  v_y <- axis(2, cex.axis = size)
  abline(v = v_x, h = v_y, lty = "dotted", col = "lightgray", lwd = 1)
  # mtext("Time (hours)", side=1, line=3, cex = size)
  
  ## plot unknowns
  matlines(fracInds, t(profs),
           col = tcol(col, 85), lty = 1, lwd = 1, type = "l")
  
  
  polygon(c(fracInds, rev(fracInds)), 
          c(quants[2, as.vector(sapply(fracInds, as.character))], 
            rev(quants[1, as.vector(sapply(fracInds, as.character))])),
          col = tcol(col, 90), border = FALSE)
  
  matlines(fracInds, meanProf,
           col = col,
           lty = 1, lwd = 2,
           type = "l")
  
  
}

## =====================PROFILES plot========================
## ==========================================================




plotRibbons2 <- function(msnset,
                         fcol = "loc.pred",
                         mrkClasses,
                         grid = FALSE,
                         unknowns = FALSE,
                         mycols = mycol) {
  
  profs <- exprs(msnset)
  if(missing(mrkClasses)) {
    mrkClasses <- getMarkerClasses(msnset)
  }
  indMrk <- lapply(mrkClasses, function(x) which(fData(msnset)[, "loc.pred"] == x))
  
  # browser()
  par(mar = c(13, 4, 1, 1), oma = c(0, 0, 0, 0))
  ylim <- c(0.05,.4)
  n <- nrow(profs)
  m <- ncol(profs)
  fracs <- colnames(profs)
  ## check if there are replicates and if their are, create breaks in the plot
  # if (!is.null(pcol)) {
  #   repInfo <- unique(pd[, pcol])
  #   repNames <- vector("list", length(repInfo))
  #   ## get fraction names by replicate
  #   fracNames <- lapply(repInfo, function(z) colnames(profs)[pd$Experiment == z])
  #   fracInds <- lapply(fracNames, function(z) which(z == colnames(profs)))
  # } else {
  fracInds <- list(seq(colnames(profs)))
  # }
  ## get unknowns
  profs_un <- profs[which(fData(msnset)[, fcol] == "unknown"), ]
  ## get quantiles for each fraction in unknowns
  quants <- apply(profs_un, MARGIN = 2, function(x) 
    quantile(x, c(0, 1), na.rm = TRUE))  # max and min for unknowns
  bound_low <- quants[1, ]
  bound_high <- quants[2, ]
  ## get quantiles for subcellular classes
  mrkProfs <- lapply(indMrk, function(z) profs[z, , drop = FALSE])   # 5% and 95% quantiles for all other classes
  quants <- lapply(mrkProfs, function(z) apply(z, MARGIN = 2, function(x) 
    quantile(x, c(0.25, .75), na.rm = TRUE)))
  meanProfs <- lapply(mrkProfs, function(z) apply(z, 2, mean, na.rm = TRUE))
  ## make polygon plots
  plot(0, ylim = ylim, xlim = c(1, m),
       type = "n", xaxt = "n", yaxt = "n", xlab = "",
       ylab = "Intensities", cex.axis = 1.2,
       cex.lab = 1.2)
  v_x <- axis(1, at = 1:m, labels = fracs, las = 2, cex.axis = 1.2)
  v_y <- axis(2)
  if (grid) abline(v = v_x, h = v_y, lty = "dotted", col = "lightgray", lwd = 1)
  mtext("Fractions", side=1, line=8, cex = 1.2)
  ## update lines on plot according to zoom
  # feats <<- which(brushBounds$i & brushBounds$j)
  # namFeats <- names(feats)[which(names(feats) %in% rownames(profs_un))]
  
  ## plot unknowns
  if (unknowns) {
    invisible(lapply(fracInds, function(x)     # plot all unknowns as lines here
      matlines(x, t(profs_un),
               col = "grey90", lty = 1, lwd = 1, type = "l")
    ))
  }
  ## markers
  
  for (i in 1:length(indMrk)) {
    if (!is.na(indMrk[[i]][1])) {
      invisible(lapply(fracInds, function(x)     # don't plot all lines
        polygon(c(x, rev(x)),
                c(quants[[i]][2, x], rev(quants[[i]][1, x])),
                col = paste0(mycol[i], 60), border = FALSE)
      ))
      invisible(lapply(fracInds, function(z)     # plot the mean profile
        matlines(z, meanProfs[[i]][z],
                 col = mycol[i],
                 lty = 1, lwd = 1,
                 type = "l")))
    }
  }
  
} ## ----------- end of function----------------

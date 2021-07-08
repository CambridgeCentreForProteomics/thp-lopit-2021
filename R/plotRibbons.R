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
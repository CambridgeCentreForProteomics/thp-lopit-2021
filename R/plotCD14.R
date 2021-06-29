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


plotCD14 <- function(data, 
                     fcol = "markers",
                     org = "Plasma Membrane",
                     protein = "CD14",
                     org.col = "#D55E00"
                     ) {
  profs_un <- exprs(data)[fData(data)[, fcol] == "unknown", ]
  profs_pm <- exprs(data)[fData(data)[, fcol] == org, ]
  fracInds <- seq(ncol(profs_un))
  quants <- apply(profs_pm, MARGIN = 2, function(x) quantile(x, c(0, 1)))  # max and min for PM
  meanProf <- apply(profs_pm, MARGIN = 2, function(x) mean(x))  # max and min for unknowns
  if (!any(fvarLabels(data) == "GN")) {
    ind <- grep(protein, fData(data)[[grep("Protein.Desc", fvarLabels(data))[1]]])
  } else {
    ind <- which(fData(data)$GN == protein)
  }
  cd14 <- exprs(data)[ind, ]
  
  par(oma = c(10,1,1,1))
  plot(0, ylim = c(0, 1), xlim = c(0, ncol(data)),
       type = "n", xaxt = "n", yaxt = "n",  
       ylab = "Normalised intensities", xlab = "")
  v_x <- axis(1, at = 1:ncol(data), labels = sampleNames(data), 
              cex.axis = 1, las = 2)
  v_y <- axis(2, cex.axis = 1)
  
  ## plot unknowns
  matlines(fracInds, t(profs_un),
           col = "#EBEBEB", lty = 1, lwd = 1, type = "l")
  
  ## plot PM
  matlines(fracInds, t(profs_pm),
           col = paste0(org.col, 30), lty = 1, lwd = 1, type = "l")
  
  # polygon(c(fracInds, rev(fracInds)), 
  #         c(quants[2, fracInds], 
  #           rev(quants[1, fracInds])),
  #         col = paste0("#D55E00", 30), border = FALSE)
  
  # matlines(fracInds, meanProf,
  #          col = paste0(, 90),
  #          lty = 3, lwd = 2,
  #          type = "l")
  
  matlines(fracInds, cd14,
           col = "black",
           lty = 1, lwd = 1,
           type = "l")
  
}

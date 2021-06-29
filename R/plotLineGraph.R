plotLineGraph <- function(data, fun = c("mean", "median"), 
                         showBars = FALSE, linecol = "#C50000",
                         .ylab = "Release (pg/ml)\n",
                         themeSize = 20) {
  sum_dat <- gather(data, replicate, measurement, X1:X6, factor_key = TRUE)
  if (fun == "mean") {
    sum_dat <- summarySE(sum_dat, measurevar = "measurement", groupvars = "Time")
    top <- sum_dat$measurement + sum_dat$se
    bottom <- sum_dat$measurement - sum_dat$se
  } 
  else if (fun == "median") {
    sum_dat <- data.frame(cbind(Time = data$Time, t(apply(t(data[, -1]), 2, summary))))
    colnames(sum_dat)[which(colnames(sum_dat) == "Median")] <- "measurement"
    top <- sum_dat[,6]
    bottom <- sum_dat[,3]
  } 
  else (stop(message(paste("Summary stat specified must be mean or median"))))
  
  p <- ggplot(sum_dat, aes(x=Time, y=measurement)) + 
    geom_line(size = 1.3, color = linecol) +
    geom_point(size = 4, color = linecol) +
    scale_x_discrete(limits = seq(0, 24, by = 2)) +
    xlab("\nTime (hours)") +
    ylab(.ylab) +
    theme_classic(base_size = themeSize) +
    theme(axis.text = element_text(color = "black", size = themeSize)) 
  
  if (showBars) 
    p <- p + geom_errorbar(aes(ymin=top, ymax=bottom), width=.3, size = 1, col = linecol)
  # ylim(c(0, .25))
  return(p)
}
# pdf(file = "../figures_thp_june2021/figures/cd14-unst.pdf", width = 12, height = 7)
# data <- unst
# # sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# sampleNames(data) <- paste0(pData(data)$Fraction, "_", pData(data)$Tag, "_Set", pData(data)$Set, "_Rep_", pData(data)$Replicate)
# .sel <- which(fData(data)$markers == "Plasma Membrane")
# fData(data)$foi <- "unknown"  
# fData(data)$foi[which(fData(data)$GN == "CD14")] <- "CD14"
# fData(data)$foi[.sel] <- "Plasma Membrane"
# col <- c("black", "#D55E00")
# plotProfiles(data, fcol = "foi")
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/cd14-lps.pdf", width = 12, height = 7)
# data <- lps
# # sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# sampleNames(data) <- paste0(pData(data)$Fraction, "_", pData(data)$Tag, "_Set", pData(data)$Set, "_Rep_", pData(data)$Replicate)
# .sel <- which(fData(data)$markers == "Plasma Membrane")
# fData(data)$foi <- "unknown"  
# fData(data)$foi[which(fData(data)$GN == "CD14")] <- "CD14"
# fData(data)$foi[.sel] <- "Plasma Membrane"
# col <- c("black", "#D55E00")
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotProfiles(data, fcol = "foi")
# dev.off()

plotProfiles <- function(data,
                         fcol = "foi",
                         col = c("black", "#D55E00"),
                         reorderByRep = TRUE,
                         replicate.column.name = "Replicate"
                         ) {
  # browser()
  if (reorderByRep) {
    pd <- pData(data)
    # if (missing(replicate.column.name))
    #   stop(message(paste("replicate.column.name must be provided")))
    new_order <- order(pd[[replicate.column.name]])
    data <- data[, new_order]
    message(paste("Reordering data according to replicate levels"))
  }
  
  # ## is the replicate information provided?
  # chk <- grep(pattern = "rep", colnames(pd1), ignore.case = TRUE)
  # if (length(chk) > 0) {
  #   if (length(chk) != 1) {
  #     stop(paste0("More than one column in pData containing replicate information"))
  #   } else {
  #     repInfo <- pd1[, chk]
  #   }
  # } else {
  #   repInfo <- rep(1, ncol(data))
  # }
  
  repInfo = 1
  
  ## get data
  fd <- fData(data)
  pd <- pData(data)
  names(col) <- getMarkerClasses(data, fcol)
  col <- c(col, "unknown" = "darkgrey")
  data <- exprs(data)
  
              
  
  ## prep data for ggplot 
  .rn <- rownames(data)
  .cn <- colnames(data)
  plot_data <- data.frame(id = rep(.rn, ncol(data)),  
                          fraction = rep(.cn, each = nrow(data)), # variable
                          intensities = as.vector(data),  # value
                          rep = factor(rep(repInfo, each = nrow(data))),
                          mrk = rep(fd[, fcol], ncol(data)))
  plot_data <- within(plot_data, fraction <- factor(fraction, levels = colnames(data)))
  
  ## extract colours for organelles in the data 
  if (is.factor(plot_data$mrk)) 
    col <- col[levels(plot_data$mrk)]
  else
    col <- col[unique(plot_data$mrk)]
  
  ## line plot
  p <- ggplot() + 
    geom_line(data = subset(plot_data, mrk=="unknown"), 
              aes(x = fraction, y = intensities, group = id, color = mrk),
              alpha=0.01)  +
    geom_line(data = subset(plot_data, mrk=="Plasma Membrane"), 
              aes(x = fraction, y = intensities, group = id, color = mrk),
              alpha=0.2) +
    geom_line(data = subset(plot_data, mrk=="CD14"), 
              aes(x = fraction, y = intensities, group = id, color = mrk),
              alpha=0.7,
              size = .8)
  
  p <- p + 
    # scale_x_discrete(limits=fracLev, breaks = fracLev[seq(1,length(fracLev),by=2)]) +  # show only every other label on x-axis
    ylab("Normalised intensities") + xlab("") +
    scale_fill_manual(values = col, aesthetics = c("fill","colour")) +
    scale_color_manual(values = col, aesthetics = c("fill, colour")) +
    theme_minimal() +
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none", 
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
          # strip.text.x = element_text(size = 8, face="bold"),
          # panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8 + 2, colour = rgb(0, 0, 0)),
          axis.text.y = element_text(size = 8, colour = rgb(0, 0, 0))) 
    
  p

  
}




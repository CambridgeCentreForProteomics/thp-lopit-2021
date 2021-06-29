# pdf(file = "../figures_thp_june2021/figures/classProfiles_lps_TAGM_LINES.pdf", width = 12, height = 7)
# data <- lps
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# sampleNames(data) <- paste0(pData(data)$Fraction, "_", pData(data)$Tag, "_Set", pData(data)$Set, "_Rep_", pData(data)$Replicate)
# plotFacetProfiles(data, fcol = "localisation.pred", col = mycol, ribbons = FALSE,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_unst_TAGM_LINES.pdf", width = 12, height = 7)
# data <- unst
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "localisation.pred", col = mycol,ribbons = FALSE,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_lps_MRK_LINES.pdf", width = 12, height = 7)
# data <- lps
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "markers", col = mycol,ribbons = FALSE,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_unst_MRK_LINES.pdf", width = 12, height = 7)
# data <- unst
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "markers", col = mycol,ribbons = FALSE,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_lps_TAGM.pdf", width = 12, height = 7)
# data <- lps
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "localisation.pred", col = mycol,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_unst_TAGM.pdf", width = 12, height = 7)
# data <- unst
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "localisation.pred", col = mycol,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_lps_MRK.pdf", width = 12, height = 7)
# data <- lps
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "markers", col = mycol,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()
# 
# pdf(file = "../figures_thp_june2021/figures/classProfiles_unst_MRK.pdf", width = 12, height = 7)
# data <- unst
# sampleNames(data) <- paste0(pData(data)$Tag, "_Replicate_", pData(data)$Replicate, "_Set", pData(data)$Set)
# plotFacetProfiles(data, fcol = "markers", col = mycol,
#                   replicate.column.name = "Replicate", plotReps = FALSE)
# dev.off()

# make data frame for ggplot
plotFacetProfiles <- function(data, 
                              fcol,
                              replicate.column.name, 
                              col,
                              reorderByRep = TRUE,
                              plotReps = TRUE,
                              ribbons = TRUE,
                              ...) {
  if (reorderByRep) {
    pd <- pData(data)
    if (missing(replicate.column.name))
      stop(message(paste("replicate.column.name must be provided")))
    new_order <- order(pd[[replicate.column.name]])
    data <- data[, new_order]
    message(paste("Reordering data according to replicate levels"))
  }
  fd <- fData(data)
  pd <- pData(data)
  names(col) <- getMarkerClasses(data, fcol)
  data <- exprs(data)
  intensities = NULL
  mrk = NULL
  if (missing(replicate.column.name)) {
    # message(paste("Replicate information not provided, assuming 1 replicate only"))
    repInfo <- rep(1, ncol(data))
    reps <- FALSE
  } else {
    repInfo <- pd[, replicate.column.name]
    reps <- TRUE
  }
  
  ## prep data for ggplot 
  .rn <- rownames(data)
  .cn <- colnames(data)
  plot_data <- data.frame("id" = rep(.rn, ncol(data)),  
                          "fraction" = rep(.cn, each = nrow(data)), # variable
                          "intensities" = as.vector(data),  # value
                          "rep" = factor(rep(repInfo, each = nrow(data))),
                          "mrk" = rep(fd[, fcol], ncol(data)))
  plot_data <- within(plot_data, fraction <- factor(fraction, levels = colnames(data)))
  
  df <- plot_data %>% group_by(mrk, fraction, rep) %>%
    dplyr::summarise(min = min(intensities, na.rm = TRUE),
                     quant_05 = quantile(intensities, 0.25, na.rm = TRUE),
                     mean = mean(intensities, na.rm = TRUE),
                     quant_95 = quantile(intensities, 0.75, na.rm = TRUE),
                     max = max(intensities, na.rm = TRUE), .groups = "keep",
                     na.rm = TRUE)
  
  fracLev <- levels(df$fraction)
  if (plotReps) {
    repLev <- levels(df$rep)
    p <- ggplot()
    for(i in seq(repLev)){
      if (ribbons) {
        p <- p +
          geom_ribbon(data = subset(df, rep == repLev[i]),
                      mapping = aes(fraction, ymin=min, ymax=max, group = mrk, 
                                    color = NA, fill = mrk), 
                      alpha=0.3) +
          geom_line(data = subset(df, rep == repLev[i]),
                    mapping = aes(fraction, mean, group = mrk, color = mrk))
      } else {
        p <- ggplot() + 
        geom_line(data = subset(plot_data, rep == repLev[i]), 
                  aes(x = fraction, y = intensities, group = id, color = mrk),
                  alpha=0.05) 
      }
    }
  } else {
    if (ribbons) {
      p <- 
        ggplot() + geom_ribbon(data = df,
                               mapping = aes(fraction, ymin=min, ymax=max, group = mrk, 
                                             color = NA, fill = mrk), 
                               alpha=0.3) +
        geom_line(data = df,
                  mapping = aes(fraction, mean, group = mrk, color = mrk))
    } else {
      p <- ggplot() + 
        geom_line(data = plot_data, 
                  aes(x = fraction, y = intensities, group = id, color = mrk),
                  alpha=0.05) 
        
    }

  }
  
  ## extract colours for organelles in the data 
  col <- c(col, "unknown" = "darkgrey")
  if (is.factor(df$mrk)) 
    col <- col[levels(df$mrk)]
  else
    col <- col[unique(df$mrk)]

  lab_x <- fracLev
  lab_x[seq(2,length(lab_x),by=2)] <- ""
  ## plot data
  p <- p + 
    scale_x_discrete(limits=fracLev, 
                     breaks = fracLev, 
                     labels = lab_x) +  # show only every other label on x-axis
    ylab("Normalised intensities") + xlab("") +
    scale_fill_manual(values = col, aesthetics = c("fill","colour")) +
    scale_color_manual(values = col, aesthetics = c("fill, colour")) +
    theme_light() +
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none", 
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(size = 8, face="bold"),
          panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8 + 2, colour = rgb(0, 0, 0)),
          axis.text.y = element_text(size = 8, colour = rgb(0, 0, 0))) +
    facet_wrap(~ mrk, scales = "fixed", ...) 
  return(p)
}

## add limma output to the msnset
addLimmaInfo <- function(msnset, limmadf) {
  ind <- match(fData(msnset)$GN, rownames(limmadf))
  limmadf <- limmadf[ind, ]
  fData(msnset)$logFC <- limmadf$logFC
  fData(msnset)$AveExpr <- limmadf$AveExpr
  fData(msnset)$t <- limmadf$t
  fData(msnset)$P.Value <- limmadf$P.Value
  fData(msnset)$adj.P.Val <- limmadf$adj.P.Val
  fData(msnset)$B <- limmadf$B
  fData(msnset)$Significance <- limmadf$Significance
  # fData(msnset)$GN <- limmadf$GN
  return(msnset)
}

# getDescr <- function(x, fcol) {
#   info <- as.character(x[[fcol]])
#   desc <- as.character(sapply(info, function(z) 
#     strsplit(z, split = "OS=")[[1]][1]))
#   return(desc)
# }

## add plotting info to display significantly differently expressed proteins
addPlottingInfo <- function(x, p = 0.01, lfc = 0.6, lfc_strong) {
  if (missing(lfc_strong)) {
    ll <- c(paste0("Upregulated (FDR < ", p, ", log2FC > ", lfc, ")"),
            paste0("Downregulated (FDR < ", p, ", log2FC < -", lfc, ")"),
            "Not significant")
    x$Significance <- ifelse(x$adj.P.Val < p & x$logFC > lfc, ll[1], ll[3])
    getS <- x$adj.P.Val < p & x$logFC < -lfc
    if (length(getS) > 0)
      x$Significance[getS] <- ll[2]
    x$Significance <- factor(x$Significance, levels = ll)
  } else {
    ll <- c(paste0("Upregulated (FDR < ", p, ", log2FC > ", lfc, ")"),
            paste0("Upregulated (FDR < ", p, ", log2FC > ", lfc_strong, ")"),
            paste0("Downregulated (FDR < ", p, ", log2FC < -", lfc, ")"),
            paste0("Downregulated (FDR < ", p, ", log2FC < -", lfc_strong, ")"),
            "Not significant")
    x$Significance <- ifelse(x$adj.P.Val < p & x$logFC > lfc, ll[1], ll[5])
    getS <- x$adj.P.Val < p & x$logFC < -lfc
    if (length(getS) > 0)
      x$Significance[getS] <- ll[3]
    .lfc_strong1 <- x$adj.P.Val < p & x$logFC > lfc_strong
    .lfc_strong2 <- x$adj.P.Val < p & x$logFC < -lfc_strong
    if (any(.lfc_strong1)) x$Significance[.lfc_strong1] <- rep(ll[2], length(which(.lfc_strong1)))
    if (any(.lfc_strong2)) x$Significance[.lfc_strong2] <- rep(ll[4], length(which(.lfc_strong2)))
    x$Significance <- factor(x$Significance, levels = ll)
  }
  ## make a new column for Gene Names in table
  x$GN <- rownames(x)
  return(x)
}

# addInfo <- function(x, fcol, p) {
#   x$GN <- as.character(sapply(as.character(x[[fcol]]), getGN))
#   x <- addPlottingInfo(x, p)
#   return(x)
# }


## ====== volcano plot =====
## highlight - specify highlighting only points < value e.g. highlight = 0.001
## topN - specify only highlight the n proteins with the smallest adj p value
ggVolcano <- function(x, mytitle = "", p = 0.01, highlight, N, 
                      base.text.size = 16, 
                      text.size = 34, 
                      geom.point.size = 6,
                      legend.text.size = 27,
                      label.text.size = 6,
                      lfc = 0.6, 
                      legend_pos_x = .77,
                      legend_pos_y = .065,
                      addLegend = TRUE,
                      ...) {
  
  ## add plotting information
  x <- addPlottingInfo(x, p = p, lfc = lfc, ...)
  
  ## find the line for adjusted threshold to plot on y-axis 
  ind <- x$adj.P.Val <= p
  y_threshold <- min(-log10(x$P.Value[ind]))
  
  
  ## create the volcano plot
  gg <- ggplot(x, aes(x = logFC, y = -log10(P.Value))) +  
    geom_point(aes(fill = Significance), shape = 21, size = geom.point.size) 
  ## plot with two logFC indicated on the plot
  if (length(levels(x$Significance)) > 3) {
    gg <- gg + scale_fill_manual(values = setNames(c("#ffa3a3",
                                                     "#C50000", 
                                                     "#cae2f9", 
                                                     "#1B83E9",
                                                     "gray"),
                                                   levels(x$Significance)))
  } else {
    gg <- gg + scale_fill_manual(values = setNames(c("#C50000", "#1B83E9", "gray"),
                                                   levels(x$Significance)))
  }
  gg <- gg + theme_classic(base_size = base.text.size) + 
    theme(legend.position = c(legend_pos_x, legend_pos_y), 
          legend.text=element_text(size = legend.text.size),
          legend.title=element_blank(),
          # legend.title=element_text(size = legend.text.size, face = "bold", colour = "black"), 
          legend.key.height=unit(.5,"line"),
          axis.title=element_text(size = text.size, face="bold", colour = "black"),
          axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"), colour = "black", size = text.size),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"), colour = "black", size = text.size),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text = element_text(color = "black", size = text.size),
          panel.background = element_blank(),
          panel.border = element_blank()
    ) 
  y <- x %>% filter(adj.P.Val < p) %>% arrange(adj.P.Val)
  if (!missing("highlight")) {
    y <- y %>% filter(y, adj.P.Val < highlight)
    if (!misisng("N")) stop(message(paste("Can not use both highlight and topN arguments together")))
  } 
  if (!missing("N")) {
    if (length(levels(x$Significance)) > 3) {
      up <- y %>% filter(Significance == levels(x$Significance)[2]) 
      up <- up[seq(N), ]
      down <- y %>% filter(Significance == levels(x$Significance)[4])
      down <- down[seq(N), ]
      y <- rbind(up, down)
    } else {
      up <- y %>% filter(Significance == levels(x$Significance)[1]) 
      up <- up[seq(N), ]
      down <- y %>% filter(Significance == levels(x$Significance)[2])
      down <- down[seq(N), ]
      y <- rbind(up, down)
    }
    
  }
  
  ## add gene labels to the plot and
  ## add lfc and p-value cutoff lines
  if (length(levels(x$Significance)) > 3) {
    gg <- gg + geom_vline(xintercept = c(-lfc, lfc), linetype = 2) +
      geom_vline(xintercept = c(-1, 1), linetype = 2) +
      geom_hline(yintercept = y_threshold, linetype = 2) 
    gg <- gg +
      geom_label_repel(data = y,    
                       aes(label = GN, colour = Significance),    
                       size = label.text.size,   
                       box.padding = unit(.2, "lines"),  
                       point.padding = unit(.15, "lines"),
                       label.padding = unit(.15, "lines"),
                       show.legend = FALSE,
                       max.overlaps = N + 10) + 
      ylab(bquote(-log[10]~'('*p-value*')')) +
      xlab(bquote(log[2]~"FC")) +
      scale_colour_manual(values = c("#C50000", 
                                     "#1B83E9"))
  } else {
    gg <- gg + geom_vline(xintercept = c(-lfc, lfc), linetype = 2) +
      geom_hline(yintercept = y_threshold, linetype = 2) 
    gg <- gg +
      geom_label_repel(data = y,    
                       aes(label = GN, colour = Significance),    
                       size = label.text.size,   
                       box.padding = unit(.2, "lines"),  
                       point.padding = unit(.15, "lines"),
                       label.padding = unit(.15, "lines"),
                       show.legend = FALSE,
                       max.overlaps = N + 10) + 
      ylab(bquote(-log[10]~'('*p-value*')')) +
      xlab(bquote(log[2]~"FC")) +
      scale_colour_manual(values = c("#C50000", "#1B83E9"))
  }
  
  if (!addLegend) 
    gg <- gg + theme(legend.position = "none")
  return(gg)
}

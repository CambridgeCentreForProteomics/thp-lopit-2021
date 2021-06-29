GObarplots <- function(filename, cols, N=25) {
  go_data <- read.csv(filename, header = 1)
  termInfo <- strsplit(go_data$Term, split = "~")
  go_data$ID <- sapply(seq(length(termInfo)), function(x) termInfo[[x]][1])
  go_data$Term <- sapply(seq(length(termInfo)), function(x) termInfo[[x]][2])
  df <- go_data[, c("Cluster", "Term", "Count", "Benjamini")]
  df$adjPvalue <- -log10(df$Benjamini)
  df$name <- paste0("Cluster ", df$Cluster, ": ", df$Term)
  df$name <- factor(df$name, levels = df$name)
  positions <- (df$name)
  names(cols) <- levels(df$name)
  p <- ggplot(df, aes(x=name, y=adjPvalue, fill=name)) + 
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_light(base_size = 16) +
    scale_fill_manual(values = cols) +
    scale_x_discrete(limits = rev(positions)) +
    theme(axis.title.y=element_blank(), 
          legend.position = "none",
          axis.text.y = element_text(color="black", 
                                     size=16),
          axis.text.x = element_text(color="black", 
                                     size=16)) +
    ylab(bquote(-log[10]~'('*adj.*' '*p-value*')')) +
    # xlab("GO BP annotation") +
    geom_text(aes(label = Count), 
              size = 6, 
              color = "darkgrey",
              # position = position_dodge(width = 0.9),
              vjust = 0,
              nudge_y = 1) +
    ylim(c(0, N))
  return(p)
}
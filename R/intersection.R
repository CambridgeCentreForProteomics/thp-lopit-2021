## code for calculating overlap and generating heatmaps
## to display overlap (code written by L Breckels 19/09/18)
library("reshape2")
library("dplyr")
library("ggplot2")

## combine datasets and get data.frame of assignments
compareDatasets <- function(MSnSet, 
                            fcol1, 
                            fcol2,
                            organelle_order = c(getMarkerClasses(MSnSet, fcol1), "unknown"), 
                            rm.unknowns = FALSE) {
  
  ll_hl <- intersect(organelle_order, getMarkerClasses(MSnSet, fcol1))
  ll_dc <- intersect(organelle_order, getMarkerClasses(MSnSet, fcol2))
  if (any(ll_hl == "missing")) {
    .ind <- which(ll_hl != "missing")
    ll_hl <- c(ll_hl[.ind], "unknown", "missing")
    .ind <- which(ll_dc != "missing")
    ll_dc <- c(ll_dc[.ind], "unknown", "missing")
  } else {
    ll_hl <- c(ll_hl, "unknown")
    ll_dc <- c(ll_dc, "unknown")
  }
  ll_hl <- ll_dc <- union(ll_hl, ll_dc)
  loc_hl <- factor(fData(MSnSet)[, fcol1], ll_hl)
  loc_dc <- factor(fData(MSnSet)[, fcol2], ll_dc)
  tt <- table(loc_hl, loc_dc)
  if (rm.unknowns) {
    rmy <- which(colnames(tt) == "unknown")
    rmx <- which(rownames(tt) == "unknown")
    tt <- tt[-rmx, -rmy]
  } 
  rs <- rowSums(tt)
  cs <- colSums(tt)
  df <- melt(tt, value.name = "overlap")
  ref_all <- ref_hl <- ref_dc <- numeric(length = nrow(df))
  for (i in 1:nrow(df)) {
    ref_all[i] <- df[i, "overlap"]/((rs[df$loc_hl[i]] + 
                                       cs[df$loc_dc[i]]) - 
                                      df[i, "overlap"])
    ref_hl[i] <- df[i, "overlap"]/(rs[df$loc_hl[i]])
    ref_dc[i] <- df[i, "overlap"]/(cs[df$loc_dc[i]])
  }
  df <- cbind(df, ref_all, ref_hl, ref_dc)
  df$overlap <- as.character(df$overlap)
  df$overlap <- recode(df$overlap, "0" = "")
  return(df)
}



## plot heatmap of assignments
ggheatmap <- function(df, title = "", size = 24){
  ggplot(df, aes(loc_hl, loc_dc, fill=100*ref_dc)) + 
    geom_tile(colour="grey80") +
    geom_text(aes(label=overlap), size=3) +
    # theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    scale_fill_gradient(low="white", high="steelblue",
                        name="% Intersection\ntagm", 
                        limits=c(0,100)) +
    xlab("\ntagm unstimulated") +
    ylab("tagm LPS\n") +
    theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, 
                                                      b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = -18, 
                                                      b = 0, l = 0))) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    theme(text = element_text(size=size), 
          legend.text = element_text(size=size-1),
          axis.text.x  = element_text(colour = "black"),
          axis.text.y  = element_text(colour = "black")) +
    ggtitle(title)
}

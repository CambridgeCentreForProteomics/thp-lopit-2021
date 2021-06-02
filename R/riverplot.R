library("ggalluvial")

## Create an alluvial/river plot of translocations
riverplot <- function(msnset,
                      msnset2, 
                      fcol1 = "localisation.pred", 
                      fcol2 = "localisation.pred", 
                      mrkCol1 = "markers",
                      mrkCol2 = "markers",
                      onlyMovers = TRUE,
                      labels = TRUE,
                      cols,
                      ...) {
  
  ## add colour scheme if not provided
  if (missing(cols)) {
    ll <- unique(c(levels(df[,1]), levels(df[,2])))
    grid.col <- segcols <- setNames(rainbow(length(ll)), ll)
  } else {
    grid.col <- cols
  }
  
  ## make data.frame of translocations
  df <- makedf(msnset, msnset2, 
               fcol1, fcol2, 
               mrkCol1,
               mrkCol2,
               onlyMovers)
  orgs <- union(levels(df$x), levels(df$y))
  colscheme <- cols[orgs]
  torm <- which(df$count == 0)
  df <- df[-torm, ]
  
  ## set colours for alluvial plot (this is a little tricky as a specific 
  ## ordering is required for ggalluvial)
  names(df) <- c("Condition1", "Condition2", "value")
  y <- setdiff(unique(df[, "Condition1"]), unique(df[, "Condition2"]))
  updateCols <- colscheme[-match(y, names(colscheme))]
  columnCols <- c(colscheme, updateCols)
  stratCols <- c(rev(colscheme), rev(updateCols))
  
  ## plot alluvial/river schematic
  if (labels) {
    q <- ggplot(df,
                aes(y = value, axis1 = Condition1, axis2 = Condition2)) +
      geom_alluvium(aes(fill = Condition1), width = 0) +
      scale_fill_manual(values = columnCols) +
      geom_stratum(width = 1/8, fill = paste0(stratCols), color = "white") +
      geom_label(stat = "stratum", aes(label = after_stat(stratum)),
      color = stratCols, fontface = "bold", size = 3) +
      scale_x_discrete(limits = c("Condition1", "Condition2"), expand = c(.09, .09)) +
      scale_y_continuous(breaks = NULL) +
      theme_minimal() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank()) +
      theme(legend.position = "none") +
      ylab(NULL)
  } else {
    q <- ggplot(df,
                aes(y = value, axis1 = Condition1, axis2 = Condition2)) +
      geom_alluvium(aes(fill = Condition1), width = 0) +
      scale_fill_manual(values = columnCols) +
      geom_stratum(width = 1/8, fill = paste0(stratCols), color = "white") +
      scale_x_discrete(limits = c("Condition1", "Condition2"), expand = c(.09, .09)) +
      scale_y_continuous(breaks = NULL) +
      theme_minimal() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank()) +
      theme(legend.position = "none") +
      ylab(NULL)
  }
  q
}

## get count data (for organelle transfers between conditions) 
## from msnset
makedf <- function(msnset1, msnset2, 
                   fcol1, fcol2, 
                   mrkCol1 = "markers",
                   mrkCol2 = "markers",
                   # includeMarkers = FALSE, ## logical whether include markers
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

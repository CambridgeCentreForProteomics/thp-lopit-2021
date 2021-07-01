## extract gene name from protein description
addGeneName <- function(msnset, fcol, col.name = "GN") {
  .descr <- fData(msnset)[[fcol]]
  .gn <- sapply(.descr, getGN, USE.NAMES = FALSE)
  if (any(duplicated(.gn))) {
    ## add an the uniprot id as an index for duplicated Gene Names 
    ## we require unique gene names for plotting
    .ind <- which(duplicated(.gn))
    .gn[.ind] <- paste0(.gn[.ind], "_", featureNames(msnset)[.ind])
  }
  fData(msnset)[[col.name]] <- .gn
  return(msnset)
}

## helpers to extract gene name for description in msnset
getGN <- function(x) {
  pn2 <- strsplit(x, split = "GN=")
  ind <- grep(pattern = "PE=", pn2[[1]])
  if (length(ind) > 0)  
    gn <- strsplit(pn2[[1]][ind], split = " PE=")[[1]][1]
  else
    gn <- pn2[[1]][length(pn2[[1]])]
  return(gn)
}
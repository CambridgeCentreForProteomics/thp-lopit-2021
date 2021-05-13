## ============================================================================
## ============================================================================
## Caclulcate the natural L2 distance
## ============================================================================
## ============================================================================
calc_L2 <- function(msnset1, msnset2, fcol_joint) {
  ## find columns which contain the joint posterior probability for each
  ## subcellular niche 
  if (missing(fcol_joint)) {
    fcol_joint_1 <- grep("tagm.mcmc.joint", fvarLabels(msnset1))
    fcol_joint_2 <- grep("tagm.mcmc.joint", fvarLabels(msnset2))
  }
  ## get tagm probability point estimates
  postprob_1 <- fData(msnset1)[, fcol_joint_1]
  postprob_2 <- fData(msnset2)[, fcol_joint_2]
  prob_diff <- rowSums(abs(postprob_1 - postprob_2)^2/2)
  return(as.numeric(format(round(prob_diff, 4), nsmall = 3)))
}


## ============================================================================
## ============================================================================
## Extract translocations based on different classifications between datasets
## and the L2 distance between posterior probability distributions
## ============================================================================
## ============================================================================
getTranslocations <- function(msnset1, msnset2, 
                              fcol = "tagm.mcmc.allocation.pred",
                              out = 1e-3,
                              type = c(1, 2, 3, 4),
                              subset_l2 = TRUE,
                              t = .5) {
  
  stopifnot(nrow(msnset1) == nrow(msnset2))
  
  ## calculate the natural L2 distance between proteins based on their posterior
  ## probability distribution over all organelles
  L2_distance <- calc_L2(msnset1, msnset2)
  fData(msnset1)$L2_distance <- fData(msnset2)$L2_distance <- L2_distance
  
  
  ## extract all proteins that have a changed location i.e. classification 1 != classification 2
  ind_orgs <- which(fData(msnset1)[, fcol] != fData(msnset2)[, fcol]) 
  org_assign1 <- which(fData(msnset1)[, fcol] != "unknown") # those assigned to a discrete niche in unstimulated
  org_assign2 <- which(fData(msnset2)[, fcol] != "unknown") # those assigned to a discete niche in LPS
  un_assign1 <- which(fData(msnset1)[, fcol] == "unknown") # extract unknowns in unstimulated
  un_assign2 <- which(fData(msnset2)[, fcol] == "unknown") # extract unknowns in LPS
  fData(msnset1)$translocations <- fData(msnset2)$translocations <- NA
  
  
  ## ----type 1 - organelle to organelle 
  if (1 %in% type) {
    ## get proteins that have different assignments in msnset1 and msnset2
    org_assign_both <- intersect(org_assign1, org_assign2)
    org_assign_both <- intersect(org_assign_both, ind_orgs)
    ## keep only translocations that >= t
    if (subset_l2) {
      l2_check <- which(fData(msnset1)[org_assign_both, "L2_distance"] >= t)
      if (length(l2_check)>0)
        org_assign_both <- org_assign_both[l2_check]
    }
    fData(msnset1)$type1_translocation <- FALSE
    fData(msnset2)$type1_translocation <- FALSE
    fData(msnset1)$type1_translocation[org_assign_both] <- fData(msnset2)$type1_translocation[org_assign_both] <- TRUE
    fData(msnset1)$translocations[org_assign_both] <- fData(msnset2)$translocations[org_assign_both] <- "type1"
  }
  
  ## ----type 2 - unknown to organelle
  if (2 %in% type) {
    ## find unknowns with high outlier prob
    un_outlier1 <- which(fData(msnset1)[, "tagm.mcmc.outlier"] >= (1 - out))
    uns1 <- intersect(un_outlier1, un_assign1)
    un_to_org <- intersect(uns1, org_assign2)
    ## keep only translocations that >= t
    if (subset_l2) {
      l2_check <- which(fData(msnset1)[un_to_org, "L2_distance"] >= t)
      if (length(l2_check)>0)
        un_to_org <- un_to_org[l2_check]
    }
    fData(msnset1)$type2_translocation <- fData(msnset2)$type2_translocation <- FALSE
    fData(msnset1)$type2_translocation[un_to_org] <- fData(msnset2)$type2_translocation[un_to_org] <- TRUE
    fData(msnset1)$translocations[un_to_org] <- fData(msnset2)$translocations[un_to_org] <- "type2"
  }
  
  ## ----type 3 - organelle to unknown
  if (3 %in% type) {
    ## find unknowns with high outlier prob
    un_outlier2 <- which(fData(msnset2)[, "tagm.mcmc.outlier"] >= (1 - out))
    uns2 <- intersect(un_outlier2, un_assign2)
    un_to_org2 <- intersect(uns2, org_assign1)
    ## keep only translocations that >= t
    if (subset_l2) {
      l2_check <- which(fData(msnset1)[un_to_org2, "L2_distance"] >= t)
      if (length(l2_check)>0)
        un_to_org2 <- un_to_org2[l2_check]
    }
    fData(msnset1)$type3_translocation <- fData(msnset2)$type3_translocation <- FALSE
    fData(msnset1)$type3_translocation[un_to_org2] <- fData(msnset2)$type3_translocation[un_to_org2] <- TRUE
    fData(msnset1)$translocations[un_to_org2] <- fData(msnset2)$translocations[un_to_org2] <- "type3"
  }

  ##  ---- type 4
  if (4 %in% type) {
    ind_orgs <- which(fData(msnset1)[, fcol] == fData(msnset2)[, fcol]) # same location in both experiments
    un <- intersect(ind_orgs, which(fData(msnset1)[, fcol] == "unknown"))
    ind <- which(fData(msnset1)[, "L2_distance"] == 1)
    cmn <- intersect(ind, un)
    fData(msnset1)$type4_translocation <- fData(msnset2)$type4_translocation <- FALSE  
    fData(msnset1)[cmn, "type4_translocation"] <- fData(msnset2)[cmn, "type4_translocation"] <- TRUE
    fData(msnset1)$translocations[cmn] <- fData(msnset2)$translocations[cmn] <- "type4"
  }

  
  ## return data with translocations appended to the fData 
  return(MSnSetList(c(msnset1, msnset2)))
}
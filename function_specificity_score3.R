## sgRNA specificity scores ####
specificity3 <- function(penalties = c("qi.mean.per.region", "qi", "custom"), 
                         custom.penalties = NULL, 
                         sgRNAs, 
                         return.offTargets = TRUE, 
                         ...){
  ####
  # Function
  ##
  # Quantifies the specificity of a sgRNA to its target by penalizing any potential off-target hits.
  #   Specificity is maximal (1), minus weighted penalty products.
  #   When multiple off-target sequences are found, lowest specificity is returned.
  #   An exact off-target match gives a specificity score of 0.
  #### 
  # Arguments
  ##
  # penalties         Mismatch position-specific penalties to be used. 
  #                   Defaults to means per region as defined in Qi et al., 
  #                   these values were estimated from Figure 5D. Ref:
  #                   Qi et al. (2013) Cell, http://dx.doi.org/10.1016/j.cell.2013.02.022. 
  #                   Other options are those estimates per base or custom user-defined penalties.
  # custom.penalties  If penalties argument is set to "custom", this should be a 
  #                   user-defined numeric vector of length 20.
  # sgRNAs            A data.frame object of which the first column contains the 
  #                   sgRNA tags and the second column contains the sgRNA sequences.
  # return.offTargets Whether to also return the full data frame with off-targets
  #                   and scores as identified and assigned by CRISPRseek::searchHits.
  # ...               Other arguments to be passed to CRISPRseek::searchHits.
  ####
  require(CRISPRseek)
  qi <- c(0.78, 0.41, 0.44, 0.63, 0.65, 0.68, 0.65, 0.63, # III (inverse counting Qi vs CRISPRseek)
          0.30, 0.25, 0.24, 0.22, 0.24,                   # II
          0.01, 0.12, 0.10, 0.07, 0.06, 0.09, 0.05)       # I
  qi.mean.per.region <- rep(tapply(qi, rep(c("III", "II", "I"), c(8, 5, 7)), mean), c(8, 5, 7))
  penalties <- switch(penalties[1], 
                      qi = qi, 
                      qi.mean.per.region = qi.mean.per.region, 
                      custom = custom.penalties)
  if(length(penalties) != 20 | !is.numeric(penalties)){
    stop("For custom penalties, please set custom.penalties to a numeric vector of length 20")
  }
  gRNAs <- DNAStringSet(unlist(sgRNAs[, 2]), use.names = FALSE)
  names(gRNAs) <- unlist(sgRNAs[, 1])
  SH <- searchHits(gRNAs = gRNAs, 
                   ...)
  res <- list("specificity" = tapply(1 - apply(SH[, 1:20], 1, function(x){
    ifelse(any(x != 0), prod(penalties[x == 1]), 0)
  }), SH$name, function(y){
    ifelse(sum(y == 1) > 1, 0, min(y))
  }))
  if(return.offTargets){
    res <- c(res, "offTargetHits" = list(SH))
  }
  return(res)
}
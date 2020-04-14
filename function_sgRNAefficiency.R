sgRNAefficiency <- function(sgRNAs, genes, 
                            reprAct = TRUE, dist2SC = TRUE, 
                            name_by = "locus_tag", 
                            penalties = c("qi.mean.per.region", "qi", "custom"), 
                            custom.penalties = NULL, 
                            ...){
  require(CRISPRseek)
  
  # Retrieve all gene names from GFF input
  tags <- unlist(lapply(genes$attribute, function(x){
    # add ";" at end in case locus_tag is last attribute (regular expression requires ending character)
    sub(paste0(".*?", name_by, "=(.*?);.*"), "\\1", paste0(x, ";"))
  }))
  
  # Checks
  if(any(duplicated(tags))){
    stop(paste("Duplicates of", name_by, "in genes."))
  }
  if(reprAct){ # set penalties if repression activity should be computed
    qi <- c(0.78, 0.41, 0.44, 0.63, 0.65, 0.68, 0.65, 0.63, # III (inverse counting Qi vs CRISPRseek)
            0.30, 0.25, 0.24, 0.22, 0.24,                   # II
            0.01, 0.12, 0.10, 0.07, 0.06, 0.09, 0.05)       # I
    # tapply sorts by names, so sort by 1:3 (CRISPRseek), while I:III is Qi et al regions
    qi.mean.per.region <- rep(tapply(qi, rep(c("1III", "2II", "3I"), c(8, 5, 7)), mean), c(8, 5, 7))
    penalties <- switch(penalties[1], 
                        qi = qi, 
                        qi.mean.per.region = qi.mean.per.region, 
                        custom = custom.penalties)
    if(length(penalties) != 20 | !is.numeric(penalties)){
      stop("For custom penalties, please set custom.penalties to a numeric vector of length 20")
    }
  }
  
  # Identify all binding sites
  message("Starting binding site identification...")
  lib_hits <- searchHits(gRNAs = sgRNAs, 
                         ...)
  
  # Get NT CDS per hit
  # identify per site in which gene it lies on NT strand, if any
  #    if gene on + strand, sgRNA seq should be on - strand
  #    then complementary to coding = NT strand
  message("Identifying gene hits per sgRNA...")
  NT_gene <- apply(lib_hits, 1, function(x){
    # find overlapping genes per site, if any
    tmp_i <- as.numeric(x["chromStart"]) - genes$end <= 0 & 
      as.numeric(x["chromEnd"]) - genes$start >= 0 & 
      x["strand"] != genes$strand
    # add details of gene hits, if any
    if(any(tmp_i)){
      # could be multiple genes that are overlapped by a site
      do.call(rbind, lapply(which(tmp_i), function(y){
        if(as.numeric(x["chromStart"]) < genes$start[y]){
          coverPart <- switch(genes$strand[y], 
                              `+` = "5-tail", 
                              `-` = "3-PAM")
          coverSize <- as.numeric(x["chromEnd"]) - genes$start[y] + 1
        } else{
          if(as.numeric(x["chromEnd"]) > genes$end[y]){
            coverPart <- switch(genes$strand[y], 
                                `+` = "3-PAM", 
                                `-` = "5-tail")
            coverSize <- genes$end[y] - as.numeric(x["chromStart"]) + 1
          } else{
            coverPart <- "complete"
            coverSize <- 23
          }
        }
        # return details if gene hits for site
        matrix(c(tags[y], coverPart, coverSize), ncol = 3)
      }))
    } else{
      # return NAs if no gene hits for site
      matrix(NA, ncol = 3, nrow = 1)
    }
  })
  # repeat row indices with x number of gene hits x times
  NT_gene_i <- unlist(mapply(rep, 
                             seq.int(nrow(lib_hits)), 
                             unlist(lapply(NT_gene, nrow))))
  # add details to data frame
  lib_hits <- cbind(lib_hits[NT_gene_i, ], 
                    "site" = NT_gene_i, 
                    do.call(rbind, NT_gene))
  colnames(lib_hits)[32:34] <- c("NTgene", "coverPart", "coverSize")
  # coverSize should be integer, not factor
  lib_hits$coverSize <- as.integer(levels(lib_hits$coverSize)[lib_hits$coverSize])
  
  # Get repression activity
  if(reprAct){
    message("Calculating repression activity per site...")
    # Get retained repression per hit
    repr_act <- apply(lib_hits[, 1:20], 1, function(x){
      prod(penalties[x == 1]) # product of empty set is 1 by definition (OK for 0 mm)
    })
    # add to data frame
    lib_hits$reprAct <- repr_act
  }
  
  # Get distance to start codon if sgRNA can (partially) bind within gene
  if(dist2SC){
    message("Calculating relative distances of sgRNA binding within gene hits...")
    # Get distance to start codon per hit
    #    normalize within-gene distances with feature scaling
    feat_scale <- function(x, min, max){(x - min) / (max - min)}
    dist_SC <- apply(lib_hits, 1, function(x){
      if(!is.na(x["NTgene"])){
        tmp_range <- genes[match(x["NTgene"], tags), c("start", "end")]
        # max is end CDS - 22 nt for binding first nt of PAM+spacer (chromStart, strand does not matter)
        #    also partial overlap (coverPart != "complete") dist_SC should always be set to [0,1]
        ifelse(genes$strand[match(x["NTgene"], tags)] == "+", 
               # strand matters: on "-" strand, start of gene is on 3-prime end so invert distance
               pmax(0, pmin(1, 
                            feat_scale(as.numeric(x["chromStart"]), 
                                       tmp_range[1], 
                                       tmp_range[2] - 22))), 
               1 - pmax(0, pmin(1, 
                                feat_scale(as.numeric(x["chromStart"]), 
                                           tmp_range[1], 
                                           tmp_range[2] - 22))))
      } else{
        NA
      }
    })
    # add to data frame
    lib_hits$dist2SC <- unlist(dist_SC)
  }
  
  # Return full data frame
  return(lib_hits)
}
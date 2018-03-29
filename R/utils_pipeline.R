
# CIDR --------------------------------------------------------------------

#' @importFrom cidr scDataConstructor determineDropoutCandidates wThreshold scDissim cidrPcoa nPC scCluster
.cidrfree.calc <- function(countData, NCores) {

  ncores = ifelse(is.null(NCores), 1, NCores)

  # construct input object
  sData <- cidr::scDataConstructor(tags=countData)

  # determine gene dropout candidates
  sData <- cidr::determineDropoutCandidates(object = sData,
                                            min1 = 3, min2 = 8,
                                            N = 2000, alpha = 0.1,
                                            fast = TRUE, zerosOnly = FALSE,
                                            bw_adjust = 1)


  # impute weighting threshold
  sData <- cidr::wThreshold(object = sData, cutoff = 0.5, plotTornado = FALSE)

  # dissimilarity matrix
  sData <- cidr::scDissim(object = sData,
                    correction = FALSE,
                    threads = ncores,
                    useStepFunction = TRUE)

  # PCA on dissimilarity matrix ( do not use scPCA since it has a built-in plotting function)
  y <- cidr::cidrPcoa(sData@dissim)
  variation <- y$values
  ## store all eigenvalues - neg, 0, & pos
  sData@eigenvalues <- variation
  ## for variation, only deal with positive eigenvalues
  variation <- variation[variation>0]
  sData@PC <- y$vectors[, 1:length(variation)]
  variation <- variation/sum(variation)
  sData@variation <- variation

  # determine optimal number of clusters
  sData <- cidr::nPC(object = sData)

  # clustering on PCs.
  sData <- cidr::scCluster(object = sData,
                           n = NULL,
                           nCluster = NULL,
                           nPC = NULL,
                           cMethod = "ward.D2")

  # return object
  return(sData@clusters)
}

#' @importFrom cidr scDataConstructor determineDropoutCandidates wThreshold scDissim cidrPcoa nPC scCluster
.cidrbound.calc <- function(countData, clustNumber, NCores) {

  ncores = ifelse(is.null(NCores), 1, NCores)

  # construct input object
  sData <- cidr::scDataConstructor(tags=countData)

  # determine gene dropout candidates
  sData <- cidr::determineDropoutCandidates(object = sData,
                                            min1 = 3, min2 = 8,
                                            N = 2000, alpha = 0.1,
                                            fast = TRUE, zerosOnly = FALSE,
                                            bw_adjust = 1)


  # impute weighting threshold
  sData <- cidr::wThreshold(object = sData, cutoff = 0.5, plotTornado = FALSE)

  # dissimilarity matrix
  sData <- cidr::scDissim(object = sData,
                          correction = FALSE,
                          threads = ncores,
                          useStepFunction = TRUE)

  # PCA on dissimilarity matrix ( do not use scPCA since it has a built-in plotting function)
  y <- cidr::cidrPcoa(sData@dissim)
  variation <- y$values
  ## store all eigenvalues - neg, 0, & pos
  sData@eigenvalues <- variation
  ## for variation, only deal with positive eigenvalues
  variation <- variation[variation>0]
  sData@PC <- y$vectors[, 1:length(variation)]
  variation <- variation/sum(variation)
  sData@variation <- variation

  # determine optimal number of clusters
  sData <- cidr::nPC(object = sData)

  # clustering on PCs.
  sData <- cidr::scCluster(object = sData,
                           n = NULL,
                           nCluster = clustNumber,
                           nPC = NULL,
                           cMethod = "ward.D2")

  # return object
  return(sData@clusters)
}

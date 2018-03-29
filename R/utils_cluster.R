
# WRAPPER CLUSTERING ------------------------------------------------------


.cluster.calc <- function(DimData,
                          ClustMethod,
                          clustNumber,
                          verbose) {
  if(ClustMethod=="kmeans") {ClusterData <- .kmeans.calc(DimData=DimData,
                                                         clustNumber=clustNumber,
                                                         verbose=verbose)}
  if(ClustMethod=="PAM") {ClusterData <- .pam.calc(DimData=DimData,
                                                         clustNumber=clustNumber,
                                                         verbose=verbose)}
  return(ClusterData)
}


# K-MEANS ----------------------------------------------------------------

#' @importFrom stats kmeans
.kmeans.calc <- function(DimData,
                         clustNumber,
                         verbose){
  k = clustNumber
  clust.res <- stats::kmeans(DimData, k)
  assignedClust <- clust.res$cluster
  return(assignedClust)
}

# PAM ---------------------------------------------------------------------

#' @importFrom cluster pam
.pam.calc <- function(DimData,
                      clustNumber,
                      verbose){
  k = clustNumber
  assignedClust <- cluster::pam(DimData, k, cluster.only = TRUE)
  return(assignedClust)
}




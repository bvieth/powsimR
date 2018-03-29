
# WRAPPER FOR CELL GROUP CLASSIFICATION -----------------------------------

## NOTE: minclust has adjusted rand index implemented

# TODO: put a counter of how many HVG were actually found (absolute and relative to total # genes!)

# NULL always means skip this step

.classify.calc <- function(Pipeline,
                           GeneSelect,
                           DimReduce,
                           ClustMethod,
                           clustNumber,
                           normData,
                           Lengths,
                           MeanFragLengths,
                           countData,
                           spikeData,
                           spikeInfo,
                           spikeIns,
                           NCores,
                           verbose) {


  if(is.null(Pipeline)) {
    # 1. expression matrix calculation
    exprData <- .expr.calc(countData=countData,
                           normData=normData,
                           Lengths=Lengths,
                           MeanFragLengths=MeanFragLengths,
                           group=NULL)

    exprData <- log2(exprData+1)

    # 2. feature selection
    FeatureSelect <- .feature.select(GeneSelect=GeneSelect,
                                     countData=countData,
                                     normData=normData,
                                     exprData=exprData,
                                     spikeData=spikeData,
                                     spikeInfo=spikeInfo,
                                     spikeIns=spikeIns,
                                     verbose=verbose)

    # 3. dimensionality reduction
    DimData <- .dim.reduce(DimReduce,
                           FeatureSelect)

    # 4. clustering
    ClusterData <- .cluster.calc(DimData,
                                 ClustMethod,
                                 clustNumber,
                                 verbose)
  }

  if(!is.null(Pipeline) && Pipeline=="CIDR_free") {
    ClusterData <- .cidrfree.calc(countData, NCores)
  }
  if(!is.null(Pipeline) && Pipeline=="CIDR_bound") {
    ClusterData <- .cidrbound.calc(countData, clustNumber, NCores)
  }

  return(ClusterData)
}

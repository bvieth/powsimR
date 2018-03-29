# methods: HVG: trend fit, CV2 (see scran; can be based on spike-ins or endogeneous); can be extended with correlated pairs (); BASiCS with highly variable gene estimation; Gini Index (see GiniClust, but rather cumbersome implementation);

# FEATURE SELECTION WRAPPER -----------------------------------------------

#' @importFrom scran DM
.feature.select <- function(GeneSelect,
                            countData,
                            normData,
                            exprData,
                            spikeData,
                            spikeInfo,
                            spikeIns,
                            verbose) {
  if(GeneSelect=="none") {FeatureExprData <- exprData}
  if(GeneSelect=="HVG") {FeatureExprData <- .HVG.calc(spikeIns=spikeIns,
                                                      spikeData=spikeData,
                                                      countData=countData,
                                                      normData=normData,
                                                      exprData=exprData,
                                                      verbose=verbose)}
  if(GeneSelect=="Gini") {FeatureExprData <- .Gini.calc(spikeIns=spikeIns,
                                                      spikeData=spikeData,
                                                      countData=countData,
                                                      normData=normData,
                                                      exprData=exprData,
                                                      verbose=verbose)}
  # remove zero variance genes from object! otherwise PCA will not work
  # Computing the DM.
  means <- rowMeans(FeatureExprData)
  vars <- apply(FeatureExprData, 1, var)
  cv2 <- vars/means^2
  dm.stat <- scran::DM(means, cv2)
  qt <- stats::quantile(dm.stat, probs = c(0.2, 0.8), na.rm=TRUE)
  rows <- sapply(dm.stat, function(x) any( x < qt[1] | x > qt[2] ) )
  rows[is.na(rows)] <- FALSE
  FeatureExprData.filt <- FeatureExprData[rows,]
  if(verbose) {
    message(paste0(nrow(FeatureExprData.filt),
                   " genes left after DM filtering. Started with ", nrow(FeatureExprData)))
  }

  return(FeatureExprData.filt)
}


# HVG CALCULATION ---------------------------------------------------------
#' @importFrom scater normalize calculateQCMetrics
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike
#' @importFrom BiocGenerics sizeFactors
#' @importFrom scran trendVar decomposeVar computeSpikeFactors
.HVG.calc <- function(spikeIns, spikeData, countData, normData, exprData, verbose) {
  if(!isTRUE(spikeIns)) {
    # create sce set for normalisation
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countData)))
    BiocGenerics::sizeFactors(sce) <- normData$size.factors
    sce <- scater::normalize(sce)
    # fit mean-dependent trend
    var.fit <- scran::trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
    # decompose variation and select GeneSelect based on test result
    var.out <- scran::decomposeVar(sce, var.fit)
    hvg.out <- var.out[which(var.out$FDR <= 0.1),]
    # stop and continue with all genes if not enough HVGs!
    if(nrow(hvg.out)<=10) {
      if (isTRUE(verbose)) {
        message(paste0("Less than 10 endogenous genes were defined as highly variable.
                       Continue with all genes instead."))
      }
      hvg.genes <- row.names(var.out)
    }
    if (nrow(hvg.out)>10) {
      hvg.genes <- row.names(hvg.out)
    }
  }
  if(isTRUE(spikeIns) && !is.null(spikeData)) {
    # create sce set for normalisation
    combined <- rbind(countData, spikeData)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(combined)))
    sce <- scater::calculateQCMetrics(sce,
                                      feature_controls = list(Spike = nrow(countData):nrow(combined)))
    SingleCellExperiment::isSpike(sce, "Spike") <- nrow(countData):nrow(combined)
    BiocGenerics::sizeFactors(sce) <- normData$size.factors
    sce <- scran::computeSpikeFactors(sce, general.use=FALSE)
    sce <- scater::normalize(sce)
    # fit mean-dependent trend
    var.fit <- scran::trendVar(sce, method="loess", use.spikes=TRUE, span=0.2)
    # decompose variation and select GeneSelect based on test result
    var.out <- scran::decomposeVar(sce, var.fit)
    #remove spikes from object otherwise they might end up in exprData!
    var.out <- var.out[!rownames(var.out) %in% rownames(spikeData),]
    # filter HVGs by adj p-value
    hvg.out <- var.out[which(var.out$FDR <= 0.1),]
    # stop and continue with all genes if not enough HVGs!
    if(nrow(hvg.out)<=10) {
      if (isTRUE(verbose)) {
        message(paste0("Less than 10 endogenous genes were defined as highly variable. Continue with all genes instead."))
      }
      hvg.genes <- row.names(exprData)
    }
    if (nrow(hvg.out)>10) {
      hvg.genes <- row.names(hvg.out)
    }
  }

  #TODO: use highly variable gene selection by BASiCS if BASiCS ran for normalisation
  if(isTRUE(spikeIns) && attr(normData, 'normFramework') == 'BASiCS') {
    feature.expr.dat <- exprData
  }

  # subset expr data based on gene selection
  feature.expr.dat <- exprData[hvg.genes,]
  return(feature.expr.dat)
}

# GINI CALCULATION --------------------------------------------------------

#TODO: Implement Gini Index calculation
.Gini.calc <- function(spikeIns, spikeData, countData, normData, exprData, verbose) {
  return(exprData)
}

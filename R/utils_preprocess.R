
# IMPUTATION WRAPPER ------------------------------------------------------

.impute.calc <- function(Imputation,
                         countData,
                         spikeData,
                         batchData,
                         clustNumber,
                         Lengths,
                         MeanFragLengths,
                         NCores,
                         verbose) {

  if(Imputation=='SAVER') {FilterData <- .saver.impute(countData = countData,
                                                   NCores = NCores,
                                                   verbose = verbose)}

  return(FilterData)
}

# imputation --------------------------------------------------------------

# SAVER
# this is one that apparently works but is really slow and needs cores!
# as an example: would need 30 minutes for 3500 genes in 500 cells on 10 cores
# so i changed the number of cells and genes used in predictions
#' @importFrom SAVER saver combine.saver
#' @importFrom utils capture.output
#' @importFrom doParallel registerDoParallel stopImplicitCluster
.saver.impute <- function(countData, NCores, verbose) {

  if(!is.null(NCores)) {
    doParallel::registerDoParallel(cores = NCores)
  }

  npred <- nrow(countData)/3 # number of genes for regression prediction

  # find genes with at least 50 % dropout and only do imputation for those
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highDrop <- p0 > 0.5

  d <- rownames(countData)[highDrop]
  max.g <- 2500
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/max.g))
  saver.outs <- sapply(1:length(d1), function(s) {
    genes <- d1[[s]]
    genes.ind <- which(rownames(countData) %in% genes)
    invisible(utils::capture.output(
      tmp <- suppressMessages(
        SAVER::saver(x = countData,
                     do.fast = TRUE,
                     size.factor = NULL,
                     npred = npred,
                     pred.cells = NULL,
                     pred.genes = genes.ind,
                     pred.genes.only = TRUE,
                     null.model = FALSE)
    )))
    tmp
  }, simplify = F)

  if(length(saver.outs)==1) {
    saver.out <- saver.outs[[1]]
  }

  if(length(saver.outs)>1) {
    saver.out <- SAVER::combine.saver(saver.outs)
  }

  if(!is.null(NCores)) {
    doParallel::stopImplicitCluster()
  }

  # internally he uses normalized data to calculate the fits. these are used to predict expression.
  # the scaled size factors are used for normalisation of raw data in prediction, but this is equal to 1 if left to default.
  # therefore the predicted estimates are like raw imputed counts albeit still rounding needed.

  ImputeData = saver.out$estimate
  # round to integer counts
  ImputeData = as.matrix(round(ImputeData))

  # combine ImputeData with countData since imputation
  # was only done for the genes with higher dropouts
  ImputeData <- ImputeData[order(match(rownames(ImputeData), rownames(countData))),]
  CombineData <- countData
  CombineData[rownames(countData) %in% rownames(ImputeData),] <- ImputeData

  rownames(CombineData) = rownames(countData)
  colnames(CombineData) = colnames(countData)

  invisible(gc())

  # return object
  return(CombineData)
}

# PREFILTER WRAPPER ---------------------------------------------------

.prefilter.calc <- function(Prefilter,
                            countData,
                            NCores) {
  if(Prefilter=='CountFilter') {FilterData <- .count.filter(countData = countData)}
    if(Prefilter=='FreqFilter') {FilterData <- .freq.filter(countData = countData)}

  return(FilterData)
}

# filtering ---------------------------------------------------------------

# apply stochastic filtering by Lun et al 2016 / Finak et al 2016
.count.filter <- function(countData) {
  highE<- rowMeans(countData) >= 0.2
  FilterData <- countData[highE,]
  return(FilterData)
}

# apply frequency filter: gene must have less than 80% dropouts
.freq.filter <- function(countData) {
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highE <- p0 < 0.8
  FilterData <- countData[highE,]
  return(FilterData)
}

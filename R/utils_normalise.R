
# Normalisation Wrapper ---------------------------------------------------

.norm.calc <- function(Normalisation,
                       sf,
                       countData,
                       spikeData,
                       spikeInfo,
                       batchData,
                       Lengths,
                       MeanFragLengths,
                       PreclustNumber,
                       Label,
                       Step,
                       Protocol,
                       NCores,
                       verbose) {
  if(Normalisation=='TMM') {
    NormData <- .TMM.calc(countData = countData,
                          verbose = verbose)
  }
  if(Normalisation=='MR') {
    NormData <- .MR.calc(countData = countData,
                         spikeData = spikeData,
                         verbose = verbose)
  }
  if(Normalisation=='scran' && Label == "none") {
    NormData <- .scran.calc(countData = countData,
                            spikeData = spikeData,
                            verbose = verbose)
  }
  if(Normalisation=='scran' && Label == "known") {
    if(Step=="Estimation"){
    NormData <- .scran.calc(countData = countData,
                            spikeData = spikeData,
                            verbose = verbose)}
    if(Step=="Simulation"){
      NormData <- .scrangroups.calc(countData = countData,
                                    batchData = batchData,
                                    verbose = verbose)}
  }
  if(Normalisation=='scran' && Label == "clustering") {
    NormData <- suppressWarnings(
    .scranclust.calc(countData = countData,
                     PreclustNumber = PreclustNumber,
                     verbose = verbose))
  }

  if(Normalisation=='depth') {
    NormData <- .depth.calc(countData = countData,
                            verbose = verbose)
  }
  if(Normalisation=='SF') {
    NormData <- .sf.calc(countData = countData,
                         sf = sf,
                         verbose = verbose)
  }
  if(Normalisation=='none') {
    NormData <- .none.calc(countData = countData,
                           verbose = verbose)
  }


  # inform the user of unusual / problematic normalisation results

  if(attr(NormData, "normFramework")=="scran" && any(NormData$size.factors<0)){
    if (verbose) { message(paste0("Negative size factors estimated.
                                  Apply stronger gene and sample filtering for parameter estimation.")) }
  }

  return(NormData)
}

# TMM --------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.TMM.calc <- function(countData, verbose) {
  norm.factors <- edgeR::calcNormFactors(object=countData, method='TMM')
  sf <- norm.factors * colSums(countData)
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "TMM"
  return(res)
}


# MR ---------------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
.MR.calc <- function(countData, spikeData, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)

  if(spike==TRUE) {
    message(paste0("Using controlGenes, i.e. spike-ins for normalisation!"))
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    controlgenes = grepl(pattern='ERCC', rownames(cnts))
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cnts,
                                                           colData=data.frame(group=rep('A', ncol(countData))),
                                                           design=~1))
    dds <- DESeq2::estimateSizeFactors(dds, type='ratio', controlGenes = controlgenes)
  }
  if(spike==FALSE) {
    message(paste0("Using all genes for normalisation!"))
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                                           colData=data.frame(group=rep('A', ncol(countData))),
                                                           design=~1))
    dds <- DESeq2::estimateSizeFactors(dds, type='ratio')
  }

  sf <- DESeq2::sizeFactors(dds)
  names(sf) <- colnames(countData)
  norm.counts <- DESeq2::counts(dds, normalized=TRUE)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "MR"
  return(res)
}

# scran ------------------------------------------------------

#' @importFrom scuttle computeSpikeFactors computePooledFactors
#' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors altExp
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
.scran.calc <- function(countData, spikeData, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)

  if(spike==TRUE) {
    message(paste0("Using computeSpikeFactors, i.e. spike-ins for normalisation!"))
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = countData))
    sce.spike <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = spikeData))
    SummarizedExperiment::rowData(sce.spike)$IDs <- rownames(spikeData)
    SingleCellExperiment::altExp(sce, "Spikes") <- sce.spike
    sce <- scuttle::computeSpikeFactors(x = sce, spikes = "Spikes", assay.type = "counts")
  }

  if(spike==FALSE) {
    message(paste0("Using computePooledFactors, i.e. deconvolution over all cells!"))
    cnts = countData
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
    if(ncol(countData)<=14) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1)))
      sf <- scuttle::computePooledFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(ncol(countData)>14 & ncol(countData)<=50) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6)))
      sce <- scuttle::computePooledFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(ncol(countData)>50 & ncol(countData)<=1000) {
      sizes <- c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6)))
      sce <- scuttle::computePooledFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(trunc(ncol(countData))>1000) {
      sizes <- c(round(seq(from=20, to=trunc(ncol(countData)/2), length.out=6)))
      sce <- scuttle::computePooledFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
  }

  sf <- SingleCellExperiment::sizeFactors(sce)
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}

#' @importFrom scuttle computeSpikeFactors computePooledFactors
#' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors altExp
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
.scrangroups.calc <- function(countData, batchData, verbose) {
  if (verbose) { message(paste0("Deconvolution within groups.")) }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countData)))

  if(is.vector(batchData)) {
    clusters <- as.factor(batchData)
  }
  if(!is.vector(batchData)) {
    clusters <- as.factor(batchData[,1])
  }

  if(ncol(countData)<=14) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(trunc(ncol(countData))>1000) {
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }

  sf <- SingleCellExperiment::sizeFactors(sce)
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}

#' @importFrom scran quickCluster
#' @importFrom scuttle computeSpikeFactors computePooledFactors
#' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors altExp
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
.scranclust.calc <- function(countData, PreclustNumber, verbose) {
  if (verbose) { message(paste0("Deconvolution within clusters.")) }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countData)))

  if(ncol(countData)<=14) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>1000 & ncol(countData)<=5000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sce <- scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(sce, method="igraph", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sce <-scuttle::computePooledFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }

  sf <- SingleCellExperiment::sizeFactors(sce)
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}

# Depth normalisation -----------------------------------------

.depth.calc <- function(countData, verbose) {
  sf <- colSums(countData) / mean(colSums(countData))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "depth"
  return(res)
}


# SF normalisation --------------------------------------------------------

.sf.calc <- function(countData, sf, verbose) {
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "SF"
  return(res)
}


# No normalisation --------------------------------------------------------

.none.calc <- function(countData, verbose) {
  sf <- rep(1, ncol(countData))
  names(sf) <- colnames(countData)
  norm.counts <- countData
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "none"
  return(res)
}



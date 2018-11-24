
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
                       NCores,
                       verbose) {
  if(Normalisation=='TMM') {NormData <- .TMM.calc(countData = countData,
                                                  verbose = verbose)}
  if(Normalisation=='UQ') {NormData <- .UQ.calc(countData = countData,
                                                verbose = verbose)}
  if(Normalisation=='MR') {NormData <- .MR.calc(countData = countData,
                                                spikeData = spikeData,
                                                verbose = verbose)}
  if(Normalisation=='PosCounts') {NormData <- .PosCounts.calc(countData = countData,
                                                              spikeData = spikeData,
                                                              verbose = verbose)}
  if(Normalisation=='Linnorm') {NormData <- .linnormnorm.calc(countData = countData,
                                                              spikeData = spikeData,
                                                              verbose = verbose)}
  if(Normalisation=='scran' && Label == "none") {NormData <- .scran.calc(countData = countData,
                                                                         spikeData = spikeData,
                                                                         verbose = verbose)}
  if(Normalisation=='scran' && Label == "known") {NormData <- .scrangroups.calc(countData = countData,
                                                                                batchData = batchData,
                                                                                verbose = verbose)}
  if(Normalisation=='scran' && Label == "clustering") {NormData <- suppressWarnings(
    .scranclust.calc(countData = countData,
                     PreclustNumber = PreclustNumber,
                     verbose = verbose))}
  if(Normalisation=='SCnorm' && Label == "none") {NormData <- .SCnorm.calc(countData = countData,
                                                                           spikeData = spikeData,
                                                                           batchData = NULL,
                                                                           NCores = NCores,
                                                                           verbose = verbose)}
  if(Normalisation=='SCnorm' && Label == "known") {NormData <- .SCnorm.calc(countData = countData,
                                                                           spikeData = spikeData,
                                                                           batchData = batchData,
                                                                           NCores = NCores,
                                                                           verbose = verbose)}
  if(Normalisation=='SCnorm' && Label == "clustering") {NormData <- .SCnormclust.calc(countData = countData,
                                                                                      spikeData = spikeData,
                                                                                      PreclustNumber=PreclustNumber,
                                                                                      NCores = NCores,
                                                                                      verbose = verbose)}
  if(Normalisation=='Census') {NormData <- .Census.calc(countData=countData,
                                                        batchData=batchData,
                                                        Lengths=Lengths,
                                                        MeanFragLengths=MeanFragLengths,
                                                        spikeData=spikeData,
                                                        spikeInfo = spikeInfo,
                                                        NCores=NCores,
                                                        verbose=verbose)}
  if(Normalisation=='RUV') {NormData <- .RUV.calc(countData = countData,
                                                  spikeData = spikeData,
                                                  batchData = batchData,
                                                  verbose = verbose)}
  # if(Normalisation=='BASiCS') {NormData <- .BASiCS.calc(countData = countData,
  #                                                       spikeData = spikeData,
  #                                                       spikeInfo = spikeInfo,
  #                                                       batchData = batchData,
  #                                                       verbose = verbose)}
  if(Normalisation=='depth') {NormData <- .depth.calc(countData = countData,
                                                      verbose = verbose)}
  if(Normalisation=='SF') {NormData <- .sf.calc(countData = countData,
                                                  sf=sf,
                                                  verbose = verbose)}
  if(Normalisation=='none') {NormData <- .none.calc(countData = countData,
                                                    verbose = verbose)}
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


# UQ ---------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.UQ.calc <- function(countData, verbose) {
  norm.factors <- edgeR::calcNormFactors(object=countData,
                                         method='upperquartile')
  sf <- norm.factors * colSums(countData)
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "UQ"
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

# poscounts -------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors
.PosCounts.calc <- function(countData, spikeData, verbose) {
  spike = ifelse(is.null(spikeData), FALSE, TRUE)

  if(spike==TRUE) {
    message(paste0("Using controlGenes, i.e. spike-ins for normalisation!"))
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    controlgenes = grepl(pattern='ERCC', rownames(cnts))
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cnts,
                                                           colData=data.frame(group=rep('A', ncol(countData))),
                                                           design=~1))
    dds <- DESeq2::estimateSizeFactors(dds, type='poscounts', controlGenes = controlgenes)
  }
  if(spike==FALSE) {
    message(paste0("Using all genes for normalisation!"))
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                                           colData=data.frame(group=rep('A', ncol(countData))),
                                                           design=~1))
    dds <- DESeq2::estimateSizeFactors(dds, type='poscounts')
  }

  sf <- DESeq2::sizeFactors(dds)
  names(sf) <- colnames(countData)
  norm.counts <- DESeq2::counts(dds, normalized=TRUE)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "PosCounts"
  return(res)
}


# Linnorm -----------------------------------------------------------------

#' @importFrom Linnorm Linnorm.Norm
#' @importFrom utils capture.output
.linnormnorm.calc <- function(countData, spikeData, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spikeID = rownames(spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
    spikeID = NULL
  }

  if(isTRUE(verbose)) {
    message(paste0("Linnorm messages:"))
    linnorm.out <- Linnorm::Linnorm.Norm(datamatrix = cnts, # matrix/data.frame of raw counts
                                         RowSamples = FALSE, # if I would switch gene and samples position in input
                                         spikein = spikeID, # names of the spike-ins
                                         spikein_log2FC = NULL,  # LFC of spike-ins, assume mix 1, so no
                                         showinfo = TRUE, # verbosity
                                         output = "Raw", # type of output: raw, ie  total count output will  be equal to median of input counts
                                         minNonZeroPortion = 0.5, # minimum porportion of nonzero values per gene, needed to put this to lower value otherwise 0 and NA normalized count matrix!
                                         BE_F_p = 0.3173, # p-value cutoff for standard deviation and skewness testing of batch effect normalisation
                                         BE_F_LC_Genes = "Auto", # porportion of lowly expressed genes to filter before batch effect normalisation
                                         BE_F_HC_Genes = 0.01, # proportion of highly expressed genes to filter before batch effect normalisation
                                         BE_strength = 0.5, # strength of batch effect normalisation
                                         max_F_LC = 0.75) # maximum threshold for filtering of low and high expression
  }

  if(!isTRUE(verbose)) {
    invisible(utils::capture.output(
      linnorm.out <- suppressMessages(Linnorm::Linnorm.Norm(datamatrix = cnts, # matrix/data.frame of raw counts
                                                            RowSamples = FALSE, # if I would switch gene and samples position in input
                                                            spikein = spikeID, # names of the spike-ins
                                                            spikein_log2FC = NULL,  # LFC of spike-ins, assume mix 1, so no
                                                            showinfo = TRUE, # verbosity
                                                            output = "Raw", # type of output: raw, ie  total count output will  be equal to median of input counts
                                                            minNonZeroPortion = 0.5, # minimum porportion of nonzero values per gene
                                                            BE_F_p = 0.3173, # p-value cutoff for standard deviation and skewness testing of batch effect normalisation
                                                            BE_F_LC_Genes = "Auto", # porportion of lowly expressed genes to filter before batch effect normalisation
                                                            BE_F_HC_Genes = 0.01, # proportion of highly expressed genes to filter before batch effect normalisation
                                                            BE_strength = 0.5, # strength of batch effect normalisation
                                                            max_F_LC = 0.75) # maximum threshold for filtering of low and high expression
      )
    ))

  }

  norm.counts <- linnorm.out[!grepl(pattern="ERCC", rownames(linnorm.out)),]
  gsf <-  countData / norm.counts
  wmu <- rowMeans(norm.counts, na.rm = TRUE)
  sf <- apply(gsf, 2, function(x) {
    stats::weighted.mean(x = x, w = wmu, na.rm = TRUE)
  })
  sf[is.infinite(sf)] <- mean(sf[is.finite(sf)])
  names(sf) <- colnames(countData)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "Linnorm"
  return(res)
}

# scran ------------------------------------------------------

#' @importFrom scran computeSumFactors computeSpikeFactors
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike
.scran.calc <- function(countData, spikeData, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)

  if(spike==TRUE) {
    message(paste0("Using computeSpikeFactors, i.e. spike-ins for normalisation!"))
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
    SingleCellExperiment::isSpike(sce, "Spike") <- nrow(countData) + seq_len(nrow(spikeData))
    sf <- scran::computeSpikeFactors(sce, sf.out=T)
  }

  if(spike==FALSE) {
    message(paste0("Using computeSumFactors, i.e. deconvolution over all cells!"))
    cnts = countData
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
    if(ncol(countData)<=14) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1)))
      sf <- scran::computeSumFactors(sce,
                                     sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(ncol(countData)>14 & ncol(countData)<=50) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::computeSumFactors(sce,
                                     sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(ncol(countData)>50 & ncol(countData)<=1000) {
      sizes <- c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::computeSumFactors(sce,
                                     sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(trunc(ncol(countData))>1000) {
      sizes <- c(round(seq(from=20, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::computeSumFactors(sce, positive=FALSE, sf.out=TRUE)
    }
  }

  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}

#' @importFrom scran computeSumFactors
#' @importFrom SingleCellExperiment SingleCellExperiment
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
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   clusters = clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   clusters = clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   clusters = clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(trunc(ncol(countData))>1000) {
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes = sizes,
                                   clusters = clusters,
                                   positive=FALSE,
                                   sf.out=TRUE)
  }

  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}


#' @importFrom scran computeSumFactors quickCluster
#' @importFrom SingleCellExperiment SingleCellExperiment
.scranclust.calc <- function(countData, PreclustNumber, verbose) {
  if (verbose) { message(paste0("Deconvolution within clusters.")) }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countData)))

  if(ncol(countData)<=14) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>1000 & ncol(countData)<=5000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(sce, method="igraph", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }

  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}




# SCnorm -----------------------------------------------------

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats weighted.mean
#' @importFrom utils capture.output
.SCnorm.calc <- function(countData, spikeData, batchData, NCores, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
  }
  if(!is.null(batchData)) {
    if(is.vector(batchData)) {
      cond <- batchData
    }
    if(!is.vector(batchData)) {
      cond <- batchData[,1]
    }
  }
  if(is.null(batchData)) {
    cond <- rep('a', ncol(cnts))
  }

  if (verbose) {
    if(is.null(batchData)) {
      message(paste0("No group annotation considered."))
    }
    if(!is.null(batchData)) {
      message(paste0("Using group annotation."))
    }

    }

  ncores = ifelse(is.null(NCores), 1, NCores)

  FilterCellNum = round(min(table(cond))*0.3)

  if(FilterCellNum<10) {
    FilterCellNum=10
  }

  FilterExpression = 0.5

  if(isTRUE(verbose)) {
    message(paste0("SCnorm messages:"))
    scnorm.out <- SCnorm::SCnorm(Data = cnts,
                                 Conditions = cond,
                                 PrintProgressPlots = FALSE,
                                 reportSF = TRUE,
                                 FilterCellNum = FilterCellNum,
                                 FilterExpression = FilterExpression,
                                 Thresh = 0.1,
                                 K = NULL,
                                 NCores = ncores,
                                 ditherCounts = TRUE,
                                 PropToUse = 0.1,
                                 Tau = 0.5,
                                 withinSample = NULL,
                                 useSpikes = spike,
                                 useZerosToScale = FALSE)
  }

  if(!isTRUE(verbose)) {
    invisible(utils::capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                 Conditions = cond,
                                 PrintProgressPlots = FALSE,
                                 reportSF = TRUE,
                                 FilterCellNum = FilterCellNum,
                                 FilterExpression = FilterExpression,
                                 Thresh = 0.1,
                                 K = NULL,
                                 NCores = ncores,
                                 ditherCounts = TRUE,
                                 PropToUse = 0.1,
                                 Tau = 0.5,
                                 withinSample = NULL,
                                 useSpikes = spike,
                                 useZerosToScale = FALSE))
    ))
  }

  if(class(scnorm.out) == "SingleCellExperiment") {
    norm.counts <- scnorm.out@assays$data$normcounts[!grepl(pattern="ERCC",
                                                             rownames(scnorm.out@assays$data$normcounts)),]
    scale.facts <- scnorm.out@metadata$ScaleFactors[!grepl(pattern="ERCC",
                                                           rownames(scnorm.out@metadata$ScaleFactors)),]
    gsf <- scale.facts
  }
  if(class(scnorm.out) == "SummarizedExperiment") {
    norm.counts <- scnorm.out@metadata$NormalizedData[!grepl(pattern="ERCC",
                                                             rownames(scnorm.out@metadata$NormalizedData)),]
    scale.facts <- scnorm.out@metadata$ScaleFactors[!grepl(pattern="ERCC",
                                                           rownames(scnorm.out@metadata$ScaleFactors)),]
    gsf <- scnorm.out@metadata$ScaleFactors
    rownames(gsf) <- rownames(cnts)
    colnames(gsf) <- colnames(cnts)
    gsf <- gsf[!grepl(pattern="ERCC", rownames(gsf)),]
  }

  wmu <- rowMeans(norm.counts, na.rm = TRUE)
  sf <- apply(scale.facts, 2, function(x) {
    stats::weighted.mean(x = x, w = wmu, na.rm = TRUE)
    })
  sf[is.infinite(sf)] <- mean(sf[is.finite(sf)])
  names(sf) <- colnames(cnts)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              scale.factors=gsf)


  attr(res, 'normFramework') <- "SCnorm"

  invisible(gc()) # hopefully this helps with the resource issue

  return(res)
}

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats weighted.mean
#' @importFrom scran quickCluster
.SCnormclust.calc <- function(countData, spikeData, PreclustNumber, NCores, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
  }

  if(ncol(countData)<=5000) {
    clusters <- scran::quickCluster(cnts, method="hclust", min.size=floor(PreclustNumber*0.75))
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(cnts, method="igraph", min.size=floor(PreclustNumber*0.75))
  }
  cond <- as.character(clusters)

  if (verbose) { message(paste0("SCnorm with preclustering to define groups.")) }

  ncores = ifelse(is.null(NCores), 1, NCores)

  FilterCellNum = round(min(table(cond))*0.3)

  if(FilterCellNum<10) {
    FilterCellNum=10
  }

  FilterExpression = 0.5

  if (verbose) { message(paste0("SCnorm messages:")) }

  scnorm.out <- SCnorm::SCnorm(Data = cnts,
                               Conditions = cond,
                               PrintProgressPlots = FALSE,
                               reportSF = TRUE,
                               FilterCellNum = FilterCellNum,
                               FilterExpression = FilterExpression,
                               Thresh = 0.1,
                               K = NULL,
                               NCores = ncores,
                               ditherCounts = TRUE,
                               PropToUse = 0.1,
                               Tau = 0.5,
                               withinSample = NULL,
                               useSpikes = spike,
                               useZerosToScale = FALSE)

  if(class(scnorm.out) == "SingleCellExperiment") {
    norm.counts <- scnorm.out@assays$data$normcounts[!grepl(pattern="ERCC",
                                                            rownames(scnorm.out@assays$data$normcounts)),]
    scale.facts <- scnorm.out@metadata$ScaleFactors[!grepl(pattern="ERCC",
                                                           rownames(scnorm.out@metadata$ScaleFactors)),]
    gsf <- scale.facts
  }
  if(class(scnorm.out) == "SummarizedExperiment") {
    norm.counts <- scnorm.out@metadata$NormalizedData[!grepl(pattern="ERCC",
                                                             rownames(scnorm.out@metadata$NormalizedData)),]
    scale.facts <- scnorm.out@metadata$ScaleFactors[!grepl(pattern="ERCC",
                                                           rownames(scnorm.out@metadata$ScaleFactors)),]
    gsf <- scnorm.out@metadata$ScaleFactors
    rownames(gsf) <- rownames(cnts)
    colnames(gsf) <- colnames(cnts)
    gsf <- gsf[!grepl(pattern="ERCC", rownames(gsf)),]
  }

  wmu <- rowMeans(norm.counts, na.rm = TRUE)
  sf <- apply(scale.facts, 2, function(x) {
    stats::weighted.mean(x = x, w = wmu, na.rm = TRUE)
  })
  sf[is.infinite(sf)] <- mean(sf[is.finite(sf)])
  names(sf) <- colnames(cnts)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "SCnorm"

  invisible(gc()) # hopefully this helps with the resource issue

  return(res)
}


# Census ----------------------------------------------------

#' @importFrom monocle newCellDataSet relative2abs
#' @importFrom Biobase exprs
#' @importFrom VGAM tobit negbinomial.size
#' @importFrom parallel detectCores
.Census.calc <- function(countData,
                         batchData,
                         spikeData,
                         spikeInfo,
                         Lengths,
                         MeanFragLengths,
                         NCores,
                         verbose) {

  ncores = ifelse(is.null(NCores), 1, NCores)

  # calculate TPM / CPM of genes
  if(!is.null(Lengths) && !is.null(MeanFragLengths)) {
    if(verbose) {message(paste0("Calculating TPM using Lengths and MeanFragLengths."))}
    ed <- .counts_to_tpm(countData=countData,
                         Lengths=Lengths,
                         MeanFragLengths=MeanFragLengths)
  }
  if(!is.null(Lengths) && is.null(MeanFragLengths)) {
    if(verbose) {message(paste0("Calculating TPM using Lengths."))}
    ed <- .calculateTPM(countData = countData, Lengths=Lengths)
  }
  if(is.null(Lengths) && is.null(MeanFragLengths)) { # if Length and MeanFragLengths is NULL, assume UMI method (3' / 5' prime counting)
    if(verbose) {message(paste0("Using counts for UMI data set."))}
    ed <- countData
  }

  # make annotated dataframes for monocle
  # for genes
  gene.dat <- data.frame(row.names = rownames(countData),
                         gene_short_name = rownames(countData),
                         biotype=rep("protein_coding", nrow(countData)),
                         num_cells_expressed=rowSums(countData>0))
  # for cells
  if(!is.null(batchData)) {
    if(!is.vector(batchData)) {
      cell.dat <- data.frame(row.names=colnames(countData),
                             Group=batchData[,1])
      ModelFormula <- "~Group"
    }
    if(is.vector(batchData)) {
      cell.dat <- data.frame(row.names=colnames(countData),
                             Group=batchData)
      ModelFormula <- "~Group"
    }


  }
  if(is.null(batchData)) {
    cell.dat <- data.frame(row.names=colnames(countData),
                           Group=rep("A", ncol(countData)))
    ModelFormula <- "~1"
  }

  fd <- new("AnnotatedDataFrame", data = gene.dat)
  pd <- new("AnnotatedDataFrame", data = cell.dat)


  if(!is.null(Lengths)) {
    # construct cell data set with expression values
    if(verbose) {message(paste0("Using tobit expression distribution."))}
    cds <- monocle::newCellDataSet(cellData = as.matrix(ed),
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit=0.1,
                                   expressionFamily=VGAM::tobit(Lower=0.1))
    # estimate RNA counts
    est_t <- .estimate_t(relative_expr_matrix=Biobase::exprs(cds),
                         relative_expr_thresh=0)

    if(is.null(spikeData) && is.null(spikeInfo)) {
      if(verbose) {message(paste0("Running monocle relative2abs."))}
      rpc_matrix <- monocle::relative2abs(relative_cds = cds,
                                          t_estimate = est_t,
                                          modelFormulaStr = ModelFormula,
                                          method = "num_genes",
                                          cores = ncores,
                                          verbose = verbose)

      # create cell data set wih RPC
      eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                     phenoData = pd,
                                     featureData = fd,
                                     lowerDetectionLimit=0.5,
                                     expressionFamily=VGAM::negbinomial.size())
      # apply normalisation
      sf <- .estimateSizeFactorsForDenseMatrix(counts=Biobase::exprs(eds),
                                               locfunc=median,
                                               round_exprs = T,
                                               method = "mean-geometric-mean-total")
      names(sf) <- colnames(Biobase::exprs(eds))

      norm.counts <- t(t(countData)/sf)
    }

    if(!is.null(spikeData) && !is.null(spikeInfo)) {
      if(verbose) {message(paste0("Running own relative2abs."))}
      # relative expression of spike-ins
      if(!is.null(spikeInfo$Lengths)) {
        ed.spike <- .calculateTPM(countData = spikeData,
                                  Lengths=spikeInfo$Lengths)
      }
      if(is.null(spikeInfo$Lengths)) {
        ed.spike <- .calculateCPM(countData = spikeData)
      }

      rel.results <- .relative2abs(relative_cds = cds,
                                   relative_spike = ed.spike,
                                   spikeInfo = spikeInfo,
                                   verbose = verbose,
                                   cores = ncores)
      rpc_matrix <- rel.results$norm_cds

      #create cell data set wih RPC
      eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                     phenoData = pd,
                                     featureData = fd,
                                     lowerDetectionLimit=0.5,
                                     expressionFamily=VGAM::negbinomial.size())

      # apply normalisation
      sf <- .estimateSizeFactorsForDenseMatrix(counts=Biobase::exprs(eds),
                                               locfunc=median,
                                               round_exprs = T,
                                               method = "mean-geometric-mean-total")
      #  It can be either "mean-geometric-mean-total" (default), "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
      names(sf) <- colnames(Biobase::exprs(eds))

      norm.counts <- t(t(countData)/sf)
    }
  }

  if(is.null(Lengths) && is.null(MeanFragLengths)) {
    if(verbose) {message(paste0("Using negbinomial.size expression distribution."))}
    # construct cell data set with expression values
    eds <- monocle::newCellDataSet(cellData = as.matrix(ed),
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit=0.5,
                                   expressionFamily=VGAM::negbinomial.size())
    # apply normalisation
    sf <- .estimateSizeFactorsForDenseMatrix(counts=Biobase::exprs(eds),
                                             locfunc=median,
                                             round_exprs = T,
                                             method = "mean-geometric-mean-total")
    names(sf) <- colnames(Biobase::exprs(eds))
    norm.counts <- t(t(countData)/sf)
    }

  # return object
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "Census"
  return(res)
}

# census estimatesizefactors function
.estimateSizeFactorsForDenseMatrix <- function(counts,
                                               locfunc,
                                               round_exprs,
                                               method){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(t(t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- .estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  return(sfs)
}

# census estimate_t function
.estimate_t <- function(relative_expr_matrix, relative_expr_thresh) {
  #apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr)
    10^mean(.dmode(log10(relative_expr[relative_expr > relative_expr_thresh])))))
}

# census mode function
#' @importFrom stats density
.dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- stats::density(x, kernel=c("gaussian"), na.rm = T)
  ( den$x[den$y==max(den$y)] )
}

#' @importFrom Biobase exprs
#' @importFrom MASS rlm
#' @importFrom stats predict
#' @importFrom parallel mcmapply
.relative2abs <- function(relative_cds,
                          relative_spike, # relative exprs matrix of spike-ins
                          spikeInfo, # spike input
                          verbose,
                          cores) {
  FPKM <- NULL
  # relative expression matrix of genes
  relative_expr_matrix <- Biobase::exprs(relative_cds)

  # relative expression of spike-ins
  ERCC_controls <- relative_spike

  # spike input information
  ERCC_annotation <- spikeInfo
  valid_ids <- which(ERCC_annotation[, "SpikeInput"] >= 0)

  # robust linear regression
  if (verbose) {message("Performing robust linear regression for each cell based
                        on the spike-in data")}
  molModels <- apply(ERCC_controls, 2, function(cell_exprs, input.ERCC.annotation, valid_ids) {
    spike_df <- input.ERCC.annotation
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$SpikeInput
    spike_df$rounded_numMolecules <- round(spike_df$numMolecules)
    if (is.null(valid_ids))
      spike_df <- subset(spike_df, FPKM >= 1e-10)
    else {
      spike_df <- spike_df[valid_ids, ]
      spike_df <- subset(spike_df, FPKM >= 1e-10)
    }
    spike_df$log_fpkm <- log10(spike_df$FPKM)
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    molModel <- tryCatch({
      molModel <- MASS::rlm(log_numMolecules ~ log_fpkm,
                            data = spike_df)
      molModel
    }, error = function(e) {
      print(e)
      NULL
    })
    molModel
  }, ERCC_annotation, valid_ids)

  if (verbose) {message("Apply the fitted robust linear regression model
                        to recover the absolute copy number for all transcripts in each cell")}
  norm_fpkms <- parallel::mcmapply(function(cell_exprs, molModel) {
    tryCatch({
      norm_df <- data.frame(log_fpkm = log10(cell_exprs))
      res <- 10^stats::predict(molModel, type = "response",
                        newdata = norm_df)
    }, error = function(e) {
      rep(NA, length(cell_exprs))
    })
  }, split(as.matrix(relative_expr_matrix), rep(1:ncol(relative_expr_matrix),
                                                each = nrow(relative_expr_matrix))), molModels, mc.cores = cores)
  k_b_solution <- data.frame(b = unlist(lapply(molModels,
                                               FUN = function(x) {
                                                 intercept = x$coefficients[1]
                                               })), k = unlist(lapply(molModels, FUN = function(x) {
                                                 slope = x$coefficients[2]
                                               })))
  kb_model <- MASS::rlm(b ~ k, data = k_b_solution)
  kb_slope <- kb_model$coefficients[2]
  kb_intercept <- kb_model$coefficients[1]

  rownames(norm_fpkms) <- rownames(relative_expr_matrix)
  colnames(norm_fpkms) <- colnames(relative_expr_matrix)

  # return object
  return(list(norm_cds = norm_fpkms,
              kb_slope = kb_slope,
              kb_intercept = kb_intercept,
              k_b_solution = k_b_solution))
}

# RUV ---------------------------------------------

#' @importFrom RUVSeq RUVg
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr %>%
#' @importFrom dplyr rename_
.RUV.calc <- function(countData, spikeData, batchData, verbose) {

  if(is.null(batchData) && !is.null(spikeData)) {
    # annotate spike-ins and genes
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spike = grepl(pattern='ERCC', rownames(cnts))

    #RUVg calculation
    RUV.out = RUVSeq::RUVg(x=as.matrix(cnts),
                            cIdx=spike,
                            k=1,
                            drop=0,
                            center=TRUE,
                            round=FALSE,
                            epsilon=1,
                            tolerance=1e-8,
                            isLog=FALSE)
  }

  if(!is.null(batchData) && !is.null(spikeData)) {
    # prepare replicate sample annotation matrix
    tmp <- batchData %>%
      dplyr::rename_(Batch = names(.)[1]) %>%
      tibble::rownames_to_column(var="SampleID")  %>%
      tibble::rownames_to_column(var="SampleNumber")
    tmp2 <- split(x = tmp, f = as.factor(tmp$Batch))
    longestbatch <- max(sapply(tmp2, nrow))
    tmp3 <- sapply(tmp2, function(x) {
      tmp <- as.numeric(x$SampleNumber)
      tmp1 <- c(tmp, rep(-1, longestbatch-length(tmp)))
      tmp1
    })
    differences = t(tmp3)

    # annotate spike-ins and genes
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spike = grepl(pattern='ERCC', rownames(cnts))

    #RUVs calculation
    RUV.out = RUVSeq::RUVs(x=as.matrix(cnts),
                           cIdx=spike,
                           scIdx = differences,
                           k=1,
                           round=FALSE,
                           epsilon=1,
                           tolerance=1e-8,
                           isLog=FALSE)
  }

  if(!is.null(batchData) && is.null(spikeData)) {
    # prepare replicate sample annotation matrix
    tmp <- batchData %>%
      dplyr::rename_(Batch = names(.)[1]) %>%
      tibble::rownames_to_column(var="SampleID")  %>%
      tibble::rownames_to_column(var="SampleNumber")
    tmp2 <- split(x = tmp, f = as.factor(tmp$Batch))
    longestbatch <- max(sapply(tmp2, nrow))
    tmp3 <- sapply(tmp2, function(x) {
      tmp <- as.numeric(x$SampleNumber)
      tmp1 <- c(tmp, rep(-1, longestbatch-length(tmp)))
      tmp1
    })
    differences = t(tmp3)

    # annotate genes
    cnts = countData
    controls = rownames(countData)

    #RUVs calculation
    RUV.out = RUVSeq::RUVs(x=as.matrix(cnts),
                           cIdx=controls,
                           scIdx = differences,
                           k=1,
                           round=FALSE,
                           epsilon=1,
                           tolerance=1e-8,
                           isLog=FALSE)
  }

  # return object
  # normalized counts
  normCounts <- RUV.out$normalizedCounts[!grepl(pattern = "ERCC",
                                                rownames(RUV.out$normalizedCounts)),]
  normCounts[normCounts<0] <- 0
  # RUV proxy size factors
  gsf <-  t(t(countData)/t(normCounts))
  sf <- apply(gsf, 2, stats::median, na.rm=T)
  sf <- sf*length(sf)/sum(sf)
  names(sf) <- colnames(countData)
  res <- list(NormCounts = normCounts,
              RoundNormCounts = round(normCounts),
              size.factors = sf,
              RUV.W = RUV.out$W)
  attr(res, 'normFramework') <- "RUV"
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



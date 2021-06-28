
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
  if(Normalisation=='UQ') {
    NormData <- .UQ.calc(countData = countData,
                         verbose = verbose)
  }
  if(Normalisation=='MR') {
    NormData <- .MR.calc(countData = countData,
                         spikeData = spikeData,
                         verbose = verbose)
  }
  if(Normalisation=='PosCounts') {
    NormData <- .PosCounts.calc(countData = countData,
                                spikeData = spikeData,
                                verbose = verbose)
  }
  if(Normalisation=='Linnorm') {
    NormData <- .linnormnorm.calc(countData = countData,
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
  if(Normalisation=='SCnorm' && Step == "Estimation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = NULL,
                             NCores = NCores,
                             verbose = verbose)
  }
  if(Normalisation=='SCnorm' && Label == "none" && Step == "Simulation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = NULL,
                             NCores = NCores,
                             verbose = verbose)
  }

  if(Normalisation=='SCnorm' && Label == "known" && Step == "Simulation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = batchData,
                             NCores = NCores,
                             verbose = verbose)
    }
  if(Normalisation=='SCnorm' && Label == "clustering" && Step == "Simulation"){
    NormData <- .SCnormclust.calc(countData = countData,
                                  spikeData = spikeData,
                                  PreclustNumber=PreclustNumber,
                                  NCores = NCores,
                                  verbose = verbose)
    }
  if(Normalisation=='sctransform') {
    NormData <- .sctransform.calc(countData = countData,
                                  batchData = batchData,
                                  Step = Step,
                                  NCores = NCores,
                                  verbose = verbose)
  }
  if(Normalisation=="bayNorm") {
    NormData <- .baynorm.calc(countData = countData,
                              batchData = batchData,
                              spikeData = spikeData,
                              spikeInfo = spikeInfo,
                              NCores  = NCores,
                              Step = Step,
                              Protocol = Protocol,
                              verbose = verbose)
  }
  if(Normalisation=='Census') {
    NormData <- .Census.calc(countData = countData,
                             Lengths = Lengths,
                             MeanFragLengths = MeanFragLengths,
                             spikeData = spikeData,
                             spikeInfo = spikeInfo,
                             Protocol = Protocol,
                             NCores = NCores,
                             verbose = verbose)
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
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "Linnorm"
  return(res)
}

# scran ------------------------------------------------------

#' @importFrom scran calculateSumFactors 
#' @importFrom scuttle computeSpikeFactors
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
    sf <- SingleCellExperiment::sizeFactors(sce)

  }

  if(spike==FALSE) {
    message(paste0("Using calculateSumFactors, i.e. deconvolution over all cells!"))
    cnts = countData
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
    if(ncol(countData)<=14) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1)))
      sf <- scran::calculateSumFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(ncol(countData)>14 & ncol(countData)<=50) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::calculateSumFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(ncol(countData)>50 & ncol(countData)<=1000) {
      sizes <- c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::calculateSumFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
    }
    if(trunc(ncol(countData))>1000) {
      sizes <- c(round(seq(from=20, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::calculateSumFactors(sce,
                                       sizes = sizes,
                                       positive = TRUE)
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
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(trunc(ncol(countData))>1000) {
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }

  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}


#' @importFrom scran calculateSumFactors quickCluster
#' @importFrom SingleCellExperiment SingleCellExperiment
.scranclust.calc <- function(countData, PreclustNumber, verbose) {
  if (verbose) { message(paste0("Deconvolution within clusters.")) }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countData)))

  if(ncol(countData)<=14) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>1000 & ncol(countData)<=5000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(sce, method="igraph", min.size=floor(PreclustNumber*0.75))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::calculateSumFactors(sce,
                                     sizes = sizes,
                                     clusters = clusters,
                                     positive = TRUE)
  }

  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}


# sctransform -------------------------------------------------------------

#' @importFrom sctransform vst correct
#' @importFrom future plan
#' @importFrom Matrix Matrix
.sctransform.calc <- function(countData, batchData, Step, NCores, verbose) {

  if(Step == "Estimation") {
    if(!is.null(batchData)) {
      if(is.vector(batchData)) {
        cond <- data.frame(batch = batchData, stringsAsFactors = FALSE)
        batch_var <- 'batch'
      }
      if(!is.vector(batchData)) {
        cond <- data.frame(batch = batchData[,1], stringsAsFactors = FALSE)
        batch_var <- 'batch'
      }
    }
    if(is.null(batchData)) {
      cond <- NULL
      batch_var <- NULL
    }
  }

  if(Step == "Simulation") {
    cond <- NULL
    batch_var <- NULL
  }

  if(!is.null(NCores)){
    future::plan(strategy = 'multicore', workers = NCores)
    options(future.globals.maxSize = 10 * 1024 ^ 3)
  }

  umi <- Matrix::Matrix(as.matrix(countData), sparse = TRUE)

  sctransform_out <- sctransform::vst(umi = umi,
                                      cell_attr = cond,
                                      latent_var = c("log_umi"),
                                      batch_var = batch_var,
                                      latent_var_nonreg = NULL,
                                      n_genes = NULL,
                                      n_cells = NULL,
                                      method = "poisson",
                                      do_regularize = TRUE,
                                      res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                                      bin_size = 256, min_cells = 5,
                                      residual_type = "pearson",
                                      return_cell_attr = TRUE,
                                      return_gene_attr = TRUE,
                                      return_corrected_umi = TRUE,
                                      bw_adjust = 3, gmean_eps = 1,
                                      theta_given = NULL,
                                      show_progress = verbose)
  norm.counts <- sctransform::correct(sctransform_out)

  ixx.valid <- rownames(countData)  %in% rownames(norm.counts)
  NormCounts <- countData
  NormCounts[ixx.valid, ] <- norm.counts

  gsf <-  countData / NormCounts
  gsf[is.infinite(as.matrix(gsf))] <- NA
  wmu <- sctransform_out$gene_attr$gmean
  sf <- apply(gsf[ixx.valid,], 2, function(x) {
    stats::weighted.mean(x = x, w = wmu, na.rm = TRUE)
  })
  sf[is.infinite(sf)] <- mean(sf[is.finite(sf)])
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)

  res <- list(NormCounts=NormCounts,
              RoundNormCounts=round(NormCounts),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "sctransform"

  invisible(gc())

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
    rownames(gsf) <- rownames(cnts)
    colnames(gsf) <- colnames(cnts)
    gsf <- gsf[!grepl(pattern="ERCC", rownames(gsf)),]
  }
  if(class(scnorm.out) == "SummarizedExperiment") {
    norm.counts <- scnorm.out@metadata$NormalizedData[!grepl(pattern="ERCC",
                                                             rownames(scnorm.out@metadata$NormalizedData)),]
    scale.facts <- scnorm.out@metadata$ScaleFactors[!grepl(pattern="ERCC",
                                                           rownames(scnorm.out@metadata$ScaleFactors)),]
    gsf <- scale.facts
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



# bayNorm -----------------------------------------------------------------

#' @importFrom bayNorm bayNorm BetaFun
.baynorm.calc <- function(countData,
                          batchData,
                          spikeData,
                          spikeInfo,
                          NCores,
                          Step,
                          Protocol,
                          verbose){
  # define parallel computation
  if(is.null(NCores)) {
    parallel = FALSE
    ncores = 1
  }
  if(!is.null(NCores)) {
    parallel = TRUE
    ncores = NCores
  }

  if(!is.null(batchData)) {
    if(is.vector(batchData)){
      cond <- batchData
    }
    if(!is.vector(batchData)) {
      cond <- batchData[,1]
    }
  }
  if(is.null(batchData)) {
    cond <- NULL
    ptype <- NULL
  }
  if(!is.null(cond)){
    ptype <- ifelse(Step == "Simulation", "LL", "GG")
  }

  if(Protocol == "Read"){
    norm.tmp <- .depth.calc(countData = countData, verbose = verbose)
    sffl <- norm.tmp$size.factors
  }

  if(Protocol == "UMI"){
    sffl <- NULL
  }

  if (verbose) {
    if(Protocol == "Read") {
      message(paste0("Since the provided count matrix is read-based, bayNorm needs scaling factors.
                     Using depth normalisation for calculation."))
    }
    if(is.null(cond)) {
      message(paste0("No condition annotation considered."))
    }
    if(!is.null(cond)) {
      message(paste0("Using condition annotation."))
    }
    if(any(is.null(spikeInfo), is.null(spikeData))) {
      message(paste0("No spike-in information provided to estimate capture efficiency beta needed for baynorm normalisation.
                     Be advised that the normalisation continues with default settings of bayNorm!"))
      beta_vec <- NULL
    }
    if(!is.null(spikeInfo) && !is.null(spikeData)) {
      message(paste0("Using spike-in information provided to estimate capture efficiency beta needed for baynorm normalisation."))
      ## cell capture efficiency
      total_ercc_molecules <- sum(spikeInfo$SpikeInput)
      cellcapture <- apply(spikeData, 2, function(i) {
        sum(i) / total_ercc_molecules
      })
      beta_vec <- bayNorm::BetaFun(Data = countData, MeanBETA = mean(cellcapture))$BETA
    }
  }

  norm_out <- bayNorm::bayNorm(Data = countData,
                               BETA_vec = beta_vec,
                               Conditions = cond,
                               UMI_sffl = sffl,
                               Prior_type = ptype,
                               mode_version = TRUE,
                               mean_version = FALSE,
                               S = 20,
                               parallel = parallel,
                               NCores = ncores,
                               FIX_MU = TRUE,
                               GR = FALSE,
                               BB_SIZE = TRUE,
                               verbose = verbose)

  if(is.null(cond)) {
    norm.counts <- norm_out$Bay_out
    gsf <- (countData/res) +1
  }
  if(!is.null(cond)) {
    res <- do.call("cbind", norm_out$Bay_out_list)
    norm.counts <- res[, match(colnames(countData), colnames(res))]
    gsf <- (countData/norm.counts) +1
  }

  wmu <- rowMeans(norm.counts, na.rm = TRUE)
  sf <- apply(gsf, 2, function(x) {
    stats::weighted.mean(x = x, w = wmu, na.rm = TRUE)
  })
  sf[is.infinite(sf)] <- mean(sf[is.finite(sf)])
  names(sf) <- colnames(norm.counts)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              scale.factors=gsf)

}



# Census ----------------------------------------------------

#' @importFrom parallel detectCores
.Census.calc <- function(countData,
                         Lengths,
                         MeanFragLengths,
                         spikeData,
                         spikeInfo,
                         Protocol,
                         NCores,
                         verbose) {

  # house keeping settings
  ncores = ifelse(is.null(NCores), 1, NCores)

  # 1. depending on protocol calculate relative expression values (CPM, TPM, FPKM)

  if(Protocol == "UMI") {
    if(verbose) {message(paste0("Census relative2abs normalisation should not be used in conjunction with UMI data!"))}
    if(verbose) {message(paste0("Calculating CPM as gene expression data."))}
    ed <- .calculateCPM(countData = countData)

    if(!is.null(spikeData) && !is.null(spikeInfo)) {
      if(verbose) {message(paste0("Transform relative expression values into absolute transcript counts
                                  using spike-ins provided."))}
      ed.spike <- .calculateCPM(countData = spikeData)
    }
    if(is.null(spikeInfo)) {
      if(verbose) {message(paste0("Transform relative expression values into absolute transcript counts."))}
      ed.spike <- NULL
    }
  }

  if(Protocol == "Read") {
    if(!is.null(Lengths) && !is.null(MeanFragLengths)) {
      if(verbose) {message(paste0("Calculating TPM Gene Expression using Lengths and MeanFragLengths."))}
      ed <- .counts_to_tpm(countData=countData,
                           Lengths=Lengths,
                           MeanFragLengths=MeanFragLengths)
    }
    if(!is.null(Lengths) && is.null(MeanFragLengths)) {
      if(verbose) {message(paste0("Calculating TPM Gene Expression using Lengths."))}
      ed <- .calculateTPM(countData = countData, Lengths=Lengths)
    }
    if(is.null(Lengths) && is.null(MeanFragLengths)) {
      if(verbose) {message(paste0("Calculating CPM Gene Expression."))}
      ed <- .calculateCPM(countData = countData)
    }

    if(!is.null(spikeData) && !is.null(spikeInfo)) {
      if(verbose) {message(paste0("Transform relative expression values into absolute transcript counts
                                  using spike-ins provided."))}
      # relative expression of spike-ins
      if(!is.null(spikeInfo$Lengths)) {
        if(verbose) {message(paste0("Calculating TPM Spike-In Expression."))}
        ed.spike <- .calculateTPM(countData = spikeData,
                                  Lengths=spikeInfo$Lengths)
      }
      if(is.null(spikeInfo$Lengths)) {
        if(verbose) {message(paste0("Calculating TPM Spike-In Expression."))}
        ed.spike <- .calculateCPM(countData = spikeData)
      }
    }

    if(is.null(spikeInfo)) {
      if(verbose) {message(paste0("Transform relative expression values into absolute transcript counts."))}
      ed.spike <- NULL
    }
  }

  # 2. find the most commonly occuring relative expression value per cell
  est_t <- .estimate_t(relative_expr_matrix = ed,
                       relative_expr_thresh = 0.1)

  # 3. apply relative2abs
  rel.results <- .relative2abs(relative_expr_matrix = ed,
                               t_estimate = est_t,
                               relative_spike = ed.spike,
                               spikeInfo = spikeInfo,
                               verbose = verbose,
                               cores = ncores)
  # print(str(rel.results))
  rpc_matrix <- rel.results$norm_cds
  # print(rpc_matrix[1:5, 1:5])

  # 4. apply normalisation
  sf <- .estimateSizeFactorsForDenseMatrix(CM=rpc_matrix,
                                           locfunc=median,
                                           round_exprs = T,
                                           method = "mean-geometric-mean-total")
  #  It can be either "mean-geometric-mean-total" (default), "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
  names(sf) <- colnames(rpc_matrix)
  norm.counts <- t(t(countData)/sf)

  # return object
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              RPC = rpc_matrix,
              size.factors=sf)
  attr(res, 'normFramework') <- "Census"
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



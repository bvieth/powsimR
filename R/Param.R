
# estimateParam ---------------------------------------------------------

#' @name estimateParam
#' @aliases estimateParam
#' @title Estimate simulation parameters
#' @description This function estimates and returns parameters needed for power simulations assuming a negative binomial read count distribution.
#' @usage estimateParam(countData,
#' batchData = NULL,
#' spikeData = NULL,
#' spikeInfo = NULL,
#' Lengths = NULL,
#' MeanFragLengths = NULL,
#' Distribution = c('NB', 'ZINB'),
#' RNAseq = c('bulk', 'singlecell'),
#' normalisation = c('TMM', 'MR', 'PosCounts', 'UQ', 'scran',
#'                   'SCnorm', 'RUV', 'BASiCS', 'Census', 'depth', 'none'),
#' sigma = 1.96,
#' NCores = NULL,
#' verbose=TRUE)
#' @param countData is a count \code{matrix}.
#' Rows correspond to genes, columns to samples.
#' The gene names should be given as rownames without "_" in the names. The samples names should be given as colnames.
#' @param batchData is a \code{data.frame} for batch annotation.
#' Rows correspond to samples. The first column should contain the batches, e.g. 'a', 'b', 'c', etc.
#' @param spikeData is a count \code{matrix}.
#' Rows correspond to spike-ins, columns to samples.
#' The order of columns should be the same as in the \code{countData}.
#' @param spikeInfo is a molecule count \code{matrix} of spike-ins.
#' Rows correspond to spike-ins. The order of rows should be the same as in the \code{spikeData}.
#' The column names should be 'SpikeID' and 'SpikeInput' for molecule counts of spike-ins.
#' @param Lengths is a numeric vector of transcript lengths with the same length as rows in countData.
#' This variable is only used for internal TPM calculations if Census normalization is specified.
#' @param MeanFragLengths is a numeric vector of mean fragment lengths with the same length as columns in countData.
#' This variable is only used for internal TPM calculations if Census normalization is specified.
#' @param Distribution is a character value: "NB" for negative binomial or "ZINB" for zero-inflated negative binomial.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param normalisation is a character value: 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'scranCLUST',
#' 'SCnorm', 'RUV', 'BASiCS', 'Census', 'depth', 'none'.
#' For more information, please consult the details section.
#' @param sigma The variability band width for mean-dispersion loess fit defining the prediction interval for read cound simulation. Default is 1.96, i.e. 95\% interval. For more information see \code{\link[msir]{loess.sd}}.
#' @param NCores The number of cores for normalisation method SCnorm and Census.
#' If \code{NULL}, the number of detected cores minus 1 is used.
#' @param verbose Logical value to indicate whether to print function information. Default is \code{TRUE}.
#' @return List object with the following entries
#' \item{seqDepth}{Library size, i.e. total number of reads per library}
#' \item{means}{Mean normalized read counts per gene.}
#' \item{dispersion}{Dispersion estimate per gene.}
#' \item{common.dispersion}{The common dispersion estimate over all genes.}
#' \item{size}{Size parameter of the negative binomial distribution, i.e. 1/dispersion.}
#' \item{p0}{Probability that the count will be zero per gene.}
#' \item{meansizefit}{A loess fit relating log2 mean to log2 size for use in simulating new data (\code{\link[msir]{loess.sd}}).}
#' \item{meandispfit}{A fit relating log2 mean to log2 dispersion used for visualizing mean-variance dependency (\code{\link[msir]{loess.sd}}).}
#' \item{p0.cut}{The knee point of meanp0fit. Log2 mean values above that value have virtually no dropouts.}
#' \item{grand.dropout}{The percentage of empty entries in the count matrix.}
#' \item{sf}{The estimated library size factor per sample.}
#' \item{totalS,totalG}{Number of samples and genes provided.}
#' \item{estS,estG}{Number of samples and genes for which parameters can be estimated.}
#' \item{RNAseq}{The type of RNAseq: bulk or single cell.}
#' \item{normFramework}{The normalisation method used to calculate normalized parameters and library size factors.}
#' \item{sigma}{The width of the variability band.}
#' @details
#' Normalisation methods
#' \describe{
#' \item{TMM, UQ}{employ the edgeR style normalization of weighted trimmed mean of M-values and upperquartile
#' as implemented in \code{\link[edgeR]{calcNormFactors}}, respectively.}
#' \item{MR, PosCounts}{employ the DESeq2 style normalization of median ratio method and a modified geometric mean method
#' as implemented in \code{\link[DESeq2]{estimateSizeFactors}}, respectively.}
#' \item{scran, SCnorm}{apply the deconvolution and quantile regression normalization methods developed for sparse RNA-seq data
#' as implemented in \code{\link[scran]{computeSumFactors}} and \code{\link[SCnorm]{SCnorm}}, respectively.
#' For \code{SCnorm}, the user can also supply \code{spikeData}.}
#' \item{RUV}{removes unwanted variation. There are two approaches implemented:
#' (1) utilizing negative control genes, i.e. spike-ins stored in \code{spikeData} (\code{\link[RUVSeq]{RUVg}}).
#' (2) utilizing replicate samples, i.e. samples for which the covariates of interest are considered constant.
#' This annotation is stored in \code{batchData} (\code{\link[RUVSeq]{RUVs}}).}
#' \item{BASiCS}{removes technical variation by utilizing negative control genes, i.e. spike-ins stored in \code{spikeData},
#' as implemented in \code{\link[BASiCS]{DenoisedCounts}}.
#' Furthermore, the molecule counts of spike-ins added to the cell lysate need to be supplied in \code{spikeInfo}.}
#' \item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
#' Census at least needs \code{Lengths} for single-end data and preferably \code{MeanFragLengths} for paired-end data.
#' Do not use this algorithm for UMI data!}
#' \item{depth}{Sequencing depth normalisation.}
#' \item{none}{No normalisation is applied. This approach can be used for prenormalized expression estimates, e.g. cufflinks, RSEM or salmon.}
#' }
#' @examples
#' \dontrun{
#' ## simulating single cell RNA-seq experiment
#' ngenes <- 10000
#' ncells <- 100
#' true.means <- 2^runif(ngenes, 3, 6)
#' true.dispersions <- 3 + 100/true.means
#' sf.values <- 2^rnorm(ncells, sd=0.5)
#' sf.means <- outer(true.means, sf.values, '*')
#' cnts <- matrix(rnbinom(ngenes*ncells,
#'                        mu=sf.means, size=1/true.dispersions),
#'                ncol=ncells)
#' ## estimating negative binomial parameters
#' estparam <- estimateParam(countData=cnts,
#'                           spikeData=NULL, spikeInfo = NULL,
#'                           Lengths=NULL, MeanFragLengths=NULL,
#'                           Distribution='NB',
#'                           RNAseq="singlecell",
#'                           normalisation='scran',
#'                           NCores=NULL,
#'                           sigma=1.96)
#' plotParam(estparam, annot=F)
#'
#' ## simulating bulk RNA-seq experiment
#' ngenes <- 10000
#' nsamples <- 10
#' true.means <- 2^rnorm(ngenes, mean=8, sd=2)
#' true.dispersions <- 3/true.means + 0.1
#' sf.values <- rnorm(nsamples, mean=1, sd=0.1)
#' sf.means <- outer(true.means, sf.values, '*')
#' cnts <- matrix(rnbinom(ngenes*nsamples,
#'                        mu=sf.means, size=1/true.dispersions),
#'                ncol=nsamples)
#' ## estimating negative binomial parameters
#' estparam <- estimateParam(countData=cnts,
#'                           Distribution='NB', RNAseq="bulk",
#'                           normalisation='MR', sigma=1.96)
#' plotParam(estparam, annot=F)
#' }
#' @author Beate Vieth
#' @rdname estimateParam
#' @export
estimateParam <- function(countData,
                          batchData = NULL,
                          spikeData = NULL,
                          spikeInfo = NULL,
                          Lengths = NULL,
                          MeanFragLengths = NULL,
                          Distribution = c('NB', 'ZINB'),
                          RNAseq = c('bulk', 'singlecell'),
                          normalisation = c('TMM', 'MR', 'PosCounts', 'UQ', 'scran',
                                            'SCnorm', 'RUV', 'BASiCS', 'Census', 'depth', 'none'),
                          sigma = 1.96,
                          NCores = NULL,
                          verbose=TRUE) {

  # Inform users of inappropiate choices
  if (RNAseq == 'bulk' && Distribution=='ZINB') {
    if(verbose) {message("Zero-inflated negative binomial is not implemented for bulk data. Setting it to negative binomial.")}
    Distribution="NB"
  }
  if (RNAseq=='singlecell' && normalisation=='MR') {
    if(verbose) {message(paste0(normalisation, " has been developed for bulk data and it is rather likely that it will not work.
                                \nPlease consider using a method that can handle single cell data, e.g. PosCounts, scran, SCnorm."))}
  }
  if (normalisation=='Census') {
    if(verbose) {message(paste0(normalisation, " should only be used for non-UMI methods! \nFor more information, please consult the monocle vignette."))}
    if (normalisation=='Census' && is.null(Lengths)) {
      if(verbose) {message(paste0(normalisation, " should be used in combination with transcript lengths.
                                   \nIf the library is paired-end, please also provide the mean fragment lengths which can be determined by e.g. Picard."))}
    }
  }

  # STOP if combination of options is not supported/would not work!
  if (RNAseq=='bulk' && normalisation %in% c('BASiCS', 'Census')) {
    stop(message(paste0(normalisation, " is only developed and implemented for single cell RNA-seq experiments.")))
  }
  if (normalisation %in% c('BASiCS') && is.null(spikeData)) {
    stop(message(paste0(normalisation, " needs spike-in information! \nPlease provide an additional table of spike-in read counts.")))
  }
  if (normalisation == 'BASiCS' && is.null(spikeInfo)) {
    stop(message(paste0(normalisation, " needs spike-in information! \nPlease provide an additional table of spike-in molecule counts.")))
  }
  if (normalisation == 'BASiCS' && !is.null(spikeInfo) && !is.null(colnames(spikeInfo)) && !all(colnames(spikeInfo) %in% c("SpikeID", "SpikeInput"))) {
    stop(message(paste0(normalisation, " needs spike-in information! \nPlease provide an additional table of spike-in molecule counts with correctly named columns.")))
  }

  if (normalisation=='none' && all((countData - round(countData)) == 0)) {
    stop(message(paste0("No normalisation should only be done with pre-normalized data, eg. RSEM output, but the provided input matrix contains only integer values!")))
  }

  # check the matching of names
  checkup <- .run.checkup(
    countData = countData,
    batchData = batchData,
    spikeData = spikeData,
    spikeInfo = spikeInfo,
    Lengths = Lengths,
    MeanFragLengths=MeanFragLengths,
    verbose=verbose
  )

  countData <- checkup$countData
  batchData <- checkup$batchData
  spikeData <- checkup$spikeData
  spikeInfo <- checkup$spikeInfo
  Lengths <- checkup$Lengths
  MeanFragLengths <- checkup$MeanFragLengths

  # run estimation
  res = .run.estParam(countData = countData,
                      batchData = batchData,
                      spikeData = spikeData,
                      spikeInfo = spikeInfo,
                      Lengths = Lengths,
                      MeanFragLengths = MeanFragLengths,
                      Distribution = Distribution,
                      RNAseq = RNAseq,
                      normalisation = normalisation,
                      sigma = sigma,
                      NCores = NCores,
                      verbose=verbose)

  # inform users
  if (RNAseq=='bulk' && c(res$grand.dropout>0.5 || is.null(res$p0.cut) ))
    if(verbose) {message('The provided bulk data has frequent dropouts.
             Consider using single cell estimation framework.')}

  # make output
  res2 <- c(res, list(RNAseq = RNAseq,
                      normFramework=normalisation,
                      sigma=sigma))
  attr(res2, 'param.type') <- "estimated"
  attr(res2, 'Distribution') <- Distribution

  return(res2)
}


# estimateSpike -----------------------------------------------------------

#' @name estimateSpike
#' @aliases estimateSpike
#' @title Estimate simulation parameters for spike-ins
#' @description This function estimates and returns parameters needed for spike-in read count simulations using supplementary code from Kim et al. 2016 (DOI: 10.1038/ncomms9687).
#' @usage estimateSpike(spikeData,
#' spikeInfo,
#' MeanFragLength = NULL,
#' batchData = NULL,
#' normalisation=c('depth','none'))
#' @param spikeData  is a count \code{matrix}. Rows correspond to spike-ins, columns to samples.
#' Rownames should contain the spike-in names, column names the sample names.
#' @param spikeInfo is a molecule count \code{matrix} of spike-ins. Rows correspond to spike-ins.
#' The order of rows should be the same as the rows in \code{spikeData}.
#' The rownames should be the same as the rownames of spikeData.
#' The first column should contain the molecule counts of spike-ins and be named 'SpikeInput'.
#' The second column can contain the sequence lengths of spike-ins in bases and be named 'Lengths'.
#' @param MeanFragLength is a numeric vector of the mean fragment length.
#' @param batchData is a \code{data.frame} for batch annotation. Rows correspond to samples.
#' The order of rows should be the same as the columns in \code{spikeData}.
#' The rownames should be the same as the column names of spikeData.
#' The first column should contain the batch annotation, e.g. 'a' for batch 1, 'b' for batch 2.
#' @param normalisation is a character value: 'depth' or 'none'. For more information, please consult the details section.
#' @return List object with the following entries:
#' \item{normCounts}{The normalised spike-in read counts.}
#' \item{sf}{The estimated library size factors.}
#' \item{EVGammaThetaEstimates}{Estimation of the four parameters capturing technical variability, namely E[\eqn{\gamma}], Var[\eqn{\gamma}], E[\eqn{\theta}] and Var[\eqn{\theta}]. For more details, please consult supplementary information of Kim et al. 2016 (DOI: 10.1038/ncomms9687). These estimates are needed for simulating spike-in read counts.}
#' \item{Input}{The input spike-in expression \code{matrix}, molecule counts \code{data.frame} and batch annotation \code{data.frame}, filtered so that only spike-ins with nonzero expression and samples with at least 100 reads are retained.}
#' \item{Settings}{Reporting the chosen normalisation framework.}
#' @examples
#' \dontrun{
#' ## not yet
#' }
#' @details
#' Normalisation methods
#' \describe{
#' \item{'depth'}{applies the depth normalization method as implemented in \code{\link[scran]{computeSpikeFactors}}.}
#' \item{'none'}{No normalisation is applied. This approach can be used for prenormalized expression estimates, e.g. TPM/FPKM/RPKM estimated by cufflinks, RSEM, salmon etc.}
#' }
#' @author Beate Vieth
#' @rdname estimateSpike
#' @importFrom stats setNames
#' @importFrom matrixStats rowSds
#' @export
estimateSpike <- function(spikeData,
                          spikeInfo,
                          MeanFragLength = NULL,
                          batchData = NULL,
                          normalisation=c('depth','none')) {
  # check provided input
  if(!is.null(batchData)) {
    if(!nrow(batchData) == ncol(spikeData)) {
      stop(message(paste0("batchData and spikeData is not the same length!")))
    }
    if(is.null(rownames(batchData))) {
      stop(message(paste0("batchData has no sample annotation as row names!")))
    }
    if(is.null(colnames(spikeData))) {
      stop(message(paste0("spikeData has no sample annotation as column names!")))
    }
  }
  if(is.null(rownames(spikeData)) || is.null(rownames(spikeInfo))) {
    stop(message(paste0("The input data frames of spikeInfo and spikeData need row names!")))
  }
  if(is.null(colnames(spikeInfo)) || !any(colnames(spikeInfo) %in% "SpikeInput")) {
    stop(message(paste0("The input data frame of spikeInfo has no column labelled SpikeInput.")))
  }
  if (!normalisation %in% c('depth', 'none')) {
    stop(message(paste0("For spike-in count data normalisation, only depth normalization or providing prenormalized data is implemented.")))
  }
  if (normalisation=='none' && all((spikeData - round(spikeData)) == 0)) {
    stop(message(paste0("Skipping the normalisation should only be done with pre-normalized data, eg. RSEM output, but the provided input matrix contains only integer values!")))
  }

  if(!is.null(rownames(spikeData)) && !is.null(rownames(spikeInfo))) {
    # kick out undetected spike-ins
    spikeData.red <- spikeData[rowSums(spikeData)>0, colSums(spikeData)>100]
    # sort them if needed
    spikeInfo.red <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData.red), , drop = FALSE]
    spikeInfo.red <- spikeInfo.red[match(rownames(spikeData.red), rownames(spikeInfo.red)), , drop = FALSE]
    # spikeInfo.red <- spikeInfo.red[,"SpikeInput", drop = FALSE]
    message(paste0(nrow(spikeInfo.red), " spike-ins have been detected."))
    if(nrow(spikeInfo.red)<10) {
      stop(message(paste0("Not enough spike-ins detected to allow reliable variance estimation.
                          Please proceed without spike-in estimation.
                          There is still the option to generate insilico spike-ins.")))
    }

    if(!is.null(batchData) && !is.null(rownames(batchData))) {
      # define batches
      batchData.red <- batchData[rownames(batchData) %in% colnames(spikeData.red), , drop = FALSE]
      batch <- stats::setNames(as.character(row.names(batchData.red)), batchData.red[,1])
      # calculate normalized spike-ins per batch
      normspikeData <- sapply(unique(names(batch)), function(b){
        tmpData <- spikeData.red[,grepl(pattern = b, names(batch))]
        normData <- .norm.calc(normalisation = normalisation,
                               countData = tmpData,
                               batchData = NULL,
                               spikeData = NULL,
                               spikeInfo = NULL,
                               Lengths = NULL,
                               MeanFragLengths = NULL,
                               NCores=NULL,
                               verbose=FALSE)
        normData
      }, simplify=FALSE)
      # estimate gamma theta per batch
      gammaThetaEstimate <- sapply(unique(names(batch)), function(b){
        tmpData <- normspikeData[[b]]$NormCounts
        tmpSF <- normspikeData[[b]]$size.factors
        if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLength)) {
          tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
            .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLength + 1) / 10^3), 1, ncol(tmpData))
        }
        if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLength)) {
          tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
            .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(tmpData))
        }
        if(is.null(spikeInfo.red$Lengths) && is.null(spikeInfo.red$MeanFragLength)) {
          tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1)
        }
        est <- .estimateGammaTheta(nCountSpikes = tmpData,
                                  numberSpikes = spikeInfo.red[,"SpikeInput", drop = FALSE],
                                  sizeFactorMatrix = tmpSFmat)
        est
      }, simplify=FALSE)
      # correct for batch effect by dividing norm counts with gamma theta
      normCountsL <- sapply(unique(names(batch)), function(b){
        tmpData <- normspikeData[[b]]$NormCounts
        tmpgammaTheta <- gammaThetaEstimate[[b]]$gammaTheta[[1]]
        tmpData / tmpgammaTheta
      }, simplify=FALSE)

      # return normalized, batch-corrected counts and size factors for spike-ins
      normCounts <- do.call('cbind', normCountsL)
      sizeFactor <- unlist(lapply(normspikeData,function(x) x$size.factors))
    }
    if(is.null(batchData)) {
      # calculate normalized spike-ins
      normspikeData <- .norm.calc(normalisation=normalisation,
                                  countData=spikeData.red,
                                  batchData = NULL,
                                  spikeData=NULL,
                                  spikeInfo=NULL,
                                  Lengths=NULL,
                                  MeanFragLengths=NULL,
                                  NCores=NULL)
      # return normalized, batch-corrected counts and size factors for spike-ins
      normCounts <- normspikeData$NormCounts
      sizeFactor <- normspikeData$size.factors
    }

    # technical noise fit
    if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLength)) {
      sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
        .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLength + 1) / 10^3), 1, ncol(normCounts))
    }
    if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLength)) {
      sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
        .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(normCounts))
    }
    if(is.null(spikeInfo.red$Lengths) && is.null(MeanFragLength)) {
      sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1)
    }
    EVGammaThetaEstimate <- .estimateEVGammaTheta(nCountSpikes = normCounts,
                                                 numberSpikes = spikeInfo.red[,"SpikeInput", drop = FALSE],
                                                 sizeFactorMatrix = sizeFactorMatrix)
  }

  normParams <- data.frame(Log2Means=rowMeans(log2(as.matrix(normCounts)+1)),
                           Log2Sd=matrixStats::rowSds(log2(as.matrix(normCounts)+1)),
                           row.names = row.names(normCounts))

  # return object
  res <- list("normCounts" = normCounts,
              "normParams" = normParams,
              "size.factors" = sizeFactor,
              "EVGammaThetaEstimates" = EVGammaThetaEstimate,
              "Input" = list("spikeData"=spikeData,
                             "spikeInfo"=spikeInfo,
                             "batchData"=batchData,
                             "MeanFragLength"=MeanFragLength),
              "Settings"=list("normFramework"=normalisation))
  return(res)
}

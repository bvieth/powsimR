
# insilicoNBParam ---------------------------------------------------------

#' @name insilicoNBParam
#' @aliases insilicoNBParam
#' @title In Silico Simulation Parameters
#' @description With this function, the user can define the parameters of the negative binomial distribution needed for the simulations. The mean, dispersion and dropout need to be specified for bulk RNAseq experiments. For single cell RNAseq experiments, the mean and dispersion should be chosen so that a high number of dropouts are generated (see details).
#' @usage insilicoNBParam(means, dispersion,
#' dropout=NULL, sf=NULL, RNAseq=c("bulk", "singlecell"))
#' @param means mean parameter of the NB read count distribution defined by a function sampling from a random distribution, e.g. function(x) rgamma(n=x, shape=2, rate=2).
#' @param dispersion dispersion parameter of the NB read count distribution. This can be:
#' (1) a constant, e.g. 0.2
#' (2) a function relating to the mean expression. The parametric fit employed in DESeq2 is recommended to get an estimate, e.g. function(x) 3 + 150/x.
#' @param dropout The probability of droput per gene. This can be:
#' (1) a constant
#' (2) a function sampling from a random distribution, e.g. function(x) runif(x, min=0, max=0.25).
#' The default is \code{NULL}, i.e. no dropout. Note that the dropout is not considered for single cell simulations. This should be inherently given by low mean expression coupled with high dispersion.
#' @param sf The size factor:
#' (1) a vector
#' (2) a function sampling from a random distribution, e.g. function(x) rnorm(x, mean=1, sd=0.2).
#' The default is \code{NULL}, i.e. equal size factor of 1.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @return List with the following vectors:
#' \item{means}{Mean log2 read counts per gene.}
#' \item{dispersion}{Log2 dispersion per gene.}
#' \item{p0}{Probability that the count will be zero per gene.}
#' \item{sf}{Size factor per sample.}
#' \item{estFramework}{"in silico" by default.}
#' @examples
#' \dontrun{
#' ## Bulk RNA-seq experiment in silico parameters:
#' insilico.bulk <- insilicoNBParam(means=function(x) 2^runif(x, 9, 12),
#' dispersion=0.2,
#' dropout=function(x) runif(x, min = 0, max = 0.25),
#' sf=function(x) truncnorm::rtruncnorm(x, a = 0.6, b = 1.4, mean=1, sd=0.2),
#' RNAseq='bulk')
#' ## Single cell RNA-seq experiment in silico parameters:
#' insilico.singlecell <- insilicoNBParam(means=function(x) rgamma(x, 4, 2),
#' dispersion=function(x) 2 + 150/x,
#' dropout=NULL,
#' sf=function(x) 2^rnorm(n=x, mean=0, sd=0.5),
#' RNAseq='singlecell')
#' }
#' @author Beate
#' @rdname insilicoNBParam
#' @export
insilicoNBParam <- function(means, dispersion, dropout=NULL, sf=NULL, RNAseq=c("bulk", "singlecell")) {
  if (is.function(means)) {
    means <- means
  }
  if (!is.function(means)) {stop("Please provide a function for the mean parameter!")}

  if (is.function(dispersion)) {
    dispersion <- dispersion
  }
  if (length(dispersion) == 1) {
    dispersion <- dispersion
  }
  # if (!length(dispersion) == 1 || !is.function(dispersion)) {message("Please provide either a function or constant value for dispersion parameter!")}

  if(RNAseq=='bulk') {
    if (is.function(dropout)) {
      dropout <- dropout
    }
    if (is.null(dropout)) {
      dropout <- NULL
    }
    # if (!is.function(dropout) || !is.function(dispersion)) {message("Please provide a function for the dropout parameter or leave it as default, i.e. NULL.")}
  }

  if(RNAseq=='singlecell') {
    dropout <- NULL
    message("The dropout is set to NULL for single cell experiments. For more information, please consult the details section of insilicoNBParam function.")
  }

  if (is.function(sf)) {
    sf <- sf
  }
  if (is.null(sf)) {
    sf <- NULL
  }
  if (is.vector(sf)) {
    sf <- sf
  }

  res <- list(means = means, dispersion = dispersion, p0 = dropout, sf = sf, RNAseq = RNAseq)
  attr(res, 'param.type') <- "insilico"
  attr(res, 'Distribution') <- "NB"

  return(res)
}

# estimateParam ---------------------------------------------------------

#' @name estimateParam
#' @aliases estimateParam
#' @title Estimate simulation parameters
#' @description This function estimates and returns parameters needed for power
#' simulations assuming a (zero-inflated) negative binomial read count distribution.
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
#' #' \item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
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
#' estparam <- estimateParam(countData=cnts, batchData = NULL,
#'                           spikeData=NULL, spikeInfo = NULL,
#'                           Lengths=NULL, MeanFragLengths=NULL,
#'                           Distribution='NB',
#'                           RNAseq="singlecell",
#'                           normalisation='scran',
#'                           NCores=NULL,
#'                           sigma=1.96, verbose=TRUE)
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
#' estparam <- estimateParam(countData=cnts, batchData = NULL,
#'                           spikeData=NULL, spikeInfo = NULL,
#'                           Lengths=NULL, MeanFragLengths=NULL,
#'                           Distribution='NB',
#'                           RNAseq="bulk",
#'                           normalisation='MR',
#'                           NCores=NULL,
#'                           sigma=1.96, verbose=TRUE)
#' plotParam(estparam, annot=F)
#'
#' # download count table
#' githubURL <- "https://github.com/bvieth/powsimRData/raw/master/data-raw/kolodziejczk_cnts.rda"
#' download.file(url = githubURL, destfile= "kolodziejczk_cnts.rda", method = "wget")
#' load('kolodziejczk_cnts.rda')
#' kolodziejczk_cnts <- kolodziejczk_cnts[, grep('standard', colnames(kolodziejczk_cnts))]
#' ## estimate NB parameters:
#' estparams <- estimateParam(countData=cnts
#'                           Distribution='NB',
#'                           RNAseq="singlecell",
#'                           normalisation='scran',
#'                           NCores=NULL,
#'                           sigma=1.96, verbose=TRUE)
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

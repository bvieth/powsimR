
# estimateParam ---------------------------------------------------------

#' @name estimateParam
#' @aliases estimateParam
#' @title Estimate simulation parameters
#' @description This function estimates and returns parameters needed for power simulations.\cr
#' The user needs to choose the following options at least: specify a gene expression matrix; the type of RNA-seq experiment, i.e. bulk or single cell; the recommended distribution is negative binomial (NB) except for single-cell full-length Smart-seq2 read data where we recommend zero-inflated NB (ZINB); the preferred normalisation method, we recommend scran for single cell and TMM or MR for bulk.\cr
#' The other parameters are optional (additional data) or have preset values (gene and sample filtering). Please consult the detailed arguments description.
#' @usage estimateParam(countData,
#' readData = NULL,
#' batchData = NULL,
#' spikeData = NULL,
#' spikeInfo = NULL,
#' Lengths = NULL,
#' MeanFragLengths = NULL,
#' RNAseq = c('bulk', 'singlecell'),
#' Protocol = c('UMI', 'Read'),
#' Distribution = c('NB', 'ZINB'),
#' Normalisation = c("TMM", "MR", "PosCounts", "UQ",
#' "scran", "Linnorm", "sctransform",
#' "SCnorm", "Census", "depth", "none"),
#' GeneFilter = 0.05,
#' SampleFilter = 5,
#' sigma = 1.96,
#' NCores = NULL,
#' verbose=TRUE)
#' @param countData is a (UMI) count \code{matrix} of gene expression.
#' Rows correspond to genes, columns to samples.
#' The gene names should be given as \code{rownames} without "_" in the names.
#' The samples names should be given as \code{colnames}.
#' The count matrix should only contain the expression of one group, e.g. wildtype / untreated control / one cell type population.
#' @param readData is a the matching read count \code{matrix} of gene expression if \code{countData} is UMI
#' and the same formatting should be applied. Default is \code{NULL} and
#' users only need to supply the read count matrix if they plan to apply downsampling of UMI counts for simulations, see \code{\link{Setup}}.
#' @param batchData is a \code{data.frame} for batch annotation.
#' Rows correspond to samples. The first column should contain the batches, e.g. 'a', 'b', 'c', etc.
#' @param spikeData is a count \code{matrix}.
#' Rows correspond to spike-ins, columns to samples.
#' The order of columns should be the same as in the \code{countData}.
#' This is only needed for spike-in aware normalisation methods ('MR', 'Linnorm', 'scran', 'SCnorm', 'bayNorm', 'Census'), see Details.
#' @param spikeInfo is a molecule count \code{matrix} of spike-ins.
#' Rows correspond to spike-ins. The order of rows should be the same as in the \code{spikeData}.
#' The column names should be 'SpikeID' and 'SpikeInput' for molecule counts of spike-ins.
#' This is only needed for spike-in aware normalisation methods (), see Details.
#' @param Lengths is a numeric vector of transcript lengths with the same length and order as the rows in countData.
#' This variable is only needed for internal gene length corrections (TPM), see Details.
#' @param MeanFragLengths is a numeric vector of mean fragment lengths with the same length as columns in countData.
#' This variable is only needed for internal gene length corrections (TPM), see Details.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param Protocol is a character value defining the type of counts given in \code{countData}.
#' Options are "UMI" (e.g. 10X Genomics, CEL-seq2) or "Read" (e.g. Smart-seq2).
#' @param Distribution is a character value: "NB" for negative binomial or "ZINB" for zero-inflated negative binomial distribution fitting.
#' @param Normalisation is a character value: 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm',
#' 'SCnorm', 'Census', 'depth', 'none'.
#' For more information, please consult the Details section.
#' @param GeneFilter is a numeric vector indicating the minimal proportion of nonzero expression values
#' for a gene across all samples to be considered expressed and used for normalisation and parameter estimation.
#' The default is \code{0.05}, i.e. at least 5\% of the expression values per gene need to be nonzero.
#' @param SampleFilter is a numeric vector indicating the minimal number of MADs (median absolute deviation)
#' away from the median number of features detected as well as sequencing depth across all samples
#' so that outlying samples are removed prior to normalisation and parameter estimation.
#' The default is \code{5}, i.e. at least 5 MADs away from the median.
#' Choose higher values if you want to filter out less samples.
#' This parameter is particularly important for single cells to ensure reliable parameter estimation.
#' For more information, please consult \code{\link[scater]{isOutlier}}.
#' @param sigma The variability band width for mean-dispersion loess fit defining the prediction interval for read count simulation. Default is 1.96, i.e. 95\% interval. For more information see \code{\link[msir]{loess.sd}}.
#' @param NCores The number of cores for normalisation method SCnorm and Census.
#' The default \code{NULL} means 1 core.
#' @param verbose Logical value to indicate whether to print function information.
#' Default is \code{TRUE}.
#'
#' @return List object with the following entries:
#' \item{Parameters}{A list object containing the estimated moments for the full, dropped out genes, dropped out samples and filtered normalized count matrix. For more information please consult the details section and the plot made with \code{\link{plotParam}}.}
#' \item{Fit}{A list object containing the fitting results of the mean-dispersion and mean-dropout relation as well as the estimated parameter data used for the fits. For more information please consult the details section and the plot made with \code{\link{plotParam}}.}
#' \item{totalS,totalG}{Number of samples and genes provided with at least one read count.}
#' \item{DropOuts}{List object containing logical vectors for gene and sample dropouts after applying gene frequency and sample outlier filtering.}
#' \item{sf}{The estimated library size factor per sample.}
#' \item{Lengths -..- SampleFilter}{The chosen parameters settings.}
#'
#' @details
#' \strong{Normalisation Methods}
#' \describe{
#' \item{TMM, UQ}{employ the edgeR style normalization of weighted trimmed mean of M-values and upperquartile
#' as implemented in \code{\link[edgeR]{calcNormFactors}}, respectively.}
#' \item{MR, PosCounts}{employ the DESeq2 style normalization of median ratio method and a modified geometric mean method
#' as implemented in \code{\link[DESeq2]{estimateSizeFactors}}, respectively.}
#' \item{scran, SCnorm}{apply the deconvolution and quantile regression normalization methods developed for sparse RNA-seq data
#' as implemented in \code{\link[scran]{calculateSumFactors}} and \code{\link[SCnorm]{SCnorm}}, respectively. Spike-ins can also be supplied for both methods via \code{spikeData}. Note, however that this means for scran that the normalisation as implemented in \code{\link[scuttle]{computeSpikeFactors}} is also applied to genes (\code{general.use=TRUE}).}
#' \item{Linnorm}{apply the normalization method for sparse RNA-seq data
#' as implemented in \code{\link[Linnorm]{Linnorm.Norm}}.
#' For \code{Linnorm}, the user can also supply \code{spikeData}.}
#' \item{sctransform}{apply the normalization method developed for single-cell
#' UMI RNA-seq data as implemented in \code{\link[sctransform]{vst}}. }
#' \item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
#' Census at least needs \code{Lengths} for single-end data and preferably \code{MeanFragLengths} for paired-end data.
#' The authors state that Census should not be used for UMI data.}
#' \item{depth}{Sequencing depth normalisation.}
#' \item{none}{No normalisation is applied. This approach can be used for prenormalized expression estimates, e.g. cufflinks, RSEM or salmon.}
#' }
#'
#' @examples
#' \dontrun{
#' # Single Cells
#' data("SmartSeq2_Gene_Read_Counts")
#' Batches <- data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
#'                       stringsAsFactors = FALSE, row.names = colnames(SmartSeq2_Gene_Read_Counts))
#' data("GeneLengths_mm10")
#' estparam <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
#'                           readData = NULL,
#'                           batchData = Batches,
#'                           spikeData = SmartSeq2_SpikeIns_Read_Counts,
#'                           spikeInfo = SmartSeq2_SpikeInfo,
#'                           Lengths = GeneLengths, MeanFragLengths = NULL,
#'                           RNAseq = 'singlecell', Protocol = 'Read',
#'                           Distribution = 'ZINB', Normalisation = "scran",
#'                           GeneFilter = 0.1, SampleFilter = 3,
#'                           sigma = 1.96, NCores = NULL, verbose = TRUE)
#'
#' # Bulk
#' data("Bulk_Read_Counts")
#' data("GeneLengths_hg19")
#' estparam <- estimateParam(countData = Bulk_Read_Counts,
#'                           readData = NULL,
#'                           batchData = NULL,
#'                           spikeData = NULL,
#'                           spikeInfo = NULL,
#'                          Lengths = GeneLengths_hg19,
#'                           MeanFragLengths = NULL,
#'                           RNAseq = 'bulk', Protocol = 'Read',
#'                           Distribution = 'NB', Normalisation = "MR",
#'                           GeneFilter = 0.1, SampleFilter = 3,
#'                           sigma = 1.96, NCores = NULL, verbose = TRUE)
#'
#' # plot the results of estimation
#' plotParam(estparam, Annot = FALSE)
#' }
#' @author Beate Vieth
#' @rdname estimateParam
#' @export
estimateParam <- function(countData,
                          readData = NULL,
                          batchData = NULL,
                          spikeData = NULL,
                          spikeInfo = NULL,
                          Lengths = NULL,
                          MeanFragLengths = NULL,
                          RNAseq = c('bulk', 'singlecell'),
                          Protocol = c('UMI', 'Read'),
                          Distribution = c('NB', 'ZINB'),
                          Normalisation = c("TMM", "MR", "PosCounts", "UQ",
                                            "scran", "Linnorm", "sctransform",
                                            "SCnorm", "Census", "depth", "none"),
                          GeneFilter = 0.05,
                          SampleFilter = 5,
                          sigma = 1.96,
                          NCores = NULL,
                          verbose = TRUE) {

  # Inform users of inappropiate choices
  if (RNAseq == 'bulk' && Distribution=='ZINB') {
    if(verbose) {message("Zero-inflated negative binomial is not implemented for bulk data. Setting it to negative binomial.")}
    Distribution="NB"
  }
  if (RNAseq=='singlecell' && Normalisation=='MR') {
    if(verbose) {message(paste0(Normalisation, " has been developed for bulk data and it is rather likely that it will not work.
                                \nPlease consider using a method that can handle single cell data, e.g. PosCounts, scran, SCnorm."))}
  }
  if (Normalisation=='Census' && Protocol == "UMI") {
    if(verbose) {message(paste0(Normalisation, " should only be used for non-UMI methods! \nFor more information, please consult the monocle vignette."))}
    if (Normalisation=='Census' && is.null(Lengths) && Protocol == "Read") {
      if(verbose) {message(paste0(Normalisation, " should be used in combination with transcript lengths.
                                   \nIf the library is paired-end, please also provide the mean fragment lengths which can be determined by e.g. Picard."))}
    }
  }

  # STOP if combination of options is not supported/would not work!
  if (RNAseq=='bulk' && Normalisation == 'Census') {
    stop(message(paste0(Normalisation, " is only developed and implemented for single cell RNA-seq experiments.")))
  }

  if (Normalisation=='none' && all((countData - round(countData)) == 0)) {
    stop(message(paste0("No normalisation should only be done with pre-normalized data, eg. RSEM output, but the provided input matrix contains only integer values!")))
  }

  if (Normalisation %in% c('TMM', 'UQ', 'depth', 'none') && !is.null(spikeData)) {
    message(paste0("The normalisation method ", Normalisation, " does not utilize spike-ins."))
    spikeData = NULL
    spikeInfo = NULL
  }

  if (!is.null(readData) && Protocol == "Read") {
    stop(message(paste0("If the matrix contains read-based derived expression, then there is no need to provide an extra readData object besides the countData.")))
  }

  # check the matching of names and objects
  checkup <- .run.checkup(
    countData = countData,
    readData = readData,
    batchData = batchData,
    spikeData = spikeData,
    spikeInfo = spikeInfo,
    Lengths = Lengths,
    MeanFragLengths=MeanFragLengths,
    RNAseq = RNAseq,
    verbose=verbose
  )

  # extract checked objects
  countData <- checkup$countData
  readData <- checkup$readData
  batchData <- checkup$batchData
  spikeData <- checkup$spikeData
  spikeInfo <- checkup$spikeInfo
  Lengths <- checkup$Lengths
  MeanFragLengths <- checkup$MeanFragLengths
  Label <- ifelse(is.null(batchData), "none", "known")

  # run estimation
  res = .run.estParam(countData = countData,
                      readData = readData,
                      batchData = batchData,
                      spikeData = spikeData,
                      spikeInfo = spikeInfo,
                      Lengths = Lengths,
                      MeanFragLengths = MeanFragLengths,
                      Distribution = Distribution,
                      RNAseq = RNAseq,
                      Protocol = Protocol,
                      Normalisation = Normalisation,
                      Label = Label,
                      GeneFilter = GeneFilter,
                      SampleFilter = SampleFilter,
                      sigma = sigma,
                      NCores = NCores,
                      verbose=verbose)

  # inform users
  if (RNAseq=='bulk' && c(res$Parameters$Filtered$grand.dropout>0.5 || is.null(res$Fit$Filtered$g0.cut) ))
    if(verbose) {message('The provided bulk data is very sparse.
             Consider using single cell estimation framework.')}

  # make output
  res2 <- c(res, list(normFramework = Normalisation,
                      sigma = sigma,
                      GeneFilter = GeneFilter,
                      SampleFilter = SampleFilter))
  attr(res2, 'RNAseq') <- RNAseq
  attr(res2, 'Distribution') <- Distribution
  attr(res2, 'Protocol') <- Protocol

  # return output
  return(res2)
}


# estimateSpike -----------------------------------------------------------

#' @name estimateSpike
#' @aliases estimateSpike
#' @title Estimate simulation parameters for spike-ins
#' @description This function estimates and returns parameters needed for spike-in count simulations using supplementary code from Kim et al. 2016 (DOI: 10.1038/ncomms9687).
#' @usage estimateSpike(spikeData,
#' spikeInfo,
#' MeanFragLengths = NULL,
#' batchData = NULL,
#' Normalisation=c('depth','none'),
#' SampleFilter = 3,
#' RNAseq = c("bulk", "singlecell"),
#' Protocol = c('UMI', 'Read'),
#' verbose = TRUE)
#' @param spikeData  is a count \code{matrix}. Rows correspond to spike-ins, columns to samples.
#' \code{rownames} should contain the spike-in names, \code{colnames} the sample names.
#' @param spikeInfo is a molecule count \code{matrix} of spike-ins. Rows correspond to spike-ins.
#' The order of rows should be the same as the rows in \code{spikeData}.
#' The rownames should be the same as the rownames of spikeData.
#' The first column must contain the molecule counts of spike-ins and be named 'SpikeInput'.
#' The second column can contain the sequence lengths of spike-ins in bases and be named 'Lengths'.
#' @param MeanFragLengths is a numeric vector of the mean fragment length.
#' @param batchData is a \code{data.frame} for batch annotation. Rows correspond to samples.
#' The order of rows should be the same as the columns in \code{spikeData}.
#' The rownames should be the same as the column names of spikeData.
#' The first column should contain the batch annotation, e.g. 'a' for batch 1, 'b' for batch 2.
#' @param Normalisation is a character value: 'depth' or 'none'. For more information, please consult the details section.
#' @param SampleFilter is a numeric vector indicating the minimal number of MADs (median absolute deviation)
#' away from the median number of features detected as well as sequencing depth across all samples
#' so that outlying samples are removed prior to normalisation and parameter estimation.
#' The default is \code{5}, i.e. at least 5 MADs away from the median.
#' Choose higher values if you want to filter out less samples.
#' This parameter is particularly important for single cells to ensure reliable parameter estimation.
#' For more information, please consult \code{\link[scater]{isOutlier}}.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param Protocol is a character value defining the type of counts given in \code{countData}.
#' Options are "UMI" (e.g. 10X Genomics, CEL-seq2) or "Read" (e.g. Smart-seq2).
#' @param verbose Logical value to indicate whether to print function information.
#' Default is \code{TRUE}.
#' @return List object with the following entries:
#' \item{normCounts}{The normalised spike-in read counts \code{data.frame}.}
#' \item{normParams}{The mean and standard deviation per spike-in using normalised read counts in a \code{data.frame}.}
#' \item{CaptureEfficiency}{The ad-hoc estimated as well as fitted detection probabilities with confidence interval per spike-in using normalised read counts in a \code{data.frame}.}
#' \item{seqDepth}{Library size, i.e. total number of reads per library}
#' \item{sf}{The estimated library size factors.}
#' \item{EVGammaThetaEstimates}{Estimation of the four parameters capturing technical variability,
#' namely E[\eqn{\gamma}], Var[\eqn{\gamma}], E[\eqn{\theta}] and Var[\eqn{\theta}].
#' For more details, please consult supplementary information of Kim et al. 2016 (DOI: 10.1038/ncomms9687).
#' These estimates are needed for simulating spike-in read counts.}
#' \item{Input}{The input spike-in expression \code{matrix}, molecule counts \code{data.frame} and batch annotation \code{data.frame},
#' filtered so that only spike-ins with nonzero expression and samples with at least 100 reads are retained.}
#' \item{Settings}{Reporting the chosen normalisation framework.}
#' @examples
#' \dontrun{
#' data("SmartSeq2_SpikeIns_Read_Counts")
#' data("SmartSeq2_SpikeInfo")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_SpikeIns_Read_Counts),
#'                                     "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SmartSeq2_SpikeIns_Read_Counts))
#' # estimation
#' spike_param <- estimateSpike(spikeData = SmartSeq2_SpikeIns_Read_Counts,
#' spikeInfo = SmartSeq2_SpikeInfo,
#' MeanFragLength = NULL,
#' batchData = Batches,
#' Normalisation = 'depth')
#' # plotting
#' plotSpike(estSpike = spike_param, Annot = FALSE)
#' }
#' @details
#' Normalisation methods
#' \describe{
#' \item{'depth'}{applies the depth normalization method as implemented in \code{\link[scuttle]{computeSpikeFactors}}.}
#' \item{'none'}{No normalisation is applied. This approach can be used for prenormalized expression estimates, e.g. TPM/FPKM/RPKM estimated by RSEM, salmon, cufflinks etc.}
#' }
#' @author Beate Vieth
#' @rdname estimateSpike
#' @importFrom Hmisc binconf
#' @importFrom matrixStats rowSds
#' @export
estimateSpike <- function(spikeData,
                          spikeInfo,
                          MeanFragLengths = NULL,
                          batchData = NULL,
                          Normalisation=c('depth','none'),
                          SampleFilter = 3,
                          RNAseq = c('bulk', 'singlecell'),
                          Protocol = c('UMI', 'Read'),
                          verbose = TRUE) {
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
  if (!Normalisation %in% c('depth', 'none')) {
    stop(message(paste0("For spike-in count data normalisation, only depth normalization or providing prenormalized data is implemented.")))
  }
  if (Normalisation=='none' && all((spikeData - round(spikeData)) == 0)) {
    stop(message(paste0("Skipping the normalisation should only be done with pre-normalized data, eg. RSEM output, but the provided input matrix contains only integer values!")))
  }


  # clean out extreme samples and spike-ins prior to normalisation
  if(!is.null(batchData)){
    bData <- as.vector(batchData[, 1])
  }
  if(is.null(batchData)){
    bData <- NULL
  }

  # define outliers as determined by SampleFilter using sequencing depth, detected features
  totCounts <- colSums(spikeData)
  libsize.drop <- scater::isOutlier(totCounts, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)
  totFeatures <- colSums(spikeData>0)
  feature.drop <- scater::isOutlier(totFeatures, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)
  minexpr.drop <- colSums(spikeData)<100
  # kick out spike-ins with no expression values
  spike.keep <- rowSums(spikeData )>0

  # remove outliers from spikeData
  spikeData.red <- spikeData[spike.keep, c(!libsize.drop | !feature.drop | !minexpr.drop)]

  # match spike Info to reduced spikeData
  spikeInfo.red <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData.red),]
  spikeInfo.red <- spikeInfo.red[match(rownames(spikeData.red), rownames(spikeInfo.red)),]

  # Mean fragment lengths
  if(!is.null(MeanFragLengths)) {
    cell.id <- colnames(spikeData.red)
    MeanFragLengths.red <- MeanFragLengths[match(cell.id,names(MeanFragLengths))]
  }
  if(is.null(MeanFragLengths)) {
    MeanFragLengths.red <- NULL
  }

  # batches
  # define batches
  if(!is.null(batchData)){
    batchData.red <- batchData[rownames(batchData) %in% colnames(spikeData.red), , drop = FALSE]
  }
  if(is.null(batchData)){
    batchData.red <- NULL
  }

  if(verbose){
    message(paste0(nrow(spikeData.red), " spike-ins have been detected in ",
                   ncol(spikeData.red), " cells."))
  }

  if(nrow(spikeInfo.red)<10) {
    stop(message(paste0("Not enough spike-ins detected to allow reliable variance estimation.
                          Please proceed without spike-in estimation.")))
  }

  if(!is.null(batchData) && !is.null(rownames(batchData))) {
    batch <- stats::setNames(as.character(row.names(batchData.red)), batchData.red[,1])
    # calculate normalized spike-ins per batch
    normspikeData <- sapply(unique(names(batch)), function(b){
      tmpData <- spikeData.red[,grepl(pattern = b, names(batch))]
      normData <- .norm.calc(Normalisation = Normalisation,
                             sf = NULL,
                             countData = tmpData,
                             spikeData = NULL,
                             spikeInfo = NULL,
                             batchData = NULL,
                             Lengths = NULL,
                             MeanFragLengths = NULL,
                             PreclustNumber = NULL,
                             Label = 'none',
                             Step = 'Estimation',
                             Protocol = Protocol,
                             NCores = NULL,
                             verbose = verbose)
      normData
    }, simplify=FALSE)

    # estimate gamma theta per batch
    gammaThetaEstimate <- sapply(unique(names(batch)), function(b){
      tmpData <- normspikeData[[b]]$NormCounts
      tmpSF <- normspikeData[[b]]$size.factors
      if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLengths.red)) {
        tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
          .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLengths.red + 1) / 10^3), 1, ncol(tmpData))
      }
      if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
        tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
          .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(tmpData))
      }
      if(is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
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
    sizeFactor <- as.vector(unlist(sapply(normspikeData, function(x) x$size.factors)))
    names(sizeFactor) <- colnames(spikeData.red)
    seqDepth <- colSums(spikeData.red)
    names(seqDepth) <- colnames(spikeData.red)
  }

  if(is.null(batchData)) {
    # calculate normalized spike-ins
    normspikeData <- .norm.calc(Normalisation = Normalisation,
                                sf = NULL,
                                countData = spikeData.red,
                                spikeData = NULL,
                                spikeInfo = NULL,
                                batchData = NULL,
                                Lengths = NULL,
                                MeanFragLengths = NULL,
                                PreclustNumber = NULL,
                                Label = 'none',
                                Step = 'Estimation',
                                Protocol = Protocol,
                                NCores = NULL,
                                verbose = verbose)
    # return normalized, batch-corrected counts and size factors for spike-ins
    normCounts <- normspikeData$NormCounts
    sizeFactor <- normspikeData$size.factors
    seqDepth <- colSums(spikeData.red)
    names(seqDepth) <- colnames(spikeData.red)
  }

  # technical noise fit
  if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
      .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLengths.red + 1) / 10^3), 1, ncol(normCounts))
  }
  if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
      .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(normCounts))
  }
  if(is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1)
  }
  EVGammaThetaEstimate <- .estimateEVGammaTheta(nCountSpikes = normCounts,
                                                numberSpikes = spikeInfo.red[,"SpikeInput", drop = FALSE],
                                                sizeFactorMatrix = sizeFactorMatrix)

  ## gene capture efficiency
  # sum of complete failures (Count=0) and successes= total- failures over all libraries
  spikecapture.dat <- data.frame(detect_fail=apply(spikeData.red==0,1,sum),
                                 detect_success=apply(spikeData.red>0,1,sum),
                                 row.names = rownames(spikeData.red))
  noobs <- ncol(spikeData.red)
  spikecapture.dat[,"p_fail"] <- spikecapture.dat["detect_fail"]/noobs
  spikecapture.dat[,"p_success"] <- spikecapture.dat["detect_success"]/noobs

  # confidence intervals of MLE p (applying Wilson interval which is a score based method as normal approximation can be biased
  # if n*p is not sufficiently large which is the case for some of the spike ins (very low p))
  spikecapture.dat [, "hat_p_success_cilower"] <- Hmisc::binconf(spikecapture.dat$detect_success, noobs, method="wilson")[,2]
  spikecapture.dat [, "hat_p_success_ciupper"] <- Hmisc::binconf(spikecapture.dat$detect_success, noobs, method="wilson")[,3]
  spikecapture.dat [, "hat_p_fail_cilower"] <- Hmisc::binconf(spikecapture.dat$detect_fail, noobs, method="wilson")[,2]
  spikecapture.dat [, "hat_p_fail_ciupper"] <- Hmisc::binconf(spikecapture.dat$detect_fail, noobs, method="wilson")[,3]

  ## cell capture efficiency
  total_ercc_molecules <- sum(spikeInfo.red$SpikeInput)
  cellcapture.dat <- apply(spikeData.red, 2, function(i) {
    sum(i) / sum(spikeInfo.red$SpikeInput)
  })

  ## estimated moments of normalized counts
  normParams <- data.frame(Log2Means=rowMeans(log2(as.matrix(normCounts)+1)),
                           Log2Sd=matrixStats::rowSds(log2(as.matrix(normCounts)+1)),
                           row.names = row.names(normCounts))

  # return object
  res <- list("normCounts" = normCounts,
              "normParams" = normParams,
              "CaptureEfficiency"=list("Spike-In"=spikecapture.dat, "Cell" =cellcapture.dat),
              "size.factors" = sizeFactor,
              "seqDepth" = seqDepth,
              "EVGammaThetaEstimates" = EVGammaThetaEstimate,
              "FilteredInput" = list("spikeData"=spikeData.red,
                                     "spikeInfo"=spikeInfo.red,
                                     "batchData"=batchData.red,
                                     "MeanFragLengths"=MeanFragLengths.red),
              "Settings"=list("normFramework"=Normalisation))
  attr(res, 'Protocol') <- Protocol
  attr(res, 'RNAseq') <- RNAseq

  return(res)
}



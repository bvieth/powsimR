# checkup -----------------------------------------------------------------

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scuttle isOutlier
.run.checkup <- function(countData,
                         readData,
                         batchData,
                         spikeData,
                         spikeInfo,
                         Lengths,
                         MeanFragLengths,
                         RNAseq,
                         verbose){
  if (!is.null(batchData) && !is.null(countData)){
    if(!ncol(countData) == nrow(batchData)) {
      stop(message(paste0("The batch information data frame and the count matrix have not the same number of samples.")))
    }
    if(ncol(countData) == nrow(batchData)){
      if (!is.null(batchData) && is.null(rownames(batchData)) && is.null(colnames(countData)) ) {
        rownames(batchData) <- paste0("S", 1:nrow(batchData))
        colnames(countData) <- paste0("S", 1:ncol(countData))
        message(paste0("No samples names were provided so that pseudo sample names are assigned."))
      }
      if (!is.null(batchData) && is.null(rownames(batchData)) && !is.null(colnames(countData)) ) {
        rownames(batchData) <- colnames(countData)
        message(paste0("No samples names were provided for batch information but countData has sample names, so that these names are assigned."))
      }
      if (!is.null(batchData) && !is.null(rownames(batchData)) && !is.null(colnames(countData)) ) {
        if (!all(rownames(batchData) %in% colnames(countData))){
          stop(message(paste0("The batch information data frame and the count matrix have not the same sample names.")))
        }
      }
    }
  }

  if (!is.null(Lengths) && is.null(names(Lengths)) && is.null(rownames(countData)) ) {
    stop(message(paste0("The gene lengths vector and the count data have no names so that correct matching is not possible. Please provide names in both input objects.")))
  }
  if (!is.null(Lengths) && !is.null(names(Lengths)) && !is.null(rownames(countData)) ) {
    if( any(duplicated(names(Lengths)), duplicated(rownames(countData))) ) {
      stop(message(paste0("The gene lengths vector and/or the count data have duplicated gene ID entries. Please provide objects with unique names only.")))
    }
    if (length(setdiff(rownames(countData), names(Lengths))) > 0) {
      stop(message(paste0("The gene lengths vector and the count data have nonoverlapping gene ID entries. Please provide matching objects.")))
    }
    if (length(setdiff(rownames(countData), names(Lengths))) == 0) {
      # match and sort Lengths based on countData
      Lengths <- Lengths[names(Lengths) %in% rownames(countData)]
      Lengths <- Lengths[match(rownames(countData), names(Lengths))]
      if (length(Lengths)==0) {
        stop(message(paste0("Please provide Lengths and countData with matching gene names!")))
      }
    }
  }

  if(!is.null(readData)){
    message(paste0("Checking countData and readData match."))
    if (is.null(rownames(readData)) ) {
      stop(message(paste0("Please provide rownames in readData input object for matching.")))
    }
    if (!is.null(rownames(readData)) && !is.null(rownames(countData)) ) {
      if (length(setdiff(rownames(countData), rownames(readData))) > 0) {
        stop(message(paste0("The read and count data have nonoverlapping gene ID entries. Please provide matching objects.")))
      }
      if (length(setdiff(rownames(countData), rownames(readData))) == 0) {
        # match and sort readData based on countData
        readData <- readData[rownames(readData) %in% rownames(countData), ]
        readData <- readData[match(rownames(countData), rownames(readData)), ]
        if (length(Lengths)==0) {
          stop(message(paste0("Please provide readData and countData with matching gene names!")))
        }
      }
    }
  }

  if (!is.null(rownames(countData)) && any(grepl(pattern = "_", rownames(countData))) ) {
    message(paste0("Some of the gene names in countData contain '_' which will interfere with simulations later on. They will be replaced with '--'."))
    rownames(countData) <- gsub(pattern = "_",
                                replacement = "--",
                                x = rownames(countData))
    names(Lengths) <- gsub(pattern = "_",
                           replacement = "--",
                           x = names(Lengths))
    if(!is.null(readData) && !is.null(rownames(readData)) && any(grepl(pattern = "_", rownames(readData))) ){
      message(paste0("Some of the gene names in readData contain '_' which will interfere with simulations later on. They will be replaced with '--'."))
      rownames(readData) <- gsub(pattern = "_",
                                  replacement = "--",
                                  x = rownames(readData))
    }
  }

  if (!is.null(MeanFragLengths) && is.null(names(MeanFragLengths)) && is.null(colnames(countData)) ) {
    stop(message(paste0("The mean fragment lengths vector and the count data samples have no names so that correct matching is not possible. Please provide names in both input objects.")))
  }
  if (!is.null(MeanFragLengths) && !is.null(names(MeanFragLengths)) && !is.null(colnames(countData)) ) {
    # match and sort mean fragment lengths based on countData
    MeanFragLengths <- MeanFragLengths[names(MeanFragLengths) %in% colnames(countData)]
    MeanFragLengths <- MeanFragLengths[match(colnames(countData), names(MeanFragLengths))]
    if (length(MeanFragLengths)==0) {
      stop(message(paste0("Please provide MeanFragLengths and countData with matching sample names!")))
    }
  }

  # stop if spikes and info do not match!
  if ( !is.null(spikeData) && is.null(rownames(spikeData)) && !is.null(spikeInfo) && is.null(rownames(spikeInfo)) ) {
    stop(message(paste0("No spike-in names were provided! Please provide spike molecule information only with matching rownames of spikeData and spikeInfo.")))
  }
  if ( !is.null(spikeData) && !is.null(rownames(spikeData)) && !is.null(spikeInfo) && !is.null(rownames(spikeInfo)) ) {
    # match and sort spikeData and spikeInfo
    spikeInfo <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData), , drop = FALSE]
    spikeInfo <- spikeInfo[match(rownames(spikeData), rownames(spikeInfo)), , drop = FALSE ]
    if (nrow(spikeInfo)==0) {
      stop(message(paste0("Please provide spike molecule information only with matching rownames of spikeData and spikeInfo.")))
    }
  }

  # check that countData and readData have the same dimensions
  if (!is.null(readData)) {
    if(!identical(dim(readData), dim(countData))) {
      message(paste0("The provided UMI and read count matrix have different dimensions."))
    }
  }

  # fill in pseudonames if missing and count data is the only input
  if(is.null(batchData) && is.null(spikeData) && is.null(spikeInfo) && is.null(Lengths)) {
    if (is.null(rownames(countData))) {
      rownames(countData) <- paste0("G", 1:nrow(countData))
      message(paste0("No gene names were provided so that pseudo gene names are assigned."))
    }
    if (is.null(colnames(countData))) {
      colnames(countData) <- paste0("S", 1:ncol(countData))
      message(paste0("No sample names were provided for counts so that pseudo sample names are assigned."))
    }
  }

  return(list(countData = countData,
              readData = readData,
              batchData = batchData,
              spikeData = spikeData,
              spikeInfo = spikeInfo,
              Lengths = Lengths,
              MeanFragLengths =MeanFragLengths))
}


# estParam ----------------------------------------------------------------

#' @importFrom parallel detectCores
#' @importFrom stats median
.run.estParam <- function(countData,
                          readData,
                          batchData,
                          spikeData,
                          spikeInfo,
                          Lengths,
                          MeanFragLengths,
                          Distribution,
                          RNAseq,
                          Protocol,
                          Normalisation,
                          Label,
                          GeneFilter,
                          SampleFilter,
                          NCores,
                          sigma,
                          verbose) {

  # name the sample unit appropiately
  samplename <- ifelse(RNAseq == "singlecell", "single cells", "bulk samples")

  # kick out empty samples and keep only expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  DetectS <- colSums(countData, na.rm = TRUE) > 0
  DetectG <- rowSums(countData, na.rm = TRUE) > 0
  detectG <- sum(DetectG)
  detectS <- sum(DetectS)

  countData <- countData[DetectG,DetectS]

  if(verbose) {
    message(paste0("The provided count matrix has ",
                   detectS, " out of ",
                   totalS," ", samplename, " and ",
                   detectG, " out of ",
                   totalG, " genes with at least 1 count."))
  }

  if(!is.null(Lengths)) {
    Lengths <- Lengths[DetectG]
  }
  if(!is.null(MeanFragLengths)) {
    MeanFragLengths <- MeanFragLengths[DetectS]
  }
  if(!is.null(batchData)) {
    batchData <- batchData[DetectS, , drop=F]
  }

  if(!is.null(spikeData) && !is.null(spikeInfo)) {
    # kick out undetected spike-ins
    spikeData <- spikeData[rowSums(spikeData)>0, DetectS]
    if(verbose) {message(paste0(nrow(spikeData), " spike-ins have been detected in ", ncol(spikeData), " ", samplename, "."))}
    # sort them if needed
    spikeInfo <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData), , drop = FALSE]
    spikeInfo <- spikeInfo[match(rownames(spikeData), rownames(spikeInfo)), , drop = FALSE ]
    if(nrow(spikeData)<10 || nrow(spikeInfo)<10) {
      stop(message(paste0("Not enough spike-ins detected to allow reliable normalisation. Please proceed with spike-in independent methods, e.g. MR, TMM, scran, SCnorm, etc.")))
      }
  }

  if(!is.null(spikeData) && is.null(spikeInfo)) {
    # kick out undetected spike-ins
    spikeData <- spikeData[rowSums(spikeData)>0, DetectS]
    if(verbose) { message(paste0(nrow(spikeData), " spike-ins have been detected in ", ncol(spikeData), " ", samplename, "."))}
    if(nrow(spikeData)<10) {
      stop(message(paste0("Not enough spike-ins detected to allow reliable normalusation.
                          Please proceed with spike-in indepdent methods, e.g. TMM, MR, scran, etc.")))
    }
  }

  # clean out extreme samples and genes prior to normalisation
  if(!is.null(batchData)){
    bData <- as.vector(batchData[, 1])
  }
  if(is.null(batchData)){
    bData <- NULL
  }
  # define outliers as determined by SampleFilter using sequencing depth, detected features and spike/gene count ratio
  totCounts <- colSums(countData)
  libsize.drop <- scuttle::isOutlier(totCounts, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)
  totFeatures <- colSums(countData>0)
  feature.drop <- scuttle::isOutlier(totFeatures, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)

  if(!is.null(spikeData)){
    totSpike <- colSums(spikeData)
    GeneSpikeRatio <- totSpike / (totSpike + totCounts) * 100
    genespike.drop <- scuttle::isOutlier(GeneSpikeRatio, nmads=SampleFilter, type="higher",
                                        log=FALSE, batch = bData)
  }
  if(is.null(spikeData)){
    genespike.drop <- rep(NA, length(totCounts))
    GeneSpikeRatio <- rep(NA, length(totCounts))
  }
  # kick out genes with no expression values as determined by GeneFilter
  gene.dropout <- rowSums(countData == 0) / ncol(countData)
  gene.drop <- gene.dropout >= 1 - GeneFilter

  # record the dropouts
  DropOuts <- list("Gene" = setNames(gene.drop, rownames(countData)),
                   "Sample" = data.frame("totCounts" = libsize.drop,
                                         "totFeatures" = feature.drop,
                                         "GeneSpike" = genespike.drop,
                                         "GeneSpikeRatio" = GeneSpikeRatio,
                                         row.names = colnames(countData)))

  sample.drop <- rowSums(DropOuts$Sample[,c(1:3)], na.rm = T) > 0L

  if(verbose) {
    message(paste0(sum(sample.drop), " out of ",
                   length(sample.drop), " ", samplename,
                   " were determined to be outliers and removed prior to normalisation."))
    message(paste0(sum(gene.drop), " genes out of ", length(gene.drop),
                   " were deemed unexpressed and removed prior to normalisation."))
  }

  # kick out features / cells that do not pass thresholds
  countData.Norm <- countData[!gene.drop, !sample.drop]

  # batch
  if(!is.null(batchData)) {
    batchData.Norm <- batchData[rownames(batchData) %in% colnames(countData.Norm), , drop = FALSE]
    batchData.Norm <- batchData.Norm[match(rownames(batchData.Norm), colnames(countData.Norm)), , drop = FALSE ]
  }
  if(is.null(batchData)) {
    batchData.Norm <- NULL
  }
  # Lengths
  if(!is.null(Lengths)) {
    gene.id <- rownames(countData.Norm)
    Lengths.Norm <- Lengths[match(gene.id,names(Lengths))]
  }
  if(is.null(Lengths)) {
    Lengths.Norm <- NULL
  }
  # Mean fragment lengths
  if(!is.null(MeanFragLengths)) {
    cell.id <- colnames(countData.Norm)
    MeanFragLengths.Norm <- MeanFragLengths[match(cell.id,names(MeanFragLengths))]
  }
  if(is.null(MeanFragLengths)) {
    MeanFragLengths.Norm <- NULL
  }
  # spike-in counts
  if(!is.null(spikeData)) {
    cell.id <- colnames(countData.Norm)
    spikeData.Norm <- spikeData[,match(cell.id, colnames(spikeData))]
  }
  if(is.null(spikeData)) {
    spikeData.Norm <- NULL
  }

  # normalisation
  NormData <- .norm.calc(Normalisation=Normalisation,
                         sf=NULL,
                         countData=countData.Norm,
                         spikeData=spikeData.Norm,
                         spikeInfo=spikeInfo,
                         batchData=batchData.Norm,
                         Lengths=Lengths.Norm,
                         MeanFragLengths=MeanFragLengths.Norm,
                         PreclustNumber=NULL,
                         Label=Label,
                         Step="Estimation",
                         Protocol=Protocol,
                         NCores=NCores,
                         verbose=verbose)

  # parameters: mean, dispersion, dropout
  if(verbose) {message('Estimating moments.')}

  # raw data
  RawData <- countData
  # normalize by sequencing depth and fill in scaling factors from normalisation
  raw.sf <- structure(colSums(RawData) / median(colSums(RawData)), names = colnames(RawData))
  raw.sf[names(NormData$size.factors)] <- NormData$size.factors
  NormRawData <- t(t(RawData)/raw.sf)
  raw.params <- .run.estparams(countData = RawData,
                                normData = NormRawData,
                                Distribution = Distribution)

  # all genes
  FullData <- countData[, !sample.drop]
  NormFullData <- t(t(FullData)/NormData$size.factors)
  full.params <- .run.estparams(countData = FullData,
                               normData = NormFullData,
                               Distribution = Distribution)
  # genes passing threshold
  FullFilterData <- countData[!gene.drop, !sample.drop]
  NormFilterData <- t(t(FullFilterData)/NormData$size.factors)
  sample.filter.params <- .run.estparams(countData = FullFilterData,
                                  normData = NormFilterData,
                                  Distribution = Distribution)

  # dropout data (genes)
  if(sum(gene.drop)/length(gene.drop) > 0.05) {
    GeneDropData <- countData[gene.drop, !sample.drop]
    NormGeneDropData <- t(t(GeneDropData)/NormData$size.factors)
    dropgene.params <- .run.estparams(countData = GeneDropData,
                                      normData = NormGeneDropData,
                                      Distribution = Distribution)
  }
  if(sum(DropOuts$Gene)/length(DropOuts$Gene) <= 0.05) {
    dropgene.params <- NA
  }
  # dropout data (samples)
  if(sum(sample.drop) > 3) {
    CellDropData <- countData[!DropOuts$Gene, sample.drop]
    sf <- colSums(CellDropData) / mean(colSums(CellDropData))
    NormCellDropData <- t(t(CellDropData)/sf)
    dropcell.params <- .run.estparams(countData = CellDropData,
                                      normData = NormCellDropData,
                                      Distribution = Distribution)
  }
  if(sum(sample.drop) <= 3) {
    dropcell.params <- NA
  }

  ParamData <- list('Raw' = raw.params,
                    'Full' = full.params,
                    'Filtered' = sample.filter.params,
                    'DropGene' = dropgene.params,
                    'DropSample' = dropcell.params)

  # fitting: mean vs dispersion, mean vs dropout
  if(verbose) {message('Fitting models.')}

  # full data
  RawFit <- .run.fits(ParamData = raw.params,
                      Distribution = Distribution,
                      RNAseq = RNAseq,
                      sigma = sigma)

  # full data
  FullFit <- .run.fits(ParamData = full.params,
                       Distribution = Distribution,
                       RNAseq = RNAseq,
                       sigma = sigma)

  # filtered data
  FilterFit <- .run.fits(ParamData = sample.filter.params,
                       Distribution = Distribution,
                       RNAseq = RNAseq,
                       sigma = sigma)

  if(sum(DropOuts$Gene)/length(DropOuts$Gene) > 0.01) {
  DropGeneFit <- .run.fits(ParamData = dropgene.params,
                           Distribution = Distribution,
                           RNAseq = RNAseq,
                           sigma = sigma)
  }
  if(sum(DropOuts$Gene)/length(DropOuts$Gene) <= 0.01) {
  DropGeneFit <- NA
  }

  # dropout data (samples)
  if(sum(sample.drop) > 3) {
  DropCellFit <- .run.fits(ParamData = dropcell.params,
                           Distribution = Distribution,
                           RNAseq = RNAseq,
                           sigma = sigma)
  }
  if(sum(sample.drop) <= 3) {
  DropCellFit <- NA
  }

  # read-umi fit
  if(!is.null(readData)) {
    UMIReadFit <- .run.fit.readumi(countData = countData,
                                   readData = readData)
  }
  if(is.null(readData)) {
    UMIReadFit <- NA
  }

  FitData <- list('Raw' = RawFit,
                  'Full' = FullFit,
                  'Filtered' = FilterFit,
                  'DropGene' = DropGeneFit,
                  'DropSample' = DropCellFit,
                  "UmiRead" = UMIReadFit)

  # inform user of final number of genes estimated and fitted for filtered data set.
  if(verbose){
    message(paste0("For ", FitData$Filtered$estG, " out of ", detectG,
                 " genes, mean, dispersion and dropout could be estimated. ",
                 FitData$Filtered$estS, " out of ", detectS, " ", samplename,
                 " were used for this."))
    }


  # result object
  res <- c(list(Parameters = ParamData),
           list(Fit = FitData),
           list(totalS = totalS,
                detectS = detectS,
                totalG = totalG,
                detectG = detectG,
                DropOuts = DropOuts,
                sf = NormData$size.factors,
                Lengths = Lengths,
                MeanFragLengths = MeanFragLengths))

  return(res)
}

# Moments Estimation ---------------------

# Estimate mean, dispersion, dropout from read counts

.run.estparams <- function(countData, normData, Distribution) {
  if(Distribution == "NB") {
    res <- .run.estparams.nb(countData, normData)
  }
  if(Distribution == "ZINB") {
    res <- .run.estparams.zinb(countData, normData)
  }
  return(res)
}

# NB
.run.estparams.nb <- function(countData, normData) {

  # sequencing depth, # features detected, # samples and # genes
  seqDepth = structure(colSums(countData), names = colnames(countData))
  nsamples = ncol(normData)
  ngenes = nrow(normData)
  totFeatures = colSums(countData>0)

  # indicator for zeroes ; dropout rates
  countData0 = countData == 0
  ng0 = rowSums(!countData0)
  nc0 = colSums(!countData0)
  g0 = (nsamples - ng0)/nsamples
  names(g0) = rownames(countData)
  c0 = (ngenes - nc0)/ngenes
  names(c0) = colnames(countData)
  grand0 = sum(countData0)/(nrow(countData)*ncol(countData))

  # calculate mean, dispersion and dropouts
  means = rowSums(normData)/nsamples
  names(means) = rownames(normData)
  s2 = rowSums((normData - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  names(size) = rownames(normData)
  dispersion = 1/size
  names(dispersion) = rownames(normData)
  common.dispersion = mean(dispersion, na.rm = TRUE)

  res = list(means = means,
             sds = s2,
             size = size,
             dispersion = dispersion,
             common.dispersion = common.dispersion,
             gene.dropout = g0,
             sample.dropout = c0,
             grand.dropout = grand0,
             seqDepth = seqDepth,
             totFeatures = totFeatures,
             nsamples = nsamples,
             ngenes = ngenes)
  attr(res, 'Distribution') <- "NB"
  return(res)
}

# ZINB
.run.estparams.zinb <- function(countData, normData) {

  # sequencing depth, # features detected, # samples and # genes
  seqDepth = structure(colSums(countData), names = colnames(countData))
  nsamples = ncol(normData)
  ngenes = nrow(normData)
  totFeatures = colSums(countData>0)

  # indicator for zeroes ; dropout rates
  countData0 = countData == 0
  ng0 = rowSums(!countData0)
  nc0 = colSums(!countData0)
  g0 = (nsamples - ng0)/nsamples
  names(g0) = rownames(countData)
  c0 = (ngenes - nc0)/ngenes
  names(c0) = colnames(countData)
  grand0 = sum(countData0)/(nrow(countData)*ncol(countData))

  # calculate mean, dispersion and dropouts including zeroes
  means = rowSums(normData)/nsamples
  names(means) = rownames(normData)
  s2 = rowSums((normData - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  names(size) = rownames(normData)
  dispersion = 1/size
  names(dispersion) = rownames(normData)
  common.dispersion = mean(dispersion, na.rm = TRUE)

  # calculate mean, dispersion and dropouts excluding zeroes
  pos.means = rowSums((!countData0)*normData)/ng0
  names(pos.means) = rownames(normData)
  pos.s2 = rowSums((!countData0)*(normData - pos.means)^2)/(ng0-1)
  pos.size = pos.means^2/(pos.s2 - pos.means + 1e-04)
  pos.size = ifelse(pos.size > 0, pos.size, NA)
  names(pos.size) = rownames(normData)
  pos.dispersion = 1/pos.size
  names(pos.dispersion) = rownames(normData)
  pos.common.dispersion = mean(pos.dispersion, na.rm = TRUE)

  res = list(means = means,
             pos.means = pos.means,
             sds = s2,
             pos.sds = pos.s2,
             size = size,
             pos.size = pos.size,
             dispersion = dispersion,
             pos.dispersion = pos.dispersion,
             common.dispersion = common.dispersion,
             pos.common.dispersion = pos.common.dispersion,
             gene.dropout = g0,
             sample.dropout = c0,
             grand.dropout = grand0,
             seqDepth = seqDepth,
             totFeatures = totFeatures,
             nsamples = nsamples,
             ngenes = ngenes)
  attr(res, 'Distribution') <- "ZINB"
  return(res)
}

.run.params <- function(countData, normData, group) {

  if (attr(normData, 'normFramework') %in% c('SCnorm', "Linnorm")) {
    # normalized count data
    norm.counts <- normData$NormCounts
  }

  if (!attr(normData, 'normFramework') %in% c('SCnorm', "Linnorm")) {
    # normalize count data
    sf <- normData$size.factors
    norm.counts <- t(t(countData)/sf)
  }

  nsamples = ncol(norm.counts)
  ngenes = nrow(norm.counts)

  # calculate mean, dispersion and dropouts (irrespective of group label)
  # means = rowMeans(norm.counts)
  # s2 = rowSums((norm.counts - means)^2)/(nsamples - 1)
  # size = means^2/(s2 - means + 1e-04)
  # size = ifelse(size > 0, size, NA)
  # dispersion = 1/size
  # counts0 = countData == 0
  # nn0 = rowSums(!counts0)
  # g0 = (nsamples - nn0)/nsamples

  # calculate mean, dispersion and dropouts (group-specific)
  norm.counts.red = norm.counts[, group == -1]
  nsamples.red = ncol(norm.counts.red)
  means = rowSums(norm.counts.red)/nsamples.red
  s2 = rowSums((norm.counts.red - means)^2)/(nsamples.red - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  g0 = (nsamples - nn0)/nsamples

  # calculate log fold changes
  lfc = .comp.FC(X = log2(norm.counts+1), L=group, is.log = T, FUN = mean)

  res = data.frame(geneIndex=rownames(countData),
                   means=means,
                   dispersion=dispersion,
                   dropout=g0,
                   lfcs=lfc,
                   stringsAsFactors = F)
  return(res)
}

.comp.FC <- function(X, L, is.log=TRUE, FUN=mean) {
    if (is.vector(X)) X <- matrix(X, byrow=TRUE)
    G1 <- X[, L == 1]
    G2 <- X[, L == -1]
    m1 <- apply(G1,1,FUN,na.rm=TRUE)
    m2 <- apply(G2,1,FUN,na.rm=TRUE)
    if (is.log) fc <- m1-m2
    else fc <- m1/m2
    return(fc)
}

# Fitting -----------------------------------------------------------------

#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @importFrom dplyr inner_join rename filter distinct
#' @importFrom magrittr %>%
.run.fit.readumi <- function(countData, readData){
  # combine umi and read count matrices, keep unique combis, log transform
  umi.dat <- reshape2::melt(as.matrix(countData))
  read.dat <- reshape2::melt(as.matrix(readData))
  vars <- c(UMI = "value.x", Read ="value.y",
            GeneID = "Var1", SampleID = 'Var2')
  UMI <- Read <- NULL
  suppressWarnings(suppressMessages(
  count.dat <- dplyr::inner_join(x = umi.dat, y = read.dat,
                    by = c(Var1 = "Var1", Var2 = "Var2")) %>%
    dplyr::rename(!!vars) %>%
    dplyr::filter(UMI > 0 & Read > 0) %>%
    dplyr::filter(UMI <= Read) %>%
    dplyr::distinct(UMI, Read, .keep_all = TRUE)
  ))
  lUMI <- log10(count.dat$UMI+1)
  lRead <- log10(count.dat$Read+1)
  lRatio <- lRead / lUMI
  keep <- complete.cases(lUMI, lRatio)
  lUMI <- lUMI[keep]
  lRatio <- lRatio[keep]
  lRead <- lRead[keep]
  ratioumi.fit <- msir::loess.sd(lRatio ~ lUMI, nsigma = 1)

  res <- list("lUMI" = lUMI,
              "lRead" = lRead,
              "lRatio" = lRatio,
              "Fit" = ratioumi.fit)

  # return object
  return(res)
}

.run.fits <- function(ParamData, sigma, Distribution, RNAseq) {
  if(Distribution == "NB") {
    res <- .run.fit.NB(ParamData, RNAseq, sigma)
  }
  if(Distribution == "ZINB") {
    res <- .run.fit.ZINB(ParamData, RNAseq, sigma)
  }
  return(res)
}

# NB Fit

#' @importFrom cobs cobs
#' @importFrom stats predict complete.cases
#' @importFrom msir loess.sd
.run.fit.NB <- function(ParamData, RNAseq, sigma) {

  mu = ParamData$means
  size = ParamData$size
  g0 = ParamData$gene.dropout
  disp =  ParamData$dispersion

  sharedgenes = Reduce(intersect, list(names(mu)[!is.na(mu)],
                                   names(g0)[!is.na(g0)],
                                   names(disp)[!is.na(disp)],
                                   names(size)[!is.na(size)]))

  mu = mu[sharedgenes]
  g0 = g0[sharedgenes]
  disp = disp[sharedgenes]
  size = size[sharedgenes]

  ldisp = log2(disp)
  lsize = log2(size)
  lmu = log2(mu+1)

  estG <- length(sharedgenes)
  estS <- ParamData$nsamples

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  if(RNAseq == "bulk"){
    cobs.fit <- cobs::cobs(x = lmu,
                           y = g0,
                           constraint = 'decrease',
                           nknots = 20,
                           print.warn = F, print.mesg = F)
    cobs.sim <- runif(1000,
                      min = min(lmu, na.rm=TRUE),
                      max = max(lmu, na.rm=TRUE))
    cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
    cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
    cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
    g0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])
  }

  # return object
  if(RNAseq == "singlecell") {
    res <- list(sharedgenes = sharedgenes,
                meansizefit = meansizefit,
                meandispfit = meandispfit,
                estS = estS,
                estG = estG)
  }

  if(RNAseq == "bulk") {
    res <- list(sharedgenes = sharedgenes,
                meansizefit = meansizefit,
                meandispfit = meandispfit,
                estS = estS,
                estG = estG,
                g0.cut = g0.cut)
  }

  return(res)
}

# ZINB fit

#' @importFrom cobs cobs
#' @importFrom stats predict complete.cases
#' @importFrom msir loess.sd
#' @importFrom matrixStats rowSds
.run.fit.ZINB <- function(ParamData, RNAseq, sigma){

  mu = ParamData$pos.means
  size = ParamData$pos.size
  g0 = ParamData$gene.dropout
  disp =  ParamData$pos.dispersion

  sharedgenes = Reduce(intersect, list(names(mu)[!is.na(mu)],
                                       names(g0)[!is.na(g0)],
                                       names(disp)[!is.na(disp)],
                                       names(size)[!is.na(size)]))

  mu = mu[sharedgenes]
  g0 = g0[sharedgenes]
  disp = disp[sharedgenes]
  size = size[sharedgenes]

  ldisp = log2(disp)
  lsize = log2(size)
  lmu = log2(mu+1)

  estG <- length(sharedgenes)
  estS <- ParamData$nsamples

  # define nonamplified and amplified genes
  mu.withzero <- ParamData$means[sharedgenes]
  sd.withzero <- ParamData$sds[sharedgenes]
  cv.withzero <- sd.withzero/mu.withzero
  nonamplified <- cv.withzero < 1.5*(sqrt(mu.withzero)/mu.withzero) & mu.withzero <=10

  # majority of genes are amplified
  if(sum(nonamplified, na.rm = T)/length(nonamplified) < 0.05){
    # mean-dispersion-size fit
    meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)
    meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

    # mean-dropout fit for all genes
    cobs.fit <- cobs::cobs(x = lmu,
                           y = g0,
                           constraint = 'decrease', nknots = 20,
                           print.warn = F, print.mesg = F)
    cobs.sim <- runif(1000,
                      min = min(lmu, na.rm=TRUE),
                      max = max(lmu, na.rm=TRUE))
    cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
    cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
    cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
    g0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

    if(length(g0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, g0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      meang0fit = loess.sd(x = all.dat$lmu,
                           y = all.dat$g0,
                           nsigma = sigma)
    }

    if(!length(g0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, g0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      all.dat <- all.dat[all.dat$lmu<g0.cut,]
      meang0fit = loess.sd(x = all.dat$lmu,
                           y = all.dat$g0,
                           nsigma = sigma)
    }

    meandispfit.amplified = meandispfit.nonamplified = meansizefit.amplified = meansizefit.nonamplified = NULL
    meang0fit.amplified = g0.cut.amplified = meang0fit.nonamplified = g0.cut.nonamplified = NULL
  }

  if(sum(nonamplified, na.rm = T)/length(nonamplified) >= 0.05){

    # split the parameters into amplified and nonamplified
    g0.amplified <- g0[!nonamplified]
    g0.nonamplified <- g0[nonamplified]

    lmu.amplified <- lmu[!nonamplified]
    lmu.nonamplified <- lmu[nonamplified]

    ldisp.amplified <- ldisp[!nonamplified]
    ldisp.nonamplified <- ldisp[nonamplified]

    lsize.amplified <- lsize[!nonamplified]
    lsize.nonamplified <- lsize[nonamplified]

    # mean-disp-size fits for amplified and nonamplified
    meandispfit.amplified =  msir::loess.sd(ldisp.amplified ~ lmu.amplified, nsigma = sigma)
    meandispfit.nonamplified =  msir::loess.sd(ldisp.nonamplified ~ lmu.nonamplified, nsigma = sigma)
    meansizefit.amplified =  msir::loess.sd(lsize.amplified ~ lmu.amplified, nsigma = sigma)
    meansizefit.nonamplified =  msir::loess.sd(lsize.nonamplified ~ lmu.nonamplified, nsigma = sigma)
    meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)
    meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

    # estimate the knee point for curves of mean versus dropout
    # amplified
    cobs.fit.amplified <- cobs::cobs(x = lmu.amplified,
                                     y = g0.amplified,
                                     constraint = 'decrease',
                                     nknots = 20,
                                     print.warn = F, print.mesg = F)
    cobs.sim.amplified <- runif(1000,
                                min = min(lmu.amplified, na.rm=T),
                                max = max(lmu.amplified, na.rm = T))
    cobs.predict.amplified <- as.data.frame(stats::predict(cobs.fit.amplified, cobs.sim.amplified))
    cobs.predict.amplified[,"fit"] <- ifelse(cobs.predict.amplified$fit < 0, 0, cobs.predict.amplified[,"fit"])
    cobs.predict.amplified <- cobs.predict.amplified[cobs.predict.amplified[,'fit'] < 0.05,]
    g0.cut.amplified <- as.numeric(cobs.predict.amplified[which.max(cobs.predict.amplified[,'fit']), 'z'])
    if(length(g0.cut.amplified)==0) {
      amplified.dat <- cbind.data.frame(lmu.amplified, g0.amplified)
      amplified.dat <- amplified.dat[stats::complete.cases(amplified.dat),]
      meang0fit.amplified = msir::loess.sd(x = amplified.dat$lmu.amplified,
                                           y = amplified.dat$g0.amplified,
                                           nsigma = sigma)
    }
    if(!length(g0.cut.amplified)==0) {
      amplified.dat <- cbind.data.frame(lmu.amplified, g0.amplified)
      amplified.dat <- amplified.dat[stats::complete.cases(amplified.dat),]
      amplified.dat <- amplified.dat[amplified.dat$lmu.amplified<g0.cut.amplified,]
      meang0fit.amplified = msir::loess.sd(x = amplified.dat$lmu.amplified,
                                           y = amplified.dat$g0.amplified,
                                           nsigma = sigma)
    }

    # nonamplified
    cobs.fit.nonamplified <- cobs::cobs(x = lmu.nonamplified,
                                        y = g0.nonamplified,
                                        constraint = 'decrease',
                                        nknots = 20,
                                        print.warn = F, print.mesg = F)
    cobs.sim.nonamplified <- runif(1000,
                                   min = min(lmu.nonamplified, na.rm=TRUE),
                                   max = max(lmu.nonamplified, na.rm=TRUE))
    cobs.predict.nonamplified <- as.data.frame(stats::predict(cobs.fit.nonamplified, cobs.sim.nonamplified))
    cobs.predict.nonamplified[,"fit"] <- ifelse(cobs.predict.nonamplified$fit < 0, 0, cobs.predict.nonamplified[,"fit"])
    cobs.predict.nonamplified <- cobs.predict.nonamplified[cobs.predict.nonamplified[,'fit'] < 0.05,]
    g0.cut.nonamplified <- log2(10)
    nonamplified.dat <- cbind.data.frame(lmu.nonamplified, g0.nonamplified)
    nonamplified.dat <- nonamplified.dat[stats::complete.cases(nonamplified.dat),]
    nonamplified.dat <- nonamplified.dat[nonamplified.dat$lmu.nonamplified<g0.cut.nonamplified,]

    if(nrow(nonamplified.dat)==0) {
      nonamplified.dat <- cbind.data.frame(lmu.nonamplified, g0.nonamplified)
      nonamplified.dat <- nonamplified.dat[stats::complete.cases(nonamplified.dat),]
      meang0fit.nonamplified = msir::loess.sd(x = nonamplified.dat$lmu.nonamplified,
                                              y = nonamplified.dat$g0.nonamplified,
                                              nsigma = sigma)
    }
    if(nrow(nonamplified.dat)>0) {
      meang0fit.nonamplified = msir::loess.sd(x = nonamplified.dat$lmu.nonamplified,
                                              y = nonamplified.dat$g0.nonamplified,
                                              nsigma = sigma)
    }

    # all
    cobs.fit <- cobs::cobs(x = lmu,
                           y = g0,
                           constraint = 'decrease', nknots = 20,
                           print.warn = F, print.mesg = F)
    cobs.sim <- runif(1000,
                      min = min(lmu, na.rm=TRUE),
                      max = max(lmu, na.rm=TRUE))
    cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
    cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
    cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
    g0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

    if(length(g0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, g0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      meang0fit = msir::loess.sd(x = all.dat$lmu,
                                 y = all.dat$g0,
                                 nsigma = sigma)
    }

    if(!length(g0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, g0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      all.dat <- all.dat[all.dat$lmu<g0.cut,]
      meang0fit = msir::loess.sd(x = all.dat$lmu,
                                 y = all.dat$g0,
                                 nsigma = sigma)
    }

  }

  nonamplified = length(which(nonamplified)) / length(nonamplified)

  # return object
  res <- list(sharedgenes = sharedgenes,
              meansizefit = meansizefit,
              meandispfit = meandispfit,
              meansizefit.amplified = meansizefit.amplified,
              meandispfit.amplified = meandispfit.amplified,
              meansizefit.nonamplified = meansizefit.nonamplified,
              meandispfit.nonamplified = meandispfit.nonamplified,
              meang0fit = meang0fit,
              g0.cut = g0.cut,
              meang0fit.amplified = meang0fit.amplified,
              g0.cut.amplified = g0.cut.amplified,
              meang0fit.nonamplified = meang0fit.nonamplified,
              g0.cut.nonamplified = g0.cut.nonamplified,
              nonamplified = nonamplified,
              estS = estS,
              estG = estG)
  return(res)
}

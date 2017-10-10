
# Estimate mean, dispersion, dropout from read counts ---------------------
#' @importFrom DEDS comp.FC
.run.params <- function(countData, normData, group) {

  # normalize count data
  sf <- normData$size.factors
  norm.counts <- t(t(countData)/sf)

  # calculate mean, dispersion and dropout
  means = rowMeans(norm.counts)
  nsamples = ncol(norm.counts)
  s2 = rowSums((norm.counts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  dropout = (nsamples - nn0)/nsamples

  #calculate log2 fold changes
  fc.foo = DEDS::comp.FC(L=group, is.log = F, FUN = mean)
  fc = fc.foo(norm.counts)
  lfc = log2(fc)

  res = data.frame(geneIndex=rownames(countData),
                    means=means,
                    dispersion=dispersion,
                    dropout=dropout,
                    lfcs=lfc,
                    stringsAsFactors = F)
  return(res)
}

# checkup -----------------------------------------------------------------

#' @importFrom scater newSCESet calculateQCMetrics isOutlier calcAverage
.run.checkup <- function(countData,
                         batchData,
                         spikeData,
                         spikeInfo,
                         Lengths,
                         MeanFragLengths,
                         RNAseq,
                         verbose) {

  if (!is.null(batchData) && is.null(rownames(batchData)) && is.null(colnames(countData)) ) {
    rownames(batchData) <- paste0("S", 1:nrow(batchData))
    message(paste0("No samples names were provided for batch information so that pseudo sample names are assigned."))
  }
  if (!is.null(batchData) && is.null(rownames(batchData)) && !is.null(colnames(countData)) ) {
    rownames(batchData) <- colnames(countData)
    message(paste0("No samples names were provided for batch information so that pseudo sample names are assigned."))
  }
  if (!is.null(Lengths) && is.null(names(Lengths)) && is.null(rownames(countData)) ) {
    stop(message(paste0("The gene lengths vector and the count data have no names so that correct matching is not possible. Please provide names in both input objects")))
  }
  if (!is.null(Lengths) && !is.null(names(Lengths)) && !is.null(rownames(countData)) ) {
    # match and sort Lengths based on countData
    Lengths <- Lengths[names(Lengths) %in% rownames(countData)]
    Lengths <- Lengths[match(rownames(countData), names(Lengths))]
    if (length(Lengths)==0) {
      stop(message(paste0("Please provide Lengths and countData with matching gene names!")))
    }
    if (!is.null(rownames(countData)) && any(grepl(pattern = "_", rownames(countData))) ) {
      message(paste0("Some of the gene names contain '_' which will interfere with simulations later on. They will be replaced with '--'."))
      rownames(countData) <- gsub(pattern = "_",
                                  replacement = "--",
                                  x = rownames(countData))
      names(Lengths) <- gsub(pattern = "_",
                             replacement = "--",
                             x = names(Lengths))
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

  if(RNAseq == 'singlecell') {
    # create SCEset
    if(is.null(batchData)) {
      sce <- scater::newSCESet(countData = countData,
                               lowerDetectionLimit = 0,
                               logExprsOffset = 1)
    }
    if(!is.null(batchData)) {
      pd <- new("AnnotatedDataFrame", data = batchData)
      sce <- scater::newSCESet(countData = countData,
                               phenoData = pd,
                               lowerDetectionLimit = 0,
                               logExprsOffset = 1)
    }
    # apply quality filters
    sce <- scater::calculateQCMetrics(sce, nmads = 3)
    # define outlier cells
    libsize.drop <- scater::isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
    feature.drop <- scater::isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
    # kick out owly expressed genes (average expression and dropout rate considered)
    ave.counts <- scater::calcAverage(sce)
    keep.ave <- ave.counts >= 0.2
    nsamples = ncol(countData)
    counts0 = countData == 0
    nn0 = rowSums(!counts0)
    p0 = (nsamples - nn0)/nsamples
    keep.p0 <- p0 < 0.9
    # kick out features / cells that do not pass thresholds
    sce <- sce[(keep.ave | keep.p0), !(libsize.drop | feature.drop)]

    # adapt the auxillary objects
    countData <- sce@assayData$counts

    # batch
    if(!is.null(batchData)) {
      batchData <- batchData[rownames(batchData) %in% colnames(countData), , drop = FALSE]
      batchData <- batchData[match(rownames(batchData), colnames(countData)), , drop = FALSE ]
    }
    # Lengths
    if(!is.null(Lengths)) {
      gene.id <- rownames(countData)
      Lengths <- Lengths[match(gene.id,names(Lengths))]
    }
    # Mean fragment lengths
    if(!is.null(MeanFragLengths)) {
      cell.id <- colnames(countData)
      MeanFragLengths <- MeanFragLengths[match(cell.id,names(MeanFragLengths))]
    }
    # spike-in counts
    if(!is.null(spikeData)) {
      cell.id <- colnames(countData)
      spikeData <- spikeData[,match(cell.id, colnames(spikeData))]
    }
  }

  return(list(countData = countData,
              batchData = batchData,
              spikeData = spikeData,
              spikeInfo = spikeInfo,
              Lengths = Lengths,
              MeanFragLengths =MeanFragLengths))
}


# estParam ----------------------------------------------------------------

#' @importFrom parallel detectCores
.run.estParam <- function(countData,
                          batchData,
                          spikeData,
                          spikeInfo,
                          Lengths,
                          MeanFragLengths,
                          Distribution,
                          RNAseq,
                          normalisation,
                          NCores,
                          sigma,
                          verbose) {
  # kick out empty samples and keep only expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData, na.rm = TRUE) > 0
  DetectG <- rowMeans(countData, na.rm = TRUE) > 0
  countData <- countData[DetectG,fullS]

  if(!is.null(Lengths)) {
    Lengths <- Lengths[DetectG]
  }
  if(!is.null(MeanFragLengths)) {
    MeanFragLengths <- MeanFragLengths[fullS]
  }
  if(!is.null(batchData)) {
    batchData <- batchData[fullS, , drop=F]
  }

  if(!is.null(spikeData) && !is.null(spikeInfo)) {
    # kick out undetected spike-ins
    spikeData <- spikeData[rowSums(spikeData)>0, colSums(spikeData)>100]
    if(verbose) {message(paste0(nrow(spikeData), " spike-ins have been detected in ", ncol(spikeData), " samples."))}
    # sort them if needed
    spikeInfo <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData), , drop = FALSE]
    spikeInfo <- spikeInfo[match(rownames(spikeData), rownames(spikeInfo)), , drop = FALSE ]
    if(nrow(spikeData)<10 || nrow(spikeInfo)<10) {
      stop(message(paste0("Not enough spike-ins detected to allow reliable normalisation. Please proceed with spike-in independent methods, e.g. TMM, scran, SCnorm, etc.")))
      }
  }

  if(!is.null(spikeData) && is.null(spikeInfo)) {
    # kick out undetected spike-ins
    spikeData <- spikeData[rowSums(spikeData)>0, colSums(spikeData)>100]
    if(verbose) { message(paste0(nrow(spikeData), " spike-ins have been detected in ", ncol(spikeData), " samples."))}
    if(nrow(spikeData)<10) {
      stop(message(paste0("Not enough spike-ins detected to allow reliable normalusation. Please proceed with spike-in indepdent methods, e.g. TMM, scran, SCnorm etc.")))
    }
  }

  countData0 <- countData == 0
  grand.dropout <- sum(countData0)/(nrow(countData)*ncol(countData))

  # normalisation
  NormData <- .norm.calc(normalisation=normalisation,
                         countData=countData,
                         spikeData=spikeData,
                         spikeInfo=spikeInfo,
                         batchData=batchData,
                         Lengths=Lengths,
                         MeanFragLengths=MeanFragLengths,
                         PreclustNumber=NULL,
                         NCores=NCores,
                         verbose=verbose)

  # parameters: mean, dispersion, dropout
  # fitting: mean vs dispersion, mean vs dropout
  if(RNAseq=='singlecell') {
    if(Distribution=='NB') {ParamData <- .est.NB.sc(countData = countData,
                                                    NormData = NormData,
                                                    sigma = sigma)}
    if(Distribution=='ZINB') {ParamData <- .est.ZINB.sc(countData = countData,
                                                        NormData = NormData,
                                                        sigma = sigma)}
  }
  if(RNAseq=='bulk') {
    if(Distribution=='NB') {ParamData <- .est.NB.bulk(countData = countData,
                                                      NormData = NormData,
                                                      sigma = sigma)}
  }

  # results return
  res <- c(ParamData, list(totalS=totalS,
                           totalG=totalG,
                           grand.dropout=grand.dropout,
                           sf=NormData$size.factors,
                           Lengths=Lengths,
                           MeanFragLengths=MeanFragLengths))
  return(res)
}

# single cell NB Parameters -----------------------------------------------

#' @importFrom cobs cobs
#' @importFrom stats predict
#' @importFrom msir loess.sd
.est.NB.sc <- function(countData, NormData, sigma) {

  # sequencing depth
  seqDepth = colSums(countData)
  names(seqDepth) = colnames(countData)

  nsamples = dim(countData)[2]
  counts0 = countData == 0
  nn0 = rowSums(!counts0)

  X2 = NormData$NormCounts

  # the negative binomial
  mu = rowSums(X2)/ncol(X2)
  s2 = rowSums((X2 - mu) ^ 2) / ncol(X2)
  size = mu ^ 2 / (s2 - mu + 1e-04)
  size = ifelse(size > 0, size, NA)
  p0 = (nsamples - nn0)/nsamples
  mu = mu[!is.na(size)]
  p0 = p0[!is.na(size)]
  size = size[!is.na(size)]
  phi.g = 1/size
  phi.c = mean(phi.g)
  ldisp = log2(phi.g)
  lsize = log2(size)
  lmu = log2(mu + 1)

  estG <- length(mu)
  estS <- length(NormData$size.factors)

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  # return object
  res <- list(seqDepth = seqDepth,
              means = mu,
              dispersion = phi.g,
              common.dispersion = phi.c,
              size = size,
              p0 = p0,
              meansizefit = meansizefit,
              meandispfit = meandispfit,
              estS = estS,
              estG = estG)
  return(res)
}

# single cell ZINB Parameters ---------------------------------------------

#' @importFrom cobs cobs
#' @importFrom stats predict complete.cases
#' @importFrom msir loess.sd
#' @importFrom matrixStats rowSds
.est.ZINB.sc <- function(countData, NormData, sigma){

  # sequencing depth
  seqDepth = colSums(countData)
  names(seqDepth) = colnames(countData)

  nsamples = dim(countData)[2]
  counts0 = countData == 0
  nn0 = rowSums(!counts0)

  X2 = NormData$NormCounts

  # the negative binomial portion
  mu = rowSums((!counts0) * countData)/nn0
  s2 = rowSums((!counts0) * (countData - mu)^2)/(nn0 - 1)
  size = mu ^ 2 / (s2 - mu + 1e-04)
  size = ifelse(size > 0, size, NA)
  p0 = (nsamples - nn0)/nsamples
  mu = mu[!is.na(size)]
  p0 = p0[!is.na(size)]
  size = size[!is.na(size)]
  phi.g = 1/size
  phi.c = mean(phi.g)
  ldisp = log2(phi.g)
  lsize = log2(size)
  lmu = log2(mu + 1)

  estG <- length(mu)
  estS <- length(NormData$size.factors)

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)
  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  mu.withzero <- rowMeans(X2)
  sd.withzero <- matrixStats::rowSds(X2)
  cv.withzero <- sd.withzero/mu.withzero
  nonamplified <- cv.withzero < 1.5*(sqrt(mu.withzero)/mu.withzero) & mu.withzero <=10

  if(sum(nonamplified, na.rm = T)/length(nonamplified) < 0.05){
    # all
    cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
    cobs.sim <- runif(1000, min(lmu, na.rm=TRUE),  max(lmu, na.rm=TRUE))
    cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
    cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
    cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
    p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

    if(length(p0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, p0)
      meanp0fit = loess.sd(x=all.dat$lmu, y=all.dat$p0, nsigma = sigma)
    }

    if(!length(p0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, p0)
      all.dat <- all.dat[all.dat$lmu<p0.cut,]
      meanp0fit = loess.sd(x=all.dat$lmu, y=all.dat$p0, nsigma = sigma)
    }

    meanp0fit.amplified = p0.cut.amplified = meanp0fit.nonamplified = p0.cut.nonamplified = NULL

  }

  if(sum(nonamplified, na.rm = T)/length(nonamplified) >= 0.05){
    p0.amplified <- p0[!nonamplified]
    p0.nonamplified <- p0[nonamplified]

    lmu.amplified <- lmu[!nonamplified]
    lmu.nonamplified <- lmu[nonamplified]

    # estimate the knee point for curves of mean versus dropout
    # amplified
    cobs.fit.amplified <- cobs::cobs(x = lmu.amplified,
                                     y = p0.amplified,
                                     constraint = 'decrease',
                                     nknots = 20,
                                     print.warn = F, print.mesg = F)
    cobs.sim.amplified <- runif(1000, min(lmu.amplified, na.rm=T),  max(lmu.amplified, na.rm = T))
    cobs.predict.amplified <- as.data.frame(stats::predict(cobs.fit.amplified, cobs.sim.amplified))
    cobs.predict.amplified[,"fit"] <- ifelse(cobs.predict.amplified$fit < 0, 0, cobs.predict.amplified[,"fit"])
    cobs.predict.amplified <- cobs.predict.amplified[cobs.predict.amplified[,'fit'] < 0.05,]
    p0.cut.amplified <- as.numeric(cobs.predict.amplified[which.max(cobs.predict.amplified[,'fit']), 'z'])
    if(length(p0.cut.amplified)==0) {
      amplified.dat <- cbind.data.frame(lmu.amplified, p0.amplified)
      amplified.dat <- amplified.dat[stats::complete.cases(amplified.dat),]
      meanp0fit.amplified = msir::loess.sd(x = amplified.dat$lmu.amplified, y=amplified.dat$p0.amplified, nsigma = sigma)
    }
    if(!length(p0.cut.amplified)==0) {
    amplified.dat <- cbind.data.frame(lmu.amplified, p0.amplified)
    amplified.dat <- amplified.dat[stats::complete.cases(amplified.dat),]
    amplified.dat <- amplified.dat[amplified.dat$lmu.amplified<p0.cut.amplified,]
    meanp0fit.amplified = msir::loess.sd(x = amplified.dat$lmu.amplified, y=amplified.dat$p0.amplified, nsigma = sigma)
    }

    # nonamplified
    cobs.fit.nonamplified <- cobs::cobs(x = lmu.nonamplified,
                                        y = p0.nonamplified,
                                        constraint = 'decrease',
                                        nknots = 20,
                                        print.warn = F, print.mesg = F)
    cobs.sim.nonamplified <- runif(1000, min(lmu.nonamplified, na.rm=TRUE),  max(lmu.nonamplified, na.rm=TRUE))
    cobs.predict.nonamplified <- as.data.frame(stats::predict(cobs.fit.nonamplified, cobs.sim.nonamplified))
    cobs.predict.nonamplified[,"fit"] <- ifelse(cobs.predict.nonamplified$fit < 0, 0, cobs.predict.nonamplified[,"fit"])
    cobs.predict.nonamplified <- cobs.predict.nonamplified[cobs.predict.nonamplified[,'fit'] < 0.05,]
    p0.cut.nonamplified <- log(10)
    nonamplified.dat <- cbind.data.frame(lmu.nonamplified, p0.nonamplified)
    nonamplified.dat <- nonamplified.dat[stats::complete.cases(nonamplified.dat),]
    nonamplified.dat <- nonamplified.dat[nonamplified.dat$lmu.nonamplified<p0.cut.nonamplified,]

    if(nrow(nonamplified.dat)==0) {
      nonamplified.dat <- cbind.data.frame(lmu.nonamplified, p0.nonamplified)
      nonamplified.dat <- nonamplified.dat[stats::complete.cases(nonamplified.dat),]
      meanp0fit.nonamplified = msir::loess.sd(x = nonamplified.dat$lmu.nonamplified, y=nonamplified.dat$p0.nonamplified, nsigma = sigma)
    }
    if(nrow(nonamplified.dat)>0) {
      meanp0fit.nonamplified = msir::loess.sd(x = nonamplified.dat$lmu.nonamplified, y=nonamplified.dat$p0.nonamplified, nsigma = sigma)
    }

    # all
    cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
    cobs.sim <- runif(1000, min(lmu, na.rm=TRUE),  max(lmu, na.rm=TRUE))
    cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
    cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
    cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
    p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

    if(length(p0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, p0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      meanp0fit = msir::loess.sd(x=all.dat$lmu, y=all.dat$p0, nsigma = sigma)
    }

    if(!length(p0.cut)==0) {
      all.dat <- cbind.data.frame(lmu, p0)
      all.dat <- all.dat[stats::complete.cases(all.dat),]
      all.dat <- all.dat[all.dat$lmu<p0.cut,]
      meanp0fit = msir::loess.sd(x=all.dat$lmu, y=all.dat$p0, nsigma = sigma)
    }

  }

  nonamplified = length(which(nonamplified)) / length(nonamplified)

  # return object
  res <- list(seqDepth = seqDepth,
              means = mu,
              dispersion = phi.g,
              common.dispersion = phi.c,
              size = size,
              p0 = p0,
              meansizefit = meansizefit,
              meandispfit = meandispfit,
              meanp0fit = meanp0fit,
              p0.cut = p0.cut,
              meanp0fit.amplified = meanp0fit.amplified,
              p0.cut.amplified = p0.cut.amplified,
              meanp0fit.nonamplified = meanp0fit.nonamplified,
              p0.cut.nonamplified = p0.cut.nonamplified,
              nonamplified = nonamplified,
              estS = estS,
              estG = estG)
  return(res)
}

# bulk NB parameters ------------------------------------------------------

#' @importFrom cobs cobs
#' @importFrom stats predict complete.cases
#' @importFrom msir loess.sd
.est.NB.bulk <- function(countData, NormData, sigma) {

  # sequencing depth
  seqDepth = colSums(countData)
  names(seqDepth) = colnames(countData)

  nsamples = dim(countData)[2]
  counts0 = countData == 0
  nn0 = rowSums(!counts0)

  X2 = NormData$NormCounts

  # the negative binomial portion
  mu = rowSums((!counts0) * countData)/nn0
  s2 = rowSums((!counts0) * (countData - mu)^2)/(nn0 - 1)
  size = mu ^ 2 / (s2 - mu + 1e-04)
  size = ifelse(size > 0, size, NA)
  p0 = (nsamples - nn0)/nsamples
  mu = mu[!is.na(size)]
  p0 = p0[!is.na(size)]
  size = size[!is.na(size)]
  phi.g = 1/size
  phi.c = mean(phi.g)
  ldisp = log2(phi.g)
  lsize = log2(size)
  lmu = log2(mu + 1)

  estG <- length(mu)
  estS <- length(NormData$size.factors)

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)
  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu, na.rm=TRUE),  max(lmu, na.rm=TRUE))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  # return object
  res <- list(seqDepth = seqDepth,
              means = mu,
              dispersion = phi.g,
              common.dispersion = phi.c,
              size = size,
              p0 = p0,
              meansizefit = meansizefit,
              meandispfit = meandispfit,
              p0.cut = p0.cut,
              estS = estS,
              estG = estG)
  return(res)
}

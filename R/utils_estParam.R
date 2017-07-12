# estParam: single cell ---------------------------------------------------

#' @importFrom parallel detectCores
.run.estParam <- function(countData,
                          spikeData,
                          spikeInfo,
                          Lengths,
                          MeanFragLengths,
                          Distribution,
                          RNAseq,
                          normalisation,
                          NCores,
                          sigma) {
  # kick out empty samples and keep only expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData, na.rm = TRUE) > 0
  DetectG <- rowMeans(countData, na.rm = TRUE) > 0
  countData <- countData[DetectG,fullS]
  countData0 <- countData == 0
  grand.dropout <- sum(countData0)/(nrow(countData)*ncol(countData))

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }

  ncores = ifelse(is.null(NCores), parallel::detectCores() - 1, NCores)

  # normalisation
  NormData <- .norm.calc(normalisation=normalisation,
                         countData=countData,
                         spikeData=spikeData,
                         spikeInfo=spikeInfo,
                         Lengths=Lengths,
                         MeanFragLengths=MeanFragLengths,
                         NCores=ncores)

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
                           sf=NormData$size.factors))
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

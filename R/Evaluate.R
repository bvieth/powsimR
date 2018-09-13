# EVALUATE DISTRIBUTIONAL FITS --------------------------------------------

#' @name evaluateDist
#' @aliases evaluateDist
#' @title Model Diagnostics for RNAseq Data
#' @description With this function, the user can determine goodness of fit for each gene.
#' @usage evaluateDist(countData, batchData =NULL,
#' spikeData = NULL, spikeInfo = NULL,
#' Lengths = NULL, MeanFragLengths = NULL,
#' RNAseq, normalisation,
#' frac.genes=1, min.meancount = 0.1,
#' max.dropout=0.7, min.libsize=1000,
#' verbose = TRUE)
#' @param countData is a count matrix (row=gene, column=sample).
#' Please provide the measurements of one group only, e.g. the control group.
#' @param batchData is a \code{data.frame} for batch annotation.
#' Rows correspond to samples. The first column should contain the batches, e.g. 'a', 'b', 'c', etc..
#' This is only used for the normalisation.
#' @param spikeData is a count \code{matrix}.
#' Rows correspond to spike-ins, columns to samples.
#' The order of columns should be the same as in the \code{countData}.
#' @param spikeInfo is a molecule count \code{matrix} of spike-ins.
#' Rows correspond to spike-ins. The order of rows should be the same as in the \code{spikeData}.
#' The column names should be 'SpikeID' and 'SpikeInput' for molecule counts of spike-ins.
#' @param Lengths is a numeric vector of transcript lengths with the same length and order as the rows in countData.
#' This variable is only used for internal TPM calculations if Census normalization is specified.
#' @param MeanFragLengths is a numeric vector of mean fragment lengths with the same length as columns in countData.
#' This variable is only used for internal TPM calculations if Census normalization is specified.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param normalisation is a character value: 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm',
#' 'SCnorm', 'RUV', 'Census', 'depth', 'none'.
#' For more information, please consult the details section of \code{\link{estimateParam}}.
#' @param frac.genes The fraction of genes to calculate goodness of fit statistics, default is 1, i.e. for all genes.
#' @param min.meancount The minimum mean normalised count per gene if fraction of genes are defined. Default is \code{0.1}.
#' @param max.dropout The maximal percentage of zero expression per gene. Genes with more than \code{max.dropout} will be remove. Default is \code{0.7}, i.e. genes with more than 70\% dropouts.
#' @param min.libsize The minimum raw read counts per sample, default is \code{1000}.
#' @param verbose Logical value to indicate whether to show progress report of simulations.
#' @return List object with the results of goodness of fit and estimated parameters:
#' \item{edgeR_res}{Goodness-of-fit statistic, degrees of freedom and associated p-value using the deviance and residual degrees of freedom from \code{\link[edgeR]{glmFit}}. Furthermore, the AIC of the edgeR model fit using the residuals of \code{\link[edgeR]{zscoreNBinom}}.}
#' \item{GOF_res}{The fitting results per distribution, including loglikelihood, goodness-of-fit statistics, AIC and predicted number of zeroes. The following distributions were considered: Poisson, negative binomial, zero-inflated poisson and negative binomial following the 'standard' (i.e. \code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}} and \code{\link[pscl]{zeroinfl}} implementation) and fitdist approach (see \code{\link[fitdistrplus]{fitdist}}) and Beta-Poisson following Marioni or Hemberg parameterisation (see \url{https://doi.org/10.1186/gb-2013-14-1-r7} and \url{https://dx.doi.org/10.1186%2Fs12859-016-0944-6}, respectively). Furthermore, model fit comparison by LRT for nested and Vuong Test for non-nested models.}
#' \item{EST_res}{The estimated parameters of distribution fitting.}
#' \item{ObservedZeros}{The number of zeroes and dropout rate per gene.}
#' @examples
#' \dontrun{
#' ## using example data set
#' data(kolodziejczk_cnts)
#' evaldist <- evaluateDist(countData = kolodziejczk_cnts,
#' RNAseq = "singlecell", normalisation="scran",
#' frac.genes=1, min.meancount = 0.1,
#' max.dropout=0.7, min.libsize=1000,
#' verbose = TRUE)
#' plotEvalDist(evaldist, annot = TRUE)
#' }
#' @author Beate Vieth
#' @rdname evaluateDist
#' @importFrom edgeR DGEList cpm.DGEList estimateDisp glmFit zscoreNBinom
#' @importFrom stats model.matrix glm logLik AIC dnbinom dpois fitted na.omit
#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl vuong
#' @importFrom fitdistrplus fitdist
#' @importFrom gamlss.dist ZIP ZINBI
#' @importFrom nonnest2 vuongtest
#' @export
evaluateDist <- function(countData, batchData =NULL,
                         spikeData = NULL, spikeInfo = NULL,
                         Lengths = NULL, MeanFragLengths = NULL,
                         RNAseq, normalisation,
                         frac.genes=1, min.meancount = 0.1,
                         max.dropout=0.7, min.libsize=1000,
                         verbose = TRUE) {

  invisible(gamlss.dist::ZINBI())
  invisible(gamlss.dist::ZIP())

  options(stringsAsFactors = F)

  # check the matching of names
  checkup <- .run.checkup(countData = countData,
                          batchData = batchData,
                          spikeData = spikeData,
                          spikeInfo = spikeInfo,
                          Lengths = Lengths,
                          MeanFragLengths = MeanFragLengths,
                          RNAseq = RNAseq,
                          verbose=verbose)

  countData <- checkup$countData
  batchData <- checkup$batchData
  spikeData <- checkup$spikeData
  spikeInfo <- checkup$spikeInfo
  Lengths <- checkup$Lengths
  MeanFragLengths <- checkup$MeanFragLengths

  # kick out dropout genes and very small libs
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highE <- p0 < max.dropout
  fullS <- colSums(countData) > min.libsize

  countData <- countData[highE,fullS]
  batchData <- batchData[fullS,]
  Lengths <- Lengths[highE]
  MeanFragLengths <- MeanFragLengths[fullS]

  # run estimation
  estParam = .run.estParam(countData = countData,
                           batchData = batchData,
                           spikeData = spikeData,
                           spikeInfo = spikeInfo,
                           Lengths = Lengths,
                           MeanFragLengths = MeanFragLengths,
                           Distribution = 'NB',
                           RNAseq = RNAseq,
                           normalisation = normalisation,
                           sigma = 1.96,
                           NCores = NULL,
                           verbose=verbose)

  countData <- countData[,colnames(countData) %in% names(estParam$seqDepth)]


  # set up dge edgeR
  sf <- estParam$sf
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/estParam$seqDepth)
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        remove.zeros = FALSE)
  # calculate normalized counts (if norm lib is true does not sum up to 1 million)
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  # estimate dispersions
  invisible(capture.output(
    dge <- suppressMessages(edgeR::estimateDisp(y=dge))
  ))
  # design.mat
  design.mat <- matrix(1,ncol(dge$counts),1)
  rownames(design.mat) <- colnames(dge$counts)
  colnames(design.mat) <- "Intercept"
  # apply edgeR glm fit
  fit.edgeR <- edgeR::glmFit(dge, design = design.mat)

  # calculate residuals
  y <- fit.edgeR$counts
  mu <- fit.edgeR$fitted.values
  phi <- fit.edgeR$dispersion
  coefs <- fit.edgeR$coefficients
  # v <- mu*(1+phi*mu)
  # d <- edgeR::nbinomUnitDeviance(y,mu,phi)
  # resid.pearson <- (y-mu) / sqrt(v)
  # resid.deviance <- sign(y-mu) * sqrt(d)
  resid.quantile <- edgeR::zscoreNBinom(y,mu=mu,size=1/phi)
  # calculate AIC
  edgeR.aic <- .myAIC(resids = resid.quantile, dat.n = ncol(y), npar = ncol(coefs), k=2)
  # calculate gof
  edgeR.gof.stats <- .myGOF(deviances = fit.edgeR$deviance, df.residuals = fit.edgeR$df.residual)

  # make output table
  genenames <- rownames(dge$counts)
  headers.1 <- c("pois_standard", "nbinom_standard", "zifpois_standard", "zifnbinom_standard")
  headers.2 <- c("loglikelihood", "loglikelihooddf", "gofstat", "gofdf", "gofpval", "predzero", "aic")
  headers_tmp1 <- paste(rep(headers.1, each=7), rep(headers.2, times=4), sep="_")
  headers.3 <- c("pois_fitdistr", "nbinom_fitdistr", "zifpois_fitdistr", "zifnbinom_fitdistr")
  headers.4 <- c("gofstat", "gofdf", "gofpval")
  headers_tmp2 <- paste(rep(headers.3, each=3), rep(headers.4, times=4), sep="_")
  headers.5 <- c("PoiBeta_Hemberg", "PoiBeta_Marioni")
  headers.6 <- c("loglikelihood", "loglikelihooddf", "predzero", 'aic')
  headers_tmp3 <- paste(rep(headers.5, each=4), rep(headers.6, times=2), sep="_")
  headers_tmp4 <- c('LRT_standard_NBPoisson', 'Vuong_standard_ZPoisson', 'Vuong_standard_ZNB', 'Vuong_standard_ZNBZPoisson')
  headersgof <- c(headers_tmp1, headers_tmp2, headers_tmp3, headers_tmp4)

  gof.res <- data.frame(matrix(vector(), length(genenames), length(headersgof),
                               dimnames=list(c(genenames), c(headersgof))))

  headersests <- c("pois_lambda", "nbinom_size", "nbinom_mu", "zifpois_mu", "zifpois_sigma", "zifnbinom_mu", "zifnbinom_sigma", "zifnbinom_nu", 'alpha_bp_hemberg', 'alpha_bp_marioni', 'beta_bp_hemberg', 'beta_bp_marioni', 'gamma_bp_hemberg', 'gamma_bp_marioni')
  estimate.res <- data.frame(matrix(vector(), length(genenames), length(headersests),
                                    dimnames=list(c(genenames), c(headersests))))

  # run model fits for the genes
  if(frac.genes==1) {
    geneid <- 1:nrow(dge$counts)
  }
  if(!frac.genes==1) { # take fraction of genes per mean bins
    meancnts = rowMeans(out.cpm)
    tmp.ecdf.mean = stats::ecdf(meancnts)
    tmp.quantile.mean = stats::quantile(tmp.ecdf.mean, probs=c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9))
    qbins.mean = c(min.meancount,unname(tmp.quantile.mean),Inf)
    bin.mean = cut(meancnts, qbins.mean)
    geneid = unname(unlist(sapply(levels(bin.mean), function(x) {
      tmp.group <- genenames[(bin.mean %in% x)]
      tmp.id <- sample(tmp.group, size= round(frac.genes*length(tmp.group)), replace = F)
    })))
  }

  for (i in geneid) {
    if(verbose) {
      message(paste0("Estimation for gene ", which(geneid == i), " out of ", length(geneid),
                     ". A total of ", nrow(dge$counts),
                     " genes are available for estimation."))
    }

    # design.mat
    design.mat <- matrix(1,ncol(dge$counts),1)
    # zero measurement indicator
    counts0 <- dge$counts[i,] == 0
    # number of samples
    nsamples=length(dge$counts[i,])
    # number of nonzero samples
    nn0 = sum(!counts0)
    # dropout
    p0 = mean(dge$counts[i,]==0)
    # mean and dispersion of nonzero portion of normalised counts
    estmu = sum((!counts0) * dge$counts[i,])/nn0
    s2 = sum((!counts0) * (dge$counts[i,] - estmu)^2)/(nn0 - 1)
    disp =  1 / (estmu^2/(s2 - estmu + 1e-04))

    # fit glm poisson
    tmp.fit.pois <- try(stats::glm(dge$counts[i,] ~ 1, family=poisson, offset = edgeR::getOffset(y=dge)), silent = TRUE)
    # fit poisson with fitdistrplus
    tmp.fit.pois.wo <- try(fitdistrplus::fitdist(data = dge$counts[i,], 'pois'), silent = TRUE)
    # fit glm negative binomial
    tmp.fit.nbinom <- try(MASS::glm.nb(dge$counts[i,] ~ 1 + offset(edgeR::getOffset(y=dge))), silent = TRUE)
    # fit negative binomial with fitdistrplus
    tmp.fit.nbinom.wo <- try(fitdistrplus::fitdist(data = dge$counts[i,], 'nbinom'), silent = TRUE)
    # fit glm zero inflated poisson
    tmp.fit.zifpois <- try(pscl::zeroinfl(dge$counts[i,] ~ 1, dist = 'poisson', offset = edgeR::getOffset(y=dge), EM = FALSE), silent = TRUE)
    # fit zero inflated poisson with fitdistrplus
    tmp.fit.zifpois.wo <- try(fitdistrplus::fitdist(dge$counts[i,], "ZIP", start = list(mu = estmu, sigma=p0), discrete = TRUE, lower = c(0, 0), upper = c(Inf, 1)), silent = TRUE)
    # fit glm zero inflated negative binomial
    tmp.fit.zifnbinom <- try(pscl::zeroinfl(dge$counts[i,] ~ 1, dist = 'negbin', offset = edgeR::getOffset(y=dge), EM = FALSE), silent = TRUE)
    # fit zero inflated negative binomial with fitdistrplus
    tmp.fit.zifnbinom.wo <- try(fitdistrplus::fitdist(dge$counts[i,], "ZINBI", start = list(mu = estmu, sigma=disp, nu=p0), discrete = TRUE, lower = c(0, 0, 0), upper = c(Inf,Inf, 1)), silent = TRUE)
    # fit poisson-beta distribution
    tmp.fit.pb <- try(.PoissonBetaFit(x.raw = dge$counts[i,], x.norm = out.cpm[i,]), silent = TRUE)

    # fill up results table
    # poisson
    if(inherits(tmp.fit.pois, 'try-error')) gof.res[i, grep("^pois_standard", colnames(gof.res), perl=TRUE)] <- NA
    else{pois.gof.stats <- .myGOF(deviances = tmp.fit.pois$deviance, df.residuals =  tmp.fit.pois$df.residual)
    pois.predzero <- round(sum(stats::dpois(0, stats::fitted(tmp.fit.pois))))
    pois.aic <- stats::AIC(tmp.fit.pois)
    pois.loglik <- as.numeric(logLik(tmp.fit.pois))
    pois.loglikdf <- as.numeric(attr(logLik(tmp.fit.pois), "df"))
    gof.res[i, grep("^pois_standard", colnames(gof.res), perl=TRUE)] <-  cbind(pois.loglik, pois.loglikdf, pois.gof.stats, pois.predzero, pois.aic)
    }
    # poisson with fitdistr
    if(inherits(tmp.fit.pois.wo, 'try-error')) gof.res[i, grep("^pois_fitdistr", colnames(gof.res), perl=TRUE)] <- NA
    else{pois.gof.fitdistr <- try(.fitdistrplusGOF(fitdistobj=tmp.fit.pois.wo), silent=T)
    if(inherits(pois.gof.fitdistr, 'try-error')) {
      gof.res[i, grep("^pois_fitdistr", colnames(gof.res), perl=TRUE)] <-  NA
    } else{
      gof.res[i, grep("^pois_fitdistr", colnames(gof.res), perl=TRUE)] <-  cbind(pois.gof.fitdistr)
      estimate.res[i, grep("^pois_", colnames(estimate.res), perl=TRUE)] <- tmp.fit.pois.wo$estimate
    }
    }
    # neg binom
    if(inherits(tmp.fit.nbinom, 'try-error')) gof.res[i, grep("^nbinom_standard", colnames(gof.res), perl=TRUE)] <- NA
    else{nbinom.gof.stats <- .myGOF(deviances = tmp.fit.nbinom$deviance, df.residuals =  tmp.fit.nbinom$df.residual)
    nbinom.predzero <- round(sum(stats::dnbinom(0, mu = stats::fitted(tmp.fit.nbinom), size = tmp.fit.nbinom$theta)))
    nbinom.aic <- stats::AIC(tmp.fit.nbinom)
    nbinom.loglik <- as.numeric(logLik(tmp.fit.nbinom))
    nbinom.loglikdf <- as.numeric(attr(logLik(tmp.fit.nbinom), "df"))
    gof.res[i, grep("^nbinom_standard", colnames(gof.res), perl=TRUE)] <-  cbind(nbinom.loglik, nbinom.loglikdf, nbinom.gof.stats, nbinom.predzero, nbinom.aic)}
    # neg binom with fitdistr
    if(inherits(tmp.fit.nbinom.wo, 'try-error')) gof.res[i, grep("^nbinom_fitdistr", colnames(gof.res), perl=TRUE)] <- NA
    else{nbinom.gof.fitdistr <- try(.fitdistrplusGOF(fitdistobj=tmp.fit.nbinom.wo), silent=T)
    if(inherits(nbinom.gof.fitdistr, 'try-error')) {
      gof.res[i, grep("^nbinom_fitdistr", colnames(gof.res), perl=TRUE)] <-  NA
    } else{
      gof.res[i, grep("^nbinom_fitdistr", colnames(gof.res), perl=TRUE)] <-  cbind(nbinom.gof.fitdistr)
      estimate.res[i, grep("^nbinom_", colnames(estimate.res), perl=TRUE)] <- tmp.fit.nbinom.wo$estimate
    }}
    # zero inflated poisson
    if(inherits(tmp.fit.zifpois, 'try-error')) gof.res[i, grep("^zifpois_standard", colnames(gof.res), perl=TRUE)] <- NA
    else{zifpois.predzero <- round(sum(predict(tmp.fit.zifpois, type = "prob")[, 1]))
    zifpois.aic <- stats::AIC(tmp.fit.zifpois)
    zifpois.loglik <- as.numeric(logLik(tmp.fit.zifpois))
    zifpois.loglikdf <- as.numeric(attr(logLik(tmp.fit.zifpois), "df"))
    zifpois.gof.stats <- .myGOF(deviances = -2*logLik(tmp.fit.zifpois), df.residuals =  length(dge$counts[i,]-2))
    gof.res[i, grep("^zifpois_standard", colnames(gof.res), perl=TRUE)] <-  cbind(zifpois.loglik, zifpois.loglikdf, zifpois.gof.stats, zifpois.predzero, zifpois.aic)}
    # zero inflated poisson with fitdistr
    if(inherits(tmp.fit.zifpois.wo, 'try-error')) gof.res[i, grep("^zifpois_fitdistr", colnames(gof.res), perl=TRUE)] <- NA
    else{zifpois.gof.fitdistr <- try(.fitdistrplusGOF(fitdistobj=tmp.fit.zifpois.wo), silent=T)
    if(inherits(zifpois.gof.fitdistr, 'try-error')) {
      gof.res[i, grep("^zifpois_fitdistr", colnames(gof.res), perl=TRUE)] <-  NA
    } else{
      gof.res[i, grep("^zifpois_fitdistr", colnames(gof.res), perl=TRUE)] <-  cbind(zifpois.gof.fitdistr)
      estimate.res[i, grep("^zifpois_", colnames(estimate.res), perl=TRUE)] <- tmp.fit.zifpois.wo$estimate
    }}
    # zero inflated neg binom
    if(inherits(tmp.fit.zifnbinom, 'try-error')) gof.res[i, grep("^zifnbinom_standard", colnames(gof.res), perl=TRUE)] <- NA
    else{zifnbinom.predzero <- round(sum(predict(tmp.fit.zifnbinom, type = "prob")[, 1]))
    zifnbinom.aic <- stats::AIC(tmp.fit.zifnbinom)
    zifnbinom.loglik <- as.numeric(logLik(tmp.fit.zifnbinom))
    zifnbinom.loglikdf <- as.numeric(attr(logLik(tmp.fit.zifnbinom), "df"))
    zifnbinom.gof.stats <- .myGOF(deviances = -2*logLik(tmp.fit.zifnbinom), df.residuals =  length(dge$counts[i,]-1))
    gof.res[i, grep("^zifnbinom_standard", colnames(gof.res), perl=TRUE)]  <-  cbind(zifnbinom.loglik, zifnbinom.loglikdf, zifnbinom.gof.stats, zifnbinom.predzero, zifnbinom.aic)}
    # zero inflated neg binom with fitdistr
    if(inherits(tmp.fit.zifnbinom.wo, 'try-error')) gof.res[i, grep("^zifnbinom_fitdistr", colnames(gof.res), perl=TRUE)] <- NA
    else{zifnbinom.gof.fitdistr <- try(.fitdistrplusGOF(fitdistobj=tmp.fit.zifnbinom.wo), silent=T)
    if(inherits(zifnbinom.gof.fitdistr, 'try-error')) {
      gof.res[i, grep("^zifnbinom_fitdistr", colnames(gof.res), perl=TRUE)] <-  NA
    } else{
      gof.res[i, grep("^zifnbinom_fitdistr", colnames(gof.res), perl=TRUE)] <-  cbind(zifnbinom.gof.fitdistr)
      estimate.res[i, grep("^zifnbinom_", colnames(estimate.res), perl=TRUE)] <- tmp.fit.zifnbinom.wo$estimate
    }}
    # LR Test for NB > P ?
    if(all(!inherits(tmp.fit.pois, 'try-error') && !inherits(tmp.fit.nbinom, 'try-error'))) {
      tmp.lrt.nbp = try(pchisq(as.numeric(2 * (logLik(tmp.fit.nbinom) - logLik(tmp.fit.pois))), df = 1, lower.tail = FALSE), silent=T)
      if(!inherits(tmp.lrt.nbp, 'try-error')) {
        gof.res[i, grep("LRT_standard_NBPoisson", colnames(gof.res), perl=TRUE)] <- tmp.lrt.nbp
      }
    }

    # Vuong Test for ZINB > NB, ZIP > P, ZNB >ZP ?
    if(all(!inherits(tmp.fit.zifpois, 'try-error') && !inherits(tmp.fit.pois, 'try-error'))) {
      tmp.vuong.p = try(nonnest2::vuongtest(tmp.fit.zifpois,tmp.fit.pois), silent=T)
      if(!inherits(tmp.vuong.p, 'try-error')) {
        gof.res[i, grep("Vuong_standard_ZPoisson", colnames(gof.res), perl=TRUE)] <- tmp.vuong.p$p_LRT$A
      }
    }
    if(all(!inherits(tmp.fit.zifnbinom, 'try-error') && !inherits(tmp.fit.nbinom, 'try-error'))) {
      tmp.vuong.nb = try(nonnest2::vuongtest(tmp.fit.zifnbinom,tmp.fit.nbinom), silent=T)
      if(!inherits(tmp.vuong.nb, 'try-error')) {
        gof.res[i, grep("Vuong_standard_ZNB", colnames(gof.res), perl=TRUE)] <- tmp.vuong.nb$p_LRT$A
      }
    }
    if(all(!inherits(tmp.fit.zifnbinom, 'try-error') && !inherits(tmp.fit.zifpois, 'try-error'))) {
      tmp.vuong.nbp = try(nonnest2::vuongtest(tmp.fit.zifnbinom,tmp.fit.zifpois), silent=T)
      if(!inherits(tmp.vuong.nbp, 'try-error')) {
        gof.res[i, grep("Vuong_standard_ZNBZPoisson", colnames(gof.res), perl=TRUE)] <- tmp.vuong.nbp$p_LRT$A
      }
    }

    # poisson beta fit
    if(inherits(tmp.fit.pb, 'try-error')) gof.res[i, grep("^PoiBeta", colnames(gof.res))] <- NA
    else{gof.res[i, grep("^PoiBeta", colnames(gof.res))] <- c(tmp.fit.pb$LogLikelihood.hemberg, 3, tmp.fit.pb$PredZero.hemberg, tmp.fit.pb$AIC.hemberg, tmp.fit.pb$LogLikelihood.marioni, 3, tmp.fit.pb$PredZero.marioni, tmp.fit.pb$AIC.marioni)
    estimate.res[i,grep('bp', colnames(estimate.res), perl=TRUE)] <- cbind(tmp.fit.pb$bp.alpha.hemberg, tmp.fit.pb$bp.alpha.marioni, tmp.fit.pb$bp.beta.hemberg, tmp.fit.pb$bp.beta.marioni,tmp.fit.pb$bp.gamma.hemberg, tmp.fit.pb$bp.gamma.marioni)}

  }

  # fill in edgeR results
  edgeR.res <- cbind(edgeR.gof.stats, edgeR.aic)
  colnames(edgeR.res) <- c('edgeRglm_gofstat', 'edgeRglm_gofdf', 'edgeRglm_gofpval', 'edgeRglm_aic')
  rownames(edgeR.res) <- rownames(dge$counts)
  # observed zeros
  ObservedZeros <- data.frame(ObsZero=rowSums(dge$counts==0),
                              Dropout=rowMeans(dge$counts==0))
  rownames(ObservedZeros) <- rownames(dge$counts)

  # list result object:
  # goodness of fit statistics per model
  # estimated parameters
  # observed zeros
  return(list(edgeR_res= edgeR.res,
              GOF_res=gof.res,
              EST_res=estimate.res,
              ObservedZeros=ObservedZeros))
}

# EVALUATE SETUP ----------------------------------------------------------
#' @name evaluateSim
#' @aliases evaluateSim
#' @title Compute the performance related metrics from simulation results.
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes several metrics that give an indication of the simulation setup performance.
#' @usage evaluateSim(simRes, timing=TRUE)
#' @param simRes The result from \code{\link{simulateDE}}.
#' @param timing A logical vector indicating whether to summarise computational time of simulation run.
#' Default is \code{TRUE}.
#' @return A list with the following entries:
#' \item{Log2FoldChange}{The absolute mean error (\code{MAE}), root mean square error (\code{RMSE})
#' and root mean square residual error of a robust linear model (\code{rRMSE}, \code{\link[MASS]{rlm}})
#' of log2 fold change differences between estimated LFC and simulated LFC.
#' Furthermore, the fraction of missing entries (\code{NAFraction}) for MAE annd RMSE.}
#' \item{SizeFactors}{The median absolute deviation (MAD) between the estimated and true size factors,
#' the root mean square residual error of a robust linear model (rRMSE, \code{\link[MASS]{rlm}}) and
#' the ratio between estimated and true size factors of the two groups (\code{GroupX}).}
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}} for negative binomial parameters,
#' \code{\link{SimSetup}} and
#' \code{\link{DESetup}} for setting up simulation parameters and
#' \code{\link{simulateDE}} for simulating differential expression and
#' @examples
#' \dontrun{
#' ## using example data set
#' data(kolodziejczk_simDE)
#' eval.sim <- evaluateSim(simRes = kolodziejczk_simDE, timing = T)
#' }
#' @rdname evaluateSim
#' @importFrom stats mad
#' @importFrom mclust adjustedRandIndex
#' @importFrom matrixStats rowSds
#' @export
evaluateSim <- function(simRes, timing=TRUE) {

  # simulation parameters
  Nreps1 = simRes$sim.settings$n1
  Nreps2 = simRes$sim.settings$n2
  ngenes = simRes$sim.settings$ngenes
  sim.opts = simRes$sim.settings
  DEids = simRes$sim.settings$DEid
  tlfcs = simRes$sim.settings$pLFC
  nsims = simRes$sim.settings$nsims

  # estimated parameters
  elfcs = simRes$elfc
  tsfs = simRes$true.sf
  esfs = simRes$est.sf

  t.designs = simRes$true.designs

  if (attr(simRes, 'Simulation') == 'Flow') {
    e.designs = simRes$def.designs
  }
  if (attr(simRes, 'Simulation') == 'DE') {
    e.designs = NULL
  }

  # create output objects
  my.names = paste0(Nreps1, " vs ", Nreps2)
  # error in log2 fold changes
  lfc.error.mat <- lapply(1:length(my.names), function(x) {
    matrix(NA, nrow =  nsims, ncol = 15,
           dimnames = list(c(paste0("Sim", 1:nsims)),
                           c(paste0(rep(x=c("ALL", "DE", "EE"),each=5),"_",
                                    c('RMSE_Value', "MAE_Value", "RMSE_NAFraction", "MAE_NAFraction", "rRMSE_Value")))))
  })
  names(lfc.error.mat) <- my.names

  # error in size factors
  sf.error.mat <- lapply(1:length(my.names), function(x) {
    matrix(NA, nrow =  nsims, ncol = 4,
           dimnames = list(c(paste0("Sim_", 1:nsims)),
                           c("MAD", "rRMSE", "Group 1","Group 2")))
  })
  names(sf.error.mat) <- my.names

  # error in clustering
  if (attr(simRes, 'Simulation') == 'Flow') {
    clust.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 1,
             dimnames = list(c(paste0("Sim_", 1:nsims)),
                             c("RandIndex")))
    })
    names(clust.error.mat) <- my.names
  }
  if (attr(simRes, 'Simulation') == 'DE') {
    clust.error.mat = NULL
  }

  ## loop over simulation and replicates
  for(i in 1:nsims) {
    # DE flag
    DEid = DEids[[i]]
    Zg = rep(0, ngenes)
    Zg[DEid] = 1
    # true log fold change of all genes
    all.tlfc = tlfcs[[i]]
    # true log fold change of DE genes
    de.tlfc = all.tlfc[which(Zg==1)]
    # true log fold change of EE genes
    ee.tlfc = all.tlfc[which(Zg==0)]

    for(j in seq(along=Nreps1)) {

      ## LOG2 FOLD CHANGES
      # estimated log fold change of all genes
      all.elfc = elfcs[, j, i]
      # estimated log fold changes of EE genes
      ix.ee.lfc = which(Zg==0)
      ee.lfc = all.elfc[ix.ee.lfc]
      # estimated log fold change of DE genes
      ix.de.lfc = which(Zg==1)
      de.lfc = all.elfc[ix.de.lfc]
      # estimate mean squared error and absolute error
      all.error <- .lfc.evaluate(truth=all.tlfc, estimated=all.elfc)
      ee.error <- .lfc.evaluate(truth=ee.tlfc, estimated=ee.lfc)
      de.error <- .lfc.evaluate(truth=de.tlfc, estimated=de.lfc)
      error.est <- c(all.error, de.error, ee.error)
      lfc.error.mat[[j]][i,] = error.est

      ## SIZE FACTORS
      # true sf over all samples, center to mean=1
      tsf = tsfs[[j]][i, ]
      n.tsf = tsf*length(tsf)/sum(tsf)

      # estimated sf over all sample, center to mean=1
      esf = esfs[[j]][i,]
      n.esf = esf*length(esf)/sum(esf)

      # MAD of log fold change difference between estimated and true size factors
      lfc.nsf = log2(esf) - log2(tsf)
      mad.nsf = stats::mad(lfc.nsf)

      # error of estimation
      error.sf <- .fiterror.sf(estimated.sf = n.esf, true.sf = n.tsf)

      # ratio of estimated and true size factors per true group assignment
      t.design = t.designs[[j]][i,]
      ratio.sf <- .ratio.sf(estimated.nsf = n.esf,
                            true.nsf = n.tsf,
                            group=t.design)
      sf.res <- unlist(c(mad.nsf, error.sf, ratio.sf))
      names(sf.res) <- NULL

      sf.error.mat[[j]][i,] = sf.res

      ## CLUSTERING
      if (attr(simRes, 'Simulation') == 'Flow') {
        e.design = e.designs[[j]][i,]
        randindex <- mclust::adjustedRandIndex(t.design, e.design)
        clust.error.mat[[j]][i,] = randindex
      }

    }
  }

  output <- list(Log2FoldChange=lfc.error.mat,
                 SizeFactors=sf.error.mat,
                 Clustering=clust.error.mat,
                 sim.settings=simRes$sim.settings)

  if(isTRUE(timing)) {
    time.taken <- simRes$time.taken
    # create output objects
    time.taken.mat <- lapply(1:length(my.names), function(x) {
      data.frame(matrix(NA, nrow = length(rownames(time.taken[,1,]))+1,
                        ncol = 3, dimnames = list(c(rownames(time.taken[,1,]), "Total"),
                                                  c("Mean", "SD", "SEM")))
      )
    })
    names(time.taken.mat) <- my.names
    for(j in seq(along=Nreps1)) {
      tmp.time <- time.taken[,j,]
      Total <- colSums(tmp.time, na.rm = T)
      tmp.time <- rbind(tmp.time, Total)
      time.taken.mat[[j]][,"Mean"] <- rowMeans(tmp.time)
      time.taken.mat[[j]][,"SD"] <- matrixStats::rowSds(tmp.time)
      time.taken.mat[[j]][,"SEM"] <- matrixStats::rowSds(tmp.time)/sqrt(nsims)
    }
    output <- c(output, list(Timing=time.taken.mat))
  }

  # return object
  attr(output, 'Simulation') = attr(simRes, 'Simulation')
  return(output)
}

# EVALUATE DIFFERENTIAL EXPRESSION ----------------------------------------

#' @name evaluateDE
#' @aliases evaluateDE
#' @title Compute the confusion matrix-related quantities from simulation results
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes quantities of the confusion matrix for statistical power evaluation.
#' @usage evaluateDE(simRes,
#' alpha.type=c("adjusted","raw"),
#' MTC=c('BY', 'BH', 'holm', 'hochberg', 'hommel', 'bonferroni', 'Storey', 'IHW'),
#' alpha.nominal=0.1,
#' stratify.by=c("mean", "dispersion", "dropout", "lfc"),
#' strata.probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#' filter.by=c("none", "mean", "dispersion", "dropout"),
#' strata.filtered=1, target.by=c("lfc", "effectsize"), delta=0)
#' @param simRes The result from \code{\link{simulateDE}}.
#' @param MTC Multiple testing correction method to use. Available options are
#' 1) see \link[stats]{p.adjust.methods},
#' 2) Storey's qvalue see \link[qvalue]{qvalue} and
#' 3) Independent Hypothesis Weighting considering mean expression as covariate (see \link[IHW]{ihw}).
#' Default is \code{BY}, i.e. Benjamini-Yekutieli FDR correction method.
#' @param alpha.type A string to represent the way to call DE genes.
#'  Available options are \code{"adjusted"} i.e. applying multiple testing correction and
#'  \code{"raw"} i.e. using p-values. Default is \code{"adjusted"}.
#' @param alpha.nominal The nomial value of significance. Default is 0.1.
#' @param stratify.by A string to represent the way to stratify genes.
#' Available options are \code{"mean"}, \code{"dispersion"}, \code{"dropout"} and \code{"lfc"},
#' for stratifying genes by average expression levels, dispersion, dropout rates or estimated log2 fold changes.
#' @param strata.probs A vector specifying the probability values for sample quantiles of the strata. See \link[qvalue]{qvalue}.
#' @param filter.by A string to represent the way to filter genes.
#' This is used in conjunction with strata.filtered for gene filtering.
#' Available options are \code{"none"}, \code{"mean"}, \code{"dispersion"} and \code{"dropout"}.
#' \code{"none"} stands for no filtering, thus all genes will be considered.
#' \code{"mean"} stands for filtering based on average gene expression levels.
#' \code{"dispersion"} stands for filtering based on gene expression dispersion.
#' \code{"dropout"} stands for filtering based on dropout rates.
#' @param strata.filtered The strata to be filtered out in computing error matrix-related quantities.
#' Genes falling into these strata will be excluded. See "Details" for more description of gene filtering.
#' @param target.by A string to specify the method to define "biologically important" DE genes.
#' Available options are (1) \code{"lfc"}: interesting genes are defined by absolute log2 fold changes.
#' (2) \code{"effectsize"}: interesting genes are defined by
#' absolute log2 fold changes divided by the square root of 1/(mean+dispersion).
#' @param delta A threshold used for defining "biologically important" genes.
#' Genes with absolute log2 fold changes (when target.by is "lfc")
#' or effect sizes (when target.by is "effectsize") greater than this value
#' are deemed DE in error rates calculations. If \code{delta=0} then no threshold is applied. See "Details" for more description.
#' @return A list with the following entries:
#' \item{TN, TP, FP, FN, TNR, TPR, FPR, FNR, FDR}{3D array representing the number of true negatives, true positives, false positives,
#' false negatives and their proportions/rates as well as false discovery rate
#' for all simulation settings. The dimension of the arrays are nstrata * N * nsims.
#' Here nstrata is number of specified strata.
#' N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TN.marginal, TP.marginal, FP.marginal, FN.marginal}{Matrix representing the number of true negatives, true positives, false positives,
#'  false negatives for all simulation settings.
#'  The dimension of the matrices are N * nsims.
#'  Here N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TNR.marginal, TPR.marginal, FPR.marginal, FNR.marginal, FDR.marginal}{Matrix representing the marginal rates for all simulation settings.
#' The dimension of the matrices are N * nsims.}
#' \item{stratagenes, stratadiffgenes}{Number of genes per stratum and number of DE genes per stratum.}
#' \item{stratify.by}{The input "stratify.by".}
#' \item{strata}{The input strata.}
#' \item{n1,n2}{Sample sizes per group.
#' This is taken from the simulation options.}
#' \item{target.by}{The input method to define "biologically important" DE genes,
#' either by log fold change or effect size.}
#' \item{delta}{The input delta for biologically important genes.
#' If delta=0, all target.by will be considered.}
#' @details This is the main function to compute various power-related quantities,
#' using stratification and filtering.
#' \describe{
#' \item{Gene stratification}{We recommend to compute and visualize error rates (especially TPR)
#' conditional on expression characteristics like mean, dispersion and/or dropout rate.
#' It is likely that the power to detect DE genes is strongly dependent on
#' mean expression levels even though the magnitude of effect sizes is the same.
#' The stratified results will provide a more comprehensive power assessment and
#' better guide the investigators in experimental designs and analysis strategies.}
#' \item{Gene filtering}{Sometimes it is advisible to filter out some genes
#' (such as the ones with very low mean expression) before DE detection.
#' The filtering option here provides an opportunity to compare the rates before and after filtering.}
#' \item{Define biologically interesting genes}{We provide two options to define biologically interesting genes:
#' by absolute values of log fold changes or effect sizes
#' (absolute values of log fold changes divided by the square root of 1/(mean+dispersions)).
#' Genes with these quantities over a threshold are deemed interesting,
#' and the rate calculations are based on these genes.}
#' }
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}} for negative binomial parameters,
#' \code{\link{SimSetup}} and
#' \code{\link{DESetup}} for setting up simulation parameters and
#' \code{\link{simulateDE}} for simulating differential expression and
#' \code{\link{plotEvalDE}} for visualisation.
#' @examples
#' \dontrun{
#' data(kolodziejczk_simDE)
#' eval.de <- evaluateDE(simRes = kolodziejczk_simDE)
#' }
#' @rdname evaluateDE
#' @importFrom stats ecdf quantile p.adjust.methods p.adjust
#' @importFrom qvalue qvalue
#' @importFrom IHW ihw adj_pvalues
#' @export
evaluateDE <- function(simRes, alpha.type=c("adjusted","raw"),
                       MTC=c('BY', 'BH', 'holm', 'hochberg', 'hommel', 'bonferroni', 'Storey', 'IHW'),
                       alpha.nominal=0.1,
                       stratify.by=c("mean", "dispersion", "dropout", "lfc"),
                       strata.probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                       filter.by=c("none", "mean", "dispersion", "dropout"),
                       strata.filtered=1,
                       target.by=c("lfc", "effectsize"), delta=0) {

  alpha.type = match.arg(alpha.type)
  MTC = match.arg(MTC)
  stratify.by = match.arg(stratify.by)
  filter.by = match.arg(filter.by)
  target.by = match.arg(target.by)

  ## some general parameters
  Nreps1 = simRes$sim.settings$n1
  Nreps2 = simRes$sim.settings$n2
  ngenes = simRes$sim.settings$ngenes
  sim.opts = simRes$sim.settings
  DEids = simRes$sim.settings$DEid
  lfcs = simRes$sim.settings$pLFC
  tlfcs = lapply(1:length(lfcs), function(i) {lfcs[[i]]})
  nsims = simRes$sim.settings$nsims
  estmeans = simRes$sim.settings$means
  estdisps = simRes$sim.settings$dispersion
  estdropout = simRes$sim.settings$p0
  mu = simRes$mu
  disp = simRes$disp
  dropout = simRes$dropout
  elfc = simRes$elfc
  DEmethod = simRes$sim.settings$DEmethod
  pvalue = simRes$pvalue
  fdr = simRes$fdr

  ## calculate strata
  tmp.ecdf.mean = stats::ecdf(log2(estmeans+1))
  tmp.quantile.mean = stats::quantile(tmp.ecdf.mean, probs=strata.probs)
  strata.mean = unique(c(0,unname(tmp.quantile.mean),Inf))
  strata.mean = unique(round(strata.mean, digits=2))
  tmp.ecdf.disps = stats::ecdf(log2(estdisps))
  tmp.quantile.disps = stats::quantile(tmp.ecdf.disps, probs=strata.probs)
  strata.disps = unique(c(0,unname(tmp.quantile.disps),Inf))
  strata.disps = unique(round(strata.disps, digits=2))
  tmp.ecdf.drop = stats::ecdf(estdropout)
  tmp.quantile.drop = stats::quantile(tmp.ecdf.drop, probs=strata.probs)
  strata.drop = unique(c(0,unname(tmp.quantile.drop),1))
  strata.drop = unique(round(strata.drop, digits=2))
  tmp.ecdf.lfc = stats::ecdf(abs(unlist(tlfcs)))
  tmp.quantile.lfc = stats::quantile(tmp.ecdf.lfc, probs=strata.probs)
  strata.lfc = unique(c(-Inf,unname(tmp.quantile.lfc),Inf))
  strata.lfc = unique(round(strata.lfc, digits=2))

  ## initialize results
  ## determine dimension of results, for filtering
  if(stratify.by=='mean') {
    nr = length(strata.mean) - 1
  }
  if(stratify.by=='dispersion') {
    nr = length(strata.disps) - 1
  }
  if(stratify.by=='dropout') {
    nr = length(strata.drop) - 1
  }
  if(stratify.by=='lfc') {
    nr = length(strata.lfc) - 1
  }
  if(filter.by %in% c("mean", "dispersion", "dropout")) {
    nr = nr - strata.filtered
  }

  TP = TN =  FP = FN = TPR = TNR = FPR = FNR = FDR = xgrl = xgrld = array(NA,dim=c(nr,length(Nreps1), nsims))
  TP.marginal = TN.marginal = FP.marginal = FN.marginal = TPR.marginal = TNR.marginal = FPR.marginal = FNR.marginal = FDR.marginal = matrix(NA,length(Nreps1), nsims)

  ## loop over simulation and replicates
  for(i in 1:nsims) {
    for(j in seq(along=Nreps1)) {
      Nrep1 = Nreps1[j]
      Nrep2 = Nreps2[j]
      ## get DE flags.
      DEid = DEids[[i]]
      lfc = lfcs[[i]]
      Zg = Zg2 = rep(0, ngenes)
      Zg[DEid] = 1
      ## find target (interesting) genes
      if(delta == 0) {
        Zg2 = Zg
      }
      if(!delta == 0) {
        if(target.by == "lfc") {
          ix = abs(lfc) > delta
        } else if (target.by == "effectsize") {
          effectsize = lfc / sqrt(1/(log2(mu[,,i])+log2(disp[,,i])))
          ix = abs(effectsize) > delta
        }
        Zg2[ix] = 1
      }

      ### STRATIFICATION
      ## calculate stratificaton
      # mean
      X.bar1 = mu[,j,i]
      ix.keep.mean = which(!is.na(X.bar1))
      xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
      xgrd.mean = cut(log2(X.bar1[DEid]+1), strata.mean)
      # dispersion
      X.disp1 = disp[,j,i]
      ix.keep.disps = which(!is.na(X.disp1))
      xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps)
      xgrd.disps = cut(log2(X.disp1[DEid]), strata.disps)
      # dropout
      X.drop1 = dropout[,j,i]
      ix.keep.drop = which(!is.na(X.drop1))
      xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
      xgrd.drop = cut(X.drop1[DEid], strata.drop)
      # lfc
      X.lfc1 = abs(elfc[,j,1])
      ix.keep.lfc = which(!is.na(X.lfc1))
      xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
      xgrd.lfc = cut(X.lfc1[DEid], strata.lfc)

      ### FILTERING
      ## calculate filtering
      ## stratify by mean
      if(stratify.by == "mean") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.mean = ix.keep.mean[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean[-strata.filt.mean])
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean[-strata.filt.mean])
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.mean = ix.keep.mean[!(xgr.disps %in% lev.disps[strata.filt.disps])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.mean = ix.keep.mean[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean)
        }
        if(filter.by == "none") {
          ix.keep.mean = ix.keep.mean
          xgr.mean = xgr.mean
          xgrd.mean = xgrd.mean
        }
      }
      ## stratify by dispersion
      if(stratify.by == "dispersion") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.disps = ix.keep.disps[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps)
          xgrd.disps = cut(log2(X.disp1[(ix.keep.disps && DEid)]), strata.disps)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.disps = ix.keep.disps[!(xgr.disps %in% lev.disps[strata.filt.disps])]
          # recut
          xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps[-strata.filt.disps])
          xgrd.disps = cut(log2(X.disp1[(ix.keep.disps && DEid)]), strata.disps[-strata.filt.disps])
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.disps = ix.keep.disps[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.disps = cut(X.disp1[ix.keep.disps], strata.disps)
          xgrd.disps = cut(X.disp1[(ix.keep.disps && DEid)], strata.disps)
        }
        if(filter.by == "none") {
          ix.keep.disps = ix.keep.disps
          xgr.disps = xgr.disps
          xgrd.disps = xgrd.disps
        }
      }
      ## stratify by dropout
      if(stratify.by == "dropout") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.drop = ix.keep.drop[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.drop = ix.keep.drop[!(xgr.disps %in% lev.mean[strata.filt.disps])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.drop = ix.keep.drop[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop[-strata.filt.drop])
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop[-strata.filt.drop])
        }
        if(filter.by == "none") {
          ix.keep.drop = ix.keep.drop
          xgr.drop = xgr.drop
          xgrd.drop = xgrd.drop
        }
      }
      ## stratify by lfc
      if(stratify.by == "lfc") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.lfc = ix.keep.lfc[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.lfc = ix.keep.lfc[!(xgr.disps %in% lev.mean[strata.filt.disps])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.lfc = ix.keep.lfc[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "none") {
          ix.keep.lfc = ix.keep.lfc
          xgr.lfc = xgr.lfc
          xgr.lfc = xgr.lfc
        }
      }

      ### SET STRATIFICATION
      if(stratify.by == "mean") {
        strata = strata.mean
        xgr = xgr.mean
        xgrd = xgrd.mean
        ix.keep = ix.keep.mean
      }
      if(stratify.by == "dispersion") {
        strata = strata.disps
        xgr = xgr.disps
        xgrd = xgrd.disps
        ix.keep = ix.keep.disps
      }
      if(stratify.by == "dropout") {
        strata = strata.drop
        xgr = xgr.drop
        xgrd = xgrd.drop
        ix.keep = ix.keep.drop
      }
      if(stratify.by == "lfc") {
        strata = strata.lfc
        xgr = xgr.lfc
        xgrd = xgrd.lfc
        ix.keep = ix.keep.lfc
      }

      ## get type I error alpha (pvalue or fdr output from testing)
      if(alpha.type == "raw") {
        if(DEmethod %in% c("edgeR-QL", "edgeR-LRT", "limma-voom", "limma-trend", "NBPSeq", "T-Test",
                           "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", "monocle", "DECENT",
                           "edgeR-zingeR", "edgeR-ZINB-WaVE", "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          x = pvalue[ix.keep,j,i]
          x[is.na(x)] = 1
        }
        if(DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod," only provides adjusted p-values."))
          x = fdr[ix.keep,j,i]
          x[is.na(x)] = 1
        }
      }
      if(alpha.type == "adjusted") {
        if(DEmethod %in% c("edgeR-QL", "edgeR-LRT", "limma-voom", "limma-trend", "NBPSeq", "T-Test",
                           "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", "monocle", "DECENT",
                           "edgeR-zingeR", "edgeR-ZINB-WaVE", "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          pval = pvalue[ix.keep,j,i]
          meanexpr = mu[ix.keep,j,i]
          if(MTC %in% stats::p.adjust.methods) {
            x = stats::p.adjust(pval, method = MTC)
            x[is.na(x)] = 1
          }
          if(MTC %in% "Storey") {
            tmp.p = pval[!is.na(pval)]
            tmp.q = qvalue::qvalue(p = tmp.p)$qvalues
            x = rep(NA, length(pval))
            x[!is.na(pval)] = tmp.q
            x[is.na(x)] = 1
          }
          if(MTC %in% "IHW") {
            in.dat = data.frame(pvalue = pval, meanexpr = meanexpr)
            tmp = IHW::ihw(pvalue ~ meanexpr, data = in.dat, alpha = alpha.nominal)
            x = IHW::adj_pvalues(tmp)
            x[is.na(x)] = 1
          }
        }
        if(DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod," only provides adjusted p-values."))
          x = fdr[ix.keep,j,i]
          x[is.na(x)] = 1
        }
      }

      ## update Zg flags after filtering
      Zg = Zg[ix.keep]
      Zg2 = Zg2[ix.keep]

      #  number of strata genes and diff strata genes in output table
      xgrl[,j,i] = table(xgr)
      xgrld[,j,i] = table(xgrd)

      ## calculate stratified power-related quantities
      error.mat = .error.matrix(p=x, p.crit=alpha.nominal, Zg=Zg, Zg2=Zg2, xgr=xgr)

      TP[,j,i] = error.mat$TP
      TN[,j,i] = error.mat$TN
      FP[,j,i] = error.mat$FP
      FN[,j,i] = error.mat$FN

      TP.marginal[j,i] = error.mat$TP.marginal
      TN.marginal[j,i] = error.mat$TN.marginal
      FP.marginal[j,i] = error.mat$FP.marginal
      FN.marginal[j,i] = error.mat$FN.marginal

      TPR[,j,i] = error.mat$TPR
      TNR[,j,i] = error.mat$TNR
      FPR[,j,i] = error.mat$FPR
      FNR[,j,i] = error.mat$FNR
      FDR[,j,i] = error.mat$FDR

      TPR.marginal[j,i] = error.mat$TPR.marginal
      TNR.marginal[j,i] = error.mat$TNR.marginal
      FPR.marginal[j,i] = error.mat$FPR.marginal
      FNR.marginal[j,i] = error.mat$FNR.marginal
      FDR.marginal[j,i] = error.mat$FDR.marginal

    }
  }

  output <- list(stratagenes=xgrl, stratadiffgenes=xgrld,
                 TN=TN, TP=TP, FP=FP, FN=FN,
                 TN.marginal=TN.marginal, TP.marginal=TP.marginal, FP.marginal=FP.marginal, FN.marginal=FN.marginal,
                 TNR=TNR, TPR=TPR, FPR=FPR, FNR=FNR,FDR=FDR,
                 TNR.marginal=TNR.marginal, TPR.marginal=TPR.marginal,
                 FPR.marginal=FPR.marginal, FNR.marginal=FNR.marginal,
                 FDR.marginal=FDR.marginal,
                 ## below are input parameters:
                 alpha.type=alpha.type, MTC=ifelse(alpha.type=="adjusted", MTC, "not applicable"), alpha.nominal=alpha.nominal,
                 stratify.by=stratify.by, strata=strata, strata.levels=levels(xgr),
                 target.by=target.by, n1=Nreps1, n2=Nreps2, delta=delta)

  return(output)
}

# EVALUATE ROCR -----------------------------------------------------------

#' @name evaluateROC
#' @aliases evaluateROC
#' @title Receiver operator characteristics from simulation results
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes receiver operator characteristics (ROC) and area under the curve values (AUC).
#' @usage evaluateROC(simRes,
#' alpha.type=c("adjusted","raw"),
#' MTC=c('BY', 'BH', 'Storey', 'IHW',
#' 'holm', 'hochberg', 'hommel', 'bonferroni'),
#' alpha.nominal = 0.1,
#' target.by=c("lfc", "effectsize"),
#' delta=0)
#' @param simRes The result from \code{\link{simulateDE}}.
#' @param alpha.type A string to represent the way to call DE genes.
#'  Available options are \code{"adjusted"} i.e. applying multiple testing correction and
#'  \code{"raw"} i.e. using p-values. Default is \code{"adjusted"}.
#' @param MTC Multiple testing correction method to use. Available options are
#' 1) see \link[stats]{p.adjust.methods},
#' 2) Storey's qvalue see \link[qvalue]{qvalue} and
#' 3) Independent Hypothesis Weighting considering mean expression as covariate (see \link[IHW]{ihw}).
#' Default is \code{BY}, i.e. Benjamini-Yekutieli FDR correction method.
#' @param alpha.nominal The nomial value of significance. Default is 0.1.
#' @param target.by A string to specify the method to define "biologically important" DE genes.
#' Available options are (1) \code{"lfc"}: interesting genes are defined by absolute log2 fold changes.
#' (2) \code{"effectsize"}: interesting genes are defined by
#' absolute log2 fold changes divided by the square root of 1/(mean+dispersion).
#' @param delta A threshold used for defining "biologically important" genes.
#' Genes with absolute log2 fold changes (when target.by is "lfc")
#' or effect sizes (when target.by is "effectsize") greater than this value
#' are deemed DE in error rates calculations.
#' If \code{delta=0} then no threshold is applied. See "Details" for more description.
#' @return A list with the following entries:
#' \item{ROCR-Prediction}{The output of \code{\link[ROCR]{prediction}} split up by group comparisons into lists.}
#' \item{ROCR-Performance}{The output of \code{\link[ROCR]{performance}} calculating TPR and FPR for ROC curve.
#' This is also split up by group comparisons into lists.}
#' \item{ROCR-AUC}{The output of \code{\link[ROCR]{performance}} calculating AUC of ROC curve.
#' split up by group comparisons into lists.}
#' \item{n1,n2}{Sample sizes per group.
#' This is taken from the simulation options.}
#' \item{target.by}{The input method to define "biologically important" DE genes,
#' either by log fold change or effect size.}
#' \item{delta}{The input delta for biologically important genes.
#' If delta=0, all target.by will be considered.}
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}} for negative binomial parameters,
#' \code{\link{SimSetup}} and
#' \code{\link{DESetup}} for setting up simulation parameters and
#' \code{\link{simulateDE}} for simulating differential expression.
#' @examples
#' \dontrun{
#' data(kolodziejczk_simDE)
#' eval.roc <- evaluateROC(simRes = kolodziejczk_simDE)
#' }
#' @rdname evaluateROC
#' @importFrom stats ecdf quantile p.adjust.methods p.adjust
#' @importFrom qvalue qvalue
#' @importFrom IHW ihw adj_pvalues
#' @importFrom ROCR prediction performance
#' @importFrom zoo rollmean
#' @importFrom iCOBRA calculate_performance COBRAData
#' @export
evaluateROC <- function(simRes, alpha.type=c("adjusted","raw"),
                        MTC=c('BY', 'BH', 'Storey', 'IHW',
                              'holm', 'hochberg', 'hommel', 'bonferroni'),
                        alpha.nominal = 0.1,
                        target.by=c("lfc", "effectsize"),
                        delta=0) {
  alpha.type = match.arg(alpha.type)
  MTC = match.arg(MTC)
  target.by = match.arg(target.by)
  Nreps1 = simRes$sim.settings$n1
  Nreps2 = simRes$sim.settings$n2
  ngenes = simRes$sim.settings$ngenes
  sim.opts = simRes$sim.settings
  DEids = simRes$sim.settings$DEid
  lfcs = simRes$sim.settings$pLFC
  elfcs = simRes$elfc
  nsims = simRes$sim.settings$nsims
  DEmethod = simRes$sim.settings$DEmethod
  pvalue = simRes$pvalue
  fdr = simRes$fdr
  mu = simRes$mu
  my.names = paste0(Nreps1, " vs ", Nreps2)
  Truths = Predictions = array(NA, dim = c(length(Nreps1),
                                           ngenes, nsims))
  predObjs = rocObjs = prcObjs = rocaucObjs = prcaucObjs = matObjs = f1Objs = stats::setNames(replicate(length(Nreps1),
                                                                                                        NULL), my.names)
  perfCOBRA = stats::setNames(replicate(length(Nreps1), NULL),
                              my.names)
  perfCOBRA <- lapply(1:length(perfCOBRA), function(x) {
    perfCOBRA[[x]] = vector("list", nsims)
  })
  for (i in 1:nsims) {
    for (j in seq(along = Nreps1)) {
      Nrep1 = Nreps1[j]
      Nrep2 = Nreps2[j]
      elfc = as.numeric(elfcs[, j, i])
      DEid = DEids[[i]]
      lfc = lfcs[[i]]
      Zg = Zg2 = rep(0, ngenes)
      Zg[DEid] = 1
      if (delta == 0) {
        Zg2 = Zg
      }
      if (!delta == 0) {
        if (target.by == "lfc") {
          ix = abs(lfc) > delta
        }
        else if (target.by == "effectsize") {
          effectsize = lfc/sqrt(1/(log2(mu[, , i]) +
                                     log2(disp[, , i])))
          ix = abs(effectsize) > delta
        }
        Zg2[ix] = 1
      }
      if (alpha.type == "raw") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", "T-Test",
                            "limma-voom", "limma-trend", "NBPSeq", "DESeq2",
                            "ROTS", "MAST", "scde", "BPSC", "scDD", "monocle",
                            "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE",
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          x = pvalue[, j, i]
          x[is.na(x)] = 1
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod,
                         " only provides adjusted p-values."))
          x = fdr[, j, i]
          x[is.na(x)] = 1
        }
      }
      if (alpha.type == "adjusted") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", "T-Test",
                            "limma-voom", "limma-trend", "NBPSeq", "DESeq2",
                            "ROTS", "MAST", "scde", "BPSC", "scDD", "monocle",
                            "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE",
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          pval = pvalue[, j, i]
          meanexpr = mu[, j, i]
          if (MTC %in% stats::p.adjust.methods) {
            x = stats::p.adjust(pval, method = MTC)
            x[is.na(x)] = 1
          }
          if (MTC %in% "Storey") {
            tmp.p = pval[!is.na(pval)]
            tmp.q = qvalue::qvalue(p = tmp.p)$qvalues
            x = rep(NA, length(pval))
            x[!is.na(pval)] = tmp.q
            x[is.na(x)] = 1
          }
          if (MTC %in% "IHW") {
            in.dat = data.frame(pvalue = pval, meanexpr = meanexpr)
            tmp = IHW::ihw(pvalue ~ meanexpr, data = in.dat,
                           alpha = alpha.nominal)
            x = IHW::adj_pvalues(tmp)
            x[is.na(x)] = 1
          }
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod,
                         " only provides adjusted p-values."))
          x = fdr[, j, i]
          x[is.na(x)] = 1
        }
      }
      Truths[j, , i] = Zg2
      Predictions[j, , i] = x
      cobradata <- iCOBRA::COBRAData(pval = data.frame(Sim = pval, row.names = paste0("G", 1:ngenes)),
                                     padj = data.frame(Sim = x, row.names = paste0("G", 1:ngenes)),
                                     score = data.frame(Sim = elfc, row.names = paste0("G", 1:ngenes)),
                                     truth = data.frame(status = Zg2, logFC = lfc, expr = meanexpr, row.names = paste0("G", 1:ngenes)))
      invisible(capture.output(perfCOBRA[[j]][[i]] <- suppressMessages(iCOBRA::calculate_performance(cobradata,
                                                                                                     binary_truth = "status", cont_truth = "logFC",
                                                                                                     thrs = seq(from = 0.01, to = 0.2, by = 0.01),
                                                                                                     svalthrs = c(0.01, 0.05, 0.1), splv = "none",
                                                                                                     maxsplit = 4, onlyshared = FALSE, thr_venn = 0.05,
                                                                                                     type_venn = "adjp", topn_venn = 100, rank_by_abs = TRUE,
                                                                                                     prefer_pval = TRUE))))
    }
  }
  for (j in seq(along = Nreps1)) {
    predObjs[[j]] <- ROCR::prediction(predictions = Predictions[j, , ],
                                      labels = Truths[j, , ],
                                      label.ordering = c(1, 0))
    rocObjs[[j]] <- ROCR::performance(predObjs[[j]], "tpr", "fpr")
    rocaucObjs[[j]] <- ROCR::performance(predObjs[[j]], "auc")
    prcObjs[[j]] <- ROCR::performance(predObjs[[j]], "prec", "rec")
    y.values <- lapply(1:nsims, function(i) {
      x <- prcObjs[[j]]@x.values[[i]]
      y <- prcObjs[[j]]@y.values[[i]]
      x1 <- x[!is.na(x) & !is.na(y)]
      y1 <- y[!is.na(x) & !is.na(y)]
      id <- order(x1)
      sum(diff(x1[id]) * zoo::rollmean(y1[id], 2))
    })
    prcaucObjs[[j]] <- rocaucObjs[[j]]
    prcaucObjs[[j]]@y.values <- y.values
    prcaucObjs[[j]]@y.name <- c("Area under the Precision-Recall ROC curve")
    matObjs[[j]] <- ROCR::performance(predObjs[[j]], "mat")
    f1Objs[[j]] <- ROCR::performance(predObjs[[j]], "f")
  }
  output <- list(`ROCR-Prediction` = predObjs,
                 `ROCR-TPRvsFPR` = rocObjs,
                 `ROCR-PrecvsRec` = prcObjs,
                 `ROCR-MAT` = matObjs,
                 `ROCR-AUC` = rocaucObjs,
                 `PRC-AUC` = prcaucObjs,
                 `PRC-Fscore` = f1Objs,
                 `iCOBRA-Performance` = perfCOBRA,
                 alpha.type = alpha.type,
                 MTC = ifelse(alpha.type == "adjusted", MTC, "not applicable"),
                 target.by = target.by,
                 n1 = Nreps1,
                 n2 = Nreps2,
                 delta = delta)
  return(output)
}





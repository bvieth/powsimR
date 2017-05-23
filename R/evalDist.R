####################################################
#### evaluateDist
####################################################
#' @name evaluateDist
#' @aliases evaluateDist
#' @title Model Diagnostics for RNAseq Data Set
#' @description With this function, the user can determine goodness of fit for each gene.
#' @usage evaluateDist(cnts, RNAseq, ncores=1, nsims=1,
#' frac.genes=1, min.meancount=1, min.libsize=1000)
#' @param cnts is a count matrix (row=gene, column=sample). Pprovide the measurements of one group only, e.g. the control group.
#' @param RNAseq Character vector for "singlecell" or "bulk".
#' @param ncores, number of cores for parallel computing, default is 1.
#' @param nsims Number of simulations for MC p-value calculations, default is 1.
#' @param frac.genes The fraction of genes to calculate goodness of fit statistics, default is 1, i.e. for all genes.
#' @param min.meancount The minimum raw mean count per gene, default is 1.
#' @param min.libsize The minimum raw read counts per sample, default is 1000.
#' @return A list object with the results of goodness of fit and estimated parameters.
#' @examples
#' \dontrun{
#' ## simulating read count matrix
#' ngenes <- 10000
#' ncells <- 100
#' ## NB genes with dropouts:
#' ZNB.genes <- 2^rgamma(ngenes*0.45, 2, 2)
#' ## NB genes with high expression:
#' NB.genes <- 2^runif(ngenes*0.45, 9, 12)
#' ## Poisson genes:
#' P.genes <- 2^runif(ngenes*0.1, 3, 6)
#' ## all means:
#' true.means <- c(ZNB.genes, NB.genes, P.genes)
#' sf.values <- rnorm(ncells, mean=1, sd=0.1)
#' sf.means <- outer(true.means, sf.values, '*')
#' true.dispersions <- 3 + 100/true.means[1:round(ngenes*0.9)]
#' ## count matrix:
#' cnts.NB <- matrix(rnbinom(ngenes*0.9*ncells,
#' mu=sf.means[1:round(ngenes*0.9)],
#' size=1/true.dispersions),
#' ncol=ncells)
#' cnts.P <- matrix(rpois(ngenes*0.1*ncells,
#' lambda=sf.means[round(ngenes*0.9)+1:ngenes]),
#' ncol=ncells)
#' cnts <- rbind(cnts.NB, cnts.P)
#' ## evaluate distribution fitting:
#' evaldistres <- evaluateDist(cnts=cnts, RNAseq='singlecell',
#' ncores=1, nsims=1, frac.genes=0.25,
#' min.meancount=0.1, min.libsize=1000)
#' ## plot the evaluation results:
#' plotEvalDist(evalDist=evaldistres, annot=TRUE)
#' }
#' @author Beate Vieth
#' @rdname evaluateDist
#' @importFrom edgeR DGEList calcNormFactors cpm.DGEList estimateDisp glmFit zscoreNBinom
#' @importFrom scater sizeFactors newSCESet
#' @importFrom stats model.matrix glm logLik AIC dnbinom dpois fitted
#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl vuong
#' @importFrom fitdistrplus fitdist
#' @importFrom NBGOF nb.gof.v
#' @importFrom gamlss.dist ZIP ZINBI
#' @importFrom nonnest2 vuongtest
#' @export
evaluateDist <- function(cnts, RNAseq, ncores=1, nsims=1, frac.genes=1, min.meancount=1, min.libsize=1000) {
  options(stringsAsFactors = F)

  # kick out very small libs and lowly expressed genes
  fullS <- colSums(cnts) > min.libsize
  DetectG <- rowMeans(cnts) >= min.meancount
  cnts <- cnts[DetectG,fullS]

  # fill in pseudonames if missing
  if(is.null(rownames(cnts))) {
    rownames(cnts) <- paste0("G", 1:nrow(cnts))
  }
  if(is.null(colnames(cnts))) {
    colnames(cnts) <- paste0("S", 1:ncol(cnts))
  }

  if (RNAseq=="bulk") {
      dge <- edgeR::DGEList(counts=cnts, group=rep(1, ncol(cnts)))
      dge <- edgeR::calcNormFactors(dge, method='TMM')
      # calculate normalized counts (if norm lib is true does not sum up to 1 million)
      out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
      # estimate dispersions
      dge <- edgeR::estimateDisp(y=dge)
      # design.mat
      design.mat <- matrix(1,ncol(dge$counts),1)
      rownames(design.mat) <- colnames(dge$counts)
      colnames(design.mat) <- "Intercept"
      # apply edgeR glm fit
      fit.edgeR <- edgeR::glmFit(dge, design = design.mat)
  }

  if (RNAseq=="singlecell") {
    # make sceset and calculate size factors
    sce <- suppressWarnings(.scran.calc(cnts = cnts))
    # kick out negative size factor samples
    sf <- scater::sizeFactors(sce)
    sce2 <- suppressWarnings(scater::newSCESet(countData=data.frame(cnts[,sf>0])))
    scater::sizeFactors(sce2) <- sf[sf>0]
    # convert to edgeR object
    dge <- .convertToedgeR(sce2)
    dds <- .convertToDESeq(sce2)
    # calculate normalised counts
    out.normcounts <- DESeq2::counts(dds, normalized=T)
    # calculate CPM values
    out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
    # estimate dispersions
    dge <- edgeR::estimateDisp(y=dge, robust=T)
    # design.mat
    design.mat <- matrix(1,ncol(dge$counts),1)
    rownames(design.mat) <- colnames(dge$counts)
    colnames(design.mat) <- "Intercept"
    # apply edgeR glm fit
    fit.edgeR <- edgeR::glmFit(dge, design = design.mat)
  }

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
  headers.7 <- c( "GOF_Poisson", "GOF_NB2")
  headers.8 <- c("pvD", "pvP")
  headers_tmp4 <- paste(rep(headers.7, each=2), rep(headers.8, times=2), sep="_")
  headers_tmp5 <- c('LRT_standard_NBPoisson', 'Vuong_standard_ZPoisson', 'Vuong_standard_ZNB', 'Vuong_standard_ZNBZPoisson')
  headersgof <- c(headers_tmp1, headers_tmp2, headers_tmp3, headers_tmp4, headers_tmp5)

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

    print(paste0("Estimation for gene ", i, " out of ", length(geneid),". A total of ", nrow(dge$counts), " genes are available for estimation."))

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
      # fit NBGOF poisson
      tmp.fit.gofpois <- try(NBGOF::nb.gof.v(y=dge$counts[i,], x=design.mat, lib.sizes=edgeR::getOffset(y=dge), sim=nsims, model = "Poisson", ncores=ncores), silent = TRUE)
      # fit glm negative binomial
      tmp.fit.nbinom <- try(MASS::glm.nb(dge$counts[i,] ~ 1 + offset(edgeR::getOffset(y=dge))), silent = TRUE)
      # fit negative binomial with fitdistrplus
      tmp.fit.nbinom.wo <- try(fitdistrplus::fitdist(data = dge$counts[i,], 'nbinom'), silent = TRUE)
      # fit NBGOF NB2
      tmp.fit.gofnb2 <- try(NBGOF::nb.gof.v(y=dge$counts[i,], x=design.mat, lib.sizes=edgeR::getOffset(y=dge), sim=nsims, model = "NB2", method="ML", ncores=ncores), silent = TRUE)
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
    # NBGOF poisson
    if(inherits(tmp.fit.gofpois, 'try-error')) gof.res[i, grep("^GOF_Poisson", colnames(gof.res))] <- NA
    else{gof.res[i, grep("^GOF_Poisson", colnames(gof.res))] <- c(tmp.fit.gofpois$new.pval, tmp.fit.gofpois$pear.pval.1)}
    # NBGOF NB2
    if(inherits(tmp.fit.gofnb2, 'try-error')) gof.res[i, grep("^GOF_NB2", colnames(gof.res))] <- NA
    else{gof.res[i, grep("^GOF_NB2", colnames(gof.res))] <- c(tmp.fit.gofnb2$new.pval, tmp.fit.gofnb2$pear.pval.1)}
  }

  # fill in edgeR results
  edgeR.res <- cbind(edgeR.gof.stats, edgeR.aic)
  colnames(edgeR.res) <- c('edgeRglm_gofstat', 'edgeRglm_gofdf', 'edgeRglm_gofpval', 'edgeRglm_aic')
  rownames(edgeR.res) <- rownames(dge$counts)
  # observed zeros
  ObservedZeros <- data.frame(ObsZero=rowSums(dge$counts==0), Dropout=rowMeans(dge$counts==0))
  rownames(ObservedZeros) <- rownames(dge$counts)

  # list result object:
  # goodness of fit statistics per model
  # estimated parameters
  # observed zeros
  # setup settings for estimation
  return(list(edgeR_res= edgeR.res, GOF_res=gof.res, EST_res=estimate.res, ObservedZeros=ObservedZeros, evalDistsetup=c(ncores=ncores, nsims=nsims, frac.genes=frac.genes, min.meancount=min.meancount, min.libsize=min.libsize)))
}

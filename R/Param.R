####################################################
#### insilicoNBParam:
####################################################
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
#' sf=function(x) rnorm(n=x, mean=1, sd=0.1), RNAseq='bulk')
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

  res <- list(means = means, dispersion = dispersion, p0 = dropout, sf = sf, RNAseq = RNAseq, estFramework = "in silico")
  attr(res, 'param.type') <- "insilico"

  return(res)
}

####################################################
#### estimateNBParam:
####################################################
#' @name estimateNBParam
#' @aliases estimateNBParam
#' @title Estimate simulation parameters
#' @description This function estimates and returns parameters needed for power simulations assuming a negative binomial read count distribution.
#' @usage estimateNBParam(countData, cData=NULL, design=NULL,
#' RNAseq=c("bulk", "singlecell"),
#' estFramework=c('edgeR', 'DESeq2', 'MatchMoments'),
#' sigma=1.96)
#' @param countData is a count \code{matrix}. Rows correspond to genes, columns to samples.
#' @param cData A \code{data.frame} with at least a single column. Rows of \code{colData} \strong{must correspond} to columns of \code{countData}. Default is \code{NULL}, i.e. the dispersion estimation is 'blind' to sample information, see \code{\link[DESeq2]{varianceStabilizingTransformation}}.
#' @param design A \code{formula} which expresses how the counts for each gene depend on the variables in colData. Designs with multiple covariates and/or interactions are supported. Default is \code{NULL}, i.e. no design information is considered.
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param estFramework is a character value: "edgeR", "DESeq2" and "MatchMoments".
#' "edgeR" or "DESeq2" employs the edgeR or DESeq2 style mean and dispersion estimation, respectively. For details, please consult \code{\link[edgeR]{estimateDisp}} and \code{\link[DESeq2]{estimateDispersions}}.
#' "MatchMoments" employs moments matching technique for of mean, dispersion and dropout estimation.
#' @param sigma The variability band width for mean-dispersion loess fit defining the prediction interval for read cound simulation. Default is 1.96, i.e. 95\% interval. For more information see \code{\link[msir]{loess.sd}}.
#' @return List with the following vectors
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
#' \item{estFramework}{The estimation framework for NB parameters.}
#' \item{sigma}{The width of the variability band.}
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
#' mu=sf.means, size=1/true.dispersions),
#' ncol=ncells)
#' ## estimating negative binomial parameters
#' estparam <- estimateNBParam(cnts, RNAseq = 'singlecell',
#' estFramework = 'MatchMoments', sigma=1.96)
#' plotNBParam(estparam)
#'
#' ## simulating bulk RNA-seq experiment
#' ngenes <- 10000
#' nsamples <- 10
#' true.means <- 2^rnorm(ngenes, mean=8, sd=2)
#' true.dispersions <- rgamma(ngenes, 2, 6)
#' sf.values <- rnorm(nsamples, mean=1, sd=0.1)
#' sf.means <- outer(true.means, sf.values, '*')
#' cnts <- matrix(rnbinom(ngenes*nsamples,
#' mu=sf.means, size=1/true.dispersions),
#' ncol=nsamples)
#' ## estimating negative binomial parameters
#' estparam <- estimateNBParam(cnts, RNAseq = 'bulk',
#' estFramework = 'MatchMoments', sigma=1.96)
#' plotNBParam(estparam)
#' }
#' @author Beate
#' @rdname estimateNBParam
#' @export
estimateNBParam <- function(countData, cData=NULL, design=NULL,
                     RNAseq=c("bulk", "singlecell"),
                     estFramework=c('edgeR', 'DESeq2', 'MatchMoments'),
                     sigma=1.96) {
  if (RNAseq == 'bulk') {
    if (estFramework == 'edgeR') { #use edgeR
      res = .getDist.edgeR.bulk(countData, cData, design, sigma)
    }
    if (estFramework == 'DESeq2') { #use DESeq2
      res = .getDist.DESeq2.bulk(countData, cData, design, sigma)
    }
    if (estFramework == 'MatchMoments') { #use MatchMoments
      message("Design information will be ignored! Please provide count measurements of one group only, e.g. the control group")
      res = .getDist.MatchMoments.bulk(countData, sigma)
    }
  }

  if (RNAseq == 'singlecell') {
    if (estFramework == 'edgeR') { #use edgeR
      res = .getDist.edgeR.sc(countData, cData, design, sigma)
    }
    if (estFramework == 'DESeq2') { #use DESeq2
      res = .getDist.DESeq2.sc(countData, cData, design, sigma)
    }
    if (estFramework == 'MatchMoments') { #use MatchMoments method
      message("Design information will be ignored! Please provide count measurements of one group only, e.g. the control group.")
      res = .getDist.MatchMoments.sc(countData, sigma)
    }
  }
  res2 <- c(res, list(RNAseq = RNAseq, estFramework = estFramework, sigma=sigma))
  attr(res2, 'param.type') <- "estimated"
  # inform users
  if (RNAseq=='bulk' && c(res$grand.dropout>0.5 || is.null(res$p0.cut) )) {message('The provided bulk data has frequent dropouts. Consider using single cell estimation framework and testing.')}

  return(res2)
}

####################################################
#### estParam:   edgeR
####################################################

#' @importFrom edgeR DGEList calcNormFactors cpm.default estimateDisp cpm.default
#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict
.getDist.edgeR.bulk <- function(countData=countData, cData=cData, design=design, sigma=sigma){

  # kick out empty samples
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  estS <- sum(fullS, na.rm = T)
  DetectG <- rowSums(countData) > 0
  countData <- countData[DetectG,fullS]
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }

  if (is.null(cData)) {
    cData <- NULL
  } else {
    cData = as.data.frame(cData, stringsAsFactors = F)
    cData = cData[fullS,]
  }

  dge = edgeR::DGEList(counts = countData,lib.size = colSums(countData), samples = cData)
  dge = edgeR::calcNormFactors(dge)

  libsize = dge$samples$lib.size
  names(libsize) = rownames(dge$samples)

  # dropout component
  if (!is.null(design)) {
    modelmat <- stats::model.matrix(design, data = cData)
    dge2 <- edgeR::estimateDisp(y = dge, design = modelmat, verbose = F)
  }
  if (is.null(design)) {
    dge2 = edgeR::estimateDisp(y = dge)
  }

  mu = base::rowMeans(dge2$counts / dge2$samples$norm.factors)
  lmu = log2(mu + 1)

  nsamples = ncol(dge2$counts)
  counts0 = dge2$counts == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples

  #meanp0fit
  # meanp0fit = suppressWarnings(drc::drm(p0 ~ lmu, fct = LL.2(), type = "binomial"))

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- stats::runif(1000, min(lmu),  max(lmu))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  # expressed component

  # kick out lowly expressed genes
  min.cpm <- as.numeric(edgeR::cpm.default(5, mean(dge$samples$lib.size)))
  keep <- rowSums(edgeR::cpm.default(countData) > min.cpm) >= 2
  estG <- sum(keep, na.rm = T)
  dge3 <- edgeR::DGEList(counts = countData[keep,],lib.size = colSums(countData[keep,]), samples = cData)
  dge3 <- edgeR::calcNormFactors(dge3)

  # mean, dispersion and size
  if (!is.null(design)) {
    modelmat = stats::model.matrix(~design, data = cData)
    dge3 = edgeR::estimateDisp(y = dge3, design = modelmat, verbose = F)
    mu2 = base::rowMeans(dge3$counts / dge3$samples$norm.factors)
    phi.g = dge3$tagwise.dispersion
    phi.c = dge3$common.dispersion
    size = 1 / phi.g
  }
  if (is.null(design)) {
    dge3 = edgeR::estimateDisp(y = dge3, verbose = F)
    mu2 = base::rowMeans(dge3$counts / dge3$samples$norm.factors)
    phi.g = dge3$tagwise.dispersion
    phi.c = dge3$common.dispersion
    size = 1 / phi.g
  }

  lmu2 = log2(mu2 + 1)
  ldisp = log2(phi.g)
  lsize = log2(size)
  sf = dge3$samples$norm.factors
  names(sf) = rownames(dge3$samples)
  p0 = p0[keep]

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu2, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu2, nsigma = sigma)

  list(seqDepth = libsize, means = mu2, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)
}

#' @importFrom edgeR DGEList calcNormFactors cpm.default estimateDisp cpm.default
#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict
#' @importFrom scater sizeFactors newSCESet
.getDist.edgeR.sc <- function(countData=countData, cData=cData, design=design, sigma=sigma){

  # kick out empty samples and lowly expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  DetectG <- rowMeans(countData) > 0
  countData <- countData[DetectG,fullS]
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }

  # scran size factors
  sce <- .scran.calc(cnts = countData)

  # kick out zero size factor samples
  sf <- scater::sizeFactors(sce)
  sce2 <- suppressWarnings(scater::newSCESet(countData = data.frame(countData[,sf > 0])))
  scater::sizeFactors(sce2) <- sf[sf > 0]
  estS <- sum(sf > 0, na.rm = T)

  # convert to edgeR object
  dge1 <- .convertToedgeR(sce2)

  if (is.null(cData)) {
    cData <- NULL
    dge <- dge1
  }else {
    cData <- as.data.frame(cData, stringsAsFactors = F)
    cData <- cData[sf > 0,]
    dge <- edgeR::DGEList(counts = dge1$counts, lib.size = dge1$samples$lib.size,
                                  norm.factors = dge1$samples$norm.factors, samples = cData)
  }

  libsize = dge$samples$lib.size
  names(libsize) = rownames(dge$samples)

  # dropout and expressed component
  if (!is.null(design)) {
    modelmat = stats::model.matrix(design, data = cData)
    dge2 = edgeR::estimateDisp(y = dge, design = modelmat, verbose = F)
    mu2 = rowMeans(dge2$counts / dge2$samples$norm.factors)
    phi.g = dge2$tagwise.dispersion
    phi.c = dge2$common.dispersion
    size = 1 / phi.g
  }
  if (is.null(design)) {
    dge2 = edgeR::estimateDisp(y = dge, verbose = F)
    mu2 = rowMeans(dge2$counts / dge2$samples$norm.factors)
    phi.g = dge2$tagwise.dispersion
    phi.c = dge2$common.dispersion
    size = 1 / phi.g
  }

  estG <- length(mu2)
  # dropout
  nsamples = ncol(dge2$counts)
  counts0 = dge2$counts == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples

  # mean, dispersion and size
  lmu2 = log2(mu2 + 1)
  ldisp = log2(phi.g)
  lsize = log2(size)
  sf = dge2$samples$norm.factors
  names(sf) = rownames(dge$samples)

  #meanp0fit
  # meanp0fit = suppressWarnings(drc::drm(p0 ~ lmu2, fct = LL.2(), type = "binomial"))

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu2, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu2),  max(lmu2))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu2, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu2, nsigma = sigma)

  list(seqDepth = libsize, means = mu2, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)
}

####################################################
#### estParam: DESeq2
####################################################

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions dispersionFunction sizeFactors counts
#' @importFrom S4Vectors mcols
#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict
.getDist.DESeq2.bulk <- function(countData=countData, cData=cData, design=design, sigma=sigma) {

  # kick out empty samples
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  estS <- sum(fullS, na.rm = T)
  DetectG <- rowSums(countData) > 0
  countData <- countData[DetectG,fullS]
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }
  if (!is.null(cData) && !rownames(cData)==colnames(countData)) {
    rownames(cData) <- colnames(countData)
  }

  if (is.null(cData)) {
    group.df <- data.frame(group = rep("A", ncol(countData)), row.names=colnames(countData), stringsAsFactors = T)
    design <- ~1
  }
  if (!is.null(cData)) {
    tmp <- cData[fullS,]
    group.df <- data.frame(tmp, stringsAsFactors = T, row.names=colnames(countData))
    colnames(group.df) <- colnames(cData)
  }

  dds <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = group.df,
                                        design = design,
                                        tidy = FALSE, ignoreRank = FALSE))
  dds  <- DESeq2::estimateSizeFactors(dds)
  dds  <- DESeq2::estimateDispersions(dds, quiet = T)
  libsize <- as.vector(colSums(DESeq2::counts(dds)))
  names(libsize) <- colnames(counts(dds))

  # dropout component
  mu <- as.vector(S4Vectors::mcols(dds)[, "baseMean"])
  lmu <- log2(mu + 1)
  nsamples = ncol(counts(dds, normalize = F))
  counts0 = counts(dds, normalize = F) == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples

  #meanp0fit
  # meanp0fit <- suppressWarnings(drc::drm(p0 ~ lmu, fct = LL.2(), type = "binomial"))

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu),  max(lmu))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  # expressed component

  # kick out lowly expressed genes
  idx <- rowSums( counts(dds, normalized = TRUE) >= 5 ) >= 1
  p0 <- p0[idx]
  dds2 <- dds[idx,]
  dds2  <- DESeq2::estimateSizeFactors(dds2)
  dds2  <- DESeq2::estimateDispersions(dds2, quiet = T)

  mindisp <- which(S4Vectors::mcols(dds2)$dispGeneEst > 1e-4)
  estG <- length(mindisp)
  mu2 <- as.vector(S4Vectors::mcols(dds2)[mindisp, "baseMean"])
  phi.g <- as.vector(S4Vectors::mcols(dds2)[mindisp, "dispGeneEst"])
  p0 <- p0[mindisp]
  phi.c <- mean(phi.g)
  size <- 1 / phi.g
  lsize = log2(size)
  lmu2 = log2(mu2 + 1)
  ldisp = log2(phi.g)

  sf <- DESeq2::sizeFactors(dds2)

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu2, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu2, nsigma = sigma)

  list(seqDepth = libsize, means = mu2, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)
}

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions dispersionFunction sizeFactors counts
#' @importFrom S4Vectors mcols
#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict
#' @importFrom scater newSCESet sizeFactors counts
.getDist.DESeq2.sc <- function(countData=countData, cData=cData, design=design, sigma=sigma) {

  # kick out empty samples and lowly expressed genes
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  DetectG <- rowMeans(countData) > 0
  countData <- countData[DetectG,fullS]
  cData <- cData[fullS,]
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }
  if (!is.null(cData) && !rownames(cData)==colnames(countData)) {
    rownames(cData) <- colnames(countData)
  }

  # scran size factors
  sce <- .scran.calc(cnts = countData)

  # kick out zero size factor samples
  sf <- scater::sizeFactors(sce)
  sce2 <- suppressWarnings(scater::newSCESet(countData = data.frame(countData[,sf > 0])))
  scater::sizeFactors(sce2) <- sf[sf > 0]
  countData <- scater::counts(sce2)

  if (is.null(cData)) {
    group.df <- data.frame(group = rep("A", ncol(countData)), row.names=colnames(countData), stringsAsFactors = T)
    design <- ~1
  }
  if (!is.null(cData)) {
    tmp <- cData[fullS,]
    group.df <- data.frame(tmp, stringsAsFactors = T, row.names=colnames(countData))
    colnames(group.df) <- colnames(cData)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = group.df,
                                        design = design,
                                        tidy = FALSE, ignoreRank = FALSE)
  DESeq2::sizeFactors(dds)  <- scater::sizeFactors(sce2)
  dds  <- DESeq2::estimateDispersions(dds, quiet = T)
  libsize <- as.vector(colSums(DESeq2::counts(dds)))
  names(libsize) <- colnames(DESeq2::counts(dds))
  sf <- DESeq2::sizeFactors(dds)

  # NB parameters and dropout
  mindisp <- which(S4Vectors::mcols(dds)$dispGeneEst > 1e-4)
  mu <- as.vector(S4Vectors::mcols(dds)[mindisp, "baseMean"])
  phi.g <- as.vector(S4Vectors::mcols(dds)[mindisp, "dispGeneEst"])
  phi.c <- mean(phi.g)
  size <- 1 / phi.g
  lsize = log2(size)
  lmu <- log2(mu + 1)
  ldisp <- log2(phi.g)
  nsamples = ncol(DESeq2::counts(dds, normalize = F))
  counts0 = DESeq2::counts(dds, normalize = F) == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples

  estS <- length(sf)
  estG <- length(mu)

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu),  max(lmu))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  #meanp0fit
  # meanp0fit = suppressWarnings(drc::drm(p0 ~ lmu, fct = LL.2(), type = "binomial"))

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  list(seqDepth = libsize, means = mu, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)
}

####################################################
#### estParam: MatchMoments
####################################################

#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict median
.getDist.MatchMoments.bulk <- function(countData = countData, sigma=sigma){

  totalS <- ncol(countData)
  totalG <- nrow(countData)
  nsamples = dim(countData)[2]

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }

  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  if (any(nn0 == 1)) {
    countData = countData[nn0 > 1, ]
    nn0 = nn0[nn0 > 1]
    counts0 = countData == 0
  }
  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  seqDepth <- colSums(countData)
  names(seqDepth) <- colnames(countData)
  sf <- seqDepth/stats::median(seqDepth)
  X2 <- sweep(countData, 2, sf, FUN = "/")
  names(sf) <- colnames(countData)

  # the negative binomial portion
  mu = rowSums((!counts0) * X2)/nn0
  s2 = rowSums((!counts0) * (X2 - mu) ^ 2)/(nn0 - 1)
  size = mu ^ 2 / (s2 - mu + 1e-04)
  size = ifelse(size > 0, size, NA)
  p0 = (nsamples - nn0)/nsamples
  mu = mu[!is.na(size)]
  p0 = p0[!is.na(size)]
  size = size[!is.na(size)]
  phi.g = 1/size
  phi.c = mean(phi.g)
  lmu = log2(mu+1)
  ldisp = log2(phi.g)
  lsize = log2(size)

  estG <- length(mu)
  estS <- length(sf)

  #meanp0fit
  # meanp0fit = suppressWarnings(drc::drm(p0 ~ lmu, fct = LL.2(), type = "binomial"))

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu),  max(lmu))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])

  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  list(seqDepth = seqDepth, means = mu, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)
}

#' @importFrom msir loess.sd
# #' @importFrom drc drm LL.2
#' @importFrom cobs cobs
#' @importFrom stats model.matrix runif predict
#' @importFrom scater newSCESet sizeFactors
.getDist.MatchMoments.sc <- function(countData = countData, sigma=sigma){

  # kick out empty samples
  totalS <- ncol(countData)
  totalG <- nrow(countData)
  fullS <- colSums(countData) > 0
  DetectG <- rowSums(countData) > 0
  countData <- countData[DetectG,fullS]

  # fill in pseudonames if missing
  if (is.null(rownames(countData))) {
    rownames(countData) <- paste0("G", 1:nrow(countData))
  }
  if (is.null(colnames(countData))) {
    colnames(countData) <- paste0("S", 1:ncol(countData))
  }

  grand.dropout <- sum(countData == 0)/(nrow(countData)*ncol(countData))

  # scran size factors
  sce <- .scran.calc(cnts = countData)

  # kick out zero size factor samples
  sf <- scater::sizeFactors(sce)
  countData <- countData[,sf > 0]
  sce2 <- suppressWarnings(newSCESet(countData  = data.frame(countData)))
  scater::sizeFactors(sce2) <- sf[sf > 0]

  sf <- scater::sizeFactors(sce2)

  nsamples = dim(countData)[2]
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  if (any(nn0 == 1)) {
    countData = countData[nn0 > 1, ]
    nn0 = nn0[nn0 > 1]
    counts0 = countData == 0
  }

  seqDepth = colSums(countData)
  names(seqDepth) = colnames(countData)
  nf <- log(sf/seqDepth)
  nf <- exp(nf - mean(nf, na.rm = T))
  sf <- nf

  X2 = sweep(countData, 2, sf, FUN = "/")


  # the negative binomial portion
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
  estS <- length(sf)

  #meanp0fit
  # meanp0fit = suppressWarnings(drc::drm(p0 ~ lmu, fct = LL.2(), type = "binomial"))

  # cut for p0
  cobs.fit <- cobs::cobs(x = lmu, y = p0, constraint = 'decrease', nknots = 20, print.warn = F, print.mesg = F)
  cobs.sim <- runif(1000, min(lmu),  max(lmu))
  cobs.predict <- as.data.frame(stats::predict(cobs.fit, cobs.sim))
  cobs.predict[,"fit"] <- ifelse(cobs.predict$fit < 0, 0, cobs.predict[,"fit"])
  cobs.predict <- cobs.predict[cobs.predict[,'fit'] < 0.05,]
  p0.cut <- as.numeric(cobs.predict[which.max(cobs.predict[,'fit']), 'z'])


  # meansizefit
  meansizefit = msir::loess.sd(lsize ~ lmu, nsigma = sigma)

  # meandispfit
  meandispfit = msir::loess.sd(ldisp ~ lmu, nsigma = sigma)

  list(seqDepth = seqDepth, means = mu, dispersion = phi.g, common.dispersion = phi.c, size = size, p0 = p0, meansizefit = meansizefit, meandispfit = meandispfit, p0.cut = p0.cut, grand.dropout = grand.dropout, sf = sf, totalS = totalS, totalG = totalG, estS = estS, estG = estG)

}

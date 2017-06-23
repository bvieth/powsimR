
# edgeR -------------------------------------------------------------------

#' @importFrom edgeR DGEList estimateGLMRobustDisp glmFit glmLRT topTags
#' @importFrom stats model.matrix
.run.edgeRLRT <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # run DE testing
  design.mat <- stats::model.matrix( ~ dat$designs)
  dge <- edgeR::estimateGLMRobustDisp(y=dge, design = design.mat)
  end.time.params <- Sys.time()
  start.time.DE <- Sys.time()
  fit.edgeR <- edgeR::glmFit(dge, design = design.mat)
  lrt.edgeR <- edgeR::glmLRT(fit.edgeR)
  res.edgeR <- edgeR::topTags(lrt.edgeR, adjust.method="BH", n=Inf, sort.by = 'none')
  end.time.DE <- Sys.time()
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()
  ## construct results
  result <- data.frame(geneIndex=rownames(res.edgeR$table),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=res.edgeR$table$PValue,
                       fdr=rep(NA, nrow(res.edgeR$table)),
                       lfc=res.edgeR$table$logFC,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}

#' @importFrom edgeR DGEList estimateGLMRobustDisp topTags glmQLFit glmQLFTest
#' @importFrom stats model.matrix
.run.edgeRQL <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # run DE testing
  design.mat <- stats::model.matrix( ~ dat$designs)
  dge <- edgeR::estimateGLMRobustDisp(y=dge, design = design.mat)
  end.time.params <- Sys.time()
  start.time.DE <- Sys.time()
  fit.edgeR <- edgeR::glmQLFit(dge, design = design.mat)
  ql.edgeR <- edgeR::glmQLFTest(fit.edgeR)
  res.edgeR <- edgeR::topTags(ql.edgeR, adjust.method="BH", n=Inf, sort.by = 'none')
  end.time.DE <- Sys.time()
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()
  ## construct results
  result <- data.frame(geneIndex=rownames(res.edgeR$table),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=res.edgeR$table$PValue,
                       fdr=rep(NA, nrow(res.edgeR$table)),
                       lfc=res.edgeR$table$logFC,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# limma -------------------------------------------------------------------

#' @importFrom limma lmFit eBayes voom topTable
#' @importFrom edgeR DGEList
#' @importFrom stats model.matrix
.run.limma.voom <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- dat$p.DE
  design.mat <- stats::model.matrix( ~ dat$designs)
  v <- limma::voom(dge, design.mat, plot=FALSE)
  end.time.params <- Sys.time()
  start.time.DE <- Sys.time()
  fit <- limma::lmFit(object = v, design = design.mat)
  fit <- limma::eBayes(fit, proportion=p.DE, robust=T)
  resT <- limma::topTable(fit=fit, coef=2, number=Inf, adjust.method = "BH", sort.by = "none")
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=rownames(resT),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=resT$P.Value,
                       fdr=rep(NA, nrow(resT)),
                       lfc=resT$logFC,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}

#' @importFrom limma lmFit eBayes voom topTable
#' @importFrom edgeR DGEList cpm
#' @importFrom stats model.matrix
.run.limma.trend <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- dat$p.DE
  design.mat <- stats::model.matrix( ~ dat$designs)
  y <- new("Elist")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  end.time.params <- Sys.time()
  start.time.DE <- Sys.time()
  fit <- limma::lmFit(object = y, design = design.mat)
  fit <- limma::eBayes(fit, trend=TRUE, proportion=p.DE, robust=TRUE)
  resT <- limma::topTable(fit=fit, coef=2, number=Inf, adjust.method = "BH", sort.by = "none")
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=rownames(resT),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=resT$P.Value,
                       fdr=rep(NA, nrow(resT)),
                       lfc=resT$logFC,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# DESeq2 ------------------------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors DESeq sizeFactors results
#' @importFrom BiocParallel MulticoreParam
#' @importFrom scater sizeFactors
#' @importFrom stats model.matrix
.run.DESeq2 <- function(dat) {
  start.time.params <- Sys.time()
  coldat <- data.frame(design=factor(dat$designs))
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  # construct input object
  dds <- DESeq2::DESeqDataSetFromMatrix(dat$counts,
                                        coldat, ~design,
                                        tidy = FALSE, ignoreRank = FALSE)
  DESeq2::sizeFactors(dds) <- sf
  end.time.params <- Sys.time()
  # run DE testing
  start.time.DE <- Sys.time()
  if (is.null(dat$ncores)) {
    fit.DeSeq <- DESeq2::DESeq(dds, test="Wald", quiet = TRUE, parallel=FALSE)
  }
  if (!is.null(dat$ncores)) {
    fit.DeSeq <- DESeq2::DESeq(dds, test="Wald", quiet = TRUE, parallel=T, BPPARAM = BiocParallel::MulticoreParam(dat$ncores))
  }
  res.DeSeq <- DESeq2::results(fit.DeSeq)
  end.time.DE <- Sys.time()
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=rownames(res.DeSeq),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=res.DeSeq$pvalue,
                       fdr=rep(NA, nrow(dat$counts)),
                       lfc=res.DeSeq$log2FoldChange,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}

# ROTS --------------------------------------------------------------------

#' @importFrom edgeR DGEList cpm.DGEList
#' @importFrom ROTS ROTS
.run.ROTS <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # size factor normalised log2(CPM+1) values. Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  out.expr <- log2(out.cpm+1)
  end.time.params <- Sys.time()

  # run DE testing
  start.time.DE <- Sys.time()
  res <-  ROTS::ROTS(data = out.expr,
                     groups = factor(dat$designs) ,
                     B = 50,
                     K = floor(nrow(out.expr)/2) , progress=F)
  end.time.DE <- Sys.time()
  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=rownames(res$data),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=res$pvalue,
                       fdr=rep(NA, nrow(dat$counts)),
                       lfc=rep(NA, nrow(dat$counts)),
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# baySeq ------------------------------------------------------------------

#' @importFrom snow makeCluster stopCluster
#' @importMethodsFrom baySeq libsizes
#' @importFrom baySeq getPriors.NB getLikelihoods topCounts
.run.baySeq <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  # set multiple cores
  if(is.null(dat$ncores)) {
    cl <- NULL
  }
  if(!is.null(dat$ncores)) {
    cl <- snow::makeCluster(dat$ncores)
  }

  # make input data sets for baySeq
  replicates <- ifelse(dat$designs==-1, "A", "B")
  groups <- list(NDE = c(rep(1, length(dat$designs))),
                 DE = c(ifelse(dat$designs==-1, 1, 2)))
  CD <- new("countData", data = dat$counts, replicates = replicates, groups = groups)
  # fill in library size factors
  CD@sampleObservables$libsizes <- sf
  CD@annotation <- data.frame(name = rownames(dat$counts),
                              stringsAsFactors = F)
  # run prior estimation
  CD <- baySeq::getPriors.NB(CD,
                             samplesize = nrow(dat$counts),
                             estimation = "QL", cl = cl,
                             equalDispersions=TRUE, verbose=F)
  end.time.params <- Sys.time()
  start.time.DE <- Sys.time()
  # run likelihood ratio test
  CD <- baySeq::getLikelihoods(CD, cl = cl,
                               bootStraps = 10,
                               verbose = FALSE)
  # get test results
  res <- baySeq::topCounts(cD=CD, group="DE",
                           decreasing = FALSE,
                           number = Inf,
                           normaliseData = FALSE)
  res <- res[match(CD@annotation$name, res$annotation),]
  end.time.DE <- Sys.time()
  # free multiple cores
  if(!is.null(dat$ncores)) {
    snow::stopCluster(cl)
  }

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  # construct result data frame
  result=data.frame(geneIndex=res$annotation,
                    means=means,
                    dispersion=dispersion,
                    dropout=p0,
                    pval=rep(NA, nrow(dat$counts)),
                    fdr=res$FDR.DE,
                    stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# NOISeq ------------------------------------------------------------------

#' @importFrom NOISeq readData noiseqbio
#' @importFrom edgeR DGEList cpm.DGEList
.run.NOISeq <- function(dat) {
  start.time.params <- Sys.time()
  groups <- data.frame(Group=factor(dat$designs))
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct dge object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # size factor normalised log2(CPM+1) values. Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  end.time.params <- Sys.time()
  # make input data set
  in.noiseq <- NOISeq::readData(data = out.cpm, factors = groups)
  end.time.params <- Sys.time()

  # run DE detection
  start.time.DE <- Sys.time()
  calc.noiseq <- NOISeq::noiseqbio(in.noiseq,
                                   k = NULL,
                                   norm = "n",
                                   nclust = 15,
                                   plot = FALSE,
                                   factor="Group",
                                   conditions = NULL,
                                   lc = 0, r = 50, adj = 1.5,
                                   a0per = 0.9, filter = 0)
  res <- calc.noiseq@results[[1]]
  res$fdr <- 1-res$prob
  end.time.DE <- Sys.time()

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  # construct result data frame
  result=data.frame(geneIndex=rownames(res),
                    means=means,
                    dispersion=dispersion,
                    dropout=p0,
                    pval=rep(NA, nrow(res)),
                    fdr=res$fdr,
                    lfc=rep(NA, nrow(res)),
                    stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# DSS ---------------------------------------------------------------------

#' @importFrom DSS newSeqCountSet estNormFactors estDispersion waldTest
#' @importFrom splines ns
#' @importFrom edgeR DGEList
.run.DSS <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  # make input data set
  designs <- ifelse(dat$designs==-1, 0, 1)
  cd <- dat$counts
  rownames(cd) <- NULL
  colnames(cd) <- NULL
  seqData <- DSS::newSeqCountSet(counts = cd, designs = designs)
  seqData@normalizationFactor <- sf
  seqData <- DSS::estDispersion(seqData)
  end.time.params <- Sys.time()

  # run DE detection
  start.time.DE <- Sys.time()
  res.dss <- suppressWarnings(DSS::waldTest(seqData = seqData,
                                            sampleA = 0, sampleB = 1))
  res.dss <- res.dss[order(res.dss$geneIndex),]
  pval <- res.dss$pval
  end.time.DE <- Sys.time()

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()


  # construct result data frame
  result=data.frame(geneIndex=rownames(dat$counts),
                    means=means,
                    dispersion=dispersion,
                    dropout=p0,
                    pval=pval,
                    fdr=rep(NA, nrow(dat$counts)),
                    lfc=res$lfc,
                    stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# EBSeq -------------------------------------------------------------------

#' @importFrom EBSeq MedianNorm EBTest
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom scater sizeFactors
.run.EBSeq <- function(dat) {
  start.time.params <- Sys.time()
  groups <- data.frame(Group=factor(dat$designs))
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  end.time.params <- Sys.time()
  # run DE detection
  start.time.DE <- Sys.time()
  calc.ebseq <- suppressMessages(
    EBSeq::EBTest(Data = dat$counts, NgVector = NULL,
                  Conditions = factor(dat$designs),
                  sizeFactors = sf,
                  maxround = 20,
                  Pool = F, NumBin = 1000,
                  ApproxVal = 10^-10, Alpha = NULL,
                  Beta = NULL, PInput = NULL,
                  RInput = NULL, PoolLower = .25,
                  PoolUpper = .75, Print = F,
                  Qtrm = 1,QtrmCut=0))
  fdr <- 1-calc.ebseq$PPDE
  end.time.DE <- Sys.time()

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=rownames(dat$counts),
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=rep(NA, nrow(dat$counts)),
                       fdr=fdr,
                       lfc=rep(NA, nrow(dat$counts)),
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# NBPSeq ------------------------------------------------------------------

#' #' @importFrom NBPSeq nbp.test
#' #' @importFrom edgeR DGEList calcNormFactors
#' .run.NBPSeq <- function(dat) {
#'
#'   dge <- edgeR::DGEList(counts=dat$counts, group=factor(dat$designs))
#'   if (dat$RNAseq=="bulk") {
#'     start.time.params <- Sys.time()
#'     dge <- edgeR::calcNormFactors(dge, method='TMM')
#'     end.time.params <- Sys.time()
#'     start.time.DE <- Sys.time()
#'     res <- NBPSeq::nbp.test(counts=dge$counts, grp.ids=dat$designs, grp1=-1, grp2=1, norm.factors = dge$samples$norm.factors, lib.sizes = colSums(dge$counts), model.disp = "NBQ", print.level = 0)
#'     end.time.DE <- Sys.time()
#'   }
#'   if (dat$RNAseq=="singlecell") {
#'     # make sceset and calculate size factors
#'     start.time.params <- Sys.time()
#'     # make sceset and calculate size factors
#'     sce <- .scran.calc(cnts = dat$counts)
#'     dge <- .convertToedgeR(sce)
#'     dge$samples$group <- factor(dat$designs)
#'     end.time.params <- Sys.time()
#'     start.time.DE <- Sys.time()
#'     res <- NBPSeq::nbp.test(counts=dge$counts, grp.ids=dat$designs, grp1=-1, grp2=1, norm.factors = dge$samples$norm.factors, lib.sizes = colSums(dge$counts), model.disp = "NBQ", print.level = 0)
#'     end.time.DE <- Sys.time()
#'   }
#'
#'   # mean, disp, dropout
#'   start.time.NB <- Sys.time()
#'   norm.counts <- dge$counts / dge$samples$norm.factors
#'   nsamples <- ncol(norm.counts)
#'   counts0 <- norm.counts == 0
#'   nn0 <- rowSums(!counts0)
#'   p0 <- (nsamples - nn0)/nsamples
#'   means = rowSums(norm.counts)/nsamples
#'   s2 = rowSums((norm.counts - means)^2)/(nsamples - 1)
#'   size = means^2/(s2 - means + 1e-04)
#'   size = ifelse(size > 0, size, NA)
#'   dispersion = 1/size
#'   end.time.NB <- Sys.time()
#'
#'   ## construct results
#'   result <- data.frame(geneIndex=rownames(dat$counts), means=means, dispersion=dispersion, dropout=p0, pval=res$pv.alues, fdr=rep(NA, nrow(dat$counts)), stringsAsFactors = F)
#'   time.taken.params <- difftime(end.time.params, start.time.params, units="mins")
#'   time.taken.DE <- difftime(end.time.DE, start.time.DE, units="mins")
#'   time.taken.NB <- difftime(end.time.NB, start.time.NB, units="mins")
#'   timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
#'   res <- list(result=result, timing=timing)
#'   return(res)
#' }


# MAST --------------------------------------------------------------------

#' @importFrom MAST FromMatrix zlm.SingleCellAssay lrTest
#' @importFrom S4Vectors mcols
#' @importFrom AnnotationDbi as.list
#' @importFrom edgeR DGEList calcNormFactors cpm.DGEList
#' @importFrom data.table data.table
#' @importFrom reshape2 melt
#' @importFrom parallel mclapply
.run.MAST <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # 1. size factor normalised log2(CPM+1) values. Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  out.expr <- log2(out.cpm+1)
  end.time.params <- Sys.time()

  # 2.: cell (sample ID, CDR, condition) and gene (gene name) annotation
  ids=colnames(out.expr)
  ngeneson=colSums(out.expr>0)
  cngeneson=ngeneson-mean(ngeneson)
  cond=factor(dat$designs)
  cdat <- data.frame(wellKey=ids,
                     ngeneson=ngeneson,
                     cngeneson=cngeneson,
                     condition=cond,
                     stringsAsFactors = F)
  fdat <- data.frame(primerid=rownames(out.expr),
                     stringsAsFactors = F)
  # 3.: construct MAST single cell assay
  sca <- MAST::FromMatrix(class = "SingleCellAssay",
                    exprsArray=out.expr,
                    cData = cdat,
                    fData = fdat)
  end.time.params <- Sys.time()
  # 4.: Model Fit
  start.time.DE <- Sys.time()
  if (!is.null(dat$ncores)) {
    options(mc.cores=dat$ncores)
  }

  zlm <- suppressMessages(
    MAST::zlm(~ condition + cngeneson,
              sca,
              method = "bayesglm",
              ebayes = TRUE,
              ebayesControl = list(method = "MLE", model = "H1"))
  )
  summaryZLM <- suppressMessages(summary(zlm))
  summaryDt <- summaryZLM$datatable
  # 5.: LRT
  lrt <- suppressMessages(MAST::lrTest(zlm, "condition"))
  # results table extraction
  res_gene <- data.table::data.table(reshape2::melt(lrt))
  res_gene_hurdle <- res_gene[metric=="Pr(>Chisq)" & test.type=="hurdle"]
  res_gene_hurdle <- data.frame(res_gene_hurdle, stringsAsFactors = F)
  res_gene_hurdle <- res_gene_hurdle[match(S4Vectors::mcols(sca)$primerid, res_gene_hurdle$primerid),]
  res_lfc_hurdle <- summaryDt[contrast=='condition' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)]
  res_lfc_hurdle <- data.frame(res_lfc_hurdle, stringsAsFactors = F)
  res_lfc_hurdle <- res_lfc_hurdle[match(S4Vectors::mcols(sca)$primerid, res_lfc_hurdle$primerid),]
  end.time.DE <- Sys.time()

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result <- data.frame(geneIndex=res$primerid,
                       means=means,
                       dispersion=dispersion,
                       dropout=p0,
                       pval=res_gene_hurdle$value,
                       fdr=rep(NA, nrow(dat$counts)),
                       lfc=res_lfc_hurdle$coef,
                       stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# scde --------------------------------------------------------------------

#' @importFrom scde scde.error.models scde.expression.prior scde.expression.difference
#' @importFrom stats pnorm
.run.scde <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))

  # make group vector
  groups <- factor(dat$designs)
  names(groups) <- colnames(counts)

  if(is.null(dat$ncores)) {
    ncores = 1
  }
  if(!is.null(dat$ncores)) {
    ncores = dat$ncores
  }

  # calculate error models
  o.ifm <- scde::scde.error.models(counts = dat$counts,
                                   groups = groups,
                                   n.cores = ncores,
                                   min.count.threshold = 1,
                                   threshold.segmentation = TRUE,
                                   save.crossfit.plots = FALSE,
                                   save.model.plots = FALSE,
                                   verbose = 0)
  # estimate gene expression prior
  o.prior <- scde::scde.expression.prior(models = o.ifm,
                                         counts = dat$counts,
                                         length.out = 400,
                                         show.plot = FALSE)
  end.time.params <- Sys.time()
  # run differential expression tests on all genes.
  start.time.DE <- Sys.time()
  ediff <- scde::scde.expression.difference(models=o.ifm,
                                            counts=dat$counts,
                                            prior=o.prior,
                                            groups = groups,
                                            n.cores = ncores,
                                            n.randomizations  =  100,
                                            verbose  =  0)
  pval <- 2 * (1 - pnorm(abs(ediff$Z)))
  end.time.DE <- Sys.time()

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  ## construct results
  result=data.frame(geneIndex=rownames(ediff),
                    means=means,
                    dispersion=dispersion,
                    dropout=p0,
                    pval=pval,
                    fdr=rep(NA, nrow(dat$counts)),
                    lfc=ediff$ce,
                    stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# BPSC --------------------------------------------------------------------

#' @importFrom BPSC BPglm
#' @importFrom gtools mixedsort
#' @importFrom edgeR DGEList cpm.DGEList
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats model.matrix
.run.BPSC <- function(dat) {
  start.time.params <- Sys.time()
  # run normalisation
  normData <- .norm.calc(normalisation=dat$normalisation,
                         countData=dat$counts,
                         spikeData=NULL,
                         spikeInfo=NULL,
                         Lengths=NULL,
                         MeanFragLengths=NULL,
                         NCores=ifelse(is.null(dat$ncores), 1, dat$ncores))
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(dat$counts))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = dat$counts,
                        lib.size = colSums(dat$counts),
                        norm.factors = nsf,
                        group = factor(dat$designs),
                        remove.zeros = FALSE)
  # size factor normalised log2(CPM+1) values. Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  exprmat <- out.cpm
  group <- dat$designs
  controlIDs <- which(group == -1)
  design.mat <- stats::model.matrix( ~ group)
  coef <- 2
  end.time.params <- Sys.time()

  if(!is.null(dat$ncores)) {
    start.time.DE <- Sys.time()
    cl <- parallel::makeCluster(dat$ncores)
    doParallel::registerDoParallel(cl)
    res <- BPglm(data = exprmat, controlIds = controlIDs, design = design.mat, coef = coef, useParallel=TRUE)
    parallel::stopCluster(cl)
    end.time.DE <- Sys.time()
  }
  if(is.null(dat$ncores)) {
    start.time.DE <- Sys.time()
    res <- BPSC::BPglm(data = exprmat, controlIds = controlIDs, design = design.mat, coef = coef, keepFit = TRUE)
    summary_BPglm <- summary(res)
    summaryDt <-  summary_BPglm$topTable
    res_lfc <- as.data.frame(summaryDt)
    res_lfc <- res_lfc[gtools::mixedsort(rownames(res_lfc)),]
    end.time.DE <- Sys.time()
  }

  # parameter estimation
  start.time.NB <- Sys.time()
  means <- rowMeans(normData$NormCounts)
  nsamples <- ncol(normData$NormCounts)
  s2 = rowSums((normData$NormCounts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  counts0 <- dat$counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  end.time.NB <- Sys.time()

  # construct result data frame
  result=data.frame(geneIndex=rownames(exprmat),
                    means=means,
                    dispersion=dispersion,
                    dropout=p0,
                    pval=res$PVAL,
                    fdr=rep(NA, nrow(dat$counts)),
                    lfc=res_lfc$Estimate,
                    stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params,
                                start.time.params,
                                units="mins")
  time.taken.DE <- difftime(end.time.DE,
                            start.time.DE,
                            units="mins")
  time.taken.NB <- difftime(end.time.NB,
                            start.time.NB,
                            units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# monocle -----------------------------------------------------------------

#' #' @importFrom monocle newCellDataSet differentialGeneTest
#' #' @importFrom VGAM tobit
#' #' @importFrom edgeR cpm.DGEList
#' #' @importFrom scater sizeFactors
#' #' @importFrom methods new
#' .run.monocle <- function(dat) {
#'   if (dat$RNAseq=="bulk") {
#'     stop("monocle is only for single cell RNAseq data analysis")
#'   }
#'   if (dat$RNAseq=="singlecell") {
#'     start.time.params <- Sys.time()
#'     # make sceset and calculate size factors
#'     sce <- .scran.calc(cnts = dat$counts)
#'     dge <- .convertToedgeR(sce)
#'     dge$samples$group <- factor(dat$designs)
#'     out.sf <- scater::sizeFactors(sce)
#'     out.sf[out.sf<0] <- min(out.sf[out.sf > 0])
#'     out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
#'   }
#'     # make annotated dataframes for monocle
#'     gene.dat <- data.frame(row.names = rownames(dge$counts), biotype=rep("protein_coding", nrow(dge$counts)), num_cells_expressed=rowSums(dge$counts>0))
#'     cell.dat <- data.frame(row.names=colnames(dge$counts), Group=dge$samples$group)
#'     fd <- new("AnnotatedDataFrame", data = gene.dat)
#'     pd <- new("AnnotatedDataFrame", data = cell.dat)
#'     ed <- out.cpm
#'     # construct cell data set
#'     cds <- monocle::newCellDataSet(cellData = ed, phenoData = pd, featureData = fd, expressionFamily = VGAM::tobit())
#'     end.time.params <- Sys.time()
#'
#'     # run the testing
#'     if(!is.null(dat$ncores)) {
#'       start.time.DE <- Sys.time()
#'       diff_test_res <- monocle::differentialGeneTest(cds, fullModelFormulaStr = "~Group", reducedModelFormulaStr = "~1", relative_expr = FALSE, cores = dat$ncores, verbose = FALSE)
#'     }
#'     if(is.null(dat$ncores)) {
#'       start.time.DE <- Sys.time()
#'       diff_test_res <- monocle::differentialGeneTest(cds, fullModelFormulaStr = "~Group", reducedModelFormulaStr = "~1", relative_expr = FALSE, cores = 1, verbose = FALSE)
#'     }
#'     res <- diff_test_res[match(rownames(dge$counts), rownames(diff_test_res)),]
#'     end.time.DE <- Sys.time()
#'
#'     # mean, disp, droput
#'     start.time.NB <- Sys.time()
#'     norm.counts <- dge$counts / dge$samples$norm.factors
#'     nsamples <- ncol(norm.counts)
#'     counts0 <- norm.counts == 0
#'     nn0 <- rowSums(!counts0)
#'     p0 <- (nsamples - nn0)/nsamples
#'     means = rowSums(norm.counts)/nsamples
#'     s2 = rowSums((norm.counts - means)^2)/(nsamples - 1)
#'     size = means^2/(s2 - means + 1e-04)
#'     size = ifelse(size > 0, size, NA)
#'     dispersion = 1/size
#'     end.time.NB <- Sys.time()
#'
#'     # construct result data frame
#'     result=data.frame(geneIndex=rownames(res), means=means, dispersion=dispersion, dropout=p0, pval=res$pval, fdr=rep(NA, nrow(res)), stringsAsFactors = F)
#'     time.taken.params <- difftime(end.time.params, start.time.params, units="mins")
#'     time.taken.DE <- difftime(end.time.DE, start.time.DE, units="mins")
#'     time.taken.NB <- difftime(end.time.NB, start.time.NB, units="mins")
#'     timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
#'     res <- list(result=result, timing=timing)
#'     return(res)
#' }


# scDD --------------------------------------------------------------------

#' @importFrom scDD scDD
#' @importFrom edgeR cpm.DGEList
#' @importFrom SummarizedExperiment SummarizedExperiment
.run.scDD <- function(dat) {
  if (dat$RNAseq=="bulk") {
    stop("scDD is only for single cell RNAseq data analysis")
  }
  if (dat$RNAseq=="singlecell") {
    start.time.params <- Sys.time()
    # make sceset and calculate size factors
    sce <- .scran.calc(cnts = dat$counts)
    dge <- .convertToedgeR(sce)
    dge$samples$group <- factor(dat$designs)
  }
  # size factor normalised CPM values.
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)

  # create input data
  exprmat <- out.cpm
  condition <- ifelse(dat$designs==-1, 1, 2)
  cell.dat <- data.frame(row.names=colnames(exprmat), condition=condition)
  SCdat <- SummarizedExperiment::SummarizedExperiment(assays=list('NormCounts'=exprmat), colData=cell.dat)
  # SCdat <- Biobase::ExpressionSet(assayData=exprmat, phenoData=as(cell.dat, "AnnotatedDataFrame"))
  end.time.params <- Sys.time()

  # DE testing
  if(!is.null(dat$ncores)) {
    start.time.DE <- Sys.time()
    res.tmp <- scDD::scDD(SCdat, prior_param = list(alpha = 0.1, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01), permutations = 0, testZeroes = FALSE, adjust.perms = FALSE, param = BiocParallel::MulticoreParam(dat$ncores), parallelBy = "Genes", condition = "condition")
    end.time.DE <- Sys.time()
  }
  if(is.null(dat$ncores)) {
    start.time.DE <- Sys.time()
    res.tmp <- scDD(SCdat, prior_param = list(alpha = 0.1, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01), permutations = 0, testZeroes = FALSE, adjust.perms = FALSE, parallelBy = "Genes", condition = "condition")
    end.time.params <- Sys.time()
  }
  res <- res.tmp$Genes

  # mean, disp, dropout
  start.time.NB <- Sys.time()
  norm.counts <- dge$counts / dge$samples$norm.factors
  nsamples <- ncol(norm.counts)
  counts0 <- norm.counts == 0
  nn0 <- rowSums(!counts0)
  p0 <- (nsamples - nn0)/nsamples
  means = rowSums(norm.counts)/nsamples
  s2 = rowSums((norm.counts - means)^2)/(nsamples - 1)
  size = means^2/(s2 - means + 1e-04)
  size = ifelse(size > 0, size, NA)
  dispersion = 1/size
  end.time.NB <- Sys.time()

  # construct result data frame
  result=data.frame(geneIndex=as.character(res$gene), means=means, dispersion=dispersion, dropout=p0, pval=res$nonzero.pvalue, fdr=rep(NA, nrow(dat$counts)), stringsAsFactors = F)
  time.taken.params <- difftime(end.time.params, start.time.params, units="mins")
  time.taken.DE <- difftime(end.time.DE, start.time.DE, units="mins")
  time.taken.NB <- difftime(end.time.NB, start.time.NB, units="mins")
  timing <- rbind(time.taken.params, time.taken.DE, time.taken.NB)
  res <- list(result=result, timing=timing)
  return(res)
}


# TODO: Do a system call since D3E is written in python

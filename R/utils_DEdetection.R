# NOTES -------------------------------------------------------------------

# Note that the function in scater gave negative values and
# when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!

# DE TOOLS WRAPPER --------------------------------------------------------

# normData is normalized counts
# countData is raw counts
# DEOpts is a list of designs, p.DE,
# spikeData is a spike-in count table
# spikeInfo is a table of molecules counts of spike-ins
# Lengths is a named vector of gene lengths
# MeanFragLengths is the average fragment length for paired end designs
# NCores is number of cores

.de.calc <- function(DEmethod,
                     normData,
                     countData,
                     Lengths,
                     MeanFragLengths,
                     DEOpts,
                     spikeData,
                     spikeInfo,
                     NCores,
                     verbose) {

  # classic testing
  if (DEmethod == "T-Test") {DERes = .run.TTest(normData=normData,
                                                countData=countData,
                                                Lengths=Lengths,
                                                MeanFragLengths=MeanFragLengths,
                                                DEOpts=DEOpts)}

  # methods developed for bulk
  if (DEmethod == "limma-voom") {DERes = .run.limma.voom(normData=normData,
                                                         countData=countData,
                                                         DEOpts=DEOpts,
                                                         verbose=verbose)}
  if (DEmethod == "limma-trend") {DERes = .run.limma.trend(normData=normData,
                                                           countData=countData,
                                                           DEOpts=DEOpts,
                                                           verbose=verbose)}
  if (DEmethod == "DESeq2") {DERes = .run.DESeq2(normData=normData,
                                                 countData=countData,
                                                 DEOpts=DEOpts,
                                                 NCores=NCores,
                                                 verbose=verbose)}

  return(DERes)
}


# limma -------------------------------------------------------------------

#' @importFrom limma lmFit eBayes voom topTable
#' @importFrom edgeR DGEList
#' @importFrom stats model.matrix
.run.limma.voom <- function(normData, countData, DEOpts, verbose) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))

  if (attr(normData, 'normFramework')  %in% c('sctransform')) {
    norm.cnts <- normData$RoundNormCounts
    ixx.valid <- rownames(countData) %in% rownames(norm.cnts)
    countData[ixx.valid, ] <- norm.cnts
  }

  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- DEOpts$p.DE
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  if (attr(normData, 'normFramework') %in% c('SCnorm', "Linnorm")) {
    scale.facts <- normData$scale.factors
    ixx.valid <- rownames(countData) %in% rownames(scale.facts)
    wgenes <- countData
    wgenes[1:nrow(wgenes), 1:ncol(wgenes)] <- NA
    wgenes[ixx.valid, ] <- scale.facts
    v <- limma::voom(dge, design.mat, plot=FALSE, weights = wgenes)
  }
  if (!attr(normData, 'normFramework')  %in% c('SCnorm', "Linnorm")) {
    v <- limma::voom(dge, design.mat, plot=FALSE)
  }
  fit <- limma::lmFit(object = v, design = design.mat)
  fit <- limma::eBayes(fit, proportion=p.DE, robust=TRUE)
  resT <- limma::topTable(fit=fit, coef=2, number=Inf, adjust.method = "BH", sort.by = "none")

  # construct results
  result <- data.frame(geneIndex=rownames(resT),
                       pval=resT$P.Value,
                       fdr=rep(NA, nrow(resT)),
                       lfc=resT$logFC,
                       stringsAsFactors = F)
  return(result)
}

#' @importFrom limma lmFit eBayes voom topTable
#' @importFrom edgeR DGEList cpm
#' @importFrom stats model.matrix
.run.limma.trend <- function(normData, countData, DEOpts, verbose) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))

  if (attr(normData, 'normFramework')  %in% c('sctransform')) {
    norm.cnts <- normData$RoundNormCounts
    ixx.valid <- rownames(countData) %in% rownames(norm.cnts)
    countData[ixx.valid, ] <- norm.cnts
  }

  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- DEOpts$p.DE
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  y <- methods::new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)

  if (attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    scale.facts <- normData$scale.factors
    ixx.valid <- rownames(countData) %in% rownames(scale.facts)
    wgenes <- countData
    wgenes[1:nrow(wgenes), 1:ncol(wgenes)] <- NA
    wgenes[ixx.valid, ] <- scale.facts
    fit <- limma::lmFit(object = y, design = design.mat, weights = wgenes)
  }
  if (!attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    fit <- limma::lmFit(object = y, design = design.mat)
  }

  fit <- limma::eBayes(fit, trend=TRUE, proportion=p.DE, robust=TRUE)
  resT <- limma::topTable(fit=fit, coef=2, number=Inf, adjust.method = "BH", sort.by = "none")

  # construct results
  result <- data.frame(geneIndex=rownames(resT),
                       pval=resT$P.Value,
                       fdr=rep(NA, nrow(resT)),
                       lfc=resT$logFC,
                       stringsAsFactors = F)
  return(result)
}

# DESeq2 ------------------------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors DESeq sizeFactors results
#' @importMethodsFrom DESeq2 sizeFactors
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats model.matrix
.run.DESeq2 <- function(normData, countData, DEOpts, NCores, verbose) {

  # construct input object
  groups <- ifelse(DEOpts$designs==min(DEOpts$designs), 1, 2)
  coldat <- data.frame(design=factor(groups))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  if (attr(normData, 'normFramework')  %in% c('sctransform')) {
    norm.cnts <- normData$RoundNormCounts
    ixx.valid <- rownames(countData) %in% rownames(norm.cnts)
    countData[ixx.valid, ] <- norm.cnts
  }

  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData,
                                                         coldat, ~design,
                                                         tidy = FALSE,
                                                         ignoreRank = FALSE))
  DESeq2::sizeFactors(dds) <- sf

  if (attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    scale.facts <- normData$scale.factors
    ixx.valid <- rownames(countData) %in% rownames(scale.facts)
    wgenes <- countData
    wgenes[1:nrow(wgenes), 1:ncol(wgenes)] <- NA
    wgenes[ixx.valid, ] <- scale.facts
    SummarizedExperiment::assays(dds)[["weights"]] <- wgenes
  }

  # run DE testing
  if (is.null(NCores)) {
    fit.DeSeq <- suppressMessages(DESeq2::DESeq(dds,
                                                test="Wald",
                                                quiet = verbose,
                                                parallel=FALSE))
  }
  if (!is.null(NCores)) {
    fit.DeSeq <- suppressMessages(DESeq2::DESeq(dds, test="Wald",
                                                quiet = verbose,
                                                parallel=TRUE,
                                                BPPARAM = BiocParallel::MulticoreParam(NCores)))
  }
  res.DeSeq <- DESeq2::results(fit.DeSeq)

  ## construct results
  result <- data.frame(geneIndex=rownames(res.DeSeq),
                       pval=res.DeSeq$pvalue,
                       fdr=rep(NA, nrow(countData)),
                       lfc=res.DeSeq$log2FoldChange,
                       stringsAsFactors = F)

  invisible(gc())

  return(result)
}

# T -TEST -----------------------------------------------------------------

#' @importFrom stats t.test
.run.TTest <- function(normData, countData, Lengths, MeanFragLengths, DEOpts) {

  if (attr(normData, 'normFramework')  %in% c('sctransform')) {
    norm.cnts <- normData$RoundNormCounts
    ixx.valid <- rownames(countData) %in% rownames(norm.cnts)
    countData[ixx.valid, ] <- norm.cnts
  }

  # 1. calculate size factor normalised expression values (ie CPM, TPM).
  out.expr <- .calculateExpr(countData=countData,
                             normData=normData,
                             Lengths=Lengths,
                             MeanFragLengths=MeanFragLengths)
  out.expr <- log2(out.expr+1)

  # 2. perform T test per gene
  cond <- factor(DEOpts$designs)
  idx <- seq_len(nrow(out.expr))
  names(idx) <- rownames(out.expr)
  ttest_res <- sapply(idx, function(i) {
    ttest <- stats::t.test(out.expr[i, ] ~ cond)
    pval <- ttest$p.value
    lfc <- as.numeric(log2(ttest$estimate[1]/ttest$estimate[2]))
    data.frame("pval"=pval, "lfc"=lfc)
  }, simplify=FALSE, USE.NAMES = TRUE)
  ttest_res <- do.call("rbind", ttest_res)
  ## construct results
  result <- data.frame(geneIndex=rownames(out.expr),
                       pval=ttest_res[,'pval'],
                       fdr=rep(NA, nrow(out.expr)),
                       lfc=ttest_res[,'lfc'],
                       stringsAsFactors = F)
  return(result)
}

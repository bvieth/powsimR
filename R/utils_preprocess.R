
# IMPUTATION WRAPPER ------------------------------------------------------

.impute.calc <- function(Impute,
                         countData,
                         spikeData,
                         batchData,
                         clustNumber,
                         NCores,
                         verbose) {
  if(Impute=='scImpute') {invisible(capture.output(
    FilterData <- suppressMessages(.scimpute.impute(countData = countData,
                                                    clustNumber = clustNumber,
                                                    NCores = NCores,
                                                    verbose = verbose))
  ))}
  if(Impute=='DrImpute') {FilterData <- .drimpute.impute(countData = countData,
                                                         verbose = verbose)}

  if(Impute=='SAVER') {FilterData <- .saver.impute(countData = countData,
                                                   NCores = NCores,
                                                   verbose = verbose)}


  if(Impute=='Seurat') {FilterData <- .seurat.impute(countData = countData,
                                                     verbose = verbose)}

  if(Impute=='scone') {FilterData <- .scone.impute(countData = countData,
                                                   spikeData = spikeData,
                                                   batchData = batchData,
                                                   NCores = NCores,
                                                   verbose = verbose)}

  if(Impute=='BISCUIT') {FilterData <- .biscuit.impute(countData = countData,
                                                       verbose = verbose)}

  ## no easy way of extracting cidr imputed values, C++ function call
  # if(Impute=='CIDR') {FilterData <- .cidr.impute(countData = countData,
  #                                                NCores = NCores)}

  # the authors do not recommend doing it and
  # it is applied after normalisation and transformation, not to help normalisaion!
  # if(Impute=='Linnorm') {FilterData <- .linnorm.impute(countData = countData)}

  # the current R implementation has a lot of coding errors,
  # and also the imputation gives weird results (all of the expr values turned into miniscule tiny values)
  # if(Impute=='MAGIC') {FilterData <- .magic.impute(countData = countData)}

  return(FilterData)
}

# imputation --------------------------------------------------------------

# adapted functions from scImpute so that input and output are R objects
.scimpute.impute <- function(countData, NCores, clustNumber, verbose) {

  Kcluster = clustNumber
  drop_thre = 0.5
  # fill in ncores
  ncores = ifelse(is.null(NCores), 1, NCores)

  # convert raw counts to log10 seq depth normalized counts
  raw_count = as.matrix(countData)
  totalCounts_by_cell = colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] = 1
  raw_count = sweep(raw_count,
                    MARGIN = 2,
                    totalCounts_by_cell/10^6,
                    FUN = "/")
  count_lnorm = log10(raw_count + 1.01)
  genenames = rownames(count_lnorm)
  cellnames = colnames(count_lnorm)

  res_imp = .imputation_model8(count = count_lnorm,
                              labeled = FALSE,
                              point = log10(1.01),
                              drop_thre = drop_thre,
                              Kcluster = Kcluster,
                              ncores = ncores)

  # backtransform to unlogged counts
  count_imp = res_imp$count_imp
  count_imp = 10^count_imp - 1.01
  rownames(count_imp) = genenames
  colnames(count_imp) = cellnames

  # remove seq depth normalization
  ImputeData = sweep(count_imp,
                     MARGIN = 2,
                     totalCounts_by_cell/10^6,
                     FUN = "*")

  # round to integer counts
  ImputeData = round(ImputeData)

  invisible(gc())

  # return imputed data
  return(ImputeData)
}

# DrImpute
#' @importFrom DrImpute preprocessSC DrImpute
#' @importFrom utils capture.output
.drimpute.impute <- function(countData, verbose) {

  if(isTRUE(verbose)) {
    # prefiltering of expression
    invisible(utils::capture.output(
      exdata <- DrImpute::preprocessSC(X=countData,
                                       min.expressed.gene = 0,
                                       min.expressed.cell = 2,
                                       max.expressed.ratio = 1,
                                       normalize.by.size.effect = FALSE)
    ))
    # read depth normalization and log transformation
    sf <- apply(exdata, 2, mean)
    npX <- t(t(exdata) / sf )
    lnpX <- log(npX+1)
    #imputation
    invisible(capture.output(
      lnpX_imp <- DrImpute::DrImpute(lnpX)
    ))
  }
  if(!isTRUE(verbose)) {
    # prefiltering of expression
    invisible(utils::capture.output(
      exdata <- suppressMessages(DrImpute::preprocessSC(X=countData,
                                                        min.expressed.gene = 0,
                                                        min.expressed.cell = 2,
                                                        max.expressed.ratio = 1,
                                                        normalize.by.size.effect = FALSE))
    ))
    # read depth normalization and log transformation
    sf <- apply(exdata, 2, mean)
    npX <- t(t(exdata) / sf )
    lnpX <- log(npX+1)
    #imputation
    invisible(capture.output(
      lnpX_imp <- suppressMessages(DrImpute::DrImpute(lnpX))
    ))
  }

  # backtransform to unlogged counts
  count_imp = exp(lnpX_imp) - 1
  # remove seq depth normalization
  ImputeData = sweep(count_imp,
                     MARGIN = 2,
                     sf,
                     FUN = "*")
  # round to integer counts
  ImputeData = round(ImputeData)
  rownames(ImputeData) = rownames(exdata)
  colnames(ImputeData) = colnames(exdata)

  invisible(gc())

  # return object
  return(ImputeData)
}

# SAVER
# this is one that apparently works but is really slow and needs cores!
# as an example: would need 30 minutes for 3500 genes in 500 cells on 10 cores
# so i changed the number of cells and genes used in predictions
#' @importFrom SAVER saver combine.saver
#' @importFrom utils capture.output
#' @importFrom doParallel registerDoParallel stopImplicitCluster
.saver.impute <- function(countData, NCores, verbose) {

  if(!is.null(NCores)) {
    doParallel::registerDoParallel(cores = NCores)
  }

  npred <- nrow(countData)/3 # number of genes for regression prediction

  # find genes with at least 50 % dropout and only do imputation for those
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highDrop <- p0 > 0.5

  d <- rownames(countData)[highDrop]
  max.g <- 2500
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/max.g))
  str(d1)
  saver.outs <- sapply(1:length(d1), function(s) {
    genes <- d1[[s]]
    genes.ind <- which(rownames(countData) %in% genes)
    invisible(utils::capture.output(
      tmp <- suppressMessages(
        SAVER::saver(x = countData,
                     do.fast = TRUE,
                     size.factor = NULL,
                     npred = npred,
                     pred.cells = NULL,
                     pred.genes = genes.ind,
                     pred.genes.only = TRUE,
                     null.model = FALSE)
    )))
    tmp
  }, simplify = F)

  if(length(saver.outs)==1) {
    saver.out <- saver.outs[[1]]
  }

  if(length(saver.outs)>1) {
    saver.out <- .combine.saver(saver.outs)
  }

  if(!is.null(NCores)) {
    doParallel::stopImplicitCluster()
  }

  # internally he uses normalized data to calculate the fits. these are used to predict expression.
  # the scaled size factors are used for normalisation of raw data in prediction, but this is equal to 1 if left to default.
  # therefore the predicted estimates are like raw imputed counts albeit still rounding needed.

  ImputeData = saver.out$estimate
  # round to integer counts
  ImputeData = as.matrix(round(ImputeData))

  dim(ImputeData)
  dim(countData)
  table(rownames(countData) == rownames(ImputeData))

  # combine ImputeData with countData since imputation
  # was only done for the genes with higher dropouts
  ImputeData <- ImputeData[order(match(rownames(ImputeData), rownames(countData))),]
  CombineData <- countData
  CombineData[rownames(countData) %in% rownames(ImputeData),] <- ImputeData

  rownames(CombineData) = rownames(countData)
  colnames(CombineData) = colnames(countData)

  invisible(gc())

  # return object
  return(CombineData)
}

# Seurat
#' @importFrom Seurat CreateSeuratObject ExpMean LogVMR
.seurat.impute <- function(countData, NCores, verbose) {

  # the original function returned Inf values for high expression genes,
  # so that I changed the function and added the variable genes to the seurat object

  # create input object
  seurat.obj <- Seurat::CreateSeuratObject(raw.data = countData,
                                           project = "SeuratProject",
                                           min.cells = 0,
                                           min.genes = 0,
                                           is.expr = 0,
                                           normalization.method = NULL,
                                           scale.factor = 10000,
                                           do.scale = FALSE,
                                           do.center = FALSE)

  # find genes with at least 50 % dropout and only do imputation for those
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highDrop <- p0 > 0.5

  # find variable genes
  seurat.obj <- .FindVarGenes(object = seurat.obj,
                              mean.function = Seurat::ExpMean,
                              dispersion.function = Seurat::LogVMR,
                              set.var.genes = TRUE,
                              x.low.cutoff = 0.1,
                              x.high.cutoff = Inf,
                              y.cutoff = 1,
                              y.high.cutoff = Inf,
                              num.bin = 20,
                              sort.results = TRUE)

  print(length(seurat.obj@var.genes))
  print(table(highDrop))

  # imputation
  seurat.obj <- .AddImputedScore(object = seurat.obj,
                                 genes.use = seurat.obj@var.genes,
                                 genes.fit = rownames(x = seurat.obj@data)[highDrop],
                                 gram = ifelse(length(seurat.obj@var.genes)<500, TRUE, FALSE))

  # Seurat uses a lasso fit and subsequent preidciton for imputation
  # resulting in negative estimates.
  # i changed these to zero.

  ImputeData = seurat.obj@imputed
  ImputeData[ImputeData<0] = 0
  # round to integer counts
  ImputeData = as.matrix(round(ImputeData))

  # combine ImputeData with countData since imputation
  # was only done for the genes with higher dropouts
  countData <- as.matrix(seurat.obj@raw.data)
  ImputeData <- ImputeData[order(match(rownames(ImputeData), rownames(countData))),]
  CombineData <- countData
  CombineData[rownames(countData) %in% rownames(ImputeData),] <- ImputeData

  rownames(CombineData) = rownames(countData)
  colnames(CombineData) = colnames(countData)

  invisible(gc())

  # return object
  return(CombineData)
}

# scone
#' @importFrom scone SconeExperiment estimate_ziber scone impute_expectation
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#' @importFrom S4Vectors SimpleList
#' @importFrom BiocParallel SerialParam MulticoreParam register bpparam
.scone.impute <- function(countData, spikeData, batchData, NCores, verbose) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    rowdata <- data.frame(negcon=c(rep(FALSE, nrow(countData)),
                                   rep(TRUE, nrow(spikeData))),
                          stringsAsFactors = F)
  }
  if(spike==FALSE) {
    cnts = countData
    rowdata <- data.frame(negcon=c(rep(FALSE, nrow(countData))),
                          stringsAsFactors = F)
  }

  if(!is.null(batchData)) {
    if(is.vector(batchData)) {
      cond <- as.factor(batchData)
    }
    if(!is.vector(batchData)) {
      cond <- as.factor(batchData[1,])
    }
    coldata <- data.frame(bio=cond)
  }
  if(is.null(batchData)) {
    cond <- as.factor(rep(1, ncol(cnts)))
    coldata <- data.frame(bio=cond)
  }
  if (!is.null(NCores)) {
    prll=BiocParallel::MulticoreParam(workers=NCores)
    BiocParallel::register(BPPARAM = prll, default=TRUE)
    bpparams = BiocParallel::bpparam("MulticoreParam")
  }
  if (is.null(NCores)) {
    bpparams = BiocParallel::SerialParam()
  }

  # create input
  se <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=cnts),
                             rowData=rowdata, colData=coldata)
  scone1 <- scone::SconeExperiment(se, which_bio=1L, which_negconeval=1L)


  # Simple FNR model estimation with SCONE::estimate_ziber
  counts0 = cnts == 0
  nn0 = rowSums(!counts0)
  nsamples = ncol(cnts)
  dropout = (nsamples - nn0)/nsamples
  posgenes = dropout < 0.15
  fnr_out = scone::estimate_ziber(x = cnts,
                                  bulk_model = TRUE,
                                  pos_controls = posgenes,
                                  maxiter = 50,
                                  em_tol = 10,
                                  verbose = verbose)

  # define scaling arguments
  scaling=list(none=identity) # Identity - do nothing

  ## Imputation List Argument
  imputation=list(expect=scone::impute_expectation) # Replace zeroes
  impute_args = list(p_nodrop = fnr_out$p_nodrop,
                     mu = exp(fnr_out$Alpha[1,]))

  scone2 <- scone::scone(x = scone1,
                         imputation = imputation,
                         impute_args = impute_args,
                         zero = "none",
                         scaling =  scaling,
                         k_ruv = 0,
                         k_qc = 0,
                         adjust_bio = "no",
                         adjust_batch = "no",
                         run = TRUE,
                         evaluate = FALSE,
                         eval_pcs = 3,
                         eval_proj = NULL,
                         eval_proj_args = NULL,
                         eval_kclust = NULL,
                         verbose = verbose,
                         stratified_pam = FALSE,
                         stratified_cor = FALSE,
                         stratified_rle = FALSE,
                         return_norm = "in_memory",
                         bpparam = bpparams)

  # extract assay
  ImputeData <- SummarizedExperiment::assays(scone2)[[1]]
  ImputeData <- ImputeData[!grepl(pattern="ERCC", rownames(ImputeData)),]

  # round to integer counts
  ImputeData = round(ImputeData)
  rownames(ImputeData) = rownames(countData)
  colnames(ImputeData) = colnames(countData)

  invisible(gc())

  # return object
  return(ImputeData)
}

# PREFILTER WRAPPER ---------------------------------------------------

.prefilter.calc <- function(Prefilter,
                            countData,
                            NCores) {
  if(Prefilter=='CountFilter') {FilterData <- .count.filter(countData = countData)}
  if(Prefilter=='FreqFilter') {FilterData <- .freq.filter(countData = countData)}

  return(FilterData)
}

# filtering ---------------------------------------------------------------

# apply stochastic filtering by Lun et al 2016 / Finak et al 2016
.count.filter <- function(countData) {
  highE<- rowMeans(countData) >= 0.2
  FilterData <- countData[highE,]
  return(FilterData)
}

# apply frequency filter: gene must have less than 80% dropouts
.freq.filter <- function(countData) {
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highE <- p0 < 0.8
  FilterData <- countData[highE,]
  return(FilterData)
}

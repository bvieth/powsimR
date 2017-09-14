
# PREPROCESSING WRAPPER ---------------------------------------------------

.preprocess.calc <- function(Preprocess,
                             countData,
                             NCores) {
  if(Preprocess=='scImpute') {invisible(capture.output(
    FilterData <- suppressMessages(.scimpute.calc(countData = countData,
                                 NCores = NCores))
  ))}
  if(Preprocess=='DrImpute') {FilterData <- .drimpute.calc(countData = countData)}
  if(Preprocess=='CountFilter') {FilterData <- .countfilter.calc(countData = countData)}
  if(Preprocess=='FreqFilter') {FilterData <- .freqfilter.calc(countData = countData)}

  return(FilterData)
}

# imputation --------------------------------------------------------------

# adapted functions from scImpute so that input and output are R objects
.scimpute.calc <- function(countData, NCores) {

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

  # estimate mixture model parameters
  estmixmodelparams = .get_mix_parameters(count = count_lnorm,
                                           point = log10(1.01),
                                           ncores = ncores)

  # impute dropouts
  count_imp = .imputation_model1(count = count_lnorm,
                                  point = log10(1.01),
                                  parslist = estmixmodelparams,
                                  drop_thre = 0.5,
                                  method = 2,
                                  ncores = ncores)

  # backtransform to unlogged counts
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

  # return imputed data
  return(ImputeData)
}

# DrImpute
#' @importFrom DrImpute preprocessSC DrImpute
#' @importFrom utils capture.output
.drimpute.calc <- function(countData) {
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

  # return object
  return(ImputeData)
}

# filtering ---------------------------------------------------------------

# apply stochastic filtering by Lun et al 2016 / Finak et al 2016
.countfilter.calc <- function(countData) {
  highE<- rowMeans(countData) >= 0.2
  FilterData <- countData[highE,]
  return(FilterData)
}

# apply frequency filter: gene must have less than 80% dropouts
.freqfilter.calc <- function(countData) {
  nsamples = ncol(countData)
  counts0 = countData == 0
  nn0 = rowSums(!counts0)
  p0 = (nsamples - nn0)/nsamples
  highE <- p0 < 0.8
  FilterData <- countData[highE,]
  return(FilterData)
}

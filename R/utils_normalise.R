
# Normalisation Wrapper ---------------------------------------------------

.norm.calc <- function(normalisation,
                       countData,
                       spikeData,
                       spikeInfo,
                       Lengths,
                       MeanFragLengths,
                       NCores) {
  if(normalisation=='TMM') {NormData <- .TMM.calc(countData = countData)}
  if(normalisation=='UQ') {NormData <- .UQ.calc(countData = countData)}
  if(normalisation=='MR') {NormData <- .MR.calc(countData = countData)}
  if(normalisation=='PosCounts') {NormData <- .PosCounts.calc(countData = countData)}
  if(normalisation=='scran') {NormData <- .scran.calc(countData = countData)}
  if(normalisation=='SCnorm') {NormData <- .SCnorm.calc(countData = countData,
                                                        spikeData = spikeData,
                                                        NCores = NCores)}
  if(normalisation=='Census') {NormData <- .Census.calc(countData=countData,
                                                        Lengths=Lengths,
                                                        MeanFragLengths=MeanFragLengths,
                                                        NCores=NCores)}
  if(normalisation=='RUVg') {NormData <- .RUVg.calc(countData = countData,
                                                    spikeData = spikeData)}
  if(normalisation=='BASiCS') {NormData <- .BASiCS.calc(countData = countData,
                                                        spikeData = spikeData,
                                                        spikeInfo = spikeInfo)}
  if(normalisation=='none') {NormData <- .none.calc(countData = countData)}
  return(NormData)
}

# TMM calculations --------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.TMM.calc <- function(countData) {
  norm.factors <- edgeR::calcNormFactors(object=countData, method='TMM')
  sf <- norm.factors * colSums(countData)
  sf <- sf / mean(sf, na.rm=T)
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}


# UQ calculations ---------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.UQ.calc <- function(countData) {
  norm.factors <- edgeR::calcNormFactors(object=countData,
                                         method='upperquartile')
  sf <- norm.factors * colSums(countData)
  sf <- sf / mean(sf, na.rm=T)
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}

# MR calculations ---------------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
.MR.calc <- function(countData) {
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData=data.frame(group=rep('A', ncol(countData))),
                                        design=~1))
  dds <- DESeq2::estimateSizeFactors(dds, type='ratio')
  sf <- DESeq2::sizeFactors(dds)
  names(sf) <- colnames(countData)
  norm.counts <- DESeq2::counts(dds, normalized=TRUE)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}

# poscounts normalization -------------------------------------------------

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors
.PosCounts.calc <- function(countData) {
  dds <- suppressMessages(
    DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                   colData=data.frame(group=rep('A', ncol(countData))),
                                   design=~1)
    )
  dds <- DESeq2::estimateSizeFactors(dds, type='poscounts')
  sf <- DESeq2::sizeFactors(dds)
  names(sf) <- colnames(countData)
  norm.counts <- DESeq2::counts(dds, normalized=TRUE)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}

# scran calculations ------------------------------------------------------

#' @importFrom scran computeSumFactors
#' @importFrom scater newSCESet
.scran.calc <- function(countData) {
  sce <- scater::newSCESet(countData=data.frame(countData))
  if(ncol(countData)<=14) {
    sf <- scran::computeSumFactors(sce,
                                   sizes=c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1))),
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    sf <- scran::computeSumFactors(sce,
                                   sizes=c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6))),
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    sf <- scran::computeSumFactors(sce,
                                   sizes=c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6))),
                                   positive=FALSE, sf.out=TRUE)
  }
  if(trunc(ncol(countData))>1000) {
    sf <- scran::computeSumFactors(sce, positive=FALSE, sf.out=TRUE)
  }
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}

# SCnorm calculations -----------------------------------------------------

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats median
.SCnorm.calc <- function(countData, spikeData, NCores) {
  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
  }

  cond <- rep("A", ncol(cnts))

  ncores = ifelse(is.null(NCores), parallel::detectCores() - 1, NCores)

  invisible(capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                                  Conditions = cond,
                                                  OutputName = "/dev/null",
                                                  SavePDF = FALSE,
                                                  PropToUse = .25,
                                                  Tau = .5,
                                                  reportSF = TRUE,
                                                  FilterCellNum = 10,
                                                  K = NULL,
                                                  NCores = ncores,
                                                  FilterExpression = 0,
                                                  Thresh = .1,
                                                  ditherCounts = FALSE,
                                                  withinSample = NULL,
                                                  useSpikes = spike))
  ))

  sf <- apply(scnorm.out$ScaleFactors, 2, stats::median)
  names(sf) <- colnames(cnts)
  res <- list(NormCounts=scnorm.out$NormalizedData,
              RoundNormCounts=round(scnorm.out$NormalizedData),
              size.factors=sf)
  return(res)
}


# Census normalisation ----------------------------------------------------

#' @importFrom monocle newCellDataSet relative2abs
# #' @importMethodsFrom monocle estimateSizeFactors sizeFactors
#' @importFrom VGAM tobit negbinomial.size
#' @importFrom parallel detectCores
.Census.calc <- function(countData, Lengths, MeanFragLengths, NCores) {

  ncores = ifelse(is.null(NCores), parallel::detectCores() - 1, NCores)

  # calculate TPM / CPM
  if(!is.null(Lengths) && !is.null(MeanFragLengths)) {
    ed <- .counts_to_tpm(countData=countData,
                         Lengths=Lengths,
                         MeanFragLengths=MeanFragLengths)
  }
  if(!is.null(Lengths) && is.null(MeanFragLengths)) {
    ed <- .calculateTPM(countData = countData, Lengths=Lengths)
  }
  if(is.null(Lengths) && is.null(MeanFragLengths)) {
    ed <- .calculateCPM(countData = countData)
  }
  # make annotated dataframes for monocle
  gene.dat <- data.frame(row.names = rownames(countData),
                         biotype=rep("protein_coding", nrow(countData)),
                         num_cells_expressed=rowSums(countData>0))
  cell.dat <- data.frame(row.names=colnames(countData),
                         Group=rep("A", ncol(countData)))
  fd <- new("AnnotatedDataFrame", data = gene.dat)
  pd <- new("AnnotatedDataFrame", data = cell.dat)

  # construct cell data set with expression values
  cds <- monocle::newCellDataSet(cellData = ed,
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit=0.1,
                                 expressionFamily=VGAM::tobit(Lower=0.1))
  # estimate RNA counts
  rpc_matrix <- monocle::relative2abs(cds, cores = ncores)
  #create cell data set wih RPC
  eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit=0.5,
                                 expressionFamily=VGAM::negbinomial.size())
  # apply normalisation
  eds.proc <- estimateSizeFactors(eds)

  norm.counts <- eds@assayData$exprs

  # extract output
  sf <- sizeFactors(eds.proc)
  names(sf) <- colnames(norm.counts)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}

.calculateTPM <- function(countData, Lengths) {
  rate <- countData / Lengths
  TPM <- rate / sum(rate) * 1e6
  return(TPM)
}

.calculateCPM <- function(countData) {
  CPM <- apply(countData,2, function(x) { (x/sum(x))*1000000 })
  return(CPM)
}

.counts_to_tpm <- function(countData, Lengths, MeanFragLengths) {

  # Ensure valid arguments.
  stopifnot(length(Lengths) == nrow(counts))
  stopifnot(length(MeanFragLengths) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    Lengths - MeanFragLengths[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  Lengths <- Lengths[idx]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

# RUVg spike-in normalisation ---------------------------------------------

#' @importFrom RUVSeq RUVg
.RUVg.calc <- function(countData, spikeData) {
  #RUVg calculation
  rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
  cnts = rbind(countData, spikeData)
  spike = grepl(pattern='ERCC', rownames(cnts))
  RUVg.out = RUVSeq::RUVg(x=as.matrix(cnts),
                          cIdx=spike,
                          k=1,
                          drop=0,
                          center=TRUE,
                          round=TRUE,
                          epsilon=1,
                          tolerance=1e-8,
                          isLog=FALSE)

  # return object
  sf <- apply(cnts/RUVg.out$normalizedCounts, 2, function(x) {median(x, na.rm = T)})
  names(sf) <- colnames(cnts)
  res <- list(NormCounts=RUVg.out$normalizedCounts,
              RoundNormCounts=round(RUVg.out$normalizedCounts),
              size.factors=sf)
  return(res)
}



# BASiCS spike-in normalisation -------------------------------------------

#' @importFrom BASiCS newBASiCS_Data BASiCS_Filter BASiCS_MCMC BASiCS_DenoisedCounts
.BASiCS.calc <- function(countData, spikeData, spikeInfo) {

  # create input
  spike = c(rep(FALSE, nrow(countData)), rep(TRUE, nrow(spikeData)))
  cnts = rbind(countData, spikeData)
  cnts = as.matrix(cnts)

  # perform filtering
  Filter = suppressMessages(BASiCS::BASiCS_Filter(Counts=cnts,
                                 Tech=spike,
                                 SpikeInput=spikeInfo$SpikeInput,
                                 BatchInfo = NULL,
                                 MinTotalCountsPerCell = 2,
                                 MinTotalCountsPerGene = 2,
                                 MinCellsWithExpression = 2,
                                 MinAvCountsPerCellsWithExpression = 2))
  SpikeInfoFilter = spikeInfo[spikeInfo$SpikeID %in% names(Filter$IncludeGenes)[Filter$IncludeGenes == TRUE],]

  invisible(capture.output(
    FilterData <- newBASiCS_Data(Counts = Filter$Counts,
                                 Tech = Filter$Tech,
                                 SpikeInfo = SpikeInfoFilter)
  ))

  # fit MCMC
  invisible(capture.output(
    MCMC_Output <- suppressMessages(BASiCS::BASiCS_MCMC(FilterData,
                                                        N = 10000,
                                                        Thin = 10,
                                                        Burn = 5000,
                                                        PrintProgress = FALSE))
  ))

  # normalised counts
  DenoisedCounts = suppressMessages(BASiCS::BASiCS_DenoisedCounts(Data = FilterData, Chain = MCMC_Output))
  DenoisedCounts = DenoisedCounts[!FilterData@Tech,]

  # size factor
  Phi = apply(MCMC_Output@phi,2,median)
  Nu = apply(MCMC_Output@nu,2,median)
  sf = Phi * Nu
  names(sf) = colnames(Filter$Counts)

  # return object
  res <- list(NormCounts=DenoisedCounts,
              RoundNormCounts=round(DenoisedCounts),
              size.factors=sf)
  return(res)

}


# No normalisation --------------------------------------------------------

.none.calc <- function(countData) {
  sf <- rep(1, ncol(countData))
  names(sf) <- colnames(countData)
  norm.counts <- countData
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}


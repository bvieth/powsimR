
# Normalisation Wrapper ---------------------------------------------------

.norm.calc <- function(normalisation,
                       countData,
                       spikeData,
                       spikeInfo,
                       batchData,
                       Lengths,
                       MeanFragLengths,
                       PreclustNumber,
                       NCores) {
  if(normalisation=='TMM') {NormData <- .TMM.calc(countData = countData)}
  if(normalisation=='UQ') {NormData <- .UQ.calc(countData = countData)}
  if(normalisation=='MR') {NormData <- .MR.calc(countData = countData)}
  if(normalisation=='PosCounts') {NormData <- .PosCounts.calc(countData = countData)}
  if(normalisation=='scran') {NormData <- .scran.calc(countData = countData)}
  if(normalisation=='scran' && !is.null(PreclustNumber)) {NormData <- suppressWarnings(
    .scranclust.calc(countData = countData,
                     PreclustNumber=PreclustNumber))}
  if(normalisation=='SCnorm') {NormData <- .SCnorm.calc(countData = countData,
                                                        spikeData = spikeData,
                                                        batchData = batchData,
                                                        NCores = NCores)}
  if(normalisation=='SCnorm' && !is.null(PreclustNumber)) {NormData <- .SCnormclust.calc(countData = countData,
                                                        spikeData = spikeData,
                                                        batchData = batchData,
                                                        PreclustNumber=PreclustNumber,
                                                        NCores = NCores)}
  if(normalisation=='Census') {NormData <- .Census.calc(countData=countData,
                                                        batchData=batchData,
                                                        Lengths=Lengths,
                                                        MeanFragLengths=MeanFragLengths,
                                                        spikeData=spikeData,
                                                        spikeInfo = spikeInfo,
                                                        NCores=NCores)}
  if(normalisation=='RUV') {NormData <- .RUV.calc(countData = countData,
                                                  spikeData = spikeData,
                                                  batchData = batchData)}
  if(normalisation=='BASiCS') {NormData <- .BASiCS.calc(countData = countData,
                                                        spikeData = spikeData,
                                                        spikeInfo = spikeInfo,
                                                        batchData = batchData)}
  if(normalisation=='depth') {NormData <- .depth.calc(countData = countData)}
  if(normalisation=='none') {NormData <- .none.calc(countData = countData)}
  return(NormData)
}

# TMM calculations --------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.TMM.calc <- function(countData) {
  norm.factors <- edgeR::calcNormFactors(object=countData, method='TMM')
  sf <- norm.factors * colSums(countData)
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "TMM"
  return(res)
}


# UQ calculations ---------------------------------------------------------

#' @importFrom edgeR calcNormFactors
.UQ.calc <- function(countData) {
  norm.factors <- edgeR::calcNormFactors(object=countData,
                                         method='upperquartile')
  sf <- norm.factors * colSums(countData)
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "UQ"
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
  attr(res, 'normFramework') <- "MR"
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
  attr(res, 'normFramework') <- "PosCounts"
  return(res)
}

# scran calculations ------------------------------------------------------

### TO DO: add clustering functionality, should help for norm calc of distinct cells!

#' @importFrom scran computeSumFactors
#' @importFrom scater newSCESet
.scran.calc <- function(countData) {
  sce <- scater::newSCESet(countData=data.frame(countData))
  if(ncol(countData)<=14) {
    sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1)))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6)))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    sizes <- c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6)))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(trunc(ncol(countData))>1000) {
    sizes <- c(round(seq(from=20, to=trunc(ncol(countData)/2), length.out=6)))
    sf <- scran::computeSumFactors(sce, positive=FALSE, sf.out=TRUE)
  }
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scran"
  return(res)
}

#' @importFrom scran computeSumFactors quickCluster
#' @importFrom scater newSCESet
.scranclust.calc <- function(countData, PreclustNumber) {
  sce <- scater::newSCESet(countData=data.frame(countData))
  if(ncol(countData)<=14) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>14 & ncol(countData)<=50) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>50 & ncol(countData)<=1000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>1000 & ncol(countData)<=5000) {
    clusters <- scran::quickCluster(sce, method="hclust", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
    sf <- scran::computeSumFactors(sce, positive=FALSE, sf.out=TRUE)
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(sce, method="igraph", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
    sf <- scran::computeSumFactors(sce, positive=FALSE, sf.out=TRUE)
  }
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scranCLUST"
  return(res)
}

# SCnorm calculations -----------------------------------------------------

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats median
.SCnorm.calc <- function(countData, spikeData, batchData, NCores) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
  }
  if(!is.null(batchData)) {
    cond <- batchData[,1]
  }
  if(is.null(batchData)) {
    cond <- rep('a', ncol(cnts))
  }

  ncores = ifelse(is.null(NCores), 1, NCores)

  invisible(capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                                  Conditions = cond,
                                                  OutputName = NULL,
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
                                                  useSpikes = spike,
                                                  useZerosToScale=FALSE))
  ))

  sf <- apply(scnorm.out@metadata$ScaleFactors, 2, stats::median)
  names(sf) <- colnames(cnts)
  gsf <- scnorm.out@metadata$ScaleFactors
  res <- list(NormCounts=scnorm.out@metadata$NormalizedData,
              RoundNormCounts=round(scnorm.out@metadata$NormalizedData),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "SCnorm"
  return(res)
}

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats median
.SCnormclust.calc <- function(countData, spikeData, batchData, PreclustNumber, NCores) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
  }
  if(!is.null(batchData)) {
    cond <- batchData[,1]
  }
  if(is.null(batchData) && !is.null(PreclustNumber)) {
    if(ncol(countData)<=5000) {
      clusters <- scran::quickCluster(cnts, method="hclust", min.size=floor(PreclustNumber/2))
      }
    if(ncol(countData)>5000) {
      clusters <- scran::quickCluster(cnts, method="igraph", min.size=floor(PreclustNumber/2))
      }
    cond <- as.character(clusters)
  }

  ncores = ifelse(is.null(NCores), 1, NCores)

  invisible(capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                                  Conditions = cond,
                                                  OutputName = NULL,
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
                                                  useSpikes = spike,
                                                  useZerosToScale=FALSE))
  ))

  sf <- apply(scnorm.out@metadata$ScaleFactors, 2, stats::median)
  names(sf) <- colnames(cnts)
  gsf <- scnorm.out@metadata$ScaleFactors
  res <- list(NormCounts=scnorm.out@metadata$NormalizedData,
              RoundNormCounts=round(scnorm.out@metadata$NormalizedData),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "SCnorm"
  return(res)
}


# Census normalisation ----------------------------------------------------

#' @importFrom monocle newCellDataSet relative2abs
#' @importFrom Biobase exprs
#' @importFrom VGAM tobit negbinomial.size
#' @importFrom parallel detectCores
.Census.calc <- function(countData, batchData, spikeData, spikeInfo, Lengths, MeanFragLengths, NCores) {

  ncores = ifelse(is.null(NCores), 1, NCores)

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
                         gene_short_name = rownames(countData),
                         biotype=rep("protein_coding", nrow(countData)),
                         num_cells_expressed=rowSums(countData>0))

  if(!is.null(batchData)) {
    cell.dat <- data.frame(row.names=colnames(countData),
                           Group=batchData[,1])
    ModelFormula <- "~Group"
  }
  if(is.null(batchData)) {
    cell.dat <- data.frame(row.names=colnames(countData),
                           Group=rep("A", ncol(countData)))
    ModelFormula <- "~1"
  }

  fd <- new("AnnotatedDataFrame", data = gene.dat)
  pd <- new("AnnotatedDataFrame", data = cell.dat)

  # construct cell data set with expression values
  cds <- monocle::newCellDataSet(cellData = as.matrix(ed),
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit=0.1,
                                 expressionFamily=VGAM::tobit(Lower=0.1))
  # estimate RNA counts
  rpc_matrix <- monocle::relative2abs(relative_cds = cds,
                                      modelFormulaStr = ModelFormula,
                                      method = "num_genes",
                                      cores = ncores)
  #create cell data set wih RPC
  eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit=0.5,
                                 expressionFamily=VGAM::negbinomial.size())
  # apply normalisation
  sf <- monocle:::estimateSizeFactorsForMatrix(counts=Biobase::exprs(eds))
  names(sf) <- colnames(Biobase::exprs(eds))

  norm.counts <- t(t(countData)/sf)

  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              RPC=as.matrix(rpc_matrix))
  attr(res, 'normFramework') <- "Census"
  return(res)
}

# RUV spike-in normalisation ---------------------------------------------

#' @importFrom RUVSeq RUVg
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr %>%
#' @importFrom dplyr rename_
.RUV.calc <- function(countData, spikeData, batchData) {

  if(is.null(batchData) && !is.null(spikeData)) {
    # annotate spike-ins and genes
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spike = grepl(pattern='ERCC', rownames(cnts))

    #RUVg calculation
    RUV.out = RUVSeq::RUVg(x=as.matrix(cnts),
                            cIdx=spike,
                            k=1,
                            drop=0,
                            center=TRUE,
                            round=FALSE,
                            epsilon=1,
                            tolerance=1e-8,
                            isLog=FALSE)
  }

  if(!is.null(batchData) && !is.null(spikeData)) {
    # prepare replicate sample annotation matrix
    tmp <- batchData %>%
      dplyr::rename_(Batch = names(.)[1]) %>%
      tibble::rownames_to_column(var="SampleID")  %>%
      tibble::rownames_to_column(var="SampleNumber")
    tmp2 <- split(x = tmp, f = as.factor(tmp$Batch))
    longestbatch <- max(sapply(tmp2, nrow))
    tmp3 <- sapply(tmp2, function(x) {
      tmp <- as.numeric(x$SampleNumber)
      tmp1 <- c(tmp, rep(-1, longestbatch-length(tmp)))
      tmp1
    })
    differences = t(tmp3)

    # annotate spike-ins and genes
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spike = grepl(pattern='ERCC', rownames(cnts))

    #RUVs calculation
    RUV.out = RUVSeq::RUVs(x=as.matrix(cnts),
                           cIdx=spike,
                           scIdx = differences,
                           k=1,
                           round=FALSE,
                           epsilon=1,
                           tolerance=1e-8,
                           isLog=FALSE)
  }

  if(!is.null(batchData) && is.null(spikeData)) {
    # prepare replicate sample annotation matrix
    tmp <- batchData %>%
      dplyr::rename_(Batch = names(.)[1]) %>%
      tibble::rownames_to_column(var="SampleID")  %>%
      tibble::rownames_to_column(var="SampleNumber")
    tmp2 <- split(x = tmp, f = as.factor(tmp$Batch))
    longestbatch <- max(sapply(tmp2, nrow))
    tmp3 <- sapply(tmp2, function(x) {
      tmp <- as.numeric(x$SampleNumber)
      tmp1 <- c(tmp, rep(-1, longestbatch-length(tmp)))
      tmp1
    })
    differences = t(tmp3)

    # annotate genes
    cnts = countData
    controls = rownames(countData)

    #RUVs calculation
    RUV.out = RUVSeq::RUVs(x=as.matrix(cnts),
                           cIdx=controls,
                           scIdx = differences,
                           k=1,
                           round=FALSE,
                           epsilon=1,
                           tolerance=1e-8,
                           isLog=FALSE)
  }

  # return object
  # normalized counts
  normCounts <- RUV.out$normalizedCounts
  normCounts[normCounts<0] <- 0
  # RUV does not return size factors, so make them constant
  sf <- rep(1, ncol(countData))
  # sf <- apply(countData/normCounts, 2, function(x) {median(x, na.rm = T)})
  names(sf) <- colnames(countData)
  res <- list(NormCounts = normCounts,
              RoundNormCounts = round(normCounts),
              size.factors = sf,
              RUV.W = RUV.out$W)
  attr(res, 'normFramework') <- "RUV"
  return(res)
}



# BASiCS spike-in normalisation -------------------------------------------

#' @importFrom BASiCS newBASiCS_Data BASiCS_Filter BASiCS_MCMC BASiCS_DenoisedCounts
.BASiCS.calc <- function(countData, spikeData, spikeInfo, batchData) {

  # create input
  spike = c(rep(FALSE, nrow(countData)), rep(TRUE, nrow(spikeData)))
  cnts = rbind(countData, spikeData)
  cnts = as.matrix(cnts)

  if(!is.null(batchData)) {
    cond <- batchData[,1]
  }
  if(is.null(batchData)) {
    cond <- NULL
  }

  # perform filtering
  Filter = suppressMessages(BASiCS::BASiCS_Filter(Counts=cnts,
                                 Tech=spike,
                                 SpikeInput=spikeInfo$SpikeInput,
                                 BatchInfo = cond,
                                 MinTotalCountsPerCell = 2,
                                 MinTotalCountsPerGene = 2,
                                 MinCellsWithExpression = 2,
                                 MinAvCountsPerCellsWithExpression = 2))
  SpikeInfoFilter = spikeInfo[rownames(spikeInfo) %in% names(Filter$IncludeGenes)[Filter$IncludeGenes == TRUE],]

  invisible(capture.output(
    FilterData <- BASiCS::newBASiCS_Data(Counts = Filter$Counts,
                                 Tech = Filter$Tech,
                                 SpikeInfo = SpikeInfoFilter,
                                 BatchInfo = Filter$BatchInfo)
  ))

  # fit MCMC
  invisible(capture.output(
    MCMC_Output <- suppressMessages(BASiCS::BASiCS_MCMC(FilterData,
                                                        N = 10000,
                                                        Thin = 10,
                                                        Burn =  5000,
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
              size.factors=sf,
              MCMC_Output=MCMC_Output,
              FilterData=FilterData)
  attr(res, 'normFramework') <- "BASiCS"
  return(res)

}


# Depth normalisation (spike-ins) -----------------------------------------

.depth.calc <- function(countData) {
  sf <- colSums(countData) / mean(colSums(countData))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "depth"
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
  attr(res, 'normFramework') <- "none"
  return(res)
}




# Normalisation Wrapper ---------------------------------------------------

.norm.calc <- function(normalisation,
                       countData,
                       spikeData,
                       spikeInfo,
                       batchData,
                       Lengths,
                       MeanFragLengths,
                       PreclustNumber,
                       NCores,
                       verbose) {
  if(normalisation=='TMM') {NormData <- .TMM.calc(countData = countData)}
  if(normalisation=='UQ') {NormData <- .UQ.calc(countData = countData)}
  if(normalisation=='MR') {NormData <- .MR.calc(countData = countData)}
  if(normalisation=='PosCounts') {NormData <- .PosCounts.calc(countData = countData)}
  if(normalisation=='Linnorm') {NormData <- .linnorm.calc(countData = countData,
                                                          spikeData = spikeData)}
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
                                                        NCores=NCores,
                                                        verbose=verbose)}
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

# TMM --------------------------------------------------------

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


# UQ ---------------------------------------------------------

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

# MR ---------------------------------------------------------

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

# poscounts -------------------------------------------------

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


# Linnorm -----------------------------------------------------------------

#' @importFrom Linnorm Linnorm.Norm
.linnorm.calc <- function(countData, spikeData) {

  spike = ifelse(is.null(spikeData), FALSE, TRUE)
  if(spike==TRUE) {
    rownames(spikeData) = paste('ERCC', 1:nrow(spikeData), sep="-")
    cnts = rbind(countData, spikeData)
    spikeID = rownames(spikeData)
  }
  if(spike==FALSE) {
    cnts = countData
    spikeID = NULL
  }

  linnorm.out <- Linnorm::Linnorm.Norm(datamatrix=cnts, # matrix/data.frame of raw counts
                                       RowSamples = FALSE, # if I would switch gene and samples position in input
                                       spikein = spikeID, # names of the spike-ins
                                       spikein_log2FC = NULL,  # LFC of spike-ins, assume mix 1, so no
                                       showinfo = FALSE, # verbosity
                                       output = "Raw", # type of output: raw, ie  total count output will  be equal to median of input counts
                                       minNonZeroPortion = 0.75, # minimum porportion of nonzero values per gene
                                       BE_F_p = 0.3173, # p-value cutoff for standard deviation and skewness testing of batch effect normalisation
                                       BE_F_LC_Genes = "Auto", # porportion of lowly expressed genes to filter before batch effect normalisation
                                       BE_F_HC_Genes = 0.01, # proportion of highly expressed genes to filter before batch effect normalisation
                                       BE_strength = 0.5, # strength of batch effect normalisation
                                       max_F_LC = 0.75) # maximum threshold for filtering of low and high expression

  norm.counts <- linnorm.out[!grepl(pattern="ERCC", rownames(linnorm.out)),]
  gsf <-  t(t(countData)/t(norm.counts))
  sf <- apply(gsf, 2, stats::median, na.rm=T)
  sf <- sf - mean(sf) + 1
  names(sf) <- colnames(countData)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "Linnorm"
  return(res)
}

# scran ------------------------------------------------------

#' @importFrom scran computeSumFactors
#' @importFrom SingleCellExperiment SingleCellExperiment
.scran.calc <- function(countData) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = countData))
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
#' @importFrom SingleCellExperiment SingleCellExperiment
.scranclust.calc <- function(countData, PreclustNumber) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = countData))
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
  }
  if(ncol(countData)>5000) {
    clusters <- scran::quickCluster(sce, method="igraph", min.size=floor(PreclustNumber/2))
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=6))))
    sf <- scran::computeSumFactors(sce,
                                   sizes=sizes,
                                   cluster=clusters,
                                   positive=FALSE, sf.out=TRUE)
  }
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "scranCLUST"
  return(res)
}

# SCnorm -----------------------------------------------------

#' @importFrom SCnorm SCnorm
#' @importFrom parallel detectCores
#' @importFrom stats median
#' @importFrom utils capture.output
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

  invisible(utils::capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                                  Conditions = cond,
                                                  PrintProgressPlots = FALSE,
                                                  reportSF = TRUE,
                                                  FilterCellNum = 10,
                                                  FilterExpression = 0,
                                                  Thresh = 0.1,
                                                  K = NULL,
                                                  NCores = ncores,
                                                  ditherCounts = TRUE,
                                                  PropToUse = 0.15,
                                                  Tau = 0.5,
                                                  withinSample = NULL,
                                                  useSpikes = spike,
                                                  useZerosToScale = FALSE))
  ))

  sf <- apply(scnorm.out@metadata$ScaleFactors, 2, stats::median)
  names(sf) <- colnames(cnts)
  gsf <- scnorm.out@metadata$ScaleFactors
  rownames(gsf) <- rownames(cnts)
  colnames(gsf) <- colnames(cnts)
  gsf <- gsf[!grepl(pattern="ERCC", rownames(gsf)),]
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

  invisible(utils::capture.output(
    scnorm.out <- suppressMessages(SCnorm::SCnorm(Data = cnts,
                                                  Conditions = cond,
                                                  PrintProgressPlots = FALSE,
                                                  reportSF = TRUE,
                                                  FilterCellNum = 10,
                                                  FilterExpression = 0,
                                                  Thresh = 0.1,
                                                  K = NULL,
                                                  NCores = ncores,
                                                  ditherCounts = TRUE,
                                                  PropToUse = 0.20,
                                                  Tau = 0.5,
                                                  withinSample = NULL,
                                                  useSpikes = spike,
                                                  useZerosToScale = FALSE))
  ))

  sf <- apply(scnorm.out@metadata$ScaleFactors, 2, stats::median)
  names(sf) <- colnames(cnts)
  gsf <- scnorm.out@metadata$ScaleFactors
  rownames(gsf) <- rownames(cnts)
  colnames(gsf) <- colnames(cnts)
  gsf <- gsf[!grepl(pattern="ERCC", rownames(gsf)),]
  res <- list(NormCounts=scnorm.out@metadata$NormalizedData,
              RoundNormCounts=round(scnorm.out@metadata$NormalizedData),
              size.factors=sf,
              scale.factors=gsf)
  attr(res, 'normFramework') <- "SCnorm"
  return(res)
}


# Census ----------------------------------------------------

#' @importFrom monocle newCellDataSet relative2abs
#' @importFrom Biobase exprs
#' @importFrom VGAM tobit negbinomial.size
#' @importFrom parallel detectCores
.Census.calc <- function(countData,
                         batchData,
                         spikeData,
                         spikeInfo,
                         Lengths,
                         MeanFragLengths,
                         NCores,
                         verbose) {

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
  est_t <- .estimate_t(relative_expr_matrix=Biobase::exprs(cds),
                       relative_expr_thresh=0)

  if(is.null(spikeData) && is.null(spikeInfo)) {
    if(verbose) {message(paste0("Running monocle relative2abs."))}
    rpc_matrix <- monocle::relative2abs(relative_cds = cds,
                                        t_estimate = est_t,
                                        modelFormulaStr = ModelFormula,
                                        method = "num_genes",
                                        cores = ncores,
                                        verbose = verbose)

    #create cell data set wih RPC
    eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit=0.5,
                                   expressionFamily=VGAM::negbinomial.size())
    # apply normalisation
    sf <- .estimateSizeFactorsForDenseMatrix(counts=Biobase::exprs(eds),
                                             locfunc=median,
                                             round_exprs = T,
                                             method = "mean-geometric-mean-total")
    names(sf) <- colnames(Biobase::exprs(eds))

    norm.counts <- t(t(countData)/sf)
  }

  if(!is.null(spikeData) && !is.null(spikeInfo)) {
    if(verbose) {message(paste0("Running own relative2abs."))}
    # relative expression of spike-ins
    if(!is.null(spikeInfo$Lengths)) {
      ed.spike <- .calculateTPM(countData = spikeData,
                                Lengths=spikeInfo$Lengths)
    }
    if(is.null(spikeInfo$Lengths)) {
      ed.spike <- .calculateCPM(countData = spikeData)
    }

    rel.results <- .relative2abs(relative_cds = cds,
                                relative_spike = ed.spike,
                                spikeInfo = spikeInfo,
                                verbose = verbose,
                                cores = ncores)
    rpc_matrix <- rel.results$norm_cds

    #create cell data set wih RPC
    eds <- monocle::newCellDataSet(cellData = as.matrix(rpc_matrix),
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit=0.5,
                                   expressionFamily=VGAM::negbinomial.size())

    # apply normalisation
    sf <- .estimateSizeFactorsForDenseMatrix(counts=Biobase::exprs(eds),
                                             locfunc=median,
                                             round_exprs = T,
                                             method = "mean-geometric-mean-total")
    #  It can be either "mean-geometric-mean-total" (default), "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
    names(sf) <- colnames(Biobase::exprs(eds))

    norm.counts <- t(t(countData)/sf)
  }

  # print(summary(sf))
  # print(table(sf == 1))

  # return object
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf,
              RPC=as.matrix(rpc_matrix))
  attr(res, 'normFramework') <- "Census"
  return(res)
}

# census estimatesizefactors function
.estimateSizeFactorsForDenseMatrix <- function(counts,
                                               locfunc,
                                               round_exprs,
                                               method){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(t(t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- .estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  return(sfs)
}

# census estimate_t function
.estimate_t <- function(relative_expr_matrix, relative_expr_thresh) {
  #apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr)
    10^mean(.dmode(log10(relative_expr[relative_expr > relative_expr_thresh])))))
}

# census mode function
#' @importFrom stats density
.dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- stats::density(x, kernel=c("gaussian"), na.rm = T)
  ( den$x[den$y==max(den$y)] )
}

#' @importFrom Biobase exprs
#' @importFrom MASS rlm
#' @importFrom stats predict
#' @importFrom parallel mcmapply
.relative2abs <- function(relative_cds,
                          relative_spike, # relative exprs matrix of spike-ins
                          spikeInfo, # spike input
                          verbose,
                          cores) {
  FPKM <- NULL
  # relative expression matrix of genes
  relative_expr_matrix <- Biobase::exprs(relative_cds)

  # relative expression of spike-ins
  ERCC_controls <- relative_spike

  # spike input information
  ERCC_annotation <- spikeInfo
  valid_ids <- which(ERCC_annotation[, "SpikeInput"] >= 0)

  # robust linear regression
  if (verbose) {message("Performing robust linear regression for each cell based
                        on the spike-in data")}
  molModels <- apply(ERCC_controls, 2, function(cell_exprs, input.ERCC.annotation, valid_ids) {
    spike_df <- input.ERCC.annotation
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$SpikeInput
    spike_df$rounded_numMolecules <- round(spike_df$numMolecules)
    if (is.null(valid_ids))
      spike_df <- subset(spike_df, FPKM >= 1e-10)
    else {
      spike_df <- spike_df[valid_ids, ]
      spike_df <- subset(spike_df, FPKM >= 1e-10)
    }
    spike_df$log_fpkm <- log10(spike_df$FPKM)
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    molModel <- tryCatch({
      molModel <- MASS::rlm(log_numMolecules ~ log_fpkm,
                            data = spike_df)
      molModel
    }, error = function(e) {
      print(e)
      NULL
    })
    molModel
  }, ERCC_annotation, valid_ids)

  if (verbose) {message("Apply the fitted robust linear regression model
                        to recover the absolute copy number for all transcripts in each cell")}
  norm_fpkms <- parallel::mcmapply(function(cell_exprs, molModel) {
    tryCatch({
      norm_df <- data.frame(log_fpkm = log10(cell_exprs))
      res <- 10^stats::predict(molModel, type = "response",
                        newdata = norm_df)
    }, error = function(e) {
      rep(NA, length(cell_exprs))
    })
  }, split(as.matrix(relative_expr_matrix), rep(1:ncol(relative_expr_matrix),
                                                each = nrow(relative_expr_matrix))), molModels, mc.cores = cores)
  k_b_solution <- data.frame(b = unlist(lapply(molModels,
                                               FUN = function(x) {
                                                 intercept = x$coefficients[1]
                                               })), k = unlist(lapply(molModels, FUN = function(x) {
                                                 slope = x$coefficients[2]
                                               })))
  kb_model <- MASS::rlm(b ~ k, data = k_b_solution)
  kb_slope <- kb_model$coefficients[2]
  kb_intercept <- kb_model$coefficients[1]

  rownames(norm_fpkms) <- rownames(relative_expr_matrix)
  colnames(norm_fpkms) <- colnames(relative_expr_matrix)

  # return object
  return(list(norm_cds = norm_fpkms,
              kb_slope = kb_slope,
              kb_intercept = kb_intercept,
              k_b_solution = k_b_solution))
}



# RUV ---------------------------------------------

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
  normCounts <- RUV.out$normalizedCounts[!grepl(pattern = "ERCC",
                                                rownames(RUV.out$normalizedCounts)),]
  normCounts[normCounts<0] <- 0
  # RUV proxy size factors
  gsf <-  t(t(countData)/t(normCounts))
  sf <- apply(gsf, 2, stats::median, na.rm=T)
  names(sf) <- colnames(countData)
  res <- list(NormCounts = normCounts,
              RoundNormCounts = round(normCounts),
              size.factors = sf,
              RUV.W = RUV.out$W)
  attr(res, 'normFramework') <- "RUV"
  return(res)
}

# BASiCS -------------------------------------------

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

  invisible(utils::capture.output(
    FilterData <- BASiCS::newBASiCS_Data(Counts = Filter$Counts,
                                 Tech = Filter$Tech,
                                 SpikeInfo = SpikeInfoFilter,
                                 BatchInfo = Filter$BatchInfo)
  ))

  # fit MCMC
  invisible(utils::capture.output(
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


# Depth normalisation -----------------------------------------

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



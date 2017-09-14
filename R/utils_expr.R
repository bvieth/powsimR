# methods:
# TPM with column sums as denominator
# CPM and RPKM with normalisation factors


# EXPRESSION VALUES WRAPPER -----------------------------------------------

.expr.calc <- function(countData, normData, Lengths, MeanFragLengths, group) {
  if(!is.null(Lengths) && !is.null(MeanFragLengths)) {
    expr.data <- .counts_to_tpm(countData=countData,
                                Lengths=Lengths,
                                MeanFragLengths=MeanFragLengths)
  }
  if(!is.null(Lengths) && is.null(MeanFragLengths)) {
    expr.data <- .rpkm.calc(countData=countData,
                            normData=normData,
                            Lengths=Lengths,
                            group=group)
  }
  if(is.null(Lengths) && is.null(MeanFragLengths)) {
    expr.data <- .cpm.calc(countData=countData,
                           normData=normData,
                           group=group)
  }
  return(expr.data)
}


# CPM  --------------------------------------------------------------------

#' @importFrom edgeR DGEList cpm.DGEList
.cpm.calc <- function(normData, countData, group) {
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = group,
                        remove.zeros = FALSE)
  # size factor normalised CPM values
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  return(out.cpm)
}

# RPKM --------------------------------------------------------------------

#' @importFrom edgeR DGEList cpm.DGEList
.rpkm.calc <- function(normData, countData, Lengths, group) {
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = group,
                        remove.zeros = FALSE)
  dge$genes$Length <- Lengths
  # size factor normalised RPKM values
  out.rpkm <- edgeR::rpkm.DGEList(dge, normalized.lib.sizes = T, log = F)
  return(out.rpkm)
}

# RELATIVE MEASURES USING SEQ DEPTH AS DENOMINATOR ------------------------

# following three functions taken and adapted from
# https://gist.github.com/slowkow/c6ab0348747f86e2748b


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





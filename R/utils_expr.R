# CPM  --------------------------------------------------------------------

.calculateCPM <- function(countData) {
  CPM <- apply(countData,2, function(x) { (x/sum(x))*1000000 })
  return(CPM)
}


# TPM ---------------------------------------------------------------------

# accounting for gene lengths
.calculateTPM <- function(countData, Lengths) {
  rate <- countData / Lengths
  TPM <- rate / sum(rate) * 1e6
  return(TPM)
}

# following three functions taken and adapted from
# https://gist.github.com/slowkow/c6ab0348747f86e2748b

# accounting for effective gene lengths
.counts_to_tpm <- function(countData, Lengths, MeanFragLengths) {

  # Ensure valid arguments.
  stopifnot(length(Lengths) == nrow(countData))
  stopifnot(length(MeanFragLengths) == ncol(countData))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(countData), function(i) {
    Lengths - MeanFragLengths[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- countData[idx,]
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


# LOG TRANSFORMED EXPRESSION ----------------------------------------------

.calculateExpr <- function(countData, normData, Lengths) {


}

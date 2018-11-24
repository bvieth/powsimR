# CPM  --------------------------------------------------------------------

.calculateCPM <- function(countData) {
  CPM <- apply(countData,2, function(x) { (x/sum(x))*1000000 })
  return(CPM)
}


# TPM ---------------------------------------------------------------------

# accounting for gene lengths
.calculateTPM <- function(countData, Lengths) {
  rate <- countData / Lengths
  TPM <- apply(rate, 2, function(x) { (x/sum(x))*1000000 })
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


# NORMALIZED EXPRESSION ----------------------------------------------

.calculateExpr <- function(countData, normData, Lengths, MeanFragLengths) {
  # normalized library sizes
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))

  # calculate TPM / CPM
  if (!attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    if(!is.null(Lengths) & !is.null(MeanFragLengths)) {
      # Compute effective lengths of features in each library.
      effLen <- do.call(cbind, lapply(1:ncol(countData), function(i) {
        tmp <- Lengths - MeanFragLengths[i] + 1
      }))
      # Exclude genes with length less than the mean fragment length.
      idx <- apply(effLen, 1, function(x) min(x) > 1)
      countData <- countData[idx,]
      effLen <- effLen[idx,]

      # TPM
      t.vals <- countData / effLen
      norm.lib.size <- colSums(t.vals) * nsf
      TPM <- t(t(t.vals) / norm.lib.size) * 1e6
    }
    if(!is.null(Lengths) & is.null(MeanFragLengths)) {
      # TPM
      t.vals <- countData / Lengths
      norm.lib.size <- colSums(t.vals) * nsf
      TPM <- t(t(t.vals) / norm.lib.size) * 1e6
    }
    if(is.null(Lengths) & is.null(MeanFragLengths)) {
      # CPM
      norm.lib.size <- colSums(countData) * nsf
      TPM <- t(t(countData) / norm.lib.size) * 1e6
    }
  }
  if (attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    # scaling factors for libs
    sf <- normData$size.factors

    # scaling factors for genes
    gsf.mat <- matrix(ncol=ncol(countData), nrow=nrow(countData),
                  dimnames = list(rownames(countData), colnames(countData)))
    gsf.vals <- normData$scale.factors
    cols <- colnames(gsf.mat)[colnames(gsf.mat) %in% colnames(gsf.vals)]
    rows <- rownames(gsf.mat)[rownames(gsf.mat) %in% rownames(gsf.vals)]
    gsf.mat[rows, cols] <- gsf.vals[rows, cols]

    # fill in missing values with sf of libs
    idx <- apply(gsf.mat, 1, function(rows) {
      any(is.na(rows))
    })
    gsf.mat[idx,] <- sf
    normData <- countData/gsf.mat

    if(!is.null(Lengths) & !is.null(MeanFragLengths)) {
      # Compute effective lengths of features in each library.
      effLen <- do.call(cbind, lapply(1:ncol(normData), function(i) {
        tmp <- Lengths - MeanFragLengths[i] + 1
      }))
      # Exclude genes with length less than the mean fragment length.
      idx <- apply(effLen, 1, function(x) min(x) > 1)
      normData <- normData[idx,]
      effLen <- effLen[idx,]
      # TPM
      t.vals <- normData / effLen
      TPM <- apply(t.vals,2, function(x) { (x/sum(x))*1000000 })
    }
    if(!is.null(Lengths) & is.null(MeanFragLengths)) {
      # TPM
      t.vals <- normData / Lengths
      TPM <- apply(t.vals,2, function(x) { (x/sum(x))*1000000 })
    }
    if(is.null(Lengths) & is.null(MeanFragLengths)) {
      # CPM
      t.vals <- normData
      TPM <- apply(t.vals,2, function(x) { (x/sum(x))*1000000 })
    }
  }

  return(TPM)
}


# WRAPPER DIMENSION REDUCTION ---------------------------------------------

.dim.reduce <- function(DimReduce,
                        FeatureSelect) {
  if(DimReduce=="Euclidean") {DimData <- .euclid.calc(FeatureSelect=FeatureSelect)}
  if(DimReduce=="Spearman") {DimData <- .spearman.calc(FeatureSelect=FeatureSelect)}
  if(DimReduce=="PCA") {DimData <- .PCA.calc(FeatureSelect=FeatureSelect)}
  if(DimReduce=="t-SNE-expr") {DimData <- .tSNE.expr.calc(FeatureSelect=FeatureSelect)}
  if(DimReduce=="t-SNE-Euclidean") {DimData <- .tSNE.euclid.calc(FeatureSelect=FeatureSelect)}
  if(DimReduce=="PCA+t-SNE") {DimData <- .PCA.tSNE.calc(FeatureSelect=FeatureSelect)}

  return(DimData)
}


# DISTANCE DIMENSION REDUCTION --------------------------------------------

#' @importFrom stats dist
.euclid.calc <- function(FeatureSelect) {
  # euclidian distance matrix
  dist.cpm <- stats::dist(t(FeatureSelect), method = 'euclidean')
  attr(dist.cpm, 'DimReduceMethod') <- "Euclidean"
  return(dist.cpm)
}

# DISTANCE DIMENSION REDUCTION --------------------------------------------

#' @importFrom stats cor as.dist
.spearman.calc <- function(FeatureSelect) {
  # spearman rank correlation
  cor.cpm <-  stats::as.dist(1 - stats::cor(FeatureSelect,method="spearman"))
  attr(cor.cpm, 'DimReduceMethod') <- "Spearman"
  return(cor.cpm)
}

# PCA DIMENSION REDUCTION -------------------------------------------------

#' @importFrom stats prcomp
.PCA.calc <- function(FeatureSelect){
  pca_data <- data.frame(t(FeatureSelect))
  pca <- stats::prcomp(~ ., data=pca_data, center = TRUE, scale. = TRUE)
  pc_names <- sprintf("PC%d", 1:ncol(pca$x))
  rownames(pca$x) <- colnames(FeatureSelect)
  colnames(pca$x) <- pc_names
  n_p <- ifelse(ncol(pca$x)>=10, 10, ncol(pca$x))
  pca_loadings <- pca$x[, 1:n_p]

  attr(pca_loadings, 'DimReduceMethod') <- "PCA"
  return(pca_loadings)
}

# t-SNE DIMENSION REDUCTION WITH RAW EXPRESSION ---------------------------

#' @importFrom Rtsne Rtsne
.tSNE.expr.calc <- function(FeatureSelect) {
  # dimension reduction on expression matrix: t-SNE only
  sample.value <- ncol(FeatureSelect) -2
  max.perplexity <- sample.value/3
  tsne.res <- Rtsne::Rtsne(X=t(FeatureSelect), pca=F, is_distance=F, perplexity=max.perplexity )
  tsne.out <- tsne.res$Y
  tsne_names <- sprintf("Dim-%d", 1:ncol(tsne.out))
  rownames(tsne.out) <- colnames(FeatureSelect)
  colnames(tsne.out) <- tsne_names

  attr(tsne.out, 'DimReduceMethod') <- "t-SNE-Expr"
  return(tsne.out)
}

# t-SNE DIMENSION REDUCTION WITH DISTANCE MEASURE -------------------------

#' @importFrom stats dist
#' @importFrom Rtsne Rtsne
.tSNE.euclid.calc <- function(FeatureSelect) {
  # dimension reduction on expression matrix: first dist then t-SNE
  dist.cpm <- stats::dist(t(FeatureSelect), method = 'euclidean')
  sample.value <- attr(dist.cpm, "Size") -2
  max.perplexity <- sample.value/3
  tsne.res <- Rtsne::Rtsne(X=dist.cpm, pca=F, is_distance=T, perplexity=max.perplexity )
  tsne.out <- tsne.res$Y
  tsne_names <- sprintf("Dim-%d", 1:ncol(tsne.out))
  rownames(tsne.out) <- colnames(FeatureSelect)
  colnames(tsne.out) <- tsne_names

  attr(tsne.out, 'DimReduceMethod') <- "t-SNE-Euclidean"
  return(tsne.out)
}

# PCA + t-SNE DIMENSION REDUCTION -----------------------------------------

#' @importFrom Rtsne Rtsne
.PCA.tSNE.calc <- function(FeatureSelect) {
  # dimension reduction on expression matrix: PCA + t-SNE
  sample.value <- ncol(FeatureSelect) -2
  max.perplexity <- sample.value/3
  tsne.res <- Rtsne::Rtsne(X=t(FeatureSelect), pca=T, is_distance=F, perplexity=max.perplexity )
  tsne.out <- tsne.res$Y
  tsne_names <- sprintf("Dim-%d", 1:ncol(tsne.out))
  rownames(tsne.out) <- colnames(FeatureSelect)
  colnames(tsne.out) <- tsne_names

  attr(tsne.out, 'DimReduceMethod') <- "PCA+t-SNE"
  return(tsne.out)
}

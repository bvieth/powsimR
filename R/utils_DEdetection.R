# NOTES -------------------------------------------------------------------

# Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!

# ROTS, NOISeq, EBSeq, monocle, scDD have no log fold changes internally calculated?

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
                     DEOpts,
                     spikeData,
                     spikeInfo,
                     Lengths,
                     MeanFragLengths,
                     NCores) {
  # methods developed for bulk
  if (DEmethod == "edgeR-LRT") {DERes = .run.edgeRLRT(normData=normData,
                                                      countData=countData,
                                                      DEOpts=DEOpts)}
  if (DEmethod == "edgeR-QL") {DERes = .run.edgeRQL(normData=normData,
                                                    countData=countData,
                                                    DEOpts=DEOpts)}
  if (DEmethod == "limma-voom") {DERes = .run.limma.voom(normData=normData,
                                                         countData=countData,
                                                         DEOpts=DEOpts)}
  if (DEmethod == "limma-trend") {DERes = .run.limma.trend(normData=normData,
                                                           countData=countData,
                                                           DEOpts=DEOpts)}
  if (DEmethod == "DESeq2") {DERes = .run.DESeq2(normData=normData,
                                                 countData=countData,
                                                 DEOpts=DEOpts,
                                                 NCores=NCores)}
  if (DEmethod == "ROTS") {DERes = .run.ROTS(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts)}
  if (DEmethod == "baySeq") {DERes = .run.baySeq(normData=normData,
                                                 countData=countData,
                                                 DEOpts=DEOpts,
                                                 NCores=NCores)}
  if (DEmethod == "NOISeq") {DERes = .run.NOISeq(normData=normData,
                                                 countData=countData,
                                                 DEOpts=DEOpts)}
  if (DEmethod == "DSS") {DERes = .run.DSS(normData=normData,
                                           countData=countData,
                                           DEOpts=DEOpts)}
  if (DEmethod=='EBSeq') {DERes = .run.EBSeq(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts)}
  # methods developed for single cell
  if (DEmethod == "MAST") {DERes = .run.MAST(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts,
                                             NCores=NCores)}
  if (DEmethod == "scde") {DERes = .run.scde(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts,
                                             NCores=NCores)}
  if (DEmethod == "BPSC") {DERes = .run.BPSC(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts,
                                             NCores)}
  if (DEmethod == "scDD") {DERes = .run.scDD(normData=normData,
                                             countData=countData,
                                             DEOpts=DEOpts,
                                             NCores=NCores)}
  if (DEmethod == "monocle") {DERes = .run.monocle(normData=normData,
                                                   countData=countData,
                                                   DEOpts=DEOpts,
                                                   NCores=NCores)}
  # if (DEmethod == "BASiCS") {DERes = .run.BASiCS(normData=normData,
  #                                                countData=countData,
  #                                                DEOpts=DEOpts,
  #                                                NCores=NCores)}
  # if (DEmethod == "D3E") {DERes = .run.D3E(normData=normData,
  #                                          countData=countData,
  #                                          DEOpts=DEOpts,
  #                                          NCores=NCores)}
  return(DERes)
}

# edgeR -------------------------------------------------------------------

#' @importFrom edgeR DGEList estimateGLMRobustDisp glmFit glmLRT topTags
#' @importFrom stats model.matrix
.run.edgeRLRT <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  dge <- edgeR::estimateGLMRobustDisp(y=dge, design = design.mat)
  fit.edgeR <- edgeR::glmFit(dge, design = design.mat)
  lrt.edgeR <- edgeR::glmLRT(fit.edgeR)
  res.edgeR <- edgeR::topTags(lrt.edgeR, adjust.method="BH", n=Inf, sort.by = 'none')

  ## construct results
  result <- data.frame(geneIndex=rownames(res.edgeR$table),
                       pval=res.edgeR$table$PValue,
                       fdr=rep(NA, nrow(res.edgeR$table)),
                       lfc=res.edgeR$table$logFC,
                       stringsAsFactors = F)
  return(result)
}

#' @importFrom edgeR DGEList estimateGLMRobustDisp topTags glmQLFit glmQLFTest
#' @importFrom stats model.matrix
.run.edgeRQL <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  design.mat <- stats::model.matrix( ~DEOpts$designs)
  dge <- edgeR::estimateDisp(y=dge, design = design.mat)
  fit.edgeR <- edgeR::glmQLFit(dge, design = design.mat, robust=TRUE)
  ql.edgeR <- edgeR::glmQLFTest(fit.edgeR)
  res.edgeR <- edgeR::topTags(ql.edgeR, adjust.method="BH", n=Inf, sort.by = 'none')

  ## construct results
  result <- data.frame(geneIndex=rownames(res.edgeR$table),
                       pval=res.edgeR$table$PValue,
                       fdr=rep(NA, nrow(res.edgeR$table)),
                       lfc=res.edgeR$table$logFC,
                       stringsAsFactors = F)
  return(result)
}


# limma -------------------------------------------------------------------

#' @importFrom limma lmFit eBayes voom topTable
#' @importFrom edgeR DGEList
#' @importFrom stats model.matrix
.run.limma.voom <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- DEOpts$p.DE
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  v <- limma::voom(dge, design.mat, plot=FALSE)
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
.run.limma.trend <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- DEOpts$p.DE
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  fit <- limma::lmFit(object = y, design = design.mat)
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
.run.DESeq2 <- function(normData, countData, DEOpts, NCores) {

  # construct input object
  coldat <- data.frame(design=factor(DEOpts$designs))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  dds <- DESeq2::DESeqDataSetFromMatrix(countData,
                                        coldat, ~design,
                                        tidy = FALSE, ignoreRank = FALSE)
  DESeq2::sizeFactors(dds) <- sf

  # run DE testing
  if (is.null(NCores)) {
    fit.DeSeq <- DESeq2::DESeq(dds,
                               test="Wald",
                               quiet = TRUE,
                               parallel=FALSE)
  }
  if (!is.null(NCores)) {
    fit.DeSeq <- DESeq2::DESeq(dds, test="Wald",
                               quiet = TRUE,
                               parallel=T,
                               BPPARAM = BiocParallel::MulticoreParam(NCores))
  }
  res.DeSeq <- DESeq2::results(fit.DeSeq)

  ## construct results
  result <- data.frame(geneIndex=rownames(res.DeSeq),
                       pval=res.DeSeq$pvalue,
                       fdr=rep(NA, nrow(countData)),
                       lfc=res.DeSeq$log2FoldChange,
                       stringsAsFactors = F)
  return(result)
}

# ROTS --------------------------------------------------------------------

#' @importFrom edgeR DGEList cpm.DGEList
#' @importFrom ROTS ROTS
#' @importFrom limma voom
#' @importFrom stats model.matrix
.run.ROTS <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # apply voom to get log2 expression values
  design.mat <- stats::model.matrix( ~ factor(DEOpts$designs))
  v <- limma::voom(dge, design.mat, plot=FALSE)
  out.expr <- v$E

  # run DE testing
  res <-  ROTS::ROTS(data = out.expr,
                     groups = factor(DEOpts$designs),
                     B = 100,
                     K = floor(nrow(out.expr)/2) , progress=F)

  # construct results
  result <- data.frame(geneIndex=rownames(res$data),
                       pval=res$pvalue,
                       fdr=rep(NA, nrow(res$data)),
                       lfc=res$logfc,
                       stringsAsFactors = F)
  return(result)
}


# baySeq ------------------------------------------------------------------

#' @importFrom snow makeCluster stopCluster
#' @importMethodsFrom baySeq libsizes
#' @importFrom baySeq getPriors.NB getLikelihoods topCounts
.run.baySeq <- function(normData, countData, DEOpts, NCores) {

  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  # set multiple cores
  if(is.null(NCores)) {
    cl <- NULL
  }
  if(!is.null(NCores)) {
    cl <- snow::makeCluster(NCores)
  }

  # make input data sets for baySeq
  replicates <- ifelse(DEOpts$designs==min(DEOpts$designs), "A", "B")
  groups <- list(NDE = c(rep(1, length(DEOpts$designs))),
                 DE = c(ifelse(DEOpts$designs==min(DEOpts$designs), 1, 2)))
  CD <- new("countData", data = countData,
            replicates = replicates,
            groups = groups)
  # fill in library size factors
  CD@sampleObservables$libsizes <- sf
  CD@annotation <- data.frame(name = rownames(countData),
                              stringsAsFactors = F)
  # run prior estimation
  CD <- baySeq::getPriors.NB(CD,
                             samplesize = ifelse(nrow(countData)<1e5, nrow(countData), 1e5),
                             estimation = "QL", cl = cl,
                             equalDispersions=TRUE, verbose=F)
  # run likelihood ratio test
  CD <- baySeq::getLikelihoods(CD, cl = cl,
                               bootStraps = 10,
                               verbose = FALSE)
  # get test results
  res <- baySeq::topCounts(cD=CD, group="DE",
                           decreasing = FALSE,
                           number = Inf,
                           normaliseData = FALSE)
  res <- res[match(CD@annotation$name, res$annotation),]

  # free multiple cores
  if(!is.null(NCores)) {
    snow::stopCluster(cl)
  }

  # construct result data frame
  result = data.frame(geneIndex=res$annotation,
                      pval=rep(NA, nrow(countData)),
                      fdr=res$FDR.DE,
                      stringsAsFactors = F)
  return(result)
}


# NOISeq ------------------------------------------------------------------

#' @importFrom NOISeq readData noiseqbio
#' @importFrom edgeR DGEList cpm.DGEList
.run.NOISeq <- function(normData, countData, DEOpts) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct dge object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # size factor normalised log2(CPM+1) values
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)

  # make input data set
  groups <- data.frame(Group=factor(DEOpts$designs))
  in.noiseq <- NOISeq::readData(data = out.cpm, factors = groups)

  # run DE detection
  calc.noiseq <- NOISeq::noiseqbio(in.noiseq,
                                   k = NULL,
                                   norm = "n",
                                   nclust = 15,
                                   plot = FALSE,
                                   factor="Group",
                                   conditions = NULL,
                                   lc = 0, r = 50, adj = 1.5,
                                   a0per = 0.9, filter = 0)
  res <- calc.noiseq@results[[1]]
  res$fdr <- 1-res$prob

  # construct result data frame
  result = data.frame(geneIndex=rownames(res),
                      pval=rep(NA, nrow(res)),
                      fdr=res$fdr,
                      lfc=rep(NA, nrow(res)),
                      stringsAsFactors = F)

  return(result)
}


# DSS ---------------------------------------------------------------------

#' @importFrom DSS newSeqCountSet estNormFactors estDispersion waldTest
#' @importFrom splines ns
#' @importFrom edgeR DGEList
.run.DSS <- function(normData, countData, DEOpts) {

  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  # make input data set
  designs <- ifelse(DEOpts$designs==min(DEOpts$designs), 0, 1)
  cd <- countData
  rownames(cd) <- NULL
  colnames(cd) <- NULL
  seqData <- DSS::newSeqCountSet(counts = cd, designs = designs)
  seqData@normalizationFactor <- sf
  seqData <- DSS::estDispersion(seqData)

  # run DE detection
  res.dss <- suppressWarnings(DSS::waldTest(seqData = seqData,
                                            sampleA = 0, sampleB = 1))
  res.dss <- res.dss[order(res.dss$geneIndex),]
  pval <- res.dss$pval

  # construct result data frame
  result = data.frame(geneIndex=rownames(countData),
                      pval=pval,
                      fdr=rep(NA, nrow(countData)),
                      lfc=res.dss$lfc,
                      stringsAsFactors = F)
  return(result)
}


# EBSeq -------------------------------------------------------------------

#' @importFrom EBSeq MedianNorm EBTest
#' @importFrom edgeR DGEList calcNormFactors
.run.EBSeq <- function(normData, countData, DEOpts) {

  groups <- data.frame(Group=factor(DEOpts$designs))
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  # run DE detection
  calc.ebseq <- suppressMessages(
    EBSeq::EBTest(Data = countData, NgVector = NULL,
                  Conditions = factor(DEOpts$designs),
                  sizeFactors = sf,
                  maxround = 20,
                  Pool = F, NumBin = 1000,
                  ApproxVal = 10^-10, Alpha = NULL,
                  Beta = NULL, PInput = NULL,
                  RInput = NULL, PoolLower = .25,
                  PoolUpper = .75, Print = F,
                  Qtrm = 1,QtrmCut=0))
  fdr <- 1-calc.ebseq$PPDE

  ## construct results
  result = data.frame(geneIndex=rownames(countData),
                      pval=rep(NA, nrow(countData)),
                      fdr=fdr,
                      lfc=rep(NA, nrow(countData)),
                      stringsAsFactors = F)
  return(result)
}


# NBPSeq ------------------------------------------------------------------

#' @importFrom NBPSeq nbp.test
#' @importFrom edgeR DGEList calcNormFactors
.run.NBPSeq <- function(normData, countData, DEOpts) {

  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])

  #run DE testing
  grp.ids <- ifelse(DEOpts$designs==min(DEOpts$designs), 1, 2)
  res <- NBPSeq::nbp.test(counts=countData,
                          grp.ids=grp.ids,
                          grp1=1, grp2=2,
                          norm.factors = sf,
                          lib.sizes = colSums(countData),
                          model.disp = "NBQ", print.level = 0)

  ## construct results
  result = data.frame(geneIndex=rownames(countData),
                      pval=res$pv.alues,
                      fdr=rep(NA, nrow(countData)),
                      lfc=res$log.fc,
                      stringsAsFactors = F)
  return(result)
}


# MAST --------------------------------------------------------------------

#' @importFrom MAST FromMatrix zlm.SingleCellAssay lrTest summary
#' @importFrom S4Vectors mcols
#' @importFrom AnnotationDbi as.list
#' @importFrom edgeR DGEList calcNormFactors cpm.DGEList
#' @importFrom data.table data.table
#' @importFrom reshape2 melt
#' @importFrom parallel mclapply
.run.MAST <- function(normData, countData, DEOpts, NCores) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # 1. size factor normalised log2(CPM+1) values.
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  out.expr <- log2(out.cpm+1)

  # 2.: cell (sample ID, CDR, condition) and gene (gene name) annotation
  ids = colnames(out.expr)
  ngeneson = colSums(out.expr>0)
  cngeneson = ngeneson-mean(ngeneson)
  cond = factor(DEOpts$designs)
  cdat <- data.frame(wellKey=ids,
                     ngeneson=ngeneson,
                     cngeneson=cngeneson,
                     condition=cond,
                     stringsAsFactors = F)
  fdat <- data.frame(primerid=rownames(out.expr),
                     stringsAsFactors = F)
  # 3.: construct MAST single cell assay
  sca <- MAST::FromMatrix(class = "SingleCellAssay",
                    exprsArray=out.expr,
                    cData = cdat,
                    fData = fdat)
  # 4.: Model Fit
  if (!is.null(NCores)) {
    options(mc.cores=NCores)
  }

  zlm <- suppressMessages(
    MAST::zlm(~ condition + cngeneson,
              sca,
              method = "bayesglm",
              ebayes = TRUE,
              ebayesControl = list(method = "MLE", model = "H1"))
  )
  summaryZLM <- suppressMessages(MAST::summary(zlm))
  summaryDt <- summaryZLM$datatable
  # 5.: LRT
  lrt <- suppressMessages(MAST::lrTest(zlm, "condition"))
  # results table extraction
  res_gene <- data.table::data.table(reshape2::melt(lrt))
  res_gene_hurdle <- res_gene[metric=="Pr(>Chisq)" & test.type=="hurdle", .(primerid, value)]
  res_gene_hurdle <- data.frame(res_gene_hurdle, stringsAsFactors = F)
  res_gene_hurdle <- res_gene_hurdle[match(S4Vectors::mcols(sca)$primerid, res_gene_hurdle$primerid),]
  res_lfc_hurdle <- summaryDt[contrast=='condition1' & component=='logFC', .(primerid, coef)]
  res_lfc_hurdle <- data.frame(res_lfc_hurdle, stringsAsFactors = F)
  res_lfc_hurdle <- res_lfc_hurdle[match(S4Vectors::mcols(sca)$primerid, res_lfc_hurdle$primerid),]

  ## construct results
  result <- data.frame(geneIndex=res_gene_hurdle$primerid,
                       pval=res_gene_hurdle$value,
                       fdr=rep(NA, nrow(res_gene_hurdle)),
                       lfc=res_lfc_hurdle$coef,
                       stringsAsFactors = F)
  return(result)
}


# scde --------------------------------------------------------------------

#' @importFrom scde scde.error.models scde.expression.prior scde.expression.difference
#' @importFrom stats pnorm
.run.scde <- function(normData, countData, DEOpts, NCores) {

  # make group vector
  groups <- factor(DEOpts$designs)
  names(groups) <- colnames(countData)

  if(is.null(NCores)) {
    ncores = 1
  }
  if(!is.null(NCores)) {
    ncores = NCores
  }

  # calculate error models
  o.ifm <- scde::scde.error.models(counts = countData,
                                   groups = groups,
                                   n.cores = ncores,
                                   min.count.threshold = 1,
                                   threshold.segmentation = TRUE,
                                   save.crossfit.plots = FALSE,
                                   save.model.plots = FALSE,
                                   verbose = 0)
  # estimate gene expression prior
  o.prior <- scde::scde.expression.prior(models = o.ifm,
                                         counts = countData,
                                         length.out = 400,
                                         show.plot = FALSE)
  # run differential expression tests on all genes.
  ediff <- scde::scde.expression.difference(models=o.ifm,
                                            counts=countData,
                                            prior=o.prior,
                                            groups = groups,
                                            n.cores = ncores,
                                            n.randomizations  =  100,
                                            verbose  =  0)
  pval <- 2 * (1 - pnorm(abs(ediff$Z)))

  ## construct results
  result = data.frame(geneIndex=rownames(ediff),
                      pval=pval,
                      fdr=rep(NA, nrow(countData)),
                      lfc=ediff$ce,
                      stringsAsFactors = F)
  return(result)
}


# BPSC --------------------------------------------------------------------

#' @importFrom BPSC BPglm
#' @importFrom gtools mixedsort
#' @importFrom edgeR DGEList cpm.DGEList
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom stats model.matrix
.run.BPSC <- function(normData, countData, DEOpts, NCores) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # size factor normalised log2(CPM+1) values. Note that the function in scater gave negative values and when cpm.DGEList was allowed to take the log itself all CPMs were nonzero!
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  exprmat <- out.cpm
  group <- DEOpts$designs
  controlIDs <- which(group == min(group))
  design.mat <- stats::model.matrix( ~ group)
  coef <- 2

  if(!is.null(NCores)) {
    doParallel::registerDoParallel(cores=NCores)
    res <- BPglm(data = exprmat,
                 controlIds = controlIDs,
                 design = design.mat,
                 coef = coef,
                 keepFit = TRUE,
                 useParallel=TRUE)
    invisible(capture.output(
      summary_BPglm <- suppressMessages(summary(res))
      ))
    summaryDt <- summary_BPglm$topTable
    res_lfc <- as.data.frame(summaryDt)
    res_lfc <- res_lfc[gtools::mixedsort(rownames(res_lfc)),]
    doParallel::stopImplicitCluster()
  }
  if(is.null(NCores)) {
    res <- BPSC::BPglm(data = exprmat,
                       controlIds = controlIDs,
                       design = design.mat,
                       coef = coef,
                       keepFit = TRUE,
                       useParallel = FALSE)
    invisible(capture.output(summary_BPglm <- suppressMessages(summary(res))))
    summaryDt <- summary_BPglm$topTable
    res <- as.data.frame(summaryDt)
    res <- res[gtools::mixedsort(rownames(res)),]
  }

  # construct result data frame
  result = data.frame(geneIndex=rownames(res),
                      pval=res$`Pr(>|t|)`,
                      fdr=rep(NA, length(rownames(res))),
                      lfc=res$Estimate,
                      stringsAsFactors = F)
  return(result)
}


# monocle -----------------------------------------------------------------

#' @importFrom monocle newCellDataSet differentialGeneTest
#' @importFrom BiocGenerics sizeFactors estimateDispersions "sizeFactors<-"
#' @importFrom VGAM negbinomial.size
#' @importFrom methods new
.run.monocle <- function(normData, countData, DEOpts, NCores) {

  if(!attr(normData, 'normFramework') == 'Census') {
    sf <- normData$size.factors
    sf[sf<0] <- min(sf[sf > 0])
    # make annotated dataframes for monocle
    gene.dat <- data.frame(row.names = rownames(countData),
                           biotype=rep("protein_coding", nrow(countData)),
                           num_cells_expressed=rowSums(countData>0))
    cell.dat <- data.frame(row.names=colnames(countData),
                           Group=DEOpts$designs)
    fd <- new("AnnotatedDataFrame", data = gene.dat)
    pd <- new("AnnotatedDataFrame", data = cell.dat)
    ed <- countData
    # construct cell data set
    cds <- monocle::newCellDataSet(cellData = ed,
                                   phenoData = pd,
                                   featureData = fd,
                                   expressionFamily = VGAM::negbinomial.size())
    sizeFactors(cds) <- sf
    cds <- estimateDispersions(cds,
                               cores = ifelse(is.null(NCores), 1, NCores))
    # run the testing
    diff_test_res <- suppressMessages(
      monocle::differentialGeneTest(cds,
                                    fullModelFormulaStr = "~Group",
                                    reducedModelFormulaStr = "~1",
                                    relative_expr = FALSE,
                                    cores = ifelse(is.null(NCores), 1, NCores),
                                    verbose = FALSE)
    )
    res <- diff_test_res[match(rownames(countData), rownames(diff_test_res)),]
  }

  if(attr(normData, 'normFramework') == 'Census') {
    sf <- normData$size.factors
    # make annotated dataframes for monocle
    gene.dat <- data.frame(row.names = rownames(countData),
                           biotype=rep("protein_coding", nrow(countData)),
                           num_cells_expressed=rowSums(countData>0))
    cell.dat <- data.frame(row.names=colnames(countData),
                           Group=DEOpts$designs)
    fd <- new("AnnotatedDataFrame", data = gene.dat)
    pd <- new("AnnotatedDataFrame", data = cell.dat)
    ed <- normData$RPC
    # construct cell data set
    cds <- monocle::newCellDataSet(cellData = ed,
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit=0.5,
                                   expressionFamily = VGAM::negbinomial.size())
    sizeFactors(cds) <- sf
    cds <- estimateDispersions(cds,
                               relative_expr = TRUE,
                               cores = ifelse(is.null(NCores), 1, NCores))
    # run the testing
    diff_test_res <- suppressMessages(
      monocle::differentialGeneTest(cds,
                                    fullModelFormulaStr = "~Group",
                                    reducedModelFormulaStr = "~1",
                                    relative_expr = TRUE,
                                    cores = ifelse(is.null(NCores), 1, NCores),
                                    verbose = FALSE)
    )
    res <- diff_test_res[match(rownames(countData), rownames(diff_test_res)),]
  }

  ## construct results
  result <- data.frame(geneIndex=rownames(res),
                       pval=res$pval,
                       fdr=rep(NA, nrow(countData)),
                       lfc=rep(NA, nrow(countData)),
                       stringsAsFactors = F)
  return(result)
}


# scDD --------------------------------------------------------------------

#' @importFrom scDD scDD
#' @importFrom edgeR cpm.DGEList
#' @importFrom SummarizedExperiment SummarizedExperiment
.run.scDD <- function(normData, countData, DEOpts, NCores) {
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # 1. size factor normalised log2(CPM+1) values.
  out.cpm <- edgeR::cpm.DGEList(dge, normalized.lib.sizes = T, log = F)
  out.expr <- out.cpm

    # create input data
  exprmat <- out.cpm
  condition <- ifelse(DEOpts$designs==min(DEOpts$designs), 1, 2)
  cell.dat <- data.frame(row.names=colnames(exprmat), condition=condition)
  SCdat <- SummarizedExperiment::SummarizedExperiment(assays=list('NormCounts'=exprmat), colData=cell.dat)

  # DE testing
  if(!is.null(NCores)) {
    res.tmp <- suppressMessages(
      scDD::scDD(SCdat,
                 prior_param = list(alpha = 0.1, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01),
                 permutations = 0,
                 testZeroes = FALSE,
                 adjust.perms = FALSE,
                 param = BiocParallel::MulticoreParam(NCores),
                 parallelBy = "Genes",
                 condition = "condition")
      )
  }
  if(is.null(NCores)) {
    res.tmp <- suppressMessages(
      scDD::scDD(SCdat,
                 prior_param = list(alpha = 0.1, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01),
                 permutations = 0,
                 testZeroes = FALSE,
                 adjust.perms = FALSE,
                 parallelBy = "Genes",
                 condition = "condition")
    )
  }
  res <- scDD::results(res.tmp)

  # construct result data frame
  result = data.frame(geneIndex=as.character(res$gene),
                    pval=res$nonzero.pvalue,
                    fdr=rep(NA, nrow(res)),
                    lfc=rep(NA, nrow(res)),
                    stringsAsFactors = F)
  return(result)
}


# BASiCS ------------------------------------------------------------------

#TODO: Implement testing for BASiCS normalized data

# monocle -----------------------------------------------------------------

#TODO: Implement testing of monocle


# D3E ---------------------------------------------------------------------

# TODO: Do a system call since D3E is written in python

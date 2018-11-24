
# simulateDE --------------------------------------------------------------

#' @name simulateDE
#' @aliases simulateDE
#' @title Simulate Differential Expression
#' @description simulateDE is the main function to simulate differential expression for RNA-seq experiments.
#' The simulation parameters are specified with \code{\link{SimSetup}}.
#' The user needs to specify furthermore
#' the number of samples per group, preprocessing, normalisation and differential testing method.
#' There is also the option to consider spike-ins. \cr
#' The return object contains DE test results from all simulations as well as descriptive statistics.
#' The error matrix calculations will be conducted with \code{\link{evaluateDE}}.\cr
#' @usage simulateDE(n1=c(20,50,100), n2=c(30,60,120),
#' sim.settings,
#' DEmethod,
#' Normalisation,
#' Label = "none",
#' Prefilter = NULL,
#' Impute = NULL,
#' DEFilter = FALSE,
#' spikeIns = FALSE,
#' NCores = NULL,
#' verbose = TRUE)
#' @param n1,n2 Integer vectors specifying the number of biological replicates in each group.
#' Default values are n1=c(20,50,100) and n2=c(30,60,120).
#' @param sim.settings This object specifies the simulation setup. This must be the return object from \code{\link{SimSetup}}.
#' @param DEmethod A character vector specifying the DE detection method to be used.
#' Please consult the Details section for available options.
#' @param Normalisation Normalisation method to use.
#' Please consult the Details section for available options.
#' @param Label A character vector to define whether information about group labels should be used for normalisation.
#' This is only implemented for scran and SCnorm. Possible options include the default \code{"none"} which means that no sample group information is considered for normalisation; \code{"known"} means that the simulated group labels are used and \code{"clustering"} which applies an unsupervised hierarchical clustering to determine the group labels (for details, see \code{\link[scran]{quickCluster}}).
#' @param Prefilter A character vector specifying the gene expression filtering method
#' to be used prior to normalisation (and possibly imputation).
#' Default is \code{NULL}, i.e. no filtering.
#' Please consult the Details section for available options.
#' @param Impute A character vector specifying the gene expression imputation method
#' to be used prior to normalisation.
#' Default is \code{NULL}, i.e. no imputation.
#' Please consult the Details section for available options.
#' @param DEFilter A logical vector indicating whether to run DE testing on filtered and/or imputed count data.
#' Default is \code{FALSE}.
#' @param spikeIns Logical vector to indicate whether to simulate spike-ins.
#' Default is \code{FALSE}.
#' @param NCores integer positive number of cores for parallel processing.
#' Default is \code{NULL}, i.e. 1 core.
#' @param verbose Logical vector to indicate whether to show progress report of simulations.
#' Default is \code{TRUE}.
#' @return A list with the following fields:
#' \item{pvalue, fdr}{3D array (ngenes * N * nsims) for p-values and FDR from each simulation.
#' Note that FDR values will be empty and the calculation will be done by \code{\link{evaluateDE}} whenever applicable.}
#' \item{mu,disp,dropout}{3D (ngenes * N * nsims) array for mean, dispersion and dropout of library size factor normalized read counts.}
#' \item{true.mu,true.disp,true.dropout}{3D (ngenes * N * nsims) array for true mean, dispersion and dropout of simulated read counts.}
#' \item{true.depth}{True simulated sequencing depth per sample.}
#' \item{est.sf}{Global library size factor estimates per sample.}
#' \item{est.gsf}{3D array (ngenes * N * nsims) for size factor estimates. These are gene- and sample-wise estimates and only for SCnorm and Linnorm normalisation.}
#' \item{elfc,rlfc}{3D array (ngenes * N * nsims) for log2 fold changes (LFC):
#' elfc is for the DE tool estimated LFC; rlfc is for the LFC estimated from the normalised read counts.}
#' \item{sim.settings}{The input sim.settings to which the specifications of \code{simulateDE} is added.}
#' \item{time.taken}{The time taken for each simulation, given for preprocessing, normalisation, differential expression testing and moment estimation.}
#' @seealso \code{\link{estimateParam}},  \code{\link{insilicoNBParam}} for negative binomial parameter specifications;\cr
#'  \code{\link{DESetup}}, \code{\link{SimSetup}} for simulation setup;\cr
#'  \code{\link{evaluateDE}} for DE evaluation.
#' @details
#' Here you can find detailed information about preprocessing, imputation, normalisation and differential testing choices.
#' @section Prefiltering prior to imputation/normalisation:
#' \describe{
#' \item{CountFilter}{removes genes that have a mean expression below 0.2.}
#' \item{FreqFilter}{removes genes that have more than 80 percent dropouts.}
#' }
#' @section Imputation prior to normalisation:
#' \describe{
#' \item{scImpute}{employs scImpute method of imputing dropouts. Imputation is only carried out for genes with more than 50 percent dropout. Please consider multiple cores to speed up computation.}
#' \item{DrImpute}{employs DrImpute method of imputing dropouts as implemented
#' in \code{\link[DrImpute]{DrImpute}}.}
#' \item{SAVER}{employs SAVER method of imputing dropouts as implemented
#' in \code{\link[SAVER]{saver}}. Imputation is only carried out for genes with more than 50 percent dropout.}
#' \item{Seurat}{employs Seurat method of imputing dropouts as implemented
#' in \code{\link[Seurat]{AddImputedScore}} using variable genes identified with \code{\link[Seurat]{FindVariableGenes}}. Imputation is only carried out for genes with more than 50 percent dropout.}
#' \item{scone}{employs scone method of imputing dropouts as implemented
#' in \code{\link[scone]{scone}} using estimated dropout probabilities of \code{\link[scone]{estimate_ziber}}.}
#' }
#' @section Normalisation applied to (imputed) read count matrix:
#' \describe{
#' \item{TMM, UQ}{employ the edgeR style normalization of weighted trimmed mean of M-values and upperquartile
#' as implemented in \code{\link[edgeR]{calcNormFactors}}, respectively.}
#' \item{MR, PosCounts}{employ the DESeq2 style normalization of median ratio method and a modified geometric mean method
#' as implemented in \code{\link[DESeq2]{estimateSizeFactors}}, respectively. Spike-ins can also be supplied for both methods via \code{spikeData}.}
#' \item{scran, SCnorm}{apply the deconvolution and quantile regression normalization methods developed for sparse RNA-seq data
#' as implemented in \code{\link[scran]{computeSumFactors}} and \code{\link[SCnorm]{SCnorm}}, respectively. Spike-ins can also be supplied for both methods via \code{spikeData}. Note, however that this means for scran that the normalisation as implemented in \code{\link[scran]{computeSpikeFactors}} is also applied to genes (\code{general.use=TRUE}). Please consider multiple cores to speed up computation for SCnorm.}
#' \item{Linnorm}{apply the normalization method for sparse RNA-seq data
#' as implemented in \code{\link[Linnorm]{Linnorm.Norm}}.
#' For \code{Linnorm}, the user can also supply \code{spikeData}.}
#' \item{RUV}{removes unwanted variation. There are two approaches implemented:
#' (1) utilizing negative control genes, i.e. spike-ins stored in \code{spikeData} (\code{\link[RUVSeq]{RUVg}}).
#' (2) utilizing replicate samples, i.e. samples for which the covariates of interest are considered constant.
#' This annotation is stored in \code{batchData} (\code{\link[RUVSeq]{RUVs}}).}
#' \item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
#' Census at least needs \code{Lengths} for single-end data and preferably \code{MeanFragLengths} for paired-end data.
#' Do not use this algorithm for UMI data!}
#' \item{depth}{Sequencing depth normalisation.}
#' }
#' @section Differential testing using raw read count matrix:
#' \describe{
#' \item{T-Test}{A T-Test per gene is applied using log2 transformed and normalized expression values (i.e. CPM or TPM).}
#' \item{limma-trend, limma-voom}{apply differential testing as implemented in \code{\link[limma]{lmFit}}
#' followed by \code{\link[limma]{eBayes}} on counts transformed by \code{\link[limma]{voom}} or by applying mean-variance trend on log2 CPM values in \code{\link[limma]{eBayes}}.}
#' \item{edgeR-LRT, edgeR-QL}{apply differential testing as implemented in \code{\link[edgeR]{glmFit}}, \code{\link[edgeR]{glmLRT}} and\code{\link[edgeR]{glmQLFit}}, \code{\link[edgeR]{glmQLFTest}}, respectively.}
#' \item{DESeq2}{applies differential testing as implemented in \code{\link[DESeq2]{DESeq}}.}
#' \item{ROTS}{applies differential testing as implemented in \code{\link[ROTS]{ROTS}} with 100 permutations on transformed counts (\code{\link[limma]{voom}}).}
#' \item{baySeq}{applies differential testing as implemented in \code{\link[baySeq]{getLikelihoods}} based on negative binomial prior estimates (\code{\link[baySeq]{getPriors.NB}}).}
#' \item{NOISeq}{applies differential testing as implemented in \code{\link[NOISeq]{noiseqbio}} based on CPM values.}
#' \item{EBSeq}{applies differential testing as implemented in \code{\link[EBSeq]{EBTest}}.}
#' \item{MAST}{applies differential testing as implemented in \code{\link[MAST]{zlm}} for zero-inflated model fitting followed by \code{\link[MAST]{lrTest}} on log2 CPM values.}
#' \item{scde}{applies differential testing as implemented in \code{\link[scde]{scde.expression.difference}}.}
#' \item{BPSC}{applies differential testing as implemented in \code{\link[BPSC]{BPglm}} on CPM values.}
#' \item{scDD}{applies differential testing as implemented in \code{\link[scDD]{scDD}} on CPM values.}
#' \item{DECENT}{applies differential testing as implemented in \code{\link[DECENT]{decent}}.}
#' \item{edgeR-zingeR, DESeq2-zingeR}{In a first step, the posterior probabilities of the zero-inflated negative binomial component are estimated (see \code{\link[zingeR]{zeroWeightsLS}}) and used to define a weight matrix for dispersion estimation in \code{\link[edgeR]{estimateDisp}}. For the edgeR approach, the generalized model as implemented in \code{\link[edgeR]{glmFit}} is fitted. This is followed by an adapted LRT for differential testing to account for the weighting (see \code{\link[zingeR]{glmWeightedF}}). For DESeq2, the generalized linear model coefficients are estimated using \code{\link[DESeq2]{nbinomWaldTest}} and the weighting is done by setting the degrees of freedom for the T distribution.}
#' \item{edgeR-ZINB-WaVE, DESeq2-ZINB-WaVE}{In a first step, a zero-inflated negative binomial regression model  is fitted (see \code{\link[zinbwave]{zinbFit}}) to estimate observational weights (see \code{\link[zinbwave]{computeObservationalWeights}}) used for dispersion estimation in \code{\link[edgeR]{estimateDisp}}. For the edgeR approach, the generalized model as implemented in \code{\link[edgeR]{glmFit}} is fitted. This is followed by an adapted LRT for differential testing to account for the weighting (see \code{\link[zinbwave]{glmWeightedF}}). For DESeq2, the generalized linear model coefficients are estimated using \code{\link[DESeq2]{nbinomWaldTest}} and the weighting is done by setting the degrees of freedom for the T distribution.}
#' }
#' @examples
#' \dontrun{
#' ## define DE parameters and set up simulations
#' de.opts <- DESetup(ngenes = 10000, nsims = 25,
#' p.DE = 0.2, pLFC = function(x) sample(c(-1,1), size=x,replace=TRUE)*rgamma(x, 3, 3),
#' p.B=0.1, bLFC = function(x) rnorm(x, mean=0, sd=1.5), bPattern="uncorrelated",
#' sim.seed = 43856)
#' sim.opts <- SimSetup(desetup = de.opts,
#' params = kolodziejczk_param,
#' spike=NULL, size.factors = "equal",
#' downsample = FALSE, geneset = FALSE)
#' ## run simulations
#' sim.res <- simulateDE(n1=c(50,96,384), n2=c(80,96,384),
#' sim.settings = sim.opts,
#' DEmethod = 'limma-trend',
#' Normalisation = 'scran',
#' Preclust = FALSE,
#' Prefilter = "FreqFilter",
#' Impute = NULL,
#' spikeIns = FALSE,
#' NCores = NULL,
#' verbose = TRUE)
#' }
#' @author Beate Vieth
#' @rdname simulateDE
#' @importFrom stats setNames
#' @export
simulateDE <- function(n1=c(20,50,100), n2=c(30,60,120),
                       sim.settings,
                       DEmethod,
                       Normalisation,
                       Label = "none",
                       Prefilter = NULL,
                       Impute = NULL,
                       DEFilter = FALSE,
                       spikeIns = FALSE,
                       NCores = NULL,
                       verbose = TRUE) {
  if (!length(n1) == length(n2)) { stop("n1 and n2 must have the same length!") }
  if(isTRUE(spikeIns) && is.null(sim.settings$spike)) {
    stop(message(paste0("For the simulation of  spike-ins, fitting information is needed but there is no 'spike' object in 'sim.settings'.  Please consult the function estimateSpike for spike fitting and SimSetup for creating simulation setup object!")))
  }

  if(!is.null(NCores) && DEmethod %in% c("edgeR-LRT", "edgeR-QL", "edgeR-zingeR", "DESeq2-zingeR", 'limma-voom', "limma-trend", "NOISeq", "EBSeq", "ROTS")) {
    if(verbose) {message(paste0(DEmethod, " has no parallel computation option!"))}
  }

  if(sim.settings$RNAseq == "singlecell" && DEmethod %in% c("edgeR-LRT", "edgeR-QL", "limma-voom", "limma-trend", "DESeq2", "baySeq", "NOISeq", "EBSeq")) {
    if(verbose) {message(paste0(DEmethod, " is developed for bulk RNA-seq experiments."))}
  }

  if(sim.settings$RNAseq == "bulk" && DEmethod %in% c("MAST", 'BPSC', "edgeR-zingeR", "DESeq2-zingeR", "edgeR-ZINB-WaVE", "DESeq2-ZINB-WaVE")) {
    if(verbose) {message(paste0(DEmethod, " is developed for single cell RNA-seq experiments."))}
  }

  if(sim.settings$RNAseq == "bulk" && DEmethod %in% c('scde', 'scDD', 'monocle', 'DECENT')) {
    stop(message(paste0(DEmethod, " is only developed and implemented for single cell RNA-seq experiments.")))
  }

  if(DEmethod == "DECENT") {
    if(verbose) {message(paste0(DEmethod, " does not require additional normalisation nor imputation."))}
    Normalisation = "none"
    Prefilter = NULL
    Impute = NULL
  }

  if(c(all(is.null(Impute), is.null(Prefilter)) && isTRUE(DEFilter))) {
    stop(message(paste0("You wish to use imputed/filtered gene expression values for DE testing but you did not specify the imputation/filtering method. Aborting.")))
  }

  # define the maximal count matrix for simulations
  max.n = max(n1,n2)
  min.n = min(n1, n2)

  # append additional settings of simulateDE to sim.settings
  sim.settings$n1 = n1
  sim.settings$n2 = n2
  sim.settings$Normalisation = Normalisation
  sim.settings$DEmethod = DEmethod
  sim.settings$spikeIns = spikeIns
  sim.settings$NCores = NCores
  sim.settings$Label = Label
  sim.settings$Prefilter = Prefilter
  sim.settings$Impute = Impute
  sim.settings$DEFilter = DEFilter
  sim.settings$clustNumber = ifelse(sim.settings$design=="2grp", 2, NULL)
  if(Label == "clustering") {PreclustNumber <- min.n}
  if(!Label == "clustering") {PreclustNumber <- NULL}
  sim.settings$PreclustNumber = PreclustNumber

  if (verbose) { message(paste0("Preparing output arrays.")) }

  my.names = paste0(n1,"vs",n2)

  #set up output arrays
  pvalues = fdrs = elfcs = rlfcs = mus = true.mus = disps = true.disps = dropouts = true.drops = array(NA,dim=c(sim.settings$ngenes,length(n1), sim.settings$nsims))

  time.taken = array(NA,dim = c(4,length(n1), sim.settings$nsims),
                     dimnames = list(c('Preprocess', "Normalisation", "DE", "Moments"),
                                     NULL, NULL))

  true.sf = stats::setNames(replicate(length(n1),NULL),my.names)
  est.sf = stats::setNames(replicate(length(n1),NULL),my.names)
  true.depth = stats::setNames(replicate(length(n1),NULL),my.names)
  true.sf <- lapply(1:length(true.sf), function(x) {
    true.sf[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })
  est.sf <- lapply(1:length(est.sf), function(x) {
    est.sf[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })
  true.depth <- lapply(1:length(true.depth), function(x) {
    true.depth[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })

  if(sim.settings$Normalisation %in% c("SCnorm", "Linnorm")) {
    est.gsf = stats::setNames(replicate(length(n1),NULL),my.names)
    est.gsf <- lapply(1:length(est.gsf), function(x) {
      est.gsf[[x]] = array(NA,dim=c(sim.settings$nsims, sim.settings$ngenes, n1[x] + n2[x]))
    })
  }
  if(!sim.settings$Normalisation %in% c("SCnorm", "Linnorm")) {
    est.gsf = NULL
  }

  true.designs = stats::setNames(replicate(length(n1),NULL),my.names)
  true.designs <- lapply(1:length(true.designs), function(x) {
    true.designs[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })

  ## start simulation
  for (i in 1:sim.settings$nsims) {
    if (verbose) { message(paste0("Simulation number ", i, "\n")) }
    ## update the simulation options by extracting the ith set and change sim.seed
    tmp.simOpts = sim.settings
    tmp.simOpts$DEid = tmp.simOpts$DEid[[i]]
    tmp.simOpts$pLFC = tmp.simOpts$pLFC[[i]]
    tmp.simOpts$bLFC = tmp.simOpts$bLFC[[i]]
    tmp.simOpts$sim.seed = tmp.simOpts$sim.seed[[i]]

    ## generate gene read counts
    if (verbose) {message(paste0("Generating gene read counts")) }
    gene.data = .simRNAseq.2grp(simOptions = tmp.simOpts,
                                n1 = max.n, n2 = max.n, verbose=verbose)

    ## generate spike-in read counts
    if(isTRUE(spikeIns)) {
      if (verbose) { message(paste0("Generating spike-in read counts")) }
      spike.data = .simSpike(SpikeOptions = tmp.simOpts$spike,
                             n1 = max.n, n2 = max.n, sf = gene.data$sf)
    }
    if(!isTRUE(spikeIns)) {
      spike.data = NULL
      spike.info = NULL
    }

    ## generate mean fragment lengths for samples
    if(!is.null(tmp.simOpts$MeanFragLengths)) {
      if (verbose) { message(paste0("Sampling from observed mean fragment lengths")) }
      MeanFrag.data = sample(tmp.simOpts$MeanFragLengths, max.n+max.n, replace = TRUE)
      names(MeanFrag.data) = colnames(gene.data$counts)
    }
    if(is.null(tmp.simOpts$MeanFragLengths)) {
      MeanFrag.data = NULL
    }

    ## match sampled gene names with given gene lengths
    if(!is.null(tmp.simOpts$Lengths)) {
      gene.id = sub('_([^_]*)$', '', rownames(gene.data$counts))
      Length.data = tmp.simOpts$Lengths
      Length.data = Length.data[match(gene.id,names(Length.data))]
    }
    if(is.null(tmp.simOpts$Lengths)) {
      Length.data = NULL
    }

    ##  for different sample sizes
    for (j in seq(along=n1)) {
      Nrep1 = n1[j]
      Nrep2 = n2[j]
      if (verbose) { message(paste0(Nrep1, " vs. ", Nrep2)) }

      ## take a subsample of simulated samples
      idx = c(1:Nrep1, max.n + (1:Nrep2))
      true.design = gene.data$designs[idx]

      ## take a subsample of the simulated read counts
      sim.cnts = gene.data$counts[,idx]
      ## take a subsample of the true size factors
      gene.sf = gene.data$sf[idx]
      ## filter out zero expression genes
      ix.valid = rowSums(sim.cnts) > 0
      count.data = sim.cnts[ix.valid,, drop = FALSE]

      ## match sampled gene names with given gene lengths
      if(!is.null(Length.data)) {
        if (verbose) { message(paste0("Associating gene lengths with sampled gene expression")) }
        gene.id = sub('_([^_]*)$', '', rownames(count.data))
        length.data = Length.data
        length.data = length.data[match(gene.id,names(length.data))]
      }
      if(is.null(Length.data)) {
        length.data = NULL
      }

      ## take a subsample of simulated spike-ins
      if(!is.null(spike.data)) {
        sim.spike <- spike.data$counts
        spike.valid = rowSums(sim.spike) > 0
        count.spike = sim.spike[spike.valid, idx, drop=FALSE]
        if(!is.null(sim.settings$spike)) {
          spike.info <- tmp.simOpts$spike$Input$spikeInfo[rownames(tmp.simOpts$spike$Input$spikeInfo)  %in% rownames(count.spike), , drop = FALSE]
          spike.info <- spike.info[match(rownames(count.spike), rownames(spike.info)), , drop = FALSE]
        }
      }
      if(is.null(spike.data)) {
        count.spike = NULL
      }
      ## take a subsample of mean fragment lengths
      if(!is.null(MeanFrag.data)) {
        meanfrag.data = MeanFrag.data[idx]
      }
      if(is.null(MeanFrag.data)) {
        meanfrag.data = NULL
      }

      def.design <- true.design

      ## perform filtering / imputation (OPTIONAL)
      start.time.preprocess <- Sys.time()
      if(!is.null(Prefilter)) {
        if (verbose) { message(paste0("Applying ",Prefilter," prefiltering")) }
        filter.data <- .prefilter.calc(Prefilter=tmp.simOpts$Prefilter,
                                       countData=count.data,
                                       NCores=tmp.simOpts$NCores)
        filter.count.data <- filter.data
        if(!is.null(Length.data)) {
          gene.id = sub('_([^_]*)$', '', rownames(filter.count.data))
          length.data = Length.data
          length.data = length.data[match(gene.id,names(length.data))]
        }
        if(is.null(Length.data)) {
          length.data = NULL
        }
      }
      if(is.null(Prefilter)) {
        filter.count.data <- count.data
      }
      if(!is.null(Impute)) {
        if (verbose) { message(paste0("Applying ", Impute, " imputation")) }
        impute.data <- .impute.calc(Impute=tmp.simOpts$Impute,
                                    countData=filter.count.data,
                                    spikeData=count.spike,
                                    batchData=def.design,
                                    clustNumber=tmp.simOpts$clustNumber,
                                    Lengths = length.data,
                                    MeanFragLengths = meanfrag.data,
                                    NCores=tmp.simOpts$NCores,
                                    verbose=verbose)
        fornorm.count.data <- impute.data
        if(!is.null(Length.data)) {
          gene.id = sub('_([^_]*)$', '', rownames(fornorm.count.data))
          length.data = Length.data
          length.data = length.data[match(gene.id,names(length.data))]
        }
        if(is.null(Length.data)) {
          length.data = NULL
        }
      }
      if(is.null(Impute)) {
        fornorm.count.data <- filter.count.data
      }
      end.time.preprocess <- Sys.time()

      ## perform normalisation
      if (verbose) { message(paste0("Applying ", Normalisation, " normalisation")) }
      start.time.norm <- Sys.time()
      norm.data <- .norm.calc(Normalisation=tmp.simOpts$Normalisation,
                              sf=gene.sf,
                              countData=fornorm.count.data,
                              spikeData=count.spike,
                              spikeInfo=spike.info,
                              batchData=def.design,
                              Lengths=length.data,
                              MeanFragLengths=meanfrag.data,
                              PreclustNumber=tmp.simOpts$PreclustNumber,
                              Label=tmp.simOpts$Label,
                              NCores=tmp.simOpts$NCores,
                              verbose=verbose)
      end.time.norm <- Sys.time()

      ## create an DE options object to pass into DE detection
      DEOpts <- list(designs=def.design, p.DE=tmp.simOpts$p.DE)

      ## rematch sampled gene names with given gene lengths
      if(!is.null(Length.data)) {
        if (verbose) { message(paste0("Reassociate gene lengths with sampled gene expression")) }
        gene.id = sub('_([^_]*)$', '', rownames(count.data))
        length.data = Length.data
        length.data = length.data[match(gene.id,names(length.data))]
      }
      if(is.null(Length.data)) {
        length.data = NULL
      }

      ## Run DE detection
      start.time.DE <- Sys.time()
      if(!isTRUE(tmp.simOpts$DEFilter)) {
        if (verbose) { message(paste0("Applying ", DEmethod,
                                      " for DE analysis on raw count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$DEmethod,
                          normData=norm.data,
                          countData=count.data,
                          Lengths=length.data,
                          MeanFragLengths=meanfrag.data,
                          DEOpts=DEOpts,
                          spikeData=count.spike,
                          spikeInfo=spike.info,
                          NCores=tmp.simOpts$NCores,
                          verbose=verbose)
      }
      if(isTRUE(tmp.simOpts$DEFilter)) {
        if (verbose) { message(paste0("Applying ", DEmethod,
                                      " for DE analysis on imputed/filtered count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$DEmethod,
                          normData=norm.data,
                          countData=fornorm.count.data,
                          Lengths=length.data,
                          DEOpts=DEOpts,
                          spikeData=count.spike,
                          spikeInfo=spike.info,
                          NCores=tmp.simOpts$NCores,
                          verbose=verbose)
      }
      end.time.DE <- Sys.time()

      if(tmp.simOpts$DEmethod =="DECENT") {
        res.de <- res.de[["DEresults"]]
        norm.data <- res.de[["NormData"]]
      }

      # move up to use in param estimation
      if(attr(norm.data, 'normFramework') %in% c("SCnorm", "Linnorm")) {
        allgenes <- rownames(sim.cnts)
        testedgenes <- rownames(norm.data$scale.factors)
        ixx.valid <- allgenes %in% testedgenes
        est.gsf[[j]][i, ixx.valid, ] = norm.data$scale.factors
        norm.data$scale.factors = est.gsf[[j]][i,,]
      }

      ## estimate moments of read counts simulated
      start.time.NB <- Sys.time()
      res.params <- .run.params(countData=count.data,
                                normData=norm.data,
                                group=DEOpts$designs)
      seq.depth <- colSums(count.data)
      end.time.NB <- Sys.time()

      # generate empty vectors
      pval = fdr = est.lfc = raw.lfc = mu.tmp = true.mu.tmp = disp.tmp = true.disp.tmp = p0.tmp = true.p0.tmp = rep(NA, nrow(sim.cnts))
      # indicator of tested genes
      allgenes <- rownames(sim.cnts)
      testedgenes <- res.de$geneIndex
      ixx.valid <- allgenes %in% testedgenes
      ## extract results of DE testing
      pval[ixx.valid] = res.de$pval
      fdr[ixx.valid] = res.de$fdr
      est.lfc[ixx.valid] = res.de$lfc
      # extract parameter estimates
      raw.lfc[ix.valid] = res.params$lfc
      mu.tmp[ix.valid] = res.params$means
      disp.tmp[ix.valid] = res.params$dispersion
      p0.tmp[ix.valid] = res.params$dropout
      # extract true parameters
      true.mu.tmp = gene.data$mus
      true.disp.tmp = gene.data$disps
      true.p0.tmp = gene.data$drops
      # copy it in 3D array of results
      pvalues[,j,i] = pval
      fdrs[,j,i] = fdr
      elfcs[,j,i] = est.lfc
      rlfcs[,j,i] = raw.lfc
      mus[,j,i] = mu.tmp
      true.mus[,j,i] = true.mu.tmp
      disps[,j,i] = disp.tmp
      true.disps[,j,i] = true.disp.tmp
      dropouts[,j,i] = p0.tmp
      true.drops[,j,i] = true.p0.tmp
      true.sf[[j]][i,] = gene.sf
      est.sf[[j]][i,] = norm.data$size.factors
      true.depth[[j]][i,] = seq.depth

      true.designs[[j]][i,] = true.design

      # time taken for each step
      # copy designs into list of matrices
      if(any(c(!is.null(Prefilter), !is.null(Impute)))) {
        time.taken.preprocess <-  difftime(end.time.preprocess,
                                           start.time.preprocess,
                                           units="mins")
      }
      if(all(c(is.null(Prefilter), is.null(Impute)))) {
        time.taken.preprocess = NA
      }

      time.taken.norm <- difftime(end.time.norm,
                                  start.time.norm,
                                  units="mins")
      time.taken.DE <- difftime(end.time.DE,
                                start.time.DE,
                                units="mins")
      time.taken.NB <- difftime(end.time.NB,
                                start.time.NB,
                                units="mins")
      timing <- rbind(time.taken.preprocess,
                      time.taken.norm,
                      time.taken.DE,
                      time.taken.NB)

      # copy time taken in 2D array of time taken
      time.taken[,j,i] = timing

    }
  }

  ## return
  res.out <- list(pvalue = pvalues,
                  fdr = fdrs,
                  elfc = elfcs,
                  rlfc = rlfcs,
                  mu = mus,
                  true.mu = true.mus,
                  disp = disps,
                  true.disp = true.disps,
                  dropout = dropouts,
                  true.dropout = true.drops,
                  true.sf = true.sf,
                  true.depth = true.depth,
                  est.sf = est.sf,
                  est.gsf = est.gsf,
                  true.designs = true.designs,
                  time.taken = time.taken,
                  sim.settings = sim.settings)

  attr(res.out, 'Simulation') <- "DE"
  return(res.out)
}


# simulateCounts ----------------------------------------------------------

#' @name simulateCounts
#' @aliases simulateCounts
#' @title Simulate Read Counts
#' @description With this function, the user can simulate realistic read counts for genes and spike-ins across two and multiple groups of samples (cells).
#' @usage simulateCounts(n=c(20,100,30,25,500), ngenes=10000,
#' p.DE=0.1, pLFC,
#' p.B=NULL, bLFC=NULL, bPattern="uncorrelated",
#' p.M=NULL, mLFC=NULL,
#' params, size.factors='equal',
#' spike=NULL, spikeIns=FALSE,
#' downsample=FALSE, geneset=FALSE,
#' sim.seed=NULL, verbose=TRUE)
#' @param n The vector of sample groups with n samples, \code{e.g c(10,50,20,16)}.
#' @param ngenes The total number of genes to simulate. Default is \code{10000}.
#' @param p.DE Numeric vector between 0 and 1 representing
#' the percentage of genes being differentially expressed due to phenotype,
#' i.e. biological signal. Default is \code{0.1}.
#' @param pLFC The log2 phenotypic fold changes for DE genes. (1) For two group simulations, this can be:
#' (a) a constant, e.g. 2;
#' (b) a vector of values with length being number of DE genes.
#' If the input is a vector and the length is not the number of DE genes,
#' it will be sampled with replacement to generate log-fold change;
#' (c) an univariate function that takes an integer n and generates vector(s) of length n,
#' e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' (2) For multigroup simulations, this can be:
#' (a) a list with number of elements equal to number of groups to simulate.
#' The element of the list are vectors of log2 fold changes.
#' (b) a multivariate function that takes an integer n
#' and generates a dataframe with number of columns equal to number of groups.
#' e.g. function(x) mvtnorm::rmvnorm(x, mean=c(4,2,1), sigma = matrix(c(4,2,2,2,3,2, 2, 2, 5), ncol=3)).
#' @param p.B Numeric vector between 0 and 1 representing the percentage of genes
#' being differentially expressed between batches. Default is \code{NULL}, i.e. no batch effect.
#' @param bLFC The log2 batch fold change for all genes. This can be:
#' (1) a constant, e.g. 2;
#' (2) a vector of values with length being number of all genes.
#' If the input is a vector and the length is not the number of total genes,
#' it will be sampled with replacement to generate log2 fold changes;
#' (3) an univariate function that takes an integer n, and generates a vector of length n,
#' e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' Note that only two batches will be simulated
#' irrespective of the number of phenotypic groups defined in pLFC.
#' @param bPattern Character vector for batch effect pattern. Possible options include:
#' "uncorrelated", "orthogonal" and " correlated". Default is \code{"uncorrelated"}.
#' @param p.M Numeric vector between 0 and 1 representing the percentage of genes
#' being differentially expressed exclusively in one group, i.e. marker genes. Default is \code{NULL}.
#' @param mLFC The log2 batch fold change for marker genes. This can be:
#' (1) a constant, e.g. 2;
#' (2) a vector of values with length being number of marker genes.
#' If the input is a vector and the length is not the number of marker genes,
#' it will be sampled with replacement to generate log2 fold changes;
#' (3) a function that takes an integer n, and generates a vector of length n,
#' e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' @param params The distributional parameters for simulations of genes,
#' i.e. the output of \code{\link{estimateParam}}.
#' @param size.factors Size factors representing sample-specific differences/biases in expected mean values of the NB distribution:
#' "equal" or "given". The default is \code{"equal"}, i.e. equal size factor of 1.
#' If the user defines it as given, the size factors are sampled from the size factors ("sf") provided by the output of \code{\link{estimateParam}}.
#' @param spike The distributional parameters for simulations of spike-ins,
#' i.e. the output of \code{\link{estimateSpike}}.
#' @param spikeIns Logical value to indicate whether to simulate spike-ins.
#' Default is \code{FALSE}.
#' @param downsample Drawing the associated dispersions after determining effective mean expressions by size factors.
#' Default is \code{FALSE}, i.e. using the true mean expression values.
#' @param geneset Sampling with replacement or filling count tables with low magnitude Poisson
#' when the estimated mean expression vector is shorter than the number of genes to be simulated.
#' Default is \code{FALSE}, i.e. random sampling of mean expression values with replacement.
#' @param sim.seed Simulation seed.
#' @param verbose Logical value to indicate whether to show progress report of simulation.
#' Default is \code{TRUE}.
#' @return List with the following vectors:
#' \item{GeneCounts}{The simulated read count matrix for genes with row=genes and columns=samples.}
#' \item{SpikeCounts}{The simulated read count matrix for spike-ins with row=spike-ins and columns=samples.}
#' \item{DEid}{A vector (length=ngenes*p.DE) for the IDs of phenotypic DE genes.}
#' \item{Bid}{A vector (length=ngenes*p.B) for the IDs of batch DE genes.}
#' \item{Mid}{A vector (length=ngenes*p.M) for the IDs of marker genes.}
#' \item{pLFC}{A vector / matrix (columns = length(n); rows = ngenes) for phenotypic log fold change of all genes, ie nonDE=0 and DE=plfc.}
#' \item{bLFC}{A vector / matrix (columns = length(n); rows = ngenes) for phenotypic log fold change of all genes, ie nonDE=0 and DE=plfc.}
#' \item{mLFC}{A vector / matrix (columns = length(n); rows = ngenes) for phenotypic log fold change of all genes, ie nonDE=0 and DE=plfc.}
#' \item{ngenes, nsims, p.DE, p.B, p.M, sim.seed, n, k}{Input parameters.}
#'
#' @examples
#' \dontrun{
#' ## define log2 fold changes
#' p.foo <- function(x) mvtnorm::rmvnorm(x, mean=c(4,2,1),
#' sigma = matrix(c(4,2,2,2,3,2, 2, 2, 5), ncol=3))
#' b.foo <- function(x) rnorm(x, mean=0, sd=1.5)
#' ## simulate 3 groups of cells
#' simcounts <- simulateCounts(n=c(100,110,90), ngenes=10000,
#' p.DE=0.05, pLFC = p.foo,
#' p.B=0.1, bLFC=b.foo, bPattern="uncorrelated",
#' p.M=NULL, mLFC=NULL,
#' params=kolodziejczk_param,
#' size.factors="equal",
#' spike=NULL, spikeIns=FALSE,
#'  downsample=FALSE, geneset=FALSE,
#'  sim.seed=34628, verbose=TRUE)
#' }
#' @author Beate Vieth
#' @rdname simulateCounts
#' @export
simulateCounts <- function(n=c(20,100,30,25,500),
                           ngenes=10000,
                           p.DE=0.1, pLFC,
                           p.B=NULL, bLFC=NULL, bPattern="uncorrelated",
                           p.M=NULL, mLFC=NULL,
                           params,
                           size.factors="equal",
                           spike=NULL,
                           spikeIns=FALSE,
                           downsample=FALSE,
                           geneset=FALSE,
                           sim.seed=NULL,
                           verbose=TRUE) {

  ## setup preparation
  if (verbose) { message(paste0("Preparing setup.")) }

  # set simulation seed
  if(is.null(sim.seed)) { sim.seed = sample(1:1000000, size = 1) }
  set.seed(sim.seed)
  # define the number of cell populations
  k = length(n)
  # name the groups of cells otherwise issue when groups have same n!
  names(n) <- paste0("n", 1:k)

  # set gene expr settings
  nDE = round(ngenes*p.DE)
  if(!is.null(p.B)) { nB = round(ngenes*p.B) }
  if(!is.null(p.M)) { nM = round(ngenes*p.M) }

  DEids = Bids = Mids = plfcs = blfcs = mlfcs = lfcs = ids = NULL

  # generate a random id for DE genes
  DEids <- sample(ngenes, nDE, replace = FALSE)
  ## generate lfc for all genes: 0 for nonDE and LFC for DE
  pdim <- ifelse(k==2,1,k)
  plfcs = matrix(0, nrow = ngenes, ncol = pdim)
  plfcs[DEids,] = .setFC(pLFC, nDE, k = k)
  if(!is.null(bLFC)) {
    # generate a random id for batch affected genes
    Bids <- sample(ngenes, nB, replace = FALSE)
    ## generate batch lfc for all genes: 0 for nonDE and LFC for DE, only two batches up until now
    blfcs = matrix(0, nrow = ngenes, ncol = 1)
    blfcs[Bids,] = .setFC(bLFC, nB, k=2)
  }
  if(!is.null(mLFC)) {
    leftgenes <- c(1:ngenes)[-DEids]
    ## generate lfc for marker genes: LFC for expressing group and 0 for all others
    mlfc = .setMarker(input=mLFC, nDEgenes=nM, idpool=leftgenes, ngenes=ngenes, k=k)
    mlfcs = mlfc$lfcs
    Mids = mlfc$id
    lfcs = plfcs
    lfcs[Mids,] = mlfcs[Mids,]
    ids = c(DEids, Mids)
  }
  if(is.null(mLFC)) {
    lfcs = plfcs
    ids = DEids
    Mids = NULL
    mlfcs = NULL
  }

  sim.settings = c(params,
                   list(DEid = ids,
                        pLFC = lfcs,
                        bLFC = blfcs,
                        bPattern = bPattern,
                        sim.seed = sim.seed,
                        spike = spike,
                        spikeIns = spikeIns,
                        RNAseq = params$RNAseq,
                        ngenes = ngenes,
                        downsample = downsample,
                        geneset = geneset,
                        size.factors = size.factors))
  attr(sim.settings, 'Distribution') = attr(params, 'Distribution')
  attr(sim.settings, 'param.type') = attr(params, 'param.type')

  ## set up output arrays
  if (verbose) { message(paste0("Preparing output arrays.")) }

  # gene read counts
  genes.mat = matrix(0, nrow = ngenes, ncol = sum(n))

  ## simulation
  if (verbose) { message(paste0("Simulating gene read counts \n")) }

  if(k==2) {
    gene.simdata = .simRNAseq.2grp(simOptions = sim.settings, n1=n[1], n2 = n[2], verbose=verbose)
    gene.data = gene.simdata$counts
  }

  if(k>2) {
    gene.simdata = .simRNAseq.multi(simOptions = sim.settings, n = n, verbose=verbose)
    gene.data = gene.simdata$counts
  }

  # spike-in read counts
  if(isTRUE(spikeIns)) {
    if(is.null(sim.settings$spike)) {
      stop(message(paste0("You requested the simulations of spike-in read counts \n
                          but the spike parameter is empty. Aborting.")))
    }
    if(!is.null(sim.settings$spike)) {
      spikes <- sim.settings$spike$spikeInfo$SpikeID
      spikes.mat = matrix(0, nrow = length(spikes), ncol = sum(n),
                          dimnames = list(spikes,NULL))
    }
  }
  if(!isTRUE(sim.settings$spikeIns)) {
    spikes.mat = NULL
  }

  ## generate spike-in read counts
  if(isTRUE(spikeIns)) {
    if (verbose) { message(paste0("Generating spike-in read counts")) }
    spike.simdata = .simSpike(SpikeOptions = sim.settings$spike,
                              n1 = floor(sum(n)/2),
                              n2 = ceiling(sum(n)/2),
                              sf = gene.data$sf)
    spike.data = spike.simdata$counts
  }
  if(!isTRUE(spikeIns)) {
    spike.data = NULL
  }

  # return object
  res <- c(list(GeneCounts=gene.data,
                SpikeCounts=spike.data,
                DEid = DEids,
                Bid = Bids,
                Mid = Mids,
                pLFC = plfcs,
                bLFC = blfcs,
                mLFC = mlfcs,
                sf = gene.simdata$sf,
                ngenes = ngenes,
                p.DE = p.DE,
                p.B = p.B,
                p.M = p.M,
                sim.seed = sim.seed,
                k = k,
                n = n)
           )

  return(res)
}


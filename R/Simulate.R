
# simulateDE --------------------------------------------------------------

#' @name simulateDE
#' @aliases simulateDE
#' @title Simulate Differential Expression Pipeline
#' @description simulateDE is the main function to simulate differential expression for RNA-seq experiments.
#' The simulation parameters are specified with \code{\link{Setup}}.
#' The user now needs to specify the RNA-seq Analysis Pipeline including preprocessing, normalisation and differential testing method.
#' The return object contains DE test results from all simulations as well as descriptive statistics.
#' The error matrix calculations will be conducted with \code{\link{evaluateDE}}.\cr
#' @usage simulateDE(SetupRes,
#' Prefilter = NULL,
#' Imputation = NULL,
#' Normalisation = c("TMM", "MR", "PosCounts", "UQ",
#' "scran", "Linnorm", "sctransform",
#' "SCnorm", "Census", "depth"),
#' Label = "none",
#' DEmethod = c("T-Test", "edgeR-LRT", "edgeR-QL",
#' "edgeR-zingeR", "edgeR-ZINB-WaVE",
#' "limma-voom", "limma-trend",
#' "DESeq2", "DESeq2-zingeR", "DESeq2-ZINB-WaVE",
#' "ROTS", "baySeq", "NOISeq", "EBSeq",
#' "MAST", "BPSC", "scDD", "DECENT"),
#' DEFilter = FALSE,
#' Counts = FALSE,
#' NCores = NULL,
#' verbose = TRUE)
#' @param SetupRes This object specifies the simulation setup.
#' This must be the return object from \code{\link{Setup}}.
#' @param Prefilter A character vector specifying the gene expression filtering method
#' to be used prior to normalisation (and possibly imputation).
#' Default is \code{NULL}, i.e. no filtering.
#' Please consult the Details section for available options.
#' @param Imputation A character vector specifying the gene expression imputation method
#' to be used prior to normalisation.
#' Default is \code{NULL}, i.e. no imputation.
#' Please consult the Details section for available options.
#' @param Normalisation Normalisation method to use.
#' Please consult the Details section for available options.
#' @param Label A character vector to define whether information about
#' group labels should be used for normalisation.
#' This is only implemented for scran and SCnorm.
#' Possible options include the default \code{"none"}
#' which means that no sample group information is considered for normalisation;
#' \code{"known"} means that the simulated group labels are used and \code{"clustering"}
#' which applies an unsupervised hierarchical clustering to determine the group labels.
#' For details, see \code{\link[scran]{quickCluster}}).
#' @param DEmethod A character vector specifying the DE detection method to be used.
#' Please consult the Details section for available options.
#' @param DEFilter A logical vector indicating whether to run DE testing on
#' filtered and/or imputed count data.
#' Default is \code{FALSE}.
#' @param Counts A logical vector indicating whether the simulated count matrix is also provided as output.
#' Default is \code{FALSE} since the output can be quite large. Note that if DEFilter is \code{TRUE},
#' then the returned count matrix will countain the filtered and/or imputed count data.
#' @param NCores integer positive number of cores for parallel processing.
#' Default is \code{NULL}, i.e. 1 core.
#' @param verbose Logical vector to indicate whether to show progress report of simulations.
#' Default is \code{TRUE}.
#' @return
#' \strong{SimulateRes: Results of DE simulations}
#' \describe{
#' \item{pvalue, fdr}{3D array (ngenes * N * nsims) for p-values and FDR from each simulation.
#' Note that FDR values will be empty and the calculation will be done by \code{\link{evaluateDE}} whenever applicable.}
#' \item{mu,disp,dropout}{3D (ngenes * N * nsims) array for mean, dispersion and dropout of library size factor normalized read counts.}
#' \item{true.mu,true.disp,true.dropout}{3D (ngenes * N * nsims) array for true mean, dispersion and dropout of simulated read counts.}
#' \item{true.depth,est.depth}{True simulated and processed (after prefiltering and/or imputation if defined) sequencing depth per sample.}
#' \item{est.sf}{Global library size factor estimates per sample.}
#' \item{est.gsf}{3D array (ngenes * N * nsims) for size factor estimates. These are gene- and sample-wise estimates and only for SCnorm and Linnorm normalisation.}
#' \item{elfc,rlfc}{3D array (ngenes * N * nsims) for log fold changes (LFC):
#' elfc is for the DE tool estimated LFC; rlfc is for the LFC estimated from the normalised read counts.}
#' \item{time.taken}{The time taken given by \code{\link[base]{proc.time}} for each simulation, given for preprocessing, normalisation, differential expression testing and moment estimation.}
#' }
#' \strong{SetupRes: Simulation specifications}
#' \describe{
#' \item{DESetup - ... - estSpikeRes}{Reiterating the simulated setup defined by \code{\link{Setup}}.}
#' \item{Pipeline}{A list of chosen pipeline tools defined by above arguments.}
#' }
#' \strong{Counts: Simulated Count Matrices}
#' \describe{
#' \item{Counts}{3D array (ngenes * N * nsims) containing simulated counts. Note that this will only be returned when \code{Counts} is \code{TRUE}. In addition, if \code{DEFilter} is \code{TRUE} then the filtered/imputed counts are returned.}
#' }
#'
#' @seealso \code{\link{estimateParam}} and \code{\link{estimateSpike}},  for parameter specifications;\cr
#'  \code{\link{Setup}} for simulation setup;\cr
#'  \code{\link{evaluateDE}}, \code{\link{evaluateROC}} and \code{\link{evaluateSim}} for evaluation of simulation results.
#'
#' @details
#' Here you can find detailed information about preprocessing, imputation, normalisation and differential testing choices.
#' For recommendations concerning single cell RNA-sequencing pipelines,
#' we kindly refer the user to \href{https://www.nature.com/articles/s41467-019-12266-7}{Vieth, et al (2019). A systematic evaluation of single cell RNA-seq analysis pipelines. Nature Communications, 10(1), 4667}.\cr
#' \strong{Prefiltering}\cr
#' \describe{
#' \item{CountFilter}{removes genes that have a mean expression below 0.2.}
#' \item{FreqFilter}{removes genes that have more than 80 \% dropouts.}
#' }
#' \strong{Imputation}\cr
#' \describe{
#' \item{scImpute}{apply the imputation as implemented in \code{\link[scImpute]{scimpute}}}
#' \item{DrImpute}{apply the imputation as implemented in \code{\link[DrImpute]{DrImpute}}.}
#' \item{SAVER}{apply the imputation as implemented in \code{\link[SAVER]{saver}}.}
#' \item{scone}{apply the imputation as implemented in \code{\link[scone]{scone}}, defining 'house keeping genes' for the FNR model estimation as those that have less than 20 \% dropouts and small variance (i.e. in the lower 20th quartile). If less than 25 genes could be identified, the genes with less than 5 \% dropouts are used. If spike-in data is provided then these are used for the FNR model estimation.}
#' \item{MAGIC}{apply the imputation as implemented in \code{\link[Rmagic]{magic}}. Please note that the python package MAGIC needs to be installed to use this implementation.}
#' }
#' \strong{Normalisation}\cr
#' \describe{
#' \item{TMM, UQ}{employ the edgeR style normalization of weighted trimmed mean of M-values and upperquartile
#' as implemented in \code{\link[edgeR]{calcNormFactors}}, respectively.}
#' \item{MR, PosCounts}{employ the DESeq2 style normalization of median ratio method and a modified geometric mean method
#' as implemented in \code{\link[DESeq2]{estimateSizeFactors}}, respectively.}
#' \item{scran, SCnorm}{apply the deconvolution and quantile regression normalization methods developed for sparse RNA-seq data
#' as implemented in \code{\link[scran]{computeSumFactors}} and \code{\link[SCnorm]{SCnorm}}, respectively. Spike-ins can also be supplied for both methods via \code{spikeData}. Note, however that this means for scran that the normalisation as implemented in \code{\link[scran]{calculateSumFactors}} is also applied to genes (\code{general.use=TRUE}).}
#' \item{Linnorm}{apply the normalization method for sparse RNA-seq data
#' as implemented in \code{\link[Linnorm]{Linnorm.Norm}}.
#' For \code{Linnorm}, the user can also supply \code{spikeData}.}
#' \item{sctransform}{apply the normalization method developed for single-cell
#' UMI RNA-seq data as implemented in \code{\link[sctransform]{vst}}. }
#' \item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
#' Census at least needs \code{Lengths} for single-end data and preferably \code{MeanFragLengths} for paired-end data.
#' The authors state that Census should not be used for UMI data.}
#' \item{depth}{Sequencing depth normalisation.}
#' }
#' \strong{Differential testing}\cr
#' \describe{
#' \item{T-Test}{A T-Test per gene is applied using log transformed and normalized expression values (i.e. CPM or TPM).}
#' \item{limma-trend, limma-voom}{apply differential testing as implemented in \code{\link[limma]{lmFit}}
#' followed by \code{\link[limma]{eBayes}} on counts transformed by \code{\link[limma]{voom}} or by applying mean-variance trend on log CPM values in \code{\link[limma]{eBayes}}.}
#' \item{edgeR-LRT, edgeR-QL}{apply differential testing as implemented in \code{\link[edgeR]{glmFit}}, \code{\link[edgeR]{glmLRT}} and\code{\link[edgeR]{glmQLFit}}, \code{\link[edgeR]{glmQLFTest}}, respectively.}
#' \item{DESeq2}{applies differential testing as implemented in \code{\link[DESeq2]{DESeq}}.}
#' \item{ROTS}{applies differential testing as implemented in \code{\link[ROTS]{ROTS}} with 100 permutations on transformed counts (\code{\link[limma]{voom}}).}
#' \item{baySeq}{applies differential testing as implemented in \code{\link[baySeq]{getLikelihoods}} based on negative binomial prior estimates (\code{\link[baySeq]{getPriors.NB}}).}
#' \item{NOISeq}{applies differential testing as implemented in \code{\link[NOISeq]{noiseqbio}} based on CPM values.}
#' \item{EBSeq}{applies differential testing as implemented in \code{\link[EBSeq]{EBTest}}.}
#' \item{MAST}{applies differential testing as implemented in \code{\link[MAST]{zlm}} for zero-inflated model fitting followed by \code{\link[MAST]{lrTest}} on log CPM values.}
#' \item{BPSC}{applies differential testing as implemented in \code{\link[BPSC]{BPglm}} on CPM values.}
#' \item{scDD}{applies differential testing as implemented in \code{\link[scDD]{scDD}} on CPM values.}
#' \item{DECENT}{applies differential testing as implemented in \code{\link[DECENT]{decent}}.}
#' \item{edgeR-zingeR, DESeq2-zingeR}{In a first step, the posterior probabilities of the zero-inflated negative binomial component are estimated (see \code{\link[zingeR]{zeroWeightsLS}}) and used to define a weight matrix for dispersion estimation in \code{\link[edgeR]{estimateDisp}}. For the edgeR approach, the generalized model as implemented in \code{\link[edgeR]{glmFit}} is fitted. This is followed by an adapted LRT for differential testing to account for the weighting (see \code{\link[zingeR]{glmWeightedF}}). For DESeq2, the generalized linear model coefficients are estimated using \code{\link[DESeq2]{nbinomWaldTest}} and the weighting is done by setting the degrees of freedom for the T distribution.}
#' \item{edgeR-ZINB-WaVE, DESeq2-ZINB-WaVE}{In a first step, a zero-inflated negative binomial regression model  is fitted (see \code{\link[zinbwave]{zinbFit}}) to estimate observational weights (see \code{\link[zinbwave]{computeObservationalWeights}}) used for dispersion estimation in \code{\link[edgeR]{estimateDisp}}. For the edgeR approach, the generalized model as implemented in \code{\link[edgeR]{glmFit}} is fitted. This is followed by an adapted LRT for differential testing to account for the weighting (see \code{\link[zinbwave]{glmWeightedF}}). For DESeq2, the generalized linear model coefficients are estimated using \code{\link[DESeq2]{nbinomWaldTest}} and the weighting is done by setting the degrees of freedom for the T distribution.}
#' }
#' @examples
#' \dontrun{
#' # estimate gene parameters
#' data("SmartSeq2_Gene_Read_Counts")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SmartSeq2_Gene_Read_Counts))
#' data("GeneLengths_mm10")
#' estparam_gene <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
#'                                readData = NULL,
#'                               batchData = Batches,
#'                                spikeData = NULL, spikeInfo = NULL,
#'                                Lengths = GeneLengths_mm10, MeanFragLengths = NULL,
#'                                RNAseq = 'singlecell', Protocol = 'Read',
#'                                Distribution = 'ZINB', Normalisation = "scran",
#'                                GeneFilter = 0.1, SampleFilter = 3,
#'                                sigma = 1.96, NCores = NULL, verbose = TRUE)
#' # estimate spike parameters
#' data("SmartSeq2_SpikeIns_Read_Counts")
#' data("SmartSeq2_SpikeInfo")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_SpikeIns_Read_Counts), "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SmartSeq2_SpikeIns_Read_Counts))
#' estparam_spike <- estimateSpike(spikeData = SmartSeq2_SpikeIns_Read_Counts,
#'                                 spikeInfo = SmartSeq2_SpikeInfo,
#'                                 MeanFragLength = NULL,
#'                                 batchData = Batches,
#'                                 Normalisation = 'depth')
#' # define log fold change
#' p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#' # set up simulations
#' setupres <- Setup(ngenes = 10000, nsims = 10,
#'                   p.DE = 0.1, pLFC = p.lfc,
#'                   n1 = c(20,50,100), n2 = c(30,60,120),
#'                   Thinning = c(1,0.9,0.8), LibSize = 'given',
#'                   estParamRes = estparam_gene,
#'                   estSpikeRes = estparam_spike,
#'                   DropGenes = FALSE,
#'                   sim.seed = 52679, verbose = TRUE)
#' # run simulation
#' simres <- simulateDE(SetupRes = setupres,
#'                      Prefilter = "FreqFilter", Imputation = NULL,
#'                      Normalisation = 'scran', Label = 'none',
#'                      DEmethod = "limma-trend", DEFilter = FALSE,
#'                      NCores = NULL, verbose = TRUE)
#' # quick evaluation
#' evalderes <- evaluateDE(simRes = simres)
#' # plot evaluation
#' plotEvalDE(evalRes = evalderes, rate = "marginal",
#'            quick = TRUE, Annot = FALSE)
#' }
#' @author Beate Vieth
#' @rdname simulateDE
#' @importFrom stats setNames
#' @importFrom edgeR thinCounts
#' @export
simulateDE <- function(SetupRes,
                       Prefilter = NULL,
                       Imputation = NULL,
                       Normalisation = c("TMM", "MR", "PosCounts", "UQ",
                                         "scran", "Linnorm", "sctransform",
                                         "SCnorm", "Census", "depth"),
                       Label = "none",
                       DEmethod = c("T-Test", "edgeR-LRT", "edgeR-QL",
                                    "edgeR-zingeR", "edgeR-ZINB-WaVE",
                                    "limma-voom", "limma-trend",
                                    "DESeq2", "DESeq2-zingeR", "DESeq2-ZINB-WaVE",
                                    "ROTS", "baySeq", "NOISeq", "EBSeq",
                                    "MAST", "BPSC", "scDD", "DECENT"),
                       DEFilter = FALSE,
                       Counts = FALSE,
                       NCores = NULL,
                       verbose = TRUE) {

  if(!is.null(NCores) && DEmethod %in% c("edgeR-LRT", "edgeR-QL", "edgeR-zingeR", "DESeq2-zingeR", 'limma-voom', "limma-trend", "NOISeq", "EBSeq", "ROTS")) {
    if(verbose) {message(paste0(DEmethod, " has no parallel computation option!"))}
  }

  if(attr(SetupRes, 'RNAseq') == "singlecell" &&
     DEmethod %in% c("edgeR-LRT", "edgeR-QL", "limma-voom", "limma-trend", "DESeq2", "baySeq", "NOISeq", "EBSeq")) {
    if(verbose) {message(paste0(DEmethod, " is developed for bulk RNA-seq experiments."))}
  }

  if(attr(SetupRes, 'RNAseq') == "bulk" && DEmethod %in% c("MAST", 'BPSC', "edgeR-zingeR", "DESeq2-zingeR", "edgeR-ZINB-WaVE", "DESeq2-ZINB-WaVE")) {
    if(verbose) {message(paste0(DEmethod, " is developed for single cell RNA-seq experiments."))}
  }

  if(attr(SetupRes, 'RNAseq') == "bulk" && DEmethod %in% c('scDD', 'DECENT')) {
    stop(message(paste0(DEmethod, " is only developed and implemented for single cell RNA-seq experiments.")))
  }

  if(DEmethod == "DECENT") {
    if(verbose) {message(paste0(DEmethod, " does not require additional normalisation nor imputation."))}
    Normalisation = "none"
    Prefilter = NULL
    Imputation = NULL
  }

  if(c(all(is.null(Imputation), is.null(Prefilter)) && isTRUE(DEFilter))) {
    stop(message(paste0("You wish to use imputed/filtered gene expression values for DE testing but you did not specify the imputation/filtering method. Aborting.")))
  }

  if(!is.null(Imputation) && attr(SetupRes,"RNAseq")== "bulk") {
    message(paste0("You wish to use imputation but powsimR has only methods implemented for single cell RNA-seq and in most cases imputation is not needed for bulk RNA-seq. Setting Imputation to NULL."))
    Imputation = NULL
  }

  start.time.pipe = proc.time()

  # define the maximal count matrix for simulations
  max.n = max(SetupRes$SimSetup$n1, SetupRes$SimSetup$n2)
  min.n = min(SetupRes$SimSetup$n1, SetupRes$SimSetup$n2)

  # append additional settings of simulateDE to SetupRes
  if(Label == "clustering") {PreclustNumber <- min.n}
  if(!Label == "clustering") {PreclustNumber <- NULL}
  Pipeline <- list(Prefilter = Prefilter,
                   Imputation = Imputation,
                   Normalisation=Normalisation,
                   Label = Label,
                   DEmethod=DEmethod,
                   DEFilter = DEFilter,
                   NCores = NCores,
                   clustNumber = ifelse(SetupRes$DESetup$design=="2grp", 2, NULL),
                   PreclustNumber = PreclustNumber)

  SetupRes <- c(SetupRes, Pipeline=list(Pipeline))

  if (verbose) { message(paste0("Preparing output arrays.")) }

  my.names = paste0(SetupRes$SimSetup$n1,"vs",SetupRes$SimSetup$n2)

  #set up output arrays
  pvalues = fdrs = elfcs = rlfcs = mus = true.mus = disps = true.disps = dropouts = true.drops = array(NA,dim=c(SetupRes$DESetup$ngenes,
                                                                                                                length(SetupRes$SimSetup$n1),
                                                                                                                SetupRes$DESetup$nsims))

  tstep <- c("Simulation",
             "Preprocessing",
            "Normalisation",
            "DE",
            "Moments",
            "Total")
  tmoment <- c('User', 'System', 'Elapsed')
  tname <- paste(rep(tstep, each=3), rep(tmoment, times=6), sep="_")

  true.sf <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol =SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x])
  })
  est.sf <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol =SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x])
  })
  true.depth <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol =SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x])
  })
  est.depth <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol =SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x])
  })
  time.taken <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol = length(tname),
           dimnames = list(NULL, tname))
  })

  true.designs <- lapply(1:length(my.names), function(x) {
    matrix(NA,
           nrow = SetupRes$DESetup$nsims,
           ncol = SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x])
  })

  names(true.sf) = names(est.sf) = names(true.depth) = names(est.depth) = names(time.taken) = names(true.designs) = my.names

  if(SetupRes$Pipeline$Normalisation %in% c("SCnorm", "Linnorm", 'sctransform', 'bayNorm')) {
    est.gsf <- lapply(1:length(my.names), function(x) {
      array(NA,dim=c(SetupRes$DESetup$nsims,
                     SetupRes$DESetup$ngenes,
                     SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x]))
    })
    names(est.gsf) = my.names
  }
  if(!SetupRes$Pipeline$Normalisation %in% c("SCnorm", "Linnorm", 'sctransform', 'bayNorm')) {
    est.gsf = NULL
  }

  if(isTRUE(Counts)){
    cnts <- lapply(1:length(my.names), function(x) {
      array(0,dim=c(SetupRes$DESetup$nsims,
                     SetupRes$DESetup$ngenes,
                     SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x]))
    })
    names(cnts) = my.names
  }

  if(isFALSE(Counts)){
    cnts = NULL
  }

  ## start simulation
  for (i in 1:SetupRes$DESetup$nsims) {
    if (verbose) { message(paste0("\n  SIMULATION   NUMBER   ", i, "\n")) }
    start.time.sim1 <- proc.time()
    ## update the simulation options by extracting the ith set and change sim.seed
    tmp.simOpts = SetupRes
    tmp.simOpts$DESetup$DEid = SetupRes$DESetup$DEid[[i]]
    tmp.simOpts$DESetup$pLFC = SetupRes$DESetup$pLFC[[i]]
    tmp.simOpts$DESetup$Bid = SetupRes$DESetup$Bid[[i]]
    tmp.simOpts$DESetup$bLFC = SetupRes$DESetup$bLFC[[i]]
    tmp.simOpts$DESetup$sim.seed = SetupRes$DESetup$sim.seed[[i]]

    ## generate gene read counts
    if (verbose) {message(paste0("Generating gene expression.")) }
    gene.data = .simRNAseq.2grp(simOptions = tmp.simOpts,
                                n1 = max.n, n2 = max.n, verbose=verbose)

    ## apply gene dropouts
    if(isTRUE(tmp.simOpts$SimSetup$DropGenes)){
      gene.data = .dropGene(simOptions = tmp.simOpts, simData = gene.data)
    }

    ## generate spike-in read counts
    if(isTRUE(tmp.simOpts$SimSetup$spikeIns)) {
      if (verbose) { message(paste0("Generating spike-in expression.")) }
      spike.data = .simSpike(SpikeOptions = tmp.simOpts$estSpikeRes,
                             n1 = max.n, n2 = max.n, sf = gene.data$sf)
      spike.info = tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
    }
    if(!isTRUE(tmp.simOpts$SimSetup$spikeIns)) {
      spike.data = NULL
      spike.info = NULL
    }

    ## generate mean fragment lengths for samples
    if(!is.null(tmp.simOpts$estParamRes$MeanFragLengths)) {
      if (verbose) { message(paste0("Sampling from observed mean fragment lengths")) }
      MeanFrag.data = sample(tmp.simOpts$estParamRes$MeanFragLengths,
                             max.n+max.n, replace = TRUE)
      names(MeanFrag.data) = colnames(gene.data$counts)
    }
    if(is.null(tmp.simOpts$estParamRes$MeanFragLengths)) {
      MeanFrag.data = NULL
    }

    ## match sampled gene names with given gene lengths
    if(!is.null(tmp.simOpts$estParamRes$Lengths)) {
      gene.id = sub('_([^_]*)$', '', rownames(gene.data$counts))
      Length.data = tmp.simOpts$estParamRes$Lengths
      Length.data = Length.data[match(gene.id,names(Length.data))]
    }
    if(is.null(tmp.simOpts$estParamRes$Lengths)) {
      Length.data = NULL
    }
    end.time.sim1 <- proc.time()

    ##  for different sample sizes
    for (j in seq(along=tmp.simOpts$SimSetup$n1)) {
      start.time.sim2 <- proc.time()
      Nrep1 = tmp.simOpts$SimSetup$n1[j]
      Nrep2 = tmp.simOpts$SimSetup$n2[j]
      Thin = tmp.simOpts$SimSetup$Thinning[j]
      if (verbose) { message(paste0(Nrep1, " vs. ", Nrep2)) }

      ## take a subsample of simulated samples
      idx = c(1:Nrep1, max.n + (1:Nrep2))
      true.design = gene.data$designs[idx]

      ## take a subsample of the simulated read counts
      sim.cnts = gene.data$counts[,idx]
      ## take a subsample of the true size factors
      gene.sf = gene.data$sf[idx]

      ## apply thinning
      if(!is.null(Thin) && !Thin == 1) {
        if (verbose) { message(paste0("Reduce the size of gene counts to ",
                                      Thin, " using binomial thinning.")) }
        sim.cnts = .run.thin(countData = sim.cnts,
                               Thin = Thin,
                               simOptions = tmp.simOpts)
      }

      ## filter out zero expression genes
      ix.valid = rowSums(sim.cnts) > 0
      count.data = sim.cnts[ix.valid,, drop = FALSE]
      ## record simulated seq depth
      sim.depth = colSums(count.data)

      ## match sampled gene names with given gene lengths
      if(!is.null(Length.data)) {
        if (verbose) { message(paste0("Associating gene lengths with gene expression")) }
        gene.id = sub('_([^_]*)$', '', rownames(count.data))
        length.data = Length.data
        length.data = length.data[match(gene.id,names(length.data))]
      }
      if(is.null(Length.data)) {
        length.data = NULL
      }

      ## take a subsample of simulated spike-ins
      if(!is.null(spike.data) && !is.null(spike.info)) {
        # counts
        sim.spike <- spike.data$counts
        spike.valid = rowSums(sim.spike) > 0
        count.spike = sim.spike[spike.valid, idx, drop=FALSE]
        # spike info table
        info.spike <- tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
        info.spike <- info.spike[rownames(info.spike)  %in% rownames(count.spike), , drop = FALSE]
        info.spike <- info.spike[match(rownames(count.spike), rownames(info.spike)), , drop = FALSE]
      }

      ## apply thinning to spike-ins
      if(c(!is.null(spike.data) && !is.null(spike.info) && !is.null(Thin) && !Thin == 1 && isTRUE(tmp.simOpts$SimSetup$thinSpike))) {
        if (verbose) { message(paste0("Reduce the size of spike-in counts to ",
                                      Thin, " using binomial thinning.")) }
        count.spike = .run.thin(countData = count.spike,
                                Thin = Thin,
                                simOptions = tmp.simOpts)
        count.spike = count.spike[rowSums(count.spike)>0,]
        # spike info table
        info.spike <- tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
        info.spike <- info.spike[rownames(info.spike)  %in% rownames(count.spike), , drop = FALSE]
        info.spike <- info.spike[match(rownames(count.spike), rownames(info.spike)), , drop = FALSE]

      }

      if(is.null(spike.data) && is.null(spike.info)) {
        count.spike = NULL
        info.spike = NULL
      }

      ## take a subsample of mean fragment lengths
      if(!is.null(MeanFrag.data)) {
        meanfrag.data = MeanFrag.data[idx]
      }
      if(is.null(MeanFrag.data)) {
        meanfrag.data = NULL
      }

      def.design <- true.design

      end.time.sim2 <- proc.time()

      ## perform filtering / imputation (OPTIONAL)
      start.time.preprocess <- proc.time()
      if(!is.null(tmp.simOpts$Pipeline$Prefilter)) {
        if (verbose) { message(paste0("Applying ",tmp.simOpts$Pipeline$Prefilter," prefiltering")) }
        filter.data <- .prefilter.calc(Prefilter=tmp.simOpts$Pipeline$Prefilter,
                                       countData=count.data,
                                       NCores=tmp.simOpts$Pipeline$NCores)
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
      if(is.null(tmp.simOpts$Pipeline$Prefilter)) {
        filter.count.data <- count.data
      }
      if(!is.null(tmp.simOpts$Pipeline$Imputation)) {
        if (verbose) { message(paste0("Applying ", tmp.simOpts$Pipeline$Imputation, " imputation")) }
        impute.data <- .impute.calc(Imputation=tmp.simOpts$Pipeline$Imputation,
                                    countData=filter.count.data,
                                    spikeData=count.spike,
                                    batchData=def.design,
                                    clustNumber=tmp.simOpts$Pipeline$clustNumber,
                                    Lengths = length.data,
                                    MeanFragLengths = meanfrag.data,
                                    NCores=tmp.simOpts$Pipeline$NCores,
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
      if(is.null(Imputation)) {
        fornorm.count.data <- filter.count.data
      }
      end.time.preprocess <- proc.time()

      ## perform normalisation
      start.time.norm <- proc.time()
      if (verbose) { message(paste0("Applying ", tmp.simOpts$Pipeline$Normalisation, " normalisation")) }

      if (tmp.simOpts$Pipeline$Normalisation == "sctransform") {
        def.design <- NULL
      }
      norm.data <- .norm.calc(Normalisation=tmp.simOpts$Pipeline$Normalisation,
                              sf=gene.sf,
                              countData=fornorm.count.data,
                              spikeData=count.spike,
                              spikeInfo=info.spike,
                              batchData=def.design,
                              Lengths=length.data,
                              MeanFragLengths=meanfrag.data,
                              PreclustNumber=tmp.simOpts$Pipeline$PreclustNumber,
                              Step="Simulation",
                              Protocol=attr(tmp.simOpts$estParamRes, 'Protocol'),
                              Label=tmp.simOpts$Pipeline$Label,
                              NCores=tmp.simOpts$Pipeline$NCores,
                              verbose=verbose)
      end.time.norm <- proc.time()

      ## create an DE options object to pass into DE detection
      DEOpts <- list(designs=def.design, p.DE=tmp.simOpts$DESetup$p.DE)

      ## rematch sampled gene names with given gene lengths
      if(!is.null(Length.data)) {
        if (verbose) { message(paste0("Reassociate gene lengths with gene expression")) }
        gene.id = sub('_([^_]*)$', '', rownames(count.data))
        length.data = Length.data
        length.data = length.data[match(gene.id,names(length.data))]
      }
      if(is.null(Length.data)) {
        length.data = NULL
      }

      ## Run DE detection
      start.time.DE <- proc.time()
      if(!isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
        if (verbose) { message(paste0("Applying ", tmp.simOpts$Pipeline$DEmethod,
                                      " for DE analysis on raw count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$Pipeline$DEmethod,
                          normData=norm.data,
                          countData=count.data,
                          Lengths=length.data,
                          MeanFragLengths=meanfrag.data,
                          DEOpts=DEOpts,
                          spikeData=count.spike,
                          spikeInfo=info.spike,
                          NCores=tmp.simOpts$Pipeline$NCores,
                          verbose=verbose)
      }
      if(isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
        if (verbose) { message(paste0("Applying ", tmp.simOpts$Pipeline$DEmethod,
                                      " for DE analysis on imputed/filtered count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$Pipeline$DEmethod,
                          normData=norm.data,
                          countData=fornorm.count.data,
                          Lengths=length.data,
                          MeanFragLengths=meanfrag.data,
                          DEOpts=DEOpts,
                          spikeData=count.spike,
                          spikeInfo=info.spike,
                          NCores=tmp.simOpts$Pipeline$NCores,
                          verbose=verbose)
      }
      end.time.DE <- proc.time()

      if(tmp.simOpts$Pipeline$DEmethod =="DECENT") {
        res.de <- res.de[["DEresults"]]
        norm.data <- res.de[["NormData"]]
      }

      # output count matrix
      if(isTRUE(Counts)){
        allgenes <- rownames(sim.cnts)
        if(isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
          if (verbose) { message(paste0("Saving simulated counts after imputation / filtering.")) }
          consideredgenes <- rownames(fornorm.count.data)
          ixx.valid <- allgenes %in% consideredgenes
          cnts[[j]][i, ixx.valid, ] <- fornorm.count.data

        }
        if(!isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
          if (verbose) { message(paste0("Saving raw simulated counts.")) }
          consideredgenes <- rownames(count.data)
          ixx.valid <- allgenes %in% consideredgenes
          cnts[[j]][i, ixx.valid, ] <- count.data
        }
        dimnames(cnts[[j]][i, , ]) <- list(allgenes, colnames(sim.cnts))
        print(cnts[[j]][i, 1:5, 1:5])
      }

      # adapt scale factor matrices to use in param estimation
      if(attr(norm.data, 'normFramework') %in% c("SCnorm", "Linnorm", "scTransform")) {
        allgenes <- rownames(sim.cnts)
        estgenes <- rownames(norm.data$scale.factors)
        ixx.valid <- allgenes %in% estgenes
        est.gsf[[j]][i, ixx.valid, ] = norm.data$scale.factors
        norm.data$scale.factors = est.gsf[[j]][i,,]
      }

      ## estimate moments of read counts after normalisation and preprocessing
      start.time.moments <- proc.time()
      if(isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
        if (verbose) { message(paste0("Estimating moments of imputed/filtered count data.")) }
        res.params <- .run.params(countData=fornorm.count.data,
                                  normData=norm.data,
                                  group=DEOpts$designs)
        est.depths <- colSums(fornorm.count.data)
      }
      if(!isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
        if (verbose) { message(paste0("Estimating moments of raw count data.")) }
        res.params <- .run.params(countData=count.data,
                                  normData=norm.data,
                                  group=DEOpts$designs)
        est.depths <- colSums(count.data)
      }
      end.time.moments <- proc.time()

      # generate empty vectors
      pval = fdr = est.lfc = raw.lfc = mu.tmp = true.mu.tmp = disp.tmp = true.disp.tmp = p0.tmp = true.p0.tmp = rep(NA, nrow(sim.cnts))
      # simulated genes
      allgenes <- rownames(sim.cnts)
      # indicator of tested genes
      testedgenes <- res.de$geneIndex
      ixx.de.valid <- allgenes %in% testedgenes
      # extract results of DE testing
      pval[ixx.de.valid] = res.de$pval
      fdr[ixx.de.valid] = res.de$fdr
      est.lfc[ixx.de.valid] = res.de$lfc
      # indicator of estimated genes
      paramgenes <- res.params$geneIndex
      ixx.param.valid <- allgenes %in% paramgenes
      # extract parameter estimates
      raw.lfc[ixx.param.valid] = res.params$lfc
      mu.tmp[ixx.param.valid] = res.params$means
      disp.tmp[ixx.param.valid] = res.params$dispersion
      p0.tmp[ixx.param.valid] = res.params$dropout
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
      true.depth[[j]][i,] = sim.depth
      est.depth[[j]][i,] = est.depths

      true.designs[[j]][i,] = true.design

      end.time.pipe = proc.time()

      # time taken for each step
      # copy designs into list of matrices
      time.taken.sim <- (end.time.sim1 - start.time.sim1) + (end.time.sim2 - start.time.sim2)
      if(any(c(!is.null(tmp.simOpts$Pipeline$Prefilter), !is.null(tmp.simOpts$Pipeline$Imputation)))) {
        time.taken.preprocess <- end.time.preprocess - start.time.preprocess
      }
      if(all(c(is.null(tmp.simOpts$Pipeline$Prefilter), is.null(tmp.simOpts$Pipeline$Imputation)))) {
        time.taken.preprocess = c(NA, NA, NA)
      }
      time.taken.norm <- end.time.norm - start.time.norm
      time.taken.DE <- end.time.DE - start.time.DE
      time.taken.moments <- end.time.moments - start.time.moments
      time.taken.total <- end.time.pipe - start.time.pipe

      timing <- c(time.taken.sim,
                  time.taken.preprocess,
                  time.taken.norm,
                  time.taken.DE,
                  time.taken.moments,
                  time.taken.total)

      timing <- timing[!grepl(pattern = "child", names(timing))]
      # copy time taken in 2D array of time taken
      time.taken[[j]][i, ] = timing

    }
  }

  ## return
  Simulate <- list(pvalue = pvalues,
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
                   est.depth = est.depth,
                   est.sf = est.sf,
                   est.gsf = est.gsf,
                   true.designs = true.designs,
                   time.taken = time.taken)
  attr(Simulate, 'Simulation') <- "DE"

  res.out <- c(SetupRes,
               SimulateRes=list(Simulate),
               Counts = list(cnts))

  return(res.out)
}



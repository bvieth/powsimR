###############################################################
## WRAPPER PowSim
###############################################################
#' @name PowSim
#' @aliases PowSim
#' @title Power Simulations for RNA-seq experiments
#' @description This function is a wrapper function including estimation, simulation and evaluation necessary for power analysis.
#' @usage PowSim(input=NULL,
#' RNAseq=c('bulk', 'singlecell'),
#' ngenes=10000, nsims=25,
#' p.DE=0.1, LFC=NULL,
#' size.factors='equal',
#' n1=NULL, n2=NULL,
#' ncores=NULL, DEmethod=NULL,
#' save.plots=TRUE, verbose=TRUE)
#' @param input This can be one of the following: 1) NULL (Default): The parameters are defined by the function \code{\link{insilicoNBParam}} (see details). 2) A matrix of read counts with rows=genes and columns=samples. In this case, the user should provide the counts of one group only, as the 'MatchMoments' estimation framework is used.  3) If a character string is defined, then the precalculated negative binomial parameters are used (see details).
#' @param RNAseq is a character value: "bulk" or "singlecell".
#' @param ngenes Number of genes to be simulated. Default is 10000.
#' @param nsims Number of simulations to run. Default is 25.
#' @param p.DE Percentage of genes being differentially expressed (DE). Default is 10\%.
#' @param LFC The log2 fold change for DE genes. This can be:
#' (1) a constant, e.g. 2;
#' (2) a vector of values with length being number of DE genes. If the input is a vector and the length is not the number of DE genes, it will be sampled with replacement to generate log-fold change;
#' (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' The default is NULL, i.e. the log2 fold changes are defined by a gamma distribution: function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, 4, 2).
#' @param size.factors Size factors representing sample-specific differences/biases in expected mean values of the NB distribution: "equal" or "given". The default is "equal", i.e. equal size factor of 1. If the user defines it as given, the size factors are sampled from the size factors provided by the output of \code{\link{estimateNBParam}} for count matrix input, precalculated estimates or defined as a sampling function (see details).
#' @param n1,n2 Integer vectors specifying the number of biological replicates in each group. Default values are n1=n2=c(2,4,6,10,15,20) and n1=n2=c(24,48,96,192,384,800) for single cell and bulk RNAseq experiments, respectively.
#' @param ncores integer positive number of cores for parallel processing. Default is NULL, i.e. 1 core.
#' @param DEmethod  String to specify the DE detection method to be used. Available options are: edgeRglm, edgeRQL, DESeq2, limma, ROTS, baySeq, NOISeq, DSS, MAST, scde, BPSC, scDD. Default is limma for bulk and MAST for single cell RNA-seq experiments.
#' @param save.plots Logical vector indicating whether plots should be saved as pdf files in the current working directory. Default is TRUE.
#' @param verbose Logical vector indicating whether progress messages and additional notifications should be printed. Default is TRUE.
#' @return
#' \item{Parameters}{The parameters of the negative binomial. This will be estimated if an input is provided, a precalculated estimates (see details) or in silico parameters. For details, see the result values of \code{\link{estimateNBParam}} and \code{\link{insilicoNBParam}}, respectively.}
#' \item{SimulationResults}{The results of differential expression simulation. For details, see the result values of \code{\link{simulateDE}}.}
#' \item{EvaluationResults}{The results of error matrices evaluation. For details, see the result values of \code{\link{evaluateSim}}.}
#' \item{SummaryTable}{The marginal TPR and FDR per sample sizes. For details, see \code{\link{printEvalRes}}.}
#' \item{MarginalPlot, ConditionalPLot}{The marginal and conditional TPR and FDR plots. See \code{\link{plotEvalRes}}.}
#' @details This function is a wrapper with a number of default settings. The insilico parameter definition depends on the RNAseq experiment.\cr
#' For single cells, the following is defined:
#' \itemize{
#' \item{Mean gene expression}{\code{function(x) rgamma(x, 4, 2)}.}
#' \item{Gene-wise dispersion}{\code{function(x) 2 + 100/x} where x is the average expression level.}
#' \item{Size factors}{\code{function(x) 2^rnorm(n=x, mean=0, sd=0.25)}}
#' }
#' For bulk, the following is defined:
#' \itemize{
#' \item{Mean gene expression}{\code{function(x) 2^rnorm(x, mean=8, sd=2).}
#' \item{Gene-wise dispersion}{\code{function(x) rgamma(x, shape = 2, rate = 6).}
#' \item{Gene-wise dropout rate}{\code{function(x) runif(x, min=0, max=0.25)}}
#' \item(Size factors){\code{function(x) rnorm(n=x, mean=1, sd=0.1)}}
#' }
#' For both types of experiment, we assume no differences in relative sequencing depth across samples.\\
#' We have included precalculated negative binomial parameters for simulations. There are parameters estimated from several real datasets distributed with the package. Available string options are: (1) "buettner". (2) "islam2011". (3) "islam2014". (4) "kolodziejczk". (5) "soumillon".\\
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}}, \code{\link{insilicoNBParam}}, \code{\link{DESetup}}, \code{\link{SimSetup}}, \code{\link{simulateDE}}, \code{\link{evaluateSim}}
#' @examples
#' \dontrun{
#' ## simulating in silico single cell RNAseq data
#' ## DE testing with  MAST
#' ## NOTE: If parallel computing is not possible, consider using DEmethod="ROTS"!
#' powsim <- PowSim(input=NULL,
#' RNAseq='singlecell',
#' ngenes=10000,
#' nsims=25,
#' p.DE=0.1,
#' LFC=function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, 3, 3),
#' size.factors='equal',
#' ncores=10,
#' DEmethod="MAST",
#' save.plots=TRUE, verbose=TRUE)
#' }
#' @rdname PowSim
#' @importFrom cowplot ggsave
#' @export
PowSim <- function(input=NULL, RNAseq=c('bulk', 'singlecell'), ngenes=10000, nsims=25, p.DE=0.1, LFC=NULL, size.factors='equal', n1=NULL, n2=NULL, ncores=NULL, DEmethod=NULL, save.plots=TRUE, verbose=TRUE) {

  # ESTIMATION
  # insilico
  if(is.null(input)) {
    if(verbose){message("No input is provided. Defining in silico negative binomial parameters!")}
    if(RNAseq=='bulk') {
      if(verbose){message("Assuming a bulk RNA-seq experiment.")}
      mu <- function(x) 2^rnorm(x, mean=8, sd=2)
      disp <- function(x) rgamma(x, shape = 2, rate = 6)
      dropout <- function(x) runif(x, 0, 0.25)
      sf <- function(x) rnorm(n=x, mean=1, sd=0.1)
      param <- suppressMessages(insilicoNBParam(means = mu, dispersion = disp, dropout = dropout, sf = sf, RNAseq = 'bulk'))
    }
    if(RNAseq=='singlecell') {
      if(verbose){message("Assuming a single cell RNA-seq experiment.")}
      mu <- function(x) rgamma(x, 4, 2)
      disp <- function(x) 2 + 100/x
      sf <- function(x) 2^rnorm(n=x, mean=0, sd=0.25)
      param <- suppressMessages(insilicoNBParam(means = mu, dispersion = disp, dropout = NULL, sf = sf, RNAseq = 'singlecell'))
    }
  }
  # estimate params
  if(is.matrix(input)) {
    if(verbose){message("Input count matrix provided. Estimating negative binomial parameters!")}
    if(RNAseq=='bulk') {
      param <- suppressMessages(estimateNBParam(countData = input, cData=NULL, design=NULL, RNAseq = 'bulk', estFramework = 'edgeR'))
      param.plot <- plotNBParam(estParam.out = param, annot=T)
    }
    if(RNAseq=='singlecell') {
      param <- supressMessages(estimateNBParam(countData = input, cData=NULL, design=NULL, RNAseq = 'singlecell', estFramework = 'MatchMoments'))
      param.plot <- plotNBParam(estParam.out = param, annot=T)
    }
  }
  # precalculated data sets
  if (is.list(input)) { ## a result list object by estimateNBParam()
    param <- input
    param.plot <- plotNBParam(estParam.out = param, annot=T)
  }

  # SIMULATION

  # DEsetup
  if(is.null(LFC)){
    if(verbose){message("No log2 fold change provided. Assuming gamma distributed fold changes!")}
    lfc.foo <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, 4, 2)
    desetup <- DESetup(ngenes=ngenes, nsims=nsims, p.DE=p.DE, LFC=lfc.foo)
  }
  if(!is.null(LFC)){
    desetup <- DESetup(ngenes=ngenes, nsims=nsims, p.DE=p.DE, LFC=LFC)
  }

  # Simsetup
  simsetup <- SimSetup(desetup=desetup, params=param, size.factors = size.factors)
  if(is.null(n1) || is.null(n2)) {
    if(verbose){message("No sample sizes provided!")}
    if(RNAseq=='bulk') {
      if(verbose){message("Assuming bulk RNA-seq with small number of replicates per group.")}
      n1 <- n2 <- c(2,4,6,10,15,20)
    }
    if(RNAseq=='singlecell') {
      if(verbose){message("Assuming single cell RNA-seq with large number of replicates per group.")}
      n1 <- n2 <- c(24,48,96,192,384,800)
    }
  }

  # simulateDE
  if(is.null(DEmethod)) {
    if(verbose){message("No DE method defined.")}
    if(RNAseq=='bulk') {
      if(verbose){message("Using limma for bulk RNA-seq experiment.")}
      simres <- simulateDE(n1=n1, n2=n2, sim.settings = simsetup, ncores=NULL, DEmethod='limma', verbose=F)
    }
    if(RNAseq=='singlecell') {
      if(verbose){
        message("Using MAST for single cell RNA-seq experiment")
        if(ncores==1){message('Consider using more cores for parallel computing to speed up calculations.')}
        }
      simres <- simulateDE(n1=n1, n2=n2, sim.settings = simsetup, ncores=ncores, DEmethod='MAST', verbose=F)
    }
  }
  if(!is.null(DEmethod)) {
    simres <- simulateDE(n1=n1, n2=n2, sim.settings = simsetup, ncores=ncores, DEmethod=DEmethod, verbose=F)
  }

  # Insilico PARAMS
  if (attr(simsetup, 'param.type') == 'insilico') {
      max.n <- max(n1,n2)
      tmp.simOpts <- simsetup
      tmp.simOpts$DEid <- tmp.simOpts$DEid[[1]]
      tmp.simOpts$DEid <- 0
      tmp.simOpts$lfc <- tmp.simOpts$lfc[[1]]
      tmp.simOpts$lfc <- rep(0, ngenes)
      tmp.simOpts$sim.seed <- tmp.simOpts$sim.seed[[1]]
      dat.sim <- .simRNAseq(simOptions = tmp.simOpts, n1 = max.n, n2 = max.n)
      estparam <- estimateNBParam(countData = dat.sim$counts, cData = NULL, design=NULL, RNAseq=RNAseq, estFramework = 'MatchMoments', sigma = 1.96)
      simres$sim.settings$means <- estparam$means
      simres$sim.settings$dispersion <- estparam$dispersion
      simres$sim.settings$p0 <- estparam$p0
  }

  # EVALUATION
  evalres <- evaluateSim(simRes=simres, alpha.type="adjusted", MTC="BH", alpha.nominal=0.1,
              stratify.by="mean",
              filter.by="none", strata.filtered=1,
              target.by="lfc", delta=0)

  # print marginal table
  if(verbose){message("Printing the marginal error rates in a table.")}
  evaltable <- printEvalRes(evalRes=evalres)

  # plot marginal and conditional error matrices
  marg.plot <- plotEvalRes(evalRes = evalres, rate='marginal', quick=T, annot=T)
  cond.plot <- plotEvalRes(evalRes = evalres, rate='stratified', quick=T, annot=T)

  if(save.plots) {
    if(verbose){message(cat("Saving marginal and conditional error rate plots in ", getwd(), "."))}
    # parameter plot
    if(!is.null(input)){
      if(verbose){message(cat("Saving the plot of estimated negative binomial parameters in ", getwd(), "."))}
    cowplot::ggsave(plot = param.plot, filename = 'EstimatedParameters.pdf', width = 250, height=200, units = 'mm')
    }
    # marginal plot
    if(verbose){message(cat("Saving the marginal and conditional error matrix plots in ", getwd(), "."))}
    cowplot::ggsave(plot = marg.plot, filename = 'MarginalResults.pdf', width = 250, height=200, units = 'mm')
    # conditional plot
    cowplot::ggsave(plot = cond.plot, filename = 'ConditionalResults.pdf', width = 250, height=200, units = 'mm')
  }

  # RESULTS OBJECT contains param, simres, evalres, evaltable, plots
  output <- (list('Parameters'=param,
                  "SimulationResults"=simres,
                  'EvaluationResults'=evalres,
                  "SummaryTable"=evaltable,
                  'MarginalPlot'=marg.plot,
                  'ConditionalPlot'=cond.plot))

  return(output)

}

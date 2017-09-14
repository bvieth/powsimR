
# DESetup -----------------------------------------------------------------

#' @name DESetup
#' @aliases DESetup
#' @title Setup options for RNA-seq count simulations.
#' @description This function generates a set of differential expressed gene IDs with associated fold changes for a given number of genes, simulations and fraction of DE genes.
#' @usage DESetup(ngenes, nsims=25,
#' p.DE=0.1, p.B=NULL,
#' pLFC, bLFC=NULL, sim.seed)
#' @param ngenes Number of genes to be simulated.
#' @param nsims Number of simulations to run. Default is 25.
#' @param p.DE Percentage of genes being differentially expressed (DE). Default is 10\%.
#' @param p.B Percentage of genes being differentially expressed between batches. Default is \code{NULL}.
#' @param pLFC The log2 phenotypic fold change for DE genes. This can be:
#' (1) a constant, e.g. 2;
#' (2) a vector of values with length being number of DE genes. If the input is a vector and the length is not the number of DE genes, it will be sampled with replacement to generate log-fold change;
#' (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' @param bLFC The log2 batch fold change for all genes. This can be:
#' (1) a constant, e.g. 2;
#' (2) a vector of values with length being number of all genes. If the input is a vector and the length is not the number of total genes, it will be sampled with replacement to generate log2 fold changes;
#' (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5).
#' @param sim.seed Simulation seed.
#' @return A list with the following entries:
#' \item{ngenes}{An integer for number of genes.}
#' \item{nsims}{An integer for number of simulations.}
#' \item{sim.seed}{The specified simulation seed.}
#' \item{p.DE}{Percentage of DE genes.}
#' \item{DEid}{A list (length=nsims) of vectors (length=ngenes*p.DE) for the IDs of DE genes.}
#' \item{glfc}{A list (length=nsims) of vectors (length=ngenes) for phenotypic log fold change of all genes, ie nonDE=0 and DE=lfc.}
#' \item{blfc}{A list (length=nsims) of vectors (length=ngenes) for batch log fold change of all genes.}
#' \item{design}{Two group comparison}
#' @author Beate Vieth
#' @examples
#' \dontrun{
#' desettings <- DESetup(ngenes=10000,
#' nsims=25,
#' p.DE=0.2,
#' pLFC=function(x) sample(c(-1,1), size=x,replace=TRUE)*rgamma(x, 3, 3))
#' }
#' @rdname DESetup
#' @export
DESetup <- function(ngenes, nsims=25,
                    p.DE=0.1, p.B=NULL,
                    pLFC, bLFC=NULL,
                    sim.seed) {
  if (missing(sim.seed))
    sim.seed = sample(1:1000000, size = 1)
  set.seed(sim.seed)

  nDE = round(ngenes*p.DE)
  if(!is.null(p.B)) { nB = round(ngenes*p.B) }
  DEids = Bids = plfcs = blfcs = NULL
  for (i in 1:nsims) {
    # generate a random id for DE genes
    DEid <- sample(ngenes, nDE, replace = FALSE)
    DEids[[i]] <- DEid
    ## generate lfc for all genes: 0 for nonDE and LFC for DE
    plfc = rep(0, ngenes)
    plfc[DEid] = .setFC(pLFC, nDE)
    plfcs[[i]] = plfc
    if(!is.null(bLFC)) {
      # generate a random id for batch affected genes
      Bid <- sample(ngenes, nB, replace = FALSE)
      Bids[[i]] <- Bid
      ## generate lfc for all genes: 0 for nonDE and LFC for DE
      blfc = rep(0, ngenes)
      blfc[Bid] = .setFC(bLFC, nB)
      blfcs[[i]] = blfc
    }
  }

  sim.seed = as.vector(sample(1:1000000, size = nsims, replace = F))

  set.seed(NULL)

  ## return
  res <- c(list(DEid = DEids,
                Bid = Bids,
                pLFC = plfcs,
                bLFC = blfcs,
                ngenes = ngenes,
                nsims = nsims,
                p.DE = p.DE,
                p.B = p.B,
                sim.seed.DESetting = sim.seed),
           list(sim.seed = sim.seed), design = "2grp")
  return(res)

}

# SimSetup ----------------------------------------------------------------

#' @name SimSetup
#' @aliases SimSetup
#' @title DEA options for RNA-seq count simulations in two-group comparison.
#' @description This function adds user provided options for simulating RNA-seq data to RNAseq.SimSetup object. The resulting output list object is the input for \code{\link{simulateDE}} function.
#' @usage SimSetup(desetup,
#' params,
#' spike=NULL,
#' size.factors='equal',
#' downsample=FALSE,
#' geneset=FALSE)
#' @param desetup The RNAseq simulation parameters created by \code{\link{DESetup}}.
#' @param params The negative binomial parameters for simulations. This can be:
#' (1) The output of \code{\link{estimateParam}}.
#' (2) A string specifying the name of precalculated estimates, see details.
#' @param spike The spike-in simulation parameters created by \code{\link{estimateSpike}}. Default is \code{NULL}.
#' These are needed for applying spike-in-dependent normalisation methods, i.e. 'RUV' and 'BASiCS'.
#' @param size.factors Size factors representing sample-specific differences/biases in expected mean values of the NB distribution:
#' "equal" or "given". The default is "equal", i.e. equal size factor of 1.
#' If the user defines it as given, the size factors are sampled from the size factors provided by the output of \code{\link{estimateParam}}.
#' @param downsample Drawing the associated dispersions after determining effective mean expressions by size factors. Default is \code{FALSE}.
#' @param geneset Sampling with replacement or filling count tables low magnitude Poisson
#' when the estimated mean expression vector is shorter than the number of genes to be simulated.
#' Default is \code{FALSE}, i.e. sampling means with replacement.
#' @return A list with the following entries:
#' \item{desetup}{The RNAseq simulation parameters.}
#' \item{params}{The distributional parameters for simulations of genes.}
#' \item{spike}{The distributional parameters for simulations of spike-ins.}
#' \item{size.factors}{Size factor definition: "equal" means no difference in size factors between samples.
#' "given" means that the size factors will be randomly drawn from the size factors provided by \code{params}.
#' The user can also provide a list object containing sampling distributions per group (n1 and n2).
#' Defaul is \code{"equal"}.}
#' @author Beate Vieth
#' @examples
#' \dontrun{
#' desettings <- DESetup(ngenes=10000,
#' nsims=25, p.DE=0.2,
#' pLFC=function(x) sample(c(-1,1), size=x,replace=TRUE)*rgamma(x, 3, 3))
#' simsettings <- SimSetup(desetup=desettings,
#' params=estparam, size.factors='equal')
#' }
#' @rdname SimSetup
#' @export
SimSetup <- function(desetup,
                     params,
                     spike=NULL,
                     size.factors='equal',
                     downsample=FALSE,
                     geneset=FALSE) {

  ## return
  res <- c(desetup,
           params,
           list(spike=spike),
           size.factors=list(size.factors),
           downsample=list(downsample),
           geneset=list(geneset))
  attr(res, 'param.type') <- attr(params, 'param.type')
  attr(res, 'Distribution') <- attr(params, 'Distribution')
  return(res)

}


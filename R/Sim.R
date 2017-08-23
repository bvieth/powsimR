
# simulateDE --------------------------------------------------------------

#' @name simulateDE
#' @aliases simulateDE
#' @title Simulate Differential Expression
#' @description This function simulates RNA-seq count matrices considering differential expression specifications (number of samples per group, effect size, number of differential expressed genes etc.). The return object contains DE test results from all simulations as well as descriptive statistics.
#' @details simulateDE is the main function to simulate differential expression for RNA-seq experiments.
#' The simulation parameters are specified with \code{\link{SimSetup}}.
#' The user needs to specify the number of samples per group and the differential expression analysis method. \cr
#' It only stores and returns the DE test results (i.e. p-values). The error matrix calculations will be conducted with \code{\link{evaluateDE}}.\cr
#' @usage simulateDE(n1=c(20,50,100), n2=c(30,60,120),
#' sim.settings,
#' DEmethod,
#' Preclust=FALSE,
#' normalisation,
#' NCores=NULL,
#' verbose=TRUE)
#' @param n1,n2 Integer vectors specifying the number of biological replicates in each group. Default values are n1=c(20,50,100) and n2=c(30,60,120).
#' @param sim.settings This object specifies the simulation setup. This must be the return object from \code{\link{SimSetup}}.
#' @param DEmethod A character vector specifying the DE detection method to be used.
#' Available options are: limma-trend, limma-voom, edgeR-LRT, edgeR-QL, DESeq2,
#' ROTS, baySeq, NOISeq, EBSeq, DSS, MAST, scde, BPSC, scDD, monocle.
#' @param Preclust Whether to run a  hierarchical clustering prior to normalisation. Default is \code{FALSE}. This is implemented for scran and SCnorm.
#' For details, see \code{\link[scran]{quickCluster}}.
#' @param normalisation Normalisation method to use.
#' @param NCores integer positive number of cores for parallel processing, default is \code{NULL}, ie 1 core.
#' @param verbose Logical value to indicate whether to show progress report of simulations. Default is \code{TRUE}.
#' @return A list with the following fields.
#' \item{pvalue, fdr}{3D array (ngenes * N * nsims) for p-values and FDR from each simulation.
#' Note that FDR values will be empty and the calculation will be done by \code{\link{evaluateDE}} whenever applicable.}
#' \item{mu,disp,dropout}{3D (ngenes * N * nsims) array for mean, dispersion and dropout of library size factor normalized read counts.}
#' \item{elfc,rlfc}{3D array (ngenes * N * nsims) for log2 fold changes (LFC):
#' elfc is for the DE tool estimated LFC; rlfc is for the LFC estimated from the normalised read counts.}
#' \item{sf.values,gsf.values}{3D array (ngenes * N * nsims) for size factor estimates.
#' Global estimates per sample in sf.values; Gene- and sample-wise estimates in gsf.values only for SCnorm normalisation.}
#' \item{true.designs}{3D array (ngenes * N * nsims) for group assignment specifications.}
#' \item{sim.settings}{The input sim.settings to which the specifications of \code{simulateDE} is added.}
#' \item{time.taken}{The time taken for each simulation, given for preprocessing, normalisation, differential expression testing and moment estimation.}
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}},  \code{\link{insilicoNBParam}} for negative binomial parameter specifications;\cr
#'  \code{\link{DESetup}}, \code{\link{SimSetup}} for simulation setup
#' @examples
#' \dontrun{
#' # download count table
#' githubURL <- "https://github.com/bvieth/powsimRData/raw/master/data-raw/kolodziejczk_cnts.rda"
#' download.file(url = githubURL, destfile= "kolodziejczk_cnts.rda", method = "wget")
#' load('kolodziejczk_cnts.rda')
#' kolodziejczk_cnts <- kolodziejczk_cnts[, grep('standard', colnames(kolodziejczk_cnts))]
#' ## estimate NB parameters:
#' TwoiLIF.params = estimateParam(countData=cnts
#'                           Distribution='NB',
#'                           RNAseq="singlecell",
#'                           normalisation='scran',
#'                           NCores=NULL,
#'                           sigma=1.96, verbose=TRUE)
#' ## define DE settings:
#' desettings <- DESetup(ngenes=10000,
#' nsims=25, p.DE=0.2,
#' pLFC=function(x) sample(c(-1,1), size=x,replace=TRUE)*rgamma(x, 3, 3))
#' ## define simulation settings for Kolodziejczk:
#' simsettings <- SimSetup(desetup=desettings, params=TwoiLIF.params, size.factors='given')
#' ## run simulations:
#' simres <- simulateDE(n1=c(50,100,300), n2=c(50,120,500),
#' sim.settings,
#' DEmethod='MAST',
#' Preclust=T,
#' normalisation='scran',
#' NCores=10,
#' verbose=TRUE)
#' ## if parallel computation unavailable, consider ROTS as DEmethod
#' }
#' @rdname simulateDE
#' @importFrom stats setNames
#' @export
simulateDE <- function(n1=c(20,50,100), n2=c(30,60,120),
                       sim.settings,
                       DEmethod,
                       Preclust=FALSE,
                       normalisation,
                       NCores=NULL,
                       verbose=TRUE) {

  Preprocess = NULL
  GeneSelect=NULL
  DimReduce=NULL
  ClustMethod=NULL
  spikeIns=FALSE

  if (!length(n1) == length(n2)) { stop("n1 and n2 must have the same length!") }
  if(isTRUE(spikeIns) && is.null(sim.settings$spike)) {
    stop(message(paste0("For the simulation of  spike-ins, fitting information is needed but there is no 'spike' object in 'sim.settings'.  Please consult the function estimateSpike for spike fitting and SimSetup for creating simulation setup object!")))
  }

  if(!is.null(NCores) && DEmethod %in% c("edgeR-LRT","edgeR-QLRT", 'limma-voom', "limma-trend", "NOISeq", "DSS", "EBSeq", "ROTS")) {
    message(paste0(DEmethod, " has no parallel computation option!"))
  }

  if(sim.settings$RNAseq == "singlecell" && DEmethod %in% c("edgeR-LRT", "edgeR-QL", "limma-voom", "limma-trend", "DESeq2", "baySeq", "NOISeq", "DSS", "EBSeq")) {
    message(paste0(DEmethod, " is developed for bulk RNA-seq experiments."))
  }

  if(sim.settings$RNAseq == "bulk" && DEmethod %in% c("MAST", 'BPSC')) {
    message(paste0(DEmethod, " is developed for single cell RNA-seq experiments."))
  }

  if(sim.settings$RNAseq == "bulk" && DEmethod %in% c('scde', 'scDD', 'monocle')) {
    stop(message(paste0(DEmethod, " is only developed and implemented for single cell RNA-seq experiments.")))
  }

  # define the maximal count matrix for simulations
  max.n = max(n1,n2)
  min.n = min(n1, n2)

  # append additional settings of simulateDE to sim.settings
  sim.settings$n1 = n1
  sim.settings$n2 = n2
  sim.settings$normalisation = normalisation
  sim.settings$DEmethod = DEmethod
  sim.settings$spikeIns = spikeIns
  sim.settings$NCores = NCores
  sim.settings$Preclust = Preclust
  sim.settings$Preprocess = Preprocess
  sim.settings$DimReduce = DimReduce
  sim.settings$GeneSelect = GeneSelect
  sim.settings$ClustMethod = ClustMethod
  sim.settings$clustNumber = ifelse(sim.settings$design=="2grp", 2, NULL)
  sim.settings$PreclustNumber = ifelse(isTRUE(Preclust), min.n, NULL)

  my.names = paste0(n1,"vs",n2)

  #set up output arrays
  pvalues = fdrs = elfcs = rlfcs = mus = disps = dropouts = array(NA,dim=c(sim.settings$ngenes,length(n1), sim.settings$nsims))
  time.taken = array(NA,dim = c(5,length(n1), sim.settings$nsims),
                     dimnames = list(c('Preprocess', "Normalisation", "Clustering", "DE", "Moments"),
                                     NULL, NULL))
  true.sf = stats::setNames(replicate(length(n1),NULL),my.names)
  est.sf = stats::setNames(replicate(length(n1),NULL),my.names)
  true.sf <- lapply(1:length(true.sf), function(x) {
    true.sf[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })
  est.sf <- lapply(1:length(est.sf), function(x) {
    est.sf[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })

  if(sim.settings$normalisation=="SCnorm") {
    est.gsf = stats::setNames(replicate(length(n1),NULL),my.names)
    est.gsf <- lapply(1:length(est.gsf), function(x) {
      est.gsf[[x]] = array(NA,dim=c(sim.settings$nsims, sim.settings$ngenes, n1[x] + n2[x]))
    })
  }
  if(!sim.settings$normalisation=="SCnorm") {
    est.gsf = NULL
  }
  # est.gsf = NULL

  true.designs = stats::setNames(replicate(length(n1),NULL),my.names)
  true.designs <- lapply(1:length(true.designs), function(x) {
    true.designs[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })

  if(!is.null(sim.settings$ClustMethod)) {
    def.designs = stats::setNames(replicate(length(n1),NULL),my.names)
    def.designs <- lapply(1:length(def.designs), function(x) {
      def.designs[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
    })
  }

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
    if (verbose) {message(paste0("Generating RNA seq read counts")) }
    gene.data = .simRNAseq.2grp(simOptions = tmp.simOpts,
                                n1 = max.n, n2 = max.n)

    ## generate spike-in read counts
    if(isTRUE(spikeIns)) {
    if (verbose) { message(paste0("Generating spike-in read counts")) }
      spike.data = .simSpike(SpikeOptions = tmp.simOpts$spike, n1 = max.n, n2 = max.n)
    }
    if(!isTRUE(spikeIns)) {
      spike.data = NULL
      spike.info = NULL
    }

    ## generate mean fragment lengths for samples
    if(!is.null(tmp.simOpts$MeanFragLengths)) {
      if (verbose) { message(paste0("Sampling from observed mean fragment lengths")) }
      meanfrag.data = sample(tmp.simOpts$MeanFragLengths, max.n+max.n, replace = TRUE)
      names(meanfrag.data) = colnames(gene.data$counts)
    }
    if(is.null(tmp.simOpts$MeanFragLengths)) {
      meanfrag.data = NULL
    }

    ## match sampled gene names with given gene lengths
    if(!is.null(tmp.simOpts$Lengths)) {
      if (verbose) { message(paste0("Associating gene lengths with sampled gene expression")) }
      gene.id = sub('_([^_]*)$', '', rownames(gene.data$counts))
      length.data = tmp.simOpts$Lengths
      length.data = length.data[match(gene.id,names(length.data))]
    }
    if(is.null(tmp.simOpts$Lengths)) {
      length.data = NULL
    }


    ##  for different sample sizes
    for (j in seq(along=n1)) {
      Nrep1 = n1[j]
      Nrep2 = n2[j]

      tmp.simOpts$PreclustNumber = min(Nrep1, Nrep2)

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
      ## filter out gene lengths belonging to zero expression
      if(!is.null(length.data)) {
       length.data = length.data[ix.valid]
      }

      ## take a subsample of simulated spike-ins
      if(!is.null(spike.data)) {
        sim.spike <- spike.data$counts
        spike.valid = rowSums(sim.spike) > 0
        count.spike = sim.spike[spike.valid, idx, drop=FALSE]
        if(!is.null(sim.settings$spike)) {
          spike.info <- tmp.simOpts$spike$Input$spikeInfo[rownames(tmp.simOpts$spike$Input$spikeInfo)
                                              %in% rownames(count.spike), , drop = FALSE]
          spike.info <- spike.info[match(rownames(count.spike),
                                         rownames(spike.info)), , drop = FALSE]
        }
      }
      ## take a subsample of mean fragment lengths
      if(!is.null(meanfrag.data)) {
        meanfrag.data = meanfrag.data[idx]
      }

      ## perform filtering / imputation (OPTIONAL)
      start.time.preprocess <- Sys.time()
      if(!is.null(Preprocess)) {
        if (verbose) { message(paste0("Applying preprocessing")) }
        filter.data <- .preprocess.calc(Preprocess=tmp.simOpts$Preprocess,
                                        countData=count.data,
                                        NCores=tmp.simOpts$NCores)
        ixx.valid <- rownames(sim.cnts) %in% rownames(filter.data)
        ix.valid <- ixx.valid
        count.data <- filter.data
        if(!is.null(length.data)) {
          gene.id = sub('_([^_]*)$', '', rownames(count.data))
          length.data = length.data[match(gene.id,names(length.data))]
        }
      }
      end.time.preprocess <- Sys.time()

      ## perform normalisation
      if (verbose) { message(paste0("Normalizing read counts")) }
      start.time.norm <- Sys.time()
      norm.data <- .norm.calc(normalisation=tmp.simOpts$normalisation,
                             countData=count.data,
                             spikeData=count.spike,
                             spikeInfo=spike.info,
                             batchData=NULL,
                             Lengths=length.data,
                             MeanFragLengths=meanfrag.data,
                             PreclustNumber=tmp.simOpts$PreclustNumber,
                             NCores=tmp.simOpts$NCores)
      end.time.norm <- Sys.time()

      def.design <- true.design

      ## perform group assignment (OPTIONAL)
      if(!is.null(DimReduce) && !is.null(ClustMethod)) {
        if (verbose) { message(paste0("Applying dimension reduction and clustering")) }
        start.time.clust <- Sys.time()
        classify.data <- .classify.calc(GeneSelect=tmp.simOpts$GeneSelect,
                                    DimReduce=tmp.simOpts$DimReduce,
                                    ClustMethod=tmp.simOpts$ClustMethod,
                                    clustNumber=tmp.simOpts$clustNumber,
                                    normData=norm.data,
                                    Lengths=length.data,
                                    MeanFragLengths=meanfrag.data,
                                    countData=count.data,
                                    spikeData=count.spike,
                                    spikeInfo=spike.info,
                                    spikeIns=spikeIns,
                                    verbose=verbose)
        def.design <- classify.data
        end.time.clust <- Sys.time()
      }

      ## create an DE options object to pass into DE detection
      DEOpts <- list(designs=def.design, p.DE=tmp.simOpts$p.DE)

      ## Run DE detection
      if (verbose) { message(paste0("Running DE tool")) }
      start.time.DE <- Sys.time()
      res.de = .de.calc(DEmethod=tmp.simOpts$DEmethod,
                        normData=norm.data,
                        countData=count.data,
                        DEOpts=DEOpts,
                        spikeData=count.spike,
                        spikeInfo=spike.info,
                        Lengths=length.data,
                        MeanFragLengths=meanfrag.data,
                        NCores=tmp.simOpts$NCores)
      end.time.DE <- Sys.time()

      ## estimate moments of read counts simulated
      start.time.NB <- Sys.time()
      res.params <- .run.params(countData=count.data,
                                normData=norm.data,
                                group=DEOpts$def.design)
      end.time.NB <- Sys.time()

      # generate empty vectors
      pval = fdr = est.lfc = raw.lfc = mu.tmp = disp.tmp = p0.tmp = rep(NA, nrow(sim.cnts))
      ## extract results of DE testing
      pval[ix.valid] = res.de$pval
      fdr[ix.valid] = res.de$fdr
      est.lfc[ix.valid] = res.de$lfc
      raw.lfc[ix.valid] = res.params$lfc
      mu.tmp[ix.valid] = res.params$means
      disp.tmp[ix.valid] = res.params$dispersion
      p0.tmp[ix.valid] = res.params$dropout
      # copy it in 3D array of results
      pvalues[,j,i] = pval
      fdrs[,j,i] = fdr
      elfcs[,j,i] = est.lfc
      rlfcs[,j,i] = raw.lfc
      mus[,j,i] = mu.tmp
      disps[,j,i] = disp.tmp
      dropouts[,j,i] = p0.tmp
      true.sf[[j]][i,] = gene.sf
      est.sf[[j]][i,] = norm.data$size.factors

      if(attr(norm.data, 'normFramework') == 'SCnorm') {
        est.gsf[[j]][i, ix.valid, ] = norm.data$scale.factors
      }

      true.designs[[j]][i,] = true.design
      # time taken for each step
      # copy designs into list of matrices
      if(!is.null(ClustMethod)) {
        time.taken.clust <- difftime(end.time.clust,
                                     start.time.clust,
                                     units="mins")
        def.designs[[j]][i,] = def.design
      }
      if (is.null(ClustMethod)) {
        time.taken.clust = NA
        def.designs = NULL
      }
      if(!is.null(Preprocess)) {
        time.taken.preprocess <-  difftime(end.time.preprocess,
                                           start.time.preprocess,
                                           units="mins")
      }
      if(is.null(Preprocess)) {
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
                      time.taken.clust,
                      time.taken.DE,
                      time.taken.NB)

      # copy time taken in 2D array of time taken
      time.taken[,j,i] = timing

    }
  }

  ## return
  list(pvalue = pvalues,
       fdr = fdrs,
       elfc = elfcs,
       rlfc = rlfcs,
       mu = mus,
       disp = disps,
       dropout = dropouts,
       true.sf = true.sf,
       est.sf = est.sf,
       est.gsf = est.gsf,
       time.taken = time.taken,
       true.designs = true.designs,
       def.designs = def.designs,
       sim.settings = sim.settings)
}


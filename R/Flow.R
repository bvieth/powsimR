
# simulateFlow --------------------------------------------------------------

#' @name simulateFlow
#' @aliases simulateFlow
#' @title Simulate Single Cell Analysis Pipeline Expression
#' @description Blabla
#' @usage simulateFlow(n1=c(20,50,100), n2=c(30,60,120),
#' sim.settings,
#' DEmethod,
#' normalisation,
#' Preclust=FALSE,
#' Prefilter = NULL,
#' Impute = NULL,
#' DEFilter = FALSE,
#' GeneSelect=NULL,
#' DimReduce=NULL,
#' ClustMethod=NULL,
#' Pipeline=NULL,
#' spikeIns=FALSE,
#' NCores=NULL,
#' verbose=TRUE)
#' @param n1,n2 Integer vectors specifying the number of biological replicates in each group. Default values are n1=c(20,50,100) and n2=c(30,60,120).
#' @param sim.settings This object specifies the simulation setup. This must be the return object from \code{\link{SimSetup}}.
#' @param DEmethod A character vector specifying the DE detection method to be used.
#' Please consult the Details section for available options.
#' @param normalisation Normalisation method to use.
#' Please consult the Details section for available options.
#' @param Preclust A logical vector indicating whether to run a hierarchical clustering prior to normalisation.
#' This is implemented for scran only. Default is \code{FALSE}.
#' For details, see \code{\link[scran]{quickCluster}}.
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
#' @param GeneSelect  A character vector specifying the gene selection method.
#' Default is \code{NULL}, i.e. no gene selection.
#' Availabe options are: HVG, Gini-Index.
#' @param DimReduce A character vector specifying the dimension reduction method. Default is \code{NULL}, i.e. no dimension reduction.
#' Availabe options are: Euclidean, Spearman, PCA, t-SNE-expr, t-SNE-Euclidean, PCA+t-SNE, CIDR.
#' @param ClustMethod A character vector specifying the clustering method. Default is \code{NULL}, i.e. no clustering.
#' Availabe options are: kmeans, PAM.
#' @param Pipeline A character vector specifying out-of-the-box single cell clustering pipelines.
#' Available options are: CIDR_free and CIDR_bound.
#' @param spikeIns Logical value to indicate whether to simulate spike-ins. Default is \code{FALSE}.
#' @param NCores integer positive number of cores for parallel processing, default is \code{NULL}, ie 1 core.
#' @param verbose Logical value to indicate whether to show progress report of simulations. Default is \code{TRUE}.
#' @return A list with the following fields:
#' \item{pvalue, fdr}{3D array (ngenes * N * nsims) for p-values and FDR from each simulation.
#' Note that FDR values will be empty and the calculation will be done by \code{\link{evaluateDE}} whenever applicable.}
#' \item{mu,disp,dropout}{3D (ngenes * N * nsims) array for mean, dispersion and dropout of library size factor normalized read counts.}
#' \item{elfc,rlfc}{3D array (ngenes * N * nsims) for log2 fold changes (LFC):
#' elfc is for the DE tool estimated LFC; rlfc is for the LFC estimated from the normalised read counts.}
#' \item{sf.values,gsf.values}{3D array (ngenes * N * nsims) for size factor estimates.
#' Global estimates per sample in sf.values; Gene- and sample-wise estimates in gsf.values only for SCnorm normalisation.}
#' \item{true.designs,def.designs}{3D array (ngenes * N * nsims) for group assignment specifications.
#' true.designs for the simulated group assignment; def.designs for the group assignment determined after clustering.}
#' \item{sim.settings}{The input sim.settings to which the specifications of \code{simulateDE} is added.}
#' \item{time.taken}{The time taken for each simulation, given for preprocessing, normalisation, clustering, differential expression testing and moment estimation.}
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}},  \code{\link{insilicoNBParam}} for negative binomial parameter specifications;\cr
#'  \code{\link{DESetup}}, \code{\link{SimSetup}} for simulation setup
#' @examples
#' \dontrun{
#' ## not yet
#' }
#' @rdname simulateFlow
#' @importFrom stats setNames
#' @export
simulateFlow <- function(n1=c(20,50,100), n2=c(30,60,120),
                         sim.settings,
                         DEmethod,
                         normalisation,
                         Preclust = FALSE,
                         Prefilter = NULL,
                         Impute = NULL,
                         DEFilter = FALSE,
                         GeneSelect = NULL,
                         DimReduce = NULL,
                         ClustMethod = NULL,
                         Pipeline = NULL,
                         spikeIns = FALSE,
                         NCores = NULL,
                         verbose = TRUE) {
  if (!length(n1) == length(n2)) { stop("n1 and n2 must have the same length!") }
  if(isTRUE(spikeIns) && is.null(sim.settings$spike)) {
    stop(message(paste0("For the simulation of  spike-ins, fitting information is needed but there is no 'spike' object in 'sim.settings'.  Please consult the function estimateSpike for spike fitting and SimSetup for creating simulation setup object!")))
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
    normalisation = "none"
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
  sim.settings$normalisation = normalisation
  sim.settings$DEmethod = DEmethod
  sim.settings$spikeIns = spikeIns
  sim.settings$NCores = NCores
  sim.settings$Pipeline = Pipeline
  sim.settings$Preclust = Preclust
  sim.settings$Prefilter = Prefilter
  sim.settings$Impute = Impute
  sim.settings$DEFilter = DEFilter
  sim.settings$DimReduce = DimReduce
  sim.settings$GeneSelect = GeneSelect
  sim.settings$ClustMethod = ClustMethod
  sim.settings$clustNumber = ifelse(sim.settings$design=="2grp", 2, NULL)
  if(isTRUE(Preclust)) {PreclustNumber <- min.n}
  if(!isTRUE(Preclust)) {PreclustNumber <- NULL}
  sim.settings$PreclustNumber = PreclustNumber

  if (verbose) { message(paste0("Preparing output arrays.")) }

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

  true.designs = stats::setNames(replicate(length(n1),NULL),my.names)
  true.designs <- lapply(1:length(true.designs), function(x) {
    true.designs[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
  })

  def.designs = stats::setNames(replicate(length(n1),NULL),my.names)
  def.designs <- lapply(1:length(def.designs), function(x) {
    def.designs[[x]] = matrix(NA, nrow = sim.settings$nsims, ncol = n1[x] + n2[x])
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
    if (verbose) {message(paste0("Generating RNA seq read counts")) }
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
      if (verbose) { message(paste0("Applying ", normalisation, " normalisation")) }
      start.time.norm <- Sys.time()
      norm.data <- .norm.calc(normalisation=tmp.simOpts$normalisation,
                              countData=fornorm.count.data,
                              spikeData=count.spike,
                              spikeInfo=spike.info,
                              batchData=def.design,
                              Lengths=length.data,
                              MeanFragLengths=meanfrag.data,
                              PreclustNumber=tmp.simOpts$PreclustNumber,
                              NCores=tmp.simOpts$NCores,
                              verbose=verbose)
      end.time.norm <- Sys.time()

      ## perform group assignment (OPTIONAL)
      if(!is.null(DimReduce) && !is.null(ClustMethod)) {
        if (verbose) { message(paste0("Applying dimension reduction and clustering")) }
        start.time.clust <- Sys.time()
        classify.data <- .classify.calc(Pipeline=tmp.simOpts$Pipeline,
                                        GeneSelect=tmp.simOpts$GeneSelect,
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
                                        NCores=tmp.simOpts$NCores,
                                        verbose=verbose)
        def.design <- classify.data
        end.time.clust <- Sys.time()
      }

      ## create an DE options object to pass into DE detection
      DEOpts <- list(designs=def.design, p.DE=tmp.simOpts$p.DE)

      ## Run DE detection
      start.time.DE <- Sys.time()
      if(!isTRUE(tmp.simOpts$DEFilter)) {
        if (verbose) { message(paste0("Applying ", DEmethod, " for DE analysis on raw count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$DEmethod,
                          normData=norm.data,
                          countData=count.data,
                          DEOpts=DEOpts,
                          spikeData=count.spike,
                          spikeInfo=spike.info,
                          NCores=tmp.simOpts$NCores,
                          verbose=verbose)
      }
      if(isTRUE(tmp.simOpts$DEFilter)) {
        if (verbose) { message(paste0("Applying ", DEmethod, " for DE analysis on imputed/filtered count data.")) }
        res.de = .de.calc(DEmethod=tmp.simOpts$DEmethod,
                          normData=norm.data,
                          countData=fornorm.count.data,
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

      ## estimate moments of read counts simulated
      start.time.NB <- Sys.time()
      res.params <- .run.params(countData=count.data,
                                normData=norm.data,
                                group=DEOpts$designs)
      end.time.NB <- Sys.time()

      # generate empty vectors
      pval = fdr = est.lfc = raw.lfc = mu.tmp = disp.tmp = p0.tmp = rep(NA, nrow(sim.cnts))
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
      # copy them in 3D array of results
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
        allgenes <- rownames(sim.cnts)
        testedgenes <- rownames(norm.data$scale.factors)
        ixx.valid <- allgenes %in% testedgenes
        est.gsf[[j]][i, ixx.valid, ] = norm.data$scale.factors
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
        def.designs[[j]][i,] = true.design
      }
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
                      time.taken.clust,
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
                  disp = disps,
                  dropout = dropouts,
                  true.sf = true.sf,
                  est.sf = est.sf,
                  est.gsf = est.gsf,
                  time.taken = time.taken,
                  true.designs = true.designs,
                  def.designs = def.designs,
                  sim.settings = sim.settings)

  attr(res.out, 'Simulation') <- "Flow"
  return(res.out)
}




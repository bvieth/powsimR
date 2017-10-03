
# EVALUATE SETUP ----------------------------------------------------------
#' @name evaluateSim
#' @aliases evaluateSim
#' @title Compute the performance related metrics from simulation results
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes several metrics that give an indication of the simulation setup performance.
#' @usage evaluateSim(simRes)
#' @param simRes The result from \code{\link{simulateDE}} or \code{\link{simulateFlow}}.
#' @return A list with the following entries:
#' \item{Log2FoldChange}{The root mean square error and absolute mean error of
#' log2 fold change differences between estimated LFC and simulated LFC. Furthermore, the fraction of missing estimated LFC.}
#' \item{SizeFactors}{The median absolute deviation between the estimated and true size factors,
#' the residual error of a robust linear model and the ratio between estimated and true size factors of the two groups.}
#' \item{Clustering}{The adjusted rand index between the true group assignment and derived group assignment by clustering (only for \code{\link{simulateFlow}}).}
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}} for negative binomial parameters,
#' \code{\link{SimSetup}} and
#' \code{\link{DESetup}} for setting up simulation parameters and
#' \code{\link{simulateDE}} and \code{\link{simulateFlow}} for simulating differential expression.
#' @examples
#' \dontrun{
#' ## not yet
#' }
#' @rdname evaluateSim
#' @importFrom stats mad
#' @importFrom mclust adjustedRandIndex
#' @export
evaluateSim <- function(simRes) {

  if (attr(simRes, 'Simulation') == 'Flow') {

    # simulation parameters
    Nreps1 = simRes$sim.settings$n1
    Nreps2 = simRes$sim.settings$n2
    ngenes = simRes$sim.settings$ngenes
    sim.opts = simRes$sim.settings
    DEids = simRes$sim.settings$DEid
    tlfcs = simRes$sim.settings$pLFC
    nsims = simRes$sim.settings$nsims

    # estimated parameters
    elfcs = simRes$elfc
    tsfs = simRes$true.sf
    esfs = simRes$est.sf

    t.designs = simRes$true.designs
    e.designs = simRes$def.designs

    # create output objects
    my.names = paste0(Nreps1, "vs", Nreps2)
    # error in log2 fold changes
    lfc.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 15,
             dimnames = list(c(paste0("Sim", 1:nsims)),
                             c(paste0(rep(x=c("ALL", "DE", "EE"),each=5),"_",
                                      c('RMSE_Value', "MAE_Value", "RMSE_NAFraction", "MAE_NAFraction", "MSE_ErrorFit")))))
    })
    names(lfc.error.mat) <- my.names

    # error in size factors
    sf.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 4,
             dimnames = list(c(paste0("Sim_", 1:nsims)),
                             c("MAD.LFC.SF", "Error.Fit", "Ratio.1","Ratio.2")))
    })
    names(sf.error.mat) <- my.names

    # error in clustering
    clust.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 1,
             dimnames = list(c(paste0("Sim_", 1:nsims)),
                             c("RandIndex")))
    })
    names(clust.error.mat) <- my.names

    ## loop over simulation and replicates
    for(i in 1:nsims) {
      # DE flag
      DEid = DEids[[i]]
      Zg = rep(0, ngenes)
      Zg[DEid] = 1
      # true log fold change of all genes
      all.tlfc = tlfcs[[i]]
      # true log fold change of DE genes
      de.tlfc = all.tlfc[which(Zg==1)]
      # true log fold change of EE genes
      ee.tlfc = all.tlfc[which(Zg==0)]

      for(j in seq(along=Nreps1)) {

        ## LOG2 FOLD CHANGES
        # estimated log fold change of all genes
        all.elfc = elfcs[, j, i]
        # estimated log fold changes of EE genes
        ix.ee.lfc = which(Zg==0)
        ee.lfc = all.elfc[ix.ee.lfc]
        # estimated log fold change of DE genes
        ix.de.lfc = which(Zg==1)
        de.lfc = all.elfc[ix.de.lfc]
        # estimate mean squared error and absolute error
        all.error <- .lfc.evaluate(truth=all.tlfc, estimated=all.elfc)
        ee.error <- .lfc.evaluate(truth=ee.tlfc, estimated=ee.lfc)
        de.error <- .lfc.evaluate(truth=de.tlfc, estimated=de.lfc)
        error.est <- c(all.error, ee.error, de.error)
        lfc.error.mat[[j]][i,] = error.est

        ## SIZE FACTORS
        # true sf over all samples, center to mean=1
        tsf = tsfs[[j]][i, ]
        n.tsf = tsf*length(tsf)/sum(tsf)

        # estimated sf over all sample, center to mean=1
        esf = esfs[[j]][i,]
        n.esf = esf*length(esf)/sum(esf)

        # MAD of log fold change difference between estimated and true size factors
        lfc.nsf = log2(n.esf) - log2(n.tsf)
        mad.nsf = stats::mad(lfc.nsf)

        # error of estimation
        error.sf <- .fiterror.sf(estimated.sf = esf, true.sf = tsf)

        # ratio of estimated and true size factors per true group assignment
        t.design = t.designs[[j]][i,]
        e.design = e.designs[[j]][i,]
        ratio.sf <- .ratio.sf(estimated.nsf = n.esf,
                              true.nsf = n.tsf,
                              group=t.design)
        sf.res <- unlist(c(mad.nsf, error.sf, ratio.sf))
        names(sf.res) <- NULL

        sf.error.mat[[j]][i,] = sf.res

        ## CLUSTERING
        randindex <- mclust::adjustedRandIndex(t.design, e.design)
        clust.error.mat[[j]][i,] = randindex

      }
    }

    output <- list(Log2FoldChange=lfc.error.mat,
                   SizeFactors=sf.error.mat,
                   Clustering=clust.error.mat,
                   sim.settings=simRes$sim.settings)
  }

  if (attr(simRes, 'Simulation') == 'DE') {
    # simulation parameters
    Nreps1 = simRes$sim.settings$n1
    Nreps2 = simRes$sim.settings$n2
    ngenes = simRes$sim.settings$ngenes
    sim.opts = simRes$sim.settings
    DEids = simRes$sim.settings$DEid
    tlfcs = simRes$sim.settings$pLFC
    nsims = simRes$sim.settings$nsims

    # estimated parameters
    elfcs = simRes$elfc
    tsfs = simRes$true.sf
    esfs = simRes$est.sf

    t.designs = simRes$true.designs

    # create output objects
    my.names = paste0(Nreps1, "vs", Nreps2)
    # error in log2 fold changes
    lfc.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 15,
             dimnames = list(c(paste0("Sim", 1:nsims)),
                             c(paste0(rep(x=c("ALL", "DE", "EE"),each=5),"_",
                                      c('RMSE_Value', "MAE_Value", "RMSE_NAFraction", "MAE_NAFraction", "MSE_ErrorFit")))))
    })
    names(lfc.error.mat) <- my.names

    # error in size factors
    sf.error.mat <- lapply(1:length(my.names), function(x) {
      matrix(NA, nrow =  nsims, ncol = 4,
             dimnames = list(c(paste0("Sim_", 1:nsims)),
                             c("MAD.LFC.SF", "Error.Fit", "Ratio.1","Ratio.2")))
    })
    names(sf.error.mat) <- my.names

    ## loop over simulation and replicates
    for(i in 1:nsims) {
      # DE flag
      DEid = DEids[[i]]
      Zg = rep(0, ngenes)
      Zg[DEid] = 1
      # true log fold change of all genes
      all.tlfc = tlfcs[[i]]
      # true log fold change of DE genes
      de.tlfc = all.tlfc[which(Zg==1)]
      # true log fold change of EE genes
      ee.tlfc = all.tlfc[which(Zg==0)]

      for(j in seq(along=Nreps1)) {

        ## LOG2 FOLD CHANGES
        # estimated log fold change of all genes
        all.elfc = elfcs[, j, i]
        # estimated log fold changes of EE genes
        ix.ee.lfc = which(Zg==0)
        ee.lfc = all.elfc[ix.ee.lfc]
        # estimated log fold change of DE genes
        ix.de.lfc = which(Zg==1)
        de.lfc = all.elfc[ix.de.lfc]
        # estimate mean squared error and absolute error
        all.error <- .lfc.evaluate(truth=all.tlfc, estimated=all.elfc)
        ee.error <- .lfc.evaluate(truth=ee.tlfc, estimated=ee.lfc)
        de.error <- .lfc.evaluate(truth=de.tlfc, estimated=de.lfc)
        error.est <- c(all.error, ee.error, de.error)
        lfc.error.mat[[j]][i,] = error.est

        ## SIZE FACTORS
        # true sf over all samples, center to mean=1
        tsf = tsfs[[j]][i, ]
        n.tsf = tsf*length(tsf)/sum(tsf)

        # estimated sf over all sample, center to mean=1
        esf = esfs[[j]][i,]
        n.esf = esf*length(esf)/sum(esf)

        # MAD of log fold change difference between estimated and true size factors
        lfc.nsf = log2(n.esf) - log2(n.tsf)
        mad.nsf = stats::mad(lfc.nsf)

        # error of estimation
        error.sf <- .fiterror.sf(estimated.sf = esf, true.sf = tsf)

        # ratio of estimated and true size factors per true group assignment
        t.design = t.designs[[j]][i,]
        ratio.sf <- .ratio.sf(estimated.nsf = n.esf,
                              true.nsf = n.tsf,
                              group=t.design)
        sf.res <- unlist(c(mad.nsf, error.sf, ratio.sf))
        names(sf.res) <- NULL

        sf.error.mat[[j]][i,] = sf.res
      }
    }

    output <- list(Log2FoldChange=lfc.error.mat,
                   SizeFactors=sf.error.mat,
                   sim.settings=simRes$sim.settings)
  }


  return(output)
}

# EVALUATE DIFFERENTIAL EXPRESSION ----------------------------------------

#' @name evaluateDE
#' @aliases evaluateDE
#' @title Compute the confusion matrix-related quantities from simulation results
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes quantities of the confusion matrix of classification testing
#' @usage evaluateDE(simRes,
#' alpha.type=c("adjusted","raw"),
#' MTC=c('BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni', 'Storey', 'IHW'),
#' alpha.nominal=0.1,
#' stratify.by=c("mean", "dispersion", "dropout", "lfc"),
#' filter.by=c("none", "mean", "dispersion", "dropout"),
#' strata.filtered=1, target.by=c("lfc", "effectsize"), delta=0)
#' @param simRes The result from \code{\link{simulateDE}}.
#' @param MTC Multiple testing correction method to use. Available options are
#' 1) see \link[stats]{p.adjust.methods},
#' 2) Storey's qvalue see \link[qvalue]{qvalue} and
#' 3) Independent Hypothesis Weighting considering mean expression as covariate (see \link[IHW]{ihw}).
#' Default is \code{BH}, i.e. Benjamini-Hochberg FDR correction method.
#' @param alpha.type A string to represent the way to call DE genes.
#'  Available options are \code{"adjusted"} i.e. applying multiple testing correction and
#'  \code{"raw"} i.e. using p-values. Default is \code{"adjusted"}.
#' @param alpha.nominal The nomial value of significance. Default is 0.1.
#' @param stratify.by A string to represent the way to stratify genes.
#' Available options are \code{"mean"}, \code{"dispersion"}, \code{"dropout"} and \code{"lfc"},
#' for stratifying genes by average expression levels, dispersion, dropout rates or estimated log2 fold changes.
#' @param filter.by A string to represent the way to filter genes.
#' This is used in conjunction with strata.filtered for gene filtering.
#' Available options are \code{"none"}, \code{"mean"}, \code{"dispersion"} and \code{"dropout"}.
#' \code{"none"} stands for no filtering, thus all genes will be considered.
#' \code{"mean"} stands for filtering based on average gene expression levels.
#' \code{"dispersion"} stands for filtering based on gene expression dispersion.
#' \code{"dropout"} stands for filtering based on dropout rates.
#' @param strata.filtered The strata to be filtered out in computing error matrix-related quantities.
#' Genes falling into these strata will be excluded. See "Details" for more description of gene filtering.
#' @param target.by A string to specify the method to define "biologically important" DE genes.
#' Available options are (1) \code{"lfc"}: interesting genes are defined by absolute log2 fold changes.
#' (2) \code{"effectsize"}: interesting genes are defined by
#' absolute log2 fold changes divided by the square root of 1/(mean+dispersion).
#' @param delta A threshold used for defining "biologically important" genes.
#' Genes with absolute log2 fold changes (when target.by is "lfc")
#' or effect sizes (when target.by is "effectsize") greater than this value
#' are deemed DE in error rates calculations. If \code{delta=0} then no threshold is applied. See "Details" for more description.
#' @return A list with the following entries:
#' \item{TN, TP, FP, FN, TNR, TPR, FPR, FNR, FDR}{3D array representing the number of true negatives, true positives, false positives,
#' false negatives and their proportions/rates as well as false discovery rate
#' for all simulation settings. The dimension of the arrays are nstrata * N * nsims.
#' Here nstrata is number of specified strata.
#' N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TN.marginal, TP.marginal, FP.marginal, FN.marginal}{Matrix representing the number of true negatives, true positives, false positives,
#'  false negatives for all simulation settings.
#'  The dimension of the matrices are N * nsims.
#'  Here N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TNR.marginal, TPR.marginal, FPR.marginal, FNR.marginal, FDR.marginal}{Matrix representing the marginal rates for all simulation settings.
#' The dimension of the matrices are N * nsims.}
#' \item{stratagenes, stratadiffgenes}{Number of genes per stratum and number of DE genes per stratum.}
#' \item{stratify.by}{The input stratify.by.}
#' \item{strata}{The input strata.}
#' \item{Nreps}{Sample sizes one wants to perform simulation on.
#' This is taken from the simulation options.}
#' \item{target.by}{The input method to define "biologically important" DE genes,
#' either by log fold change or effect size.}
#' \item{delta}{The input delta for biologically important genes.
#' If delta=0, all target.by will be considered.}
#' @details This is the main function to compute various power-related quantities,
#' using stratification and filtering.
#' \describe{
#' \item{Gene stratification}{We recommend to compute and visualize error rates (especially TPR)
#' conditional on expression characteristics like mean, dispersion and dropout rate.
#' It is likely that the power to detect DE genes is stronly dependent on
#' mean expression levels even though the magnitude of effect sizes is the same.
#' The stratified results will provide a more comprehensive power assessment and
#' better guide the investigators in experimental designs and analysis strategies.}
#' \item{Gene filtering}{Sometimes it is advisible to filter out some genes
#' (such as the ones with very low mean expression) before DE detection.
#' The filtering option here provides an opportunity to compare the powers before and after filtering.}
#' \item{Define biologically interesting genes}{We provide two options to define biologically interesting genes:
#' by absolute values of log fold changes or effect sizes
#' (absolute values of log fold changes divided by the square root of 1/(mean+dispersions)).
#' Genes with these quantities over a threshold are deemed interesting,
#' and the rate calculations are based on these genes.}
#' }
#' @author Beate Vieth
#' @seealso \code{\link{estimateParam}} for negative binomial parameters,
#' \code{\link{SimSetup}} and
#' \code{\link{DESetup}} for setting up simulation parameters and
#' \code{\link{simulateDE}} for simulating differential expression.
#' @examples
#' \dontrun{
#' ## for example DE simulation result see simulateDE
#' evalres <- evaluateSim(simRes=simres,
#' alpha.type="adjusted",
#' MTC="BH", alpha.nominal=0.1,
#' stratify.by="mean",
#' filter.by="none", target.by="lfc",
#' delta=0)
#' }
#' @rdname evaluateDE
#' @importFrom stats ecdf quantile p.adjust.methods p.adjust
#' @importFrom qvalue qvalue
#' @importFrom IHW ihw adj_pvalues
#' @export
evaluateDE <- function(simRes, alpha.type=c("adjusted","raw"),
                       MTC=c('BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni', 'Storey', 'IHW'),
                       alpha.nominal=0.1,
                       stratify.by=c("mean", "dispersion", "dropout", "lfc"),
                       filter.by=c("none", "mean", "dispersion", "dropout"),
                       strata.filtered=1,
                       target.by=c("lfc", "effectsize"), delta=0) {

  alpha.type = match.arg(alpha.type)
  MTC = match.arg(MTC)
  stratify.by = match.arg(stratify.by)
  filter.by = match.arg(filter.by)
  target.by = match.arg(target.by)

  ## some general parameters
  Nreps1 = simRes$sim.settings$n1
  Nreps2 = simRes$sim.settings$n2
  ngenes = simRes$sim.settings$ngenes
  sim.opts = simRes$sim.settings
  DEids = simRes$sim.settings$DEid
  lfcs = simRes$sim.settings$pLFC
  tlfcs = lapply(1:length(lfcs), function(i) {lfcs[[i]]})
  nsims = simRes$sim.settings$nsims
  estmeans = simRes$sim.settings$means
  estdisps = simRes$sim.settings$dispersion
  estdropout = simRes$sim.settings$p0
  mu = simRes$mu
  disp = simRes$disp
  dropout = simRes$dropout
  elfc = simRes$elfc
  DEmethod = simRes$sim.settings$DEmethod
  pvalue = simRes$pvalue
  fdr = simRes$fdr

  ## calculate strata
  tmp.ecdf.mean = stats::ecdf(log2(estmeans+1))
  tmp.quantile.mean = stats::quantile(tmp.ecdf.mean, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  strata.mean = unique(c(0,unname(tmp.quantile.mean),Inf))
  strata.mean = unique(round(strata.mean, digits=2))
  tmp.ecdf.disps = stats::ecdf(log2(estdisps))
  tmp.quantile.disps = stats::quantile(tmp.ecdf.disps, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  strata.disps = unique(c(0,unname(tmp.quantile.disps),Inf))
  strata.disps = unique(round(strata.disps, digits=2))
  tmp.ecdf.drop = stats::ecdf(estdropout)
  tmp.quantile.drop = stats::quantile(tmp.ecdf.drop, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
  strata.drop = unique(c(0,unname(tmp.quantile.drop),1))
  strata.drop = unique(round(strata.drop, digits=2))
  tmp.ecdf.lfc = stats::ecdf(unlist(tlfcs))
  tmp.quantile.lfc = stats::quantile(tmp.ecdf.lfc, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  strata.lfc = unique(c(-Inf,unname(tmp.quantile.lfc),Inf))
  strata.lfc = unique(round(strata.lfc, digits=2))

  ## initialize results
  ## determine dimension of results, for filtering
  if(stratify.by=='mean') {
    nr = length(strata.mean) - 1
  }
  if(stratify.by=='dispersion') {
    nr = length(strata.disps) - 1
  }
  if(stratify.by=='dropout') {
    nr = length(strata.drop) - 1
  }
  if(stratify.by=='lfc') {
    nr = length(strata.lfc) - 1
  }
  if(filter.by %in% c("mean", "dispersion", "dropout")) {
    nr = nr - strata.filtered
  }

  TP = TN =  FP = FN = TPR = TNR = FPR = FNR = FDR = xgrl = xgrld = array(NA,dim=c(nr,length(Nreps1), nsims))
  TP.marginal = TN.marginal = FP.marginal = FN.marginal = TPR.marginal = TNR.marginal = FPR.marginal = FNR.marginal = FDR.marginal = matrix(NA,length(Nreps1), nsims)

  ## loop over simulation and replicates
  for(i in 1:nsims) {
    for(j in seq(along=Nreps1)) {
      Nrep1 = Nreps1[j]
      Nrep2 = Nreps2[j]
      ## get DE flags.
      DEid = DEids[[i]]
      lfc = lfcs[[i]]
      Zg = Zg2 = rep(0, ngenes)
      Zg[DEid] = 1
      ## find target (interesting) genes
      if(target.by == "lfc") {
        ix = abs(lfc) > delta
      } else if (target.by == "effectsize") {
        effectsize = lfc / sqrt(1/(log2(mu[,,i])+log2(disp[,,i])))
        ix = abs(effectsize) > delta
      }
      Zg2[DEid[ix]] = 1

      ### STRATIFICATION
      ## calculate stratificaton
      # mean
      X.bar1 = mu[,j,i]
      ix.keep.mean = which(!is.na(X.bar1))
      xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
      xgrd.mean = cut(log2(X.bar1[DEid]+1), strata.mean)
      # dispersion
      X.disp1 = disp[,j,i]
      ix.keep.disps = which(!is.na(X.disp1))
      xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps)
      xgrd.disps = cut(log2(X.disp1[DEid]), strata.disps)
      # dropout
      X.drop1 = dropout[,j,i]
      ix.keep.drop = which(!is.na(X.drop1))
      xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
      xgrd.drop = cut(X.drop1[DEid], strata.drop)
      # lfc
      X.lfc1 = elfc[,j,1]
      ix.keep.lfc = which(!is.na(X.lfc1))
      xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
      xgrd.lfc = cut(X.lfc1[DEid], strata.lfc)


      ### FILTERING
      ## calculate filtering
      ## stratify by mean
      if(stratify.by == "mean") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.mean = ix.keep.mean[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean[-strata.filt.mean])
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean[-strata.filt.mean])
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.mean = ix.keep.mean[!(xgr.disps %in% lev.disps[strata.filt.disps])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.mean = ix.keep.mean[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.mean = cut(log2(X.bar1[ix.keep.mean]+1), strata.mean)
          xgrd.mean = cut(log2(X.bar1[(ix.keep.mean && DEid)]+1), strata.mean)
        }
        if(filter.by == "none") {
          ix.keep.mean = ix.keep.mean
          xgr.mean = xgr.mean
          xgrd.mean = xgrd.mean
        }
      }
      ## stratify by dispersion
      if(stratify.by == "dispersion") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.disps = ix.keep.disps[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps)
          xgrd.disps = cut(log2(X.disp1[(ix.keep.disps && DEid)]), strata.disps)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.disps = ix.keep.disps[!(xgr.disps %in% lev.disps[strata.filt.disps])]
          # recut
          xgr.disps = cut(log2(X.disp1[ix.keep.disps]), strata.disps[-strata.filt.disps])
          xgrd.disps = cut(log2(X.disp1[(ix.keep.disps && DEid)]), strata.disps[-strata.filt.disps])
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.disps = ix.keep.disps[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.disps = cut(X.disp1[ix.keep.disps], strata.disps)
          xgrd.disps = cut(X.disp1[(ix.keep.disps && DEid)], strata.disps)
        }
        if(filter.by == "none") {
          ix.keep.disps = ix.keep.disps
          xgr.disps = xgr.disps
          xgrd.disps = xgrd.disps
        }
      }
      ## stratify by dropout
      if(stratify.by == "dropout") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.drop = ix.keep.drop[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.drop = ix.keep.drop[!(xgr.disps %in% lev.mean[strata.filt.disps])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.drop = ix.keep.drop[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop[-strata.filt.drop])
          xgrd.drop = cut(X.drop1[(ix.keep.drop && DEid)], strata.drop[-strata.filt.drop])
        }
        if(filter.by == "none") {
          ix.keep.drop = ix.keep.drop
          xgr.drop = xgr.drop
          xgrd.drop = xgrd.drop
        }
      }
      ## stratify by lfc
      if(stratify.by == "lfc") {
        if(filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          ix.keep.lfc = ix.keep.lfc[!(xgr.mean %in% lev.mean[strata.filt.mean])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((max(nlevels(lev.disps))-(1-strata.filtered)):max(nlevels(lev.disps)))
          ix.keep.lfc = ix.keep.lfc[!(xgr.disps %in% lev.mean[strata.filt.disps])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((max(nlevels(lev.drop))-(1-strata.filtered)):max(nlevels(lev.drop)))
          ix.keep.lfc = ix.keep.lfc[!(xgr.drop %in% lev.drop[strata.filt.drop])]
          # recut
          xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
          xgrd.lfc = cut(X.lfc1[(ix.keep.lfc && DEid)], strata.lfc)
        }
        if(filter.by == "none") {
          ix.keep.lfc = ix.keep.lfc
          xgr.lfc = xgr.lfc
          xgr.lfc = xgr.lfc
        }
      }

      ### SET STRATIFICATION
      if(stratify.by == "mean") {
        strata = strata.mean
        xgr = xgr.mean
        xgrd = xgrd.mean
        ix.keep = ix.keep.mean
      }
      if(stratify.by == "dispersion") {
        strata = strata.disps
        xgr = xgr.disps
        xgrd = xgrd.disps
        ix.keep = ix.keep.disps
      }
      if(stratify.by == "dropout") {
        strata = strata.drop
        xgr = xgr.drop
        xgrd = xgrd.drop
        ix.keep = ix.keep.drop
      }
      if(stratify.by == "lfc") {
        strata = strata.lfc
        xgr = xgr.lfc
        xgrd = xgrd.lfc
        ix.keep = ix.keep.lfc
      }

      ## get type I error alpha (pvalue or fdr output from testing)
      if(alpha.type == "raw") {
        if(DEmethod %in% c("edgeRglm", "edgeRql", "limma-voom", "limma-trend",
                           "DESeq2", "ROTS", "DSS", "MAST", "scde", "BPSC", "scDD")) {
          x = pvalue[ix.keep,j,i]
          x[is.na(x)] = 1
        }
        if(DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod," only provides adjusted p-values."))
          x = fdr[ix.keep,j,i]
          x[is.na(x)] = 1
        }
      }
      if(alpha.type == "adjusted") {
        if(DEmethod %in% c("edgeRglm", "edgeRql", "limma-voom", "limma-trend",
                           "DESeq2", "ROTS", "DSS", "NBPSeq", "MAST", "scde", "BPSC", "monocle", "scDD")) {
          pval = pvalue[ix.keep,j,i]
          meanexpr = mu[ix.keep,j,i]
          if(MTC %in% stats::p.adjust.methods) {
            x = stats::p.adjust(pval, method = MTC)
            x[is.na(x)] = 1
          }
          if(MTC %in% "Storey") {
            x = rep(NA, length(pval))
            x[!is.na(pval)] = qvalue::qvalue(p = pval[!is.na(pval)])
            x[is.na(x)] = 1
          }
          if(MTC %in% "IHW") {
            in.dat = data.frame(pvalue = pval, meanexpr = meanexpr)
            tmp = IHW::ihw(pvalue ~ meanexpr, data = in.dat, alpha = alpha.nominal)
            x = IHW::adj_pvalues(tmp)
            x[is.na(x)] = 1
          }
        }
        if(DEmethod %in% c("baySeq", "NOISeq", "EBSeq", "TSPM")) {
          message(paste0("The DE method ", DEmethod," only provides adjusted p-values."))
          x = fdr[ix.keep,j,i]
          x[is.na(x)] = 1
        }
      }

      ## update Zg flags after filtering
      Zg = Zg[ix.keep]
      Zg2 = Zg2[ix.keep]

      #  number of strata genes and diff strata genes in output table
      xgrl[,j,i] = table(xgr)
      xgrld[,j,i] = table(xgrd)

      ## calculate stratified power-related quantities
      error.mat = .error.matrix(p=x, p.crit=alpha.nominal, Zg=Zg, Zg2=Zg2, xgr=xgr)

      TP[,j,i] = error.mat$TP
      TN[,j,i] = error.mat$TN
      FP[,j,i] = error.mat$FP
      FN[,j,i] = error.mat$FN

      TP.marginal[j,i] = error.mat$TP.marginal
      TN.marginal[j,i] = error.mat$TN.marginal
      FP.marginal[j,i] = error.mat$FP.marginal
      FN.marginal[j,i] = error.mat$FN.marginal

      TPR[,j,i] = error.mat$TPR
      TNR[,j,i] = error.mat$TNR
      FPR[,j,i] = error.mat$FPR
      FNR[,j,i] = error.mat$FNR
      FDR[,j,i] = error.mat$FDR

      TPR.marginal[j,i] = error.mat$TPR.marginal
      TNR.marginal[j,i] = error.mat$TNR.marginal
      FPR.marginal[j,i] = error.mat$FPR.marginal
      FNR.marginal[j,i] = error.mat$FNR.marginal
      FDR.marginal[j,i] = error.mat$FDR.marginal

    }
  }

  output <- list(stratagenes=xgrl, stratadiffgenes=xgrld,
                 TN=TN, TP=TP, FP=FP, FN=FN,
                 TN.marginal=TN.marginal, TP.marginal=TP.marginal, FP.marginal=FP.marginal, FN.marginal=FN.marginal,
                 TNR=TNR, TPR=TPR, FPR=FPR, FNR=FNR,FDR=FDR,
                 TNR.marginal=TNR.marginal, TPR.marginal=TPR.marginal,
                 FPR.marginal=FPR.marginal, FNR.marginal=FNR.marginal,
                 FDR.marginal=FDR.marginal,
                 ## below are input parameters:
                 alpha.type=alpha.type, MTC=ifelse(alpha.type=="adjusted", MTC, "not applicable"), alpha.nominal=alpha.nominal,
                 stratify.by=stratify.by, strata=strata, strata.levels=levels(xgr),
                 target.by=target.by, n1=Nreps1, n2=Nreps2, delta=delta)

  return(output)
}






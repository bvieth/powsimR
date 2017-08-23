#' @name evaluateDE
#' @aliases evaluateDE
#' @title Compute the confusion matrix-related quantities from simulation results
#' @description This function takes the simulation output from \code{\link{simulateDE}}
#' and computes quantities of the confusion matrix of classification testing
#' @usage evaluateDE(simRes, alpha.type=c("adjusted","raw"),
#' MTC=c('BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni', 'Storey', 'IHW'),
#' alpha.nominal=0.1,
#' stratify.by=c("mean", "dispersion", "dropout", "lfc"),
#' filter.by=c("none", "mean", "dispersion", "dropout"),
#' strata.filtered=1,
#' target.by=c("lfc", "effectsize"), delta=0)
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
#' Available options are \code{"mean"}, \code{"dispersion"} and \code{"dropout"},
#' for stratifying genes by average expression levels, dispersion or dropout rates.
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
#' \item{TN, TP, FP, FN, TNR, TPR, FPR, FNR, FDR}
#' {3D array representing the number of true negatives, true positives, false positives,
#' false negatives and their proportions/rates as well as false discovery rate
#'  for all simulation settings. The dimension of the arrays are nstrata * N * nsims.
#'   Here nstrata is number of specified strata.
#'   N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TN.marginal, TP.marginal,FP.marginal, FN.marginal}{Matrix representing the number of true negatives, true positives, false positives,
#'  false negatives for all simulation settings.
#'  The dimension of the matrices are N * nsims.
#'  Here N is number of different sample sizes settings, and nsims is number of simulations.}
#' \item{TNR.marginal, TPR.marginal, FPR.marginal, FNR.marginal, FDR.marginal}
#' {Matrix representing the marginal rates for all simulation settings. The dimension of the matrices are N * nsims.}
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
#' and the power calculation are based on these genes.}
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
  dlfcs = lapply(1:length(lfcs), function(i) {lfcs[[i]][DEids[[i]]]})
  nsims = simRes$sim.settings$nsims
  estmeans = simRes$sim.settings$means
  estdisps = simRes$sim.settings$dispersion
  estdropout = simRes$sim.settings$p0
  mu = simRes$mu
  disp = simRes$disp
  dropout = simRes$dropout
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
  tmp.ecdf.lfc = stats::ecdf(unlist(dlfcs))
  tmp.quantile.lfc = stats::quantile(tmp.ecdf.lfc, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
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
      X.lfc1 = lfcs[[i]]
      ix.keep.lfc = which(!X.lfc1==0)
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
            x[!is.na(x)] = qvalue::qvalue(p = x[!is.na(x)])
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



# PROPORTIONS AND RATES ---------------------------------------------------

## compute the proportions and rates of the confusion/error matrix
## containing classification test results (marginal and per stratum)
## TP, FP, TN, FN
## TPR, FPR, TNR, FNR, FDR

.error.matrix <- function(p, p.crit, Zg, Zg2, xgr){
  ## p is input raw p-value or  q-value.
  ## p.crit is cutoff for significance (nominal alpha level)
  ## Zg is the indicator for genes with lfc added
  ## Zg2 is the indicator for genes with lfc added and with "meaningful" effect size/ above certain delta
  ## xgr is stratum

  ##  R (the number of rejected nulll hypothesis)
  ix.R = p <= p.crit # genes that are called differentially expressed according to p-value of test and chosen cutoff
  R = sum(ix.R)  # total number of null hypothesis rejected
  R.stratified = tapply(ix.R, xgr, sum) # total number of null hypothesis rejected per stratum
  #
  ## G (the number of null hypothesis accepted)
  ix.G = p > p.crit # genes that are nondifferentially expressed according to p-value of test and chosen cutoff
  G = sum(ix.G)  # total number of null hypothesis accepted
  G.stratified = tapply(ix.G, xgr, sum) # total number of null hypothesis accepted per stratum


  ##  TP: condition is positive, H0 is rejected
  id.TP = Zg2==1
  TP = tapply(p[id.TP] <= p.crit, xgr[id.TP], sum)
  TP.mar = sum(p[id.TP] <= p.crit)
  ##  TN: condition is negative, H0 is accepted
  id.TN = Zg==0
  TN = tapply(p[id.TN] > p.crit, xgr[id.TN], sum)
  TN.mar = sum(p[id.TN] > p.crit)
  ## FN: condition is positive, H0 is accepted
  id.FN = Zg2==1
  FN = tapply(p[id.FN] > p.crit, xgr[id.FN], sum)
  FN.mar = sum(p[id.FN] > p.crit)
  ## FP: condition is negative, H0 is rejected
  id.FP = Zg==0
  FP = tapply(p[id.FP] <= p.crit, xgr[id.FP], sum)
  FP.mar = sum(p[id.FP] <= p.crit)

  ## false discovery rate (FP/(TP+FP))
  FDR = FP / R.stratified
  FDR.marginal = FP.mar / R
  ## type I error rate / FPR / alpha / fall-out (FP/(FP+TN))
  FPR = as.vector(FP/table(xgr[id.FP]))
  FPR.marginal = FP.mar / sum(id.FP)
  ## True positive rate / sensitivity / power (TP/(TP+FN))
  TPR = as.vector(TP/table(xgr[id.TP]))
  TPR.marginal =  TP.mar / sum(id.TP)
  ## True negative rate / specificity (TN/(TN+FP))
  TNR = as.vector(TN/table(xgr[id.TN]))
  TNR.marginal = TN.mar / sum(id.TN)
  ## False negative rate / miss rate (FN/(FN+TP))
  FNR = as.vector(FN/table(xgr[id.FN]))
  FNR.marginal =  FN.mar / sum(id.FN)


  # output
  list(TN=TN,TN.marginal=TN.mar, TP=TP, TP.marginal=TP.mar, FP=FP, FP.marginal=FP.mar, FN=FN, FN.marginal=FN.mar,
       TNR=TNR, TPR=TPR, FPR=FPR, FNR=FNR,FDR=FDR,
       TNR.marginal=TNR.marginal, TPR.marginal=TPR.marginal,
       FPR.marginal=FPR.marginal, FNR.marginal=FNR.marginal,
       FDR.marginal=FDR.marginal)
}



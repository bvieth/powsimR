# Taken from https://github.com/Vivianstats/scImpute

# find_hv_genes
#' @importFrom stats quantile
.find_hv_genes <- function(count, I, J){
  count_nzero = lapply(1:I, function(i) setdiff(count[i, ], log10(1.01)))
  mu = sapply(count_nzero, mean)
  mu[is.na(mu)] = 0
  sd = sapply(count_nzero, sd)
  sd[is.na(sd)] = 0
  cv = sd/mu
  cv[is.na(cv)] = 0
  # sum(mu >= 1 & cv >= quantile(cv, 0.25), na.rm = TRUE)
  high_var_genes = which(mu >= 1 & cv >= stats::quantile(cv, 0.25))
  if(length(high_var_genes) < 500){
    high_var_genes = 1:I}
  count_hv = count[high_var_genes, ]
  return(count_hv)
}

# find_neighbors
#' @importFrom stats prcomp quantile
#' @importFrom kernlab specc
#' @importFrom parallel mclapply
.find_neighbors <- function(count_hv,
                            labeled,
                            J,
                            Kcluster = NULL,
                            ncores,
                            cell_labels = NULL){
  ## dimeansion reduction
  pca = stats::prcomp(t(count_hv))
  eigs = (pca$sdev)^2
  var_cum = cumsum(eigs)/sum(eigs)
  npc = which.max(var_cum > 0.4)
  if (labeled == FALSE){ npc = max(npc, Kcluster) }
  if (npc < 3){ npc = 3 }
  mat_pcs = t(pca$x[, 1:npc]) # columns are cells

  ## detect outliers
  dist_cells_list = parallel::mclapply(1:J, function(id1){
    sapply(1:J, function(id2){
      if(id1 <= id2) return(0)
      sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
      sqrt(sse)
    })
  }, mc.cores = ncores)
  dist_cells = matrix(0, nrow = J, ncol = J)
  for(cellid in 1:J){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
  dist_cells = dist_cells + t(dist_cells)

  if (labeled == FALSE){
    min_dist = sapply(1:J, function(i){
      min(dist_cells[i, -i])
    })
    iqr = stats::quantile(min_dist, 0.75) - stats::quantile(min_dist, 0.25)
    outliers = which(min_dist > 1.5 * iqr + stats::quantile(min_dist, 0.75))

    ## clustering
    non_out = setdiff(1:J, outliers)
    spec_res = kernlab::specc(t(mat_pcs[, non_out]), centers = Kcluster, kernel = "rbfdot")
    nbs = rep(NA, J)
    nbs[non_out] = spec_res
    return(list(dist_cells = dist_cells, clust = nbs))
  }

  if(labeled == TRUE){
    if(class(cell_labels) == "character"){
      labels_uniq = unique(cell_labels)
      labels_mth = 1:length(labels_uniq)
      names(labels_mth) = labels_uniq
      clust = labels_mth[cell_labels]
    }else{
      clust = cell_labels
    }
    return(list(dist_cells = dist_cells, clust = clust))
  }
}

### root-finding equation
.fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

### update parameters in gamma distribution
#' @importFrom stats uniroot
.update_gmm_pars = function(x, wt){
  tp_s = sum(wt)
  tp_t = sum(wt * x)
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = stats::uniroot(.fn, c(0.9, 1.1) * alpha0,
                             target = tp_v, extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}

### estimate parameters in the mixture distribution
.get_mix <- function(xdata, point){
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0

  while(eps > 0.5) {
    wt = .calculate_weight(xdata, paramt)
    paramt[1] = sum(wt[, 1])/nrow(wt)
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
    paramt[2:3] = .update_gmm_pars(x=xdata, wt=wt[,1])

    loglik = sum(log10(.dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100)
      break
  }
  return(paramt)
}

#' @importFrom stats dgamma dnorm
.dmix <- function (x, pars) {
    pars[1] * stats::dgamma(x, shape = pars[2],
                     rate = pars[3]) + (1 - pars[1]) * stats::dnorm(x, mean = pars[4], sd = pars[5])
}

#' @importFrom parallel mclapply
.get_mix_parameters <- function (count,
                                 point = log10(1.01),
                                 ncores = 8) {
    count = as.matrix(count)
    null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
    parslist = parallel::mclapply(1:nrow(count), function(ii) {
      if (ii %% 2000 == 0) {
        gc()
        print(ii)
      }
      if (ii %in% null_genes) {
        return(rep(NA, 5))
      }
      xdata = count[ii, ]
      paramt = try(.get_mix(xdata, point), silent = TRUE)
      if (class(paramt) == "try-error"){
        paramt = rep(NA, 5)
      }
      return(paramt)
    }, mc.cores = ncores)
    parslist = Reduce(rbind, parslist)
    colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
    return(parslist)
  }


# find_va_genes
#' @importFrom stats complete.cases dgamma dnorm
.find_va_genes = function(parslist, subcount){
  point = log10(1.01)
  valid_genes = which( (rowSums(subcount) > point * ncol(subcount)) &
                         stats::complete.cases(parslist) )
  if(length(valid_genes) == 0) return(valid_genes)
  # find out genes that violate assumption
  mu = parslist[, "mu"]
  sgene1 = which(mu <= log10(1+1.01))
  # sgene2 = which(mu <= log10(10+1.01) & mu - parslist[,5] > log10(1.01))

  dcheck1 = stats::dgamma(mu+1, shape = parslist[, "alpha"], rate = parslist[, "beta"])
  dcheck2 = stats::dnorm(mu+1, mean = parslist[, "mu"], sd = parslist[, "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)
  return(valid_genes)
}


# calculate_weight
#' @importFrom stats dgamma dnorm
.calculate_weight <- function (x, paramt){
    pz1 = paramt[1] * stats::dgamma(x, shape = paramt[2], rate = paramt[3])
    pz2 = (1 - paramt[1]) * stats::dnorm(x, mean = paramt[4], sd = paramt[5])
    pz = pz1/(pz1 + pz2)
    pz[pz1 == 0] = 0
    return(cbind(pz, 1 - pz))
}

# impute_nnls
#' @importFrom penalized penalized predict
.impute_nnls <- function(Ic,
                         cellid,
                         subcount,
                         droprate,
                         geneid_drop,
                         geneid_obs,
                         nbs,
                         distc){
  yobs = subcount[ ,cellid]
  if (length(geneid_drop) == 0 | length(geneid_drop) == Ic) {
    return(yobs) }
  yimpute = rep(0, Ic)

  xx = subcount[geneid_obs, nbs]
  yy = subcount[geneid_obs, cellid]
  ximpute = subcount[geneid_drop, nbs]
  num_thre = 1000
  if(ncol(xx) >= min(num_thre, nrow(xx))){
    if (num_thre >= nrow(xx)){
      new_thre = round((2*nrow(xx)/3))
    }else{ new_thre = num_thre}
    filterid = order(distc[cellid, -cellid])[1: new_thre]
    xx = xx[, filterid, drop = FALSE]
    ximpute = ximpute[, filterid, drop = FALSE]
  }
  set.seed(cellid)
  nnls = penalized::penalized(yy, penalized = xx, unpenalized = ~0,
                   positive = TRUE, lambda1 = 0, lambda2 = 0,
                   maxiter = 3000, trace = FALSE)
  ynew = penalized::predict(nnls, penalized = ximpute, unpenalized = ~0)[,1]

  yimpute[geneid_drop] = ynew
  yimpute[geneid_obs] = yobs[geneid_obs]
  return(yimpute)
}

#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
.imputation_model8 = function(count,
                             labeled,
                             point,
                             drop_thre = 0.5,
                             Kcluster = 10,
                             ncores){
  count = as.matrix(count)
  I = nrow(count)
  J = ncol(count)
  count_imp = count

  if(Kcluster == 1){
    clust = rep(1, J)
  }else{
    # find highly variable genes
    count_hv = .find_hv_genes(count, I, J)

    set.seed(Kcluster)
    neighbors_res = .find_neighbors(count_hv = count_hv,
                                   labeled = FALSE,
                                   J = J,
                                   Kcluster = Kcluster,
                                   ncores = ncores)
    dist_cells = neighbors_res$dist_cells
    clust = neighbors_res$clust
  }

  # mixture model
  nclust = sum(!is.na(unique(clust)))
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  for(cc in 1:nclust){
    params <- .get_mix_parameters(count = count[, which(clust == cc), drop = FALSE],
                                  point = log10(1.01),
                                  ncores = ncores)

    cells = which(clust == cc)
    if(length(cells) <= 1) { next }
    parslist = params
    valid_genes = .find_va_genes(parslist, subcount = count[, cells])
    if(length(valid_genes) <= 10){ next }

    subcount = count[valid_genes, cells, drop = FALSE]
    Ic = length(valid_genes)
    Jc = ncol(subcount)
    parslist = parslist[valid_genes, ]

    droprate = t(sapply(1:Ic, function(i) {
      wt = .calculate_weight(subcount[i, ], parslist[i, ])
      return(wt[, 1])
    }))
    mucheck = sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
    droprate[mucheck & droprate > drop_thre] = 0
    # dropouts
    setA = lapply(1:Jc, function(cellid){
      which(droprate[, cellid] > drop_thre)
    })
    # non-dropouts
    setB = lapply(1:Jc, function(cellid){
      which(droprate[, cellid] <= drop_thre)
    })
    # imputation
    # registerDoSNOW(cl)
    subres = foreach::foreach(cellid = 1:Jc, .packages = c("penalized"),
                     .combine = cbind, .export = c(".impute_nnls")) %dopar% {
                       if (cellid %% 100 == 0) {gc()}
                       nbs = setdiff(1:Jc, cellid)
                       if (length(nbs) == 0) {return(NULL)}
                       geneid_drop = setA[[cellid]]
                       geneid_obs = setB[[cellid]]
                       y = try(.impute_nnls(Ic,
                                            cellid,
                                            subcount,
                                            droprate,
                                            geneid_drop,
                                            geneid_obs,
                                            nbs,
                                            distc = dist_cells[cells, cells]),
                               silent = TRUE)
                       if (class(y) == "try-error") {
                         # print(y)
                         y = subcount[, cellid, drop = FALSE]
                       }
                       return(y)
                     }

    count_imp[valid_genes, cells] = subres
  }
  parallel::stopCluster(cl)
  outlier = which(is.na(clust))
  count_imp[count_imp < point] = point
  return(list(count_imp = count_imp, outlier = outlier))
}


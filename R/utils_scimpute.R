
# MIXTURE MODEL PARAMETERS ------------------------------------------------

# root-finding equation
.fn <- function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

# update parameters in gamma distribution
#' @importFrom stats uniroot
.update_gmm_pars <- function(x, wt){
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
      alpha = stats::uniroot(.fn, c(0.9, 1.1) * alpha0, target = tp_v,
                      extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}

# estimate parameters in the mixture distribution
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

#' @importFrom parallel mclapply
.get_mix_parameters <- function(count, point = log10(1.01), ncores) {
  count = as.matrix(count)
  null_genes = which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
  parslist = parallel::mclapply(1:nrow(count), function(ii) {
    if (ii%%500 == 0) {
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

#' @importFrom stats dgamma dnorm
.calculate_weight <- function (x, paramt) {
  pz1 = paramt[1] * stats::dgamma(x, shape = paramt[2], rate = paramt[3])
  pz2 = (1 - paramt[1]) * stats::dnorm(x, mean = paramt[4], sd = paramt[5])
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}

#' @importFrom stats dgamma dnorm
.dmix <- function (x, pars) {
  pars[1] * dgamma(x, shape = pars[2], rate = pars[3]) + (1 - pars[1]) * dnorm(x, mean = pars[4], sd = pars[5])
}

# IMPUTATION --------------------------------------------------------------

#' @importFrom stats complete.cases
#' @importFrom glmnet cv.glmnet predict.cv.glmnet
.imputation_model1 <- function(count, point, parslist, drop_thre = 0.5, method = 2, ncores) {
  count = as.matrix(count)
  count_imp = count
  valid_genes = which( (rowSums(count) > point * ncol(count)) &
                         stats::complete.cases(parslist) )
  count = count[valid_genes, , drop = FALSE]
  parslist = parslist[valid_genes, , drop = FALSE]
  I = nrow(count)
  J = ncol(count)
  droprate = t(sapply(1:I, function(i) {
    wt = .calculate_weight(count[i, ], parslist[i, ])
    return(wt[, 1])
  }))
  impute = function(I, cellid, count, parslist, droprate, method) {
    yobs = count[, cellid]
    yimpute = rep(0, I)
    geneid_drop = which(droprate[, cellid] >= drop_thre)
    if (length(geneid_drop) == 0 | length(geneid_drop) == I) {
      yimpute = yobs
    }
    else {
      geneid_obs = setdiff(1:I, geneid_drop)
      xx = count[geneid_obs, -cellid]
      yy = count[geneid_obs, cellid]
      weight = 1 - parslist[geneid_obs, "rate"]
      set.seed(cellid)
      model_cv = glmnet::cv.glmnet(xx, yy, alpha = 1, weights = weight,
                                   nfolds = 5)
      vars_nonzero = glmnet::predict.cv.glmnet(model_cv, type = "nonzero",
                                     s = "lambda.1se")
      if (class(vars_nonzero) == "list") {
        return(yobs)
      }
      else {
        vars_nonzero = vars_nonzero[, 1]
        xselect = xx[, vars_nonzero, drop = FALSE]
        model_ols = lm(yy ~ xselect, weights = weight)
        ximpute = count[geneid_drop, -cellid]
        ximpute = ximpute[, vars_nonzero, drop = FALSE]
        xnew = cbind(rep(1, length(geneid_drop)), ximpute)
        ynew = xnew %*% matrix(coef(model_ols), ncol = 1)
        if (method == 2) {
          check = apply(ximpute, 1, function(x) {
            sum(x == log10(1.01))
          })
          ynew[check == ncol(ximpute)] = point
        }
        yimpute[geneid_drop] = ynew
        yimpute[geneid_obs] = yobs[geneid_obs]
        return(yimpute)
      }
    }
  }
  res = mclapply(1:J, function(cellid) {
    if (cellid %% 50 == 0)
      print(cellid)
    y = try(impute(I, cellid, count, parslist, droprate,
                   method), silent = TRUE)
    if (class(y) == "try-error") {
      print(y)
      y = count[, cellid]
    }
    return(y)
  }, mc.cores = ncores)
  res = Reduce(cbind, res)
  count_imp[valid_genes, ] = res
  count_imp[count_imp < point] = point
  return(count_imp)
}



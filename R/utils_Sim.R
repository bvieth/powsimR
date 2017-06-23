
# Simulate DE between 2 groups (DE of mean) -------------------------------

#' @importFrom stats model.matrix
.simRNAseq.2grp <- function(simOptions, n1, n2) {

  ## make group labels
  design <- c(rep(-1, n1), rep(1, n2))
  ## make model matrix
  mod = stats::model.matrix(~-1 + design)

  # generate read counts
  if (simOptions$RNAseq == "singlecell") {
    if (attr(sim.options, 'Distribution') == 'NB') {
      simcounts = .sc.NB.RNAseq_counts(sim.options = simOptions, modelmatrix = mod)
    }
    if (attr(sim.options, 'Distribution') == 'ZINB') {
      simcounts = .sc.ZINB.RNAseq_counts(sim.options = simOptions, modelmatrix = mod)
    }
  }
  if (simOptions$RNAseq == "bulk") {
    simcounts = .bulk.NB.RNAseq_counts(sim.options = simOptions, modelmatrix = mod)
  }

  simcounts <- apply(simcounts,2,function(x) {storage.mode(x) <- 'integer'; x})
  colnames(simcounts) <- paste0("S", 1:ncol(simcounts))
  rownames(simcounts) <- paste0("G", 1:nrow(simcounts))

  ## return
  list(counts = simcounts, designs = design, simOptions = simOptions)
}

#' @importFrom stats rnorm rnbinom
.sc.NB.RNAseq_counts <- function(sim.options, modelmatrix) {

  set.seed(sim.options$sim.seed)
  print(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = cbind(sim.options$lfc)

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)


    # capture efficiency (size factor)
    if (sim.options$size.factors == "equal") {
      all.facs <- rep(1, nsamples)
    }
    if (sim.options$size.factors == "given") {
      if (is.function(sim.options$sf)) {
        all.facs <- sim.options$sf(nsamples)
      }
      if (is.vector(sim.options$sf)) {
        all.facs <- sample(sim.options$sf, nsamples, replace = TRUE)
      }
      if (is.null(sim.options$sf)) {
        all.facs <- rep(1, nsamples)
      }
    }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = 0

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec), ncol = nsamples, nrow = ngenes)
  }

  if (attr(sim.options, 'param.type') == 'insilico') {

    # sample from the mean function provided
    true.means = sim.options$means(ngenes)

    # associate with dispersion parameter
    if (length(sim.options$dispersion) == 1) {
      disp = sim.options$dispersion
    }
    if (is.function(sim.options$dispersion)) {
      disp = sim.options$dispersion(true.means)
    }

    # capture efficiency (size factor)
    if (sim.options$size.factors == "equal") {
      all.facs <- rep(1, nsamples)
    }
    if (sim.options$size.factors == "given") {
      if (is.function(sim.options$sf)) {
        all.facs <- sim.options$sf(nsamples)
      }
      if (is.vector(sim.options$sf)) {
        all.facs <- sample(sim.options$sf, nsamples, replace = TRUE)
      }
      if (is.null(sim.options$sf)) {
        all.facs <- rep(1, nsamples)
      }
    }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = 0

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 1/disp), ncol = nsamples, nrow = ngenes)
  }

  return(counts)
}

#' @importFrom stats rnorm rnbinom runif approx
.sc.ZINB.RNAseq_counts <- function(sim.options, modelmatrix) {

  set.seed(sim.options$sim.seed)
  print(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = cbind(sim.options$lfc)

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # capture efficiency (size factor)
    if (sim.options$size.factors == "equal") {
      all.facs <- rep(1, nsamples)
    }
    if (sim.options$size.factors == "given") {
      if (is.function(sim.options$sf)) {
        all.facs <- sim.options$sf(nsamples)
      }
      if (is.vector(sim.options$sf)) {
        all.facs <- sample(sim.options$sf, nsamples, replace = TRUE)
      }
      if (is.null(sim.options$sf)) {
        all.facs <- rep(1, nsamples)
      }
    }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = 0
    muvec = as.vector(mumat)

    meanp0fit.nonamplified = sim.options$meanp0fit.nonamplified
    meanp0fit.amplified = sim.options$meanp0fit.amplified
    nonamplified = sim.options$nonamplified

    if(!is.null(meanp0fit.nonamplified) && !is.null(meanp0fit.amplified)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(muvec), prob = nonamplified, size = 1)
      muvec.dat.nonamplified = abs(muvec)[p0.split==1]
      muvec.dat.amplified = abs(muvec)[p0.split==0]

      predp0.amplified.mean = approx(meanp0fit.amplified$x,
                                     meanp0fit.amplified$y,
                                     xout=muvec.dat.amplified,
                                     rule = 2:1)$y
      predp0.amplified.mean[is.na(predp0.amplified.mean)] = 0
      predp0.amplified.sd = approx(meanp0fit.amplified$x,
                                   meanp0fit.amplified$sd,
                                   xout=muvec.dat.amplified,
                                   rule=2:1)$y
      predp0.amplified.sd[is.na(predp0.amplified.sd)] = min(meanp0fit.amplified$sd)/2
      predp0.nonamplified.mean = approx(meanp0fit.nonamplified$x,
                                        meanp0fit.nonamplified$y,
                                        xout=muvec.dat.nonamplified,
                                        rule=2:1)$y
      predp0.nonamplified.mean[is.na(predp0.nonamplified.mean)] = 0
      predp0.nonamplified.sd = approx(meanp0fit.nonamplified$x,
                                      meanp0fit.nonamplified$sd,
                                      xout=muvec.dat.nonamplified,
                                      rule=2:1)$y
      predp0.nonamplified.sd[is.na(predp0.nonamplified.sd)] = min(meanp0fit.amplified$sd)/2

      p0vec = rep(NA, length(muvec))
      p0vec[p0.split==0] = rnorm(n=length(which(p0.split==0)),
                                 mean=predp0.amplified.mean,
                                 sd=predp0.amplified.sd)
      p0vec[p0.split==1] = rnorm(n=length(which(p0.split==1)),
                                 mean=predp0.nonamplified.mean,
                                 sd=predp0.nonamplified.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply (p0vec, function(x) {
        rbinom(length(x), prob = (1 - x), size = 1)
      })

      p0mat = matrix(zero.mark, nrow = ngenes)

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes) * p0mat
      counts = counts.nb * p0mat
    }

    if(any(is.null(meanp0fit.nonamplified), is.null(meanp0fit.amplified)) && !is.null(meanp0fit)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(muvec), prob = nonamplified, size = 1)
      muvec.dat = abs(muvec)

      predp0.mean = approx(meanp0fit$x,
                           meanp0fit$y,
                           xout=muvec.dat,
                           rule = 2:1)$y
      predp0.mean[is.na(predp0.mean)] = 0
      predp0.sd = approx(meanp0fit$x,
                         meanp0fit$sd,
                         xout=muvec.dat,
                         rule=2:1)$y
      predp0.sd[is.na(predp0.sd)] = min(meanp0fit$sd)/2


      p0vec = rnorm(n=length(muvec.dat),
                    mean=predp0.mean,
                    sd=predp0.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply (p0vec, function(x) {
        rbinom(length(x), prob = (1 - x), size = 1)
      })

      p0mat = matrix(zero.mark, nrow = ngenes)

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes) * p0mat
      counts = counts.nb * p0mat
    }

    if(is.null(meanp0fit.nonamplified) && is.null(meanp0fit.amplified) && is.null(meanp0fit)) {
      # result count matrix
      counts = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
    }
  }

  return(counts)
}

#' @importFrom stats rnorm rbinom rnbinom approx
.bulk.NB.RNAseq_counts <- function(sim.options, modelmatrix) {

  set.seed(sim.options$sim.seed)
  print(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = cbind(sim.options$lfc)

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with mean values
    lmu = log2(true.means + 1)
    predsize.mean = stats::approx(meansizefit$x, meansizefit$y, xout = lmu)$y
    predsize.sd = stats::approx(meansizefit$x, meansizefit$sd, xout = lmu)$y
    sizevec = stats::rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # associate dropout probability with mean parameter
    p0vec = rep(0, ngenes)
    p0s = ifelse(lmu < sim.options$p0.cut, sample(x = sim.options$p0,size = 1), p0vec)
    p0mat = t(sapply(p0s, function(x) {
      stats::rbinom(nsamples, prob = 1 - x, size = 1)
    }, simplify = T))

    # capture efficiency (size factor)
    if (sim.options$size.factors == "equal") {
      all.facs <- rep(1, nsamples)
    }
    if (sim.options$size.factors == "given") {
      if (is.function(sim.options$sf)) {
        all.facs <- sim.options$sf(nsamples)
      }
      if (is.vector(sim.options$sf)) {
        all.facs <- sample(sim.options$sf, nsamples, replace = TRUE)
      }
      if (is.null(sim.options$sf)) {
        all.facs <- rep(1, nsamples)
      }
    }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = 0

    # result count matrix
    cnts = matrix(stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec), ncol = nsamples, nrow = ngenes)
    counts = p0mat * cnts
  }

  if (attr(sim.options, 'param.type') == 'insilico') {

    # sample from the mean function provided
    true.means = sim.options$means(ngenes)

    # sample dropout parameter
    if (is.null(sim.options$p0)) {
      p0s <- rep(0, ngenes)
    }
    if (is.function(sim.options$p0)) {
      p0s <- sim.options$p0(ngenes)
    }
    p0mat = t(sapply(p0s, function(x) {
      stats::rbinom(nsamples, prob = 1 - x, size = 1)
    }, simplify = T))

    # associate with dispersion parameter
    if (length(sim.options$dispersion) == 1) {
      disp = sim.options$dispersion
    }
    if (is.function(sim.options$dispersion)) {
      disp = sim.options$dispersion(true.means)
    }

    # capture efficiency (size factor)
    if (sim.options$size.factors == "equal") {
      all.facs <- rep(1, nsamples)
    }
    if (sim.options$size.factors == "given") {
      if (is.function(sim.options$sf)) {
        all.facs <- sim.options$sf(nsamples)
      }
      if (is.vector(sim.options$sf)) {
        all.facs <- sample(sim.options$sf, nsamples, replace = TRUE)
      }
      if (is.null(sim.options$sf)) {
        all.facs <- rep(1, nsamples)
      }
    }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = 0

    # result count matrix
    cnts = matrix(stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 1/disp), ncol = nsamples, nrow = ngenes)
    counts = p0mat * cnts
  }

  return(counts)
}

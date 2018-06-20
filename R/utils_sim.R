# TODO
# change behaviour of minus mean values
# change the poisson filling?


# Simulate DE between 2 groups (DE of mean) -------------------------------

#' @importFrom stats model.matrix coef poisson sd
.simRNAseq.2grp <- function(simOptions, n1, n2, verbose) {

  if(is.null(simOptions$bLFC)) {
    ## make group labels for phenotype LFC
    phenotype <- c(rep(-1, n1), rep(1, n2))
    batch <- NULL
    ## make model matrix
    mod = stats::model.matrix(~-1 + phenotype)
    coeffs = cbind(simOptions$pLFC)
  }

  if(!is.null(simOptions$bLFC)) {
    if(simOptions$bPattern=="uncorrelated") {
      # make group labels for phenotype LFC
      phenotype <- c(rep(-1, n1), rep(1, n2))
      # make batch labels for batch LFC
      batch <- rep_len(c(-1,1), n1+n2)
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
    if(simOptions$bPattern=="orthogonal") {
      # make group labels for phenotype LFC
      phenotype <- c(rep(-1, n1), rep(1, n2))
      # make batch labels for batch LFC
      batch <-  1 - 2*rbinom(n1+n2, size=1, prob=0.5)
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
    if(simOptions$bPattern=="correlated") {
      # make group labels for phenotype LFC
      phenotype <- c(rep(-1, n1), rep(1, n2))
      # make batch labels for batch LFC
      flip <- rbinom(n1+n2, size = 1, prob = 0.9)
      batch <- phenotype*flip + -phenotype*(1-flip)
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
  }

  # generate read counts
  if (simOptions$RNAseq == "singlecell") {
    # take mean dispersion relationship as is
    if (!isTRUE(simOptions$downsample)) {
      if (verbose) {message(paste0("Predict dispersion based on true mean values.")) }
      if (attr(simOptions, 'Distribution') == 'NB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.NB.SampleWReplace_counts(sim.options = simOptions,
                                                  phenotype = phenotype,
                                                  modelmatrix = mod,
                                                  coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.NB.FillIn_counts(sim.options = simOptions,
                                          phenotype = phenotype,
                                          modelmatrix = mod,
                                          coef.dat=coeffs)
        }
      }
      if (attr(simOptions, 'Distribution') == 'ZINB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.ZINB.SampleWReplace_counts(sim.options = simOptions,
                                                    phenotype = phenotype,
                                                    modelmatrix = mod,
                                                    coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.ZINB.FillIn_counts(sim.options = simOptions,
                                            phenotype = phenotype,
                                            modelmatrix = mod,
                                            coef.dat=coeffs)
        }
      }
    }

    # downsample mean and then draw dispersion
    if (isTRUE(simOptions$downsample)) {
      if (verbose) {message(paste0("Predict dispersion based on effective mean values.")) }
      if (attr(simOptions, 'Distribution') == 'NB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.NB.SampleWReplace_downsample_counts(sim.options = simOptions,
                                                             phenotype = phenotype,
                                                             modelmatrix = mod,
                                                             coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.NB.FillIn_downsample_counts(sim.options = simOptions,
                                                     phenotype = phenotype,
                                                     modelmatrix = mod,
                                                     coef.dat=coeffs)
        }
      }
      if (attr(simOptions, 'Distribution') == 'ZINB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.ZINB.SampleWReplace_downsample_counts(sim.options = simOptions,
                                                               phenotype = phenotype,
                                                               modelmatrix = mod,
                                                               coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.ZINB.FillIn_downsample_counts(sim.options = simOptions,
                                                       phenotype = phenotype,
                                                       modelmatrix = mod,
                                                       coef.dat=coeffs)
        }
      }
    }
  }

  if (simOptions$RNAseq == "bulk") {
    sim.data = .bulk.NB.RNAseq_counts(sim.options = simOptions,
                                      phenotype = phenotype,
                                      modelmatrix = mod,
                                      coef.dat=coeffs)
  }

  sim.counts = sim.data$counts
  sf.value = sim.data$sf
  mus = sim.data$mus
  disps = sim.data$disps
  drops = sim.data$drops

  sim.counts <- apply(sim.counts,2,function(x) {storage.mode(x) <- 'integer'; x})

  # fill in gene pseudonames if missing (only true for in silico)
  if (is.null(rownames(sim.counts))) {
    rownames(sim.counts) <- paste0("G", 1:nrow(sim.counts))
    names(mus) <- names(disps) <- names(drops) <- rownames(sim.counts)
  }
  # add human readable sample names
  if (is.null(colnames(sim.counts))) {
    if(is.null(batch)) {
      ngroup <- c(rep('n1', n1), rep('n2', n2))
      colnames(sim.counts) <- paste0("S", 1:ncol(sim.counts), "_", ngroup)
    }
    if(!is.null(batch)) {
      ngroup <- c(rep('n1', n1), rep('n2', n2))
      nbatch <- ifelse(batch==-1, "b1", "b2")
      colnames(sim.counts) <- paste0("S", 1:ncol(sim.counts), "_", ngroup, "_", nbatch)
    }
  }

  ## return
  list(counts = sim.counts,
       mus=mus,
       disps=disps,
       drops=drops,
       designs = phenotype,
       sf = sf.value,
       simOptions = simOptions)
}


# Simulate multiple groups ------------------------------------------------

#' @importFrom stats model.matrix coef poisson sd rbinom
.simRNAseq.multi <- function(simOptions, n, verbose) {

  if(is.null(simOptions$bLFC)) {
    ## make group labels for phenotype LFC
    ## solve the issue by naming n!
    phenotype <- as.factor(unlist(sapply(names(n), function(i) {rep(i, n[i])}, simplify = F, USE.NAMES = F)))
    ## make model matrix
    mod = stats::model.matrix(~-1 + phenotype)
    coeffs = cbind(simOptions$pLFC)
    batch = NULL
  }

  if(!is.null(simOptions$bLFC)) {
    if(simOptions$bPattern=="uncorrelated") {
      ## make group labels for phenotype LFC
      phenotype <- as.factor(unlist(sapply(names(n), function(i) {rep(i, n[i])}, simplify = F, USE.NAMES = F)))
      # make batch labels for batch LFC
      batch <- as.factor(unlist(sapply(names(n), function(i) {rep_len(c(-1,1), n[i])}, simplify = F, USE.NAMES = F)))
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
    if(simOptions$bPattern=="orthogonal") {
      ## make group labels for phenotype LFC
      phenotype <- as.factor(unlist(sapply(names(n), function(i) {rep(i, n[i])}, simplify = F, USE.NAMES = F)))
      # make batch labels for batch LFC
      batch <- as.factor(unlist(sapply(names(n), function(i) {
        1 - 2*stats::rbinom(n[i],size=1,prob=0.5)
        }, simplify = F)))
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
    if(simOptions$bPattern=="correlated") {
      ## make group labels for phenotype LFC
      phenotype <- as.factor(unlist(sapply(names(n), function(i) {rep(i, n[i])}, simplify = F, USE.NAMES = F)))
      # make batch labels for batch LFC
      flip <- stats::rbinom(sum(n), size = 1, prob = 0.9)
      corgroup <- sample(levels(phenotype), size = 1)
      phenogroups <- ifelse(phenotype==corgroup, -1, 1)
      batch <- as.numeric(phenogroups)*flip + -as.numeric(phenogroups)*(1-flip)
      # make model matrix
      mod = stats::model.matrix(~-1 + phenotype + batch)
      coeffs = cbind(simOptions$pLFC, simOptions$bLFC)
    }
  }

  # generate read counts
  if (simOptions$RNAseq == "singlecell") {
    # take mean dispersion relationship as is
    if (!isTRUE(simOptions$downsample)) {
      if (verbose) {message(paste0("Predict dispersion based on true mean values.")) }
      if (attr(simOptions, 'Distribution') == 'NB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.NB.SampleWReplace_counts(sim.options = simOptions,
                                                  phenotype = phenotype,
                                                  modelmatrix = mod,
                                                  coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.NB.FillIn_counts(sim.options = simOptions,
                                          phenotype = phenotype,
                                          modelmatrix = mod,
                                          coef.dat=coeffs)
        }
      }
      if (attr(simOptions, 'Distribution') == 'ZINB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.ZINB.SampleWReplace_counts(sim.options = simOptions,
                                                    phenotype = phenotype,
                                                    modelmatrix = mod,
                                                    coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.ZINB.FillIn_counts(sim.options = simOptions,
                                            phenotype = phenotype,
                                            modelmatrix = mod,
                                            coef.dat=coeffs)
        }
      }
    }

    # downsample mean and then draw dispersion
    if (isTRUE(simOptions$downsample)) {
      if (verbose) {message(paste0("Predict dispersion based on effective mean values.")) }
      if (attr(simOptions, 'Distribution') == 'NB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.NB.SampleWReplace_downsample_counts(sim.options = simOptions,
                                                             phenotype = phenotype,
                                                             modelmatrix = mod,
                                                             coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.NB.FillIn_downsample_counts(sim.options = simOptions,
                                                     phenotype = phenotype,
                                                     modelmatrix = mod,
                                                     coef.dat=coeffs)
        }
      }
      if (attr(simOptions, 'Distribution') == 'ZINB') {
        if (!isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Sampling with replace.")) }
          sim.data = .sc.ZINB.SampleWReplace_downsample_counts(sim.options = simOptions,
                                                               phenotype = phenotype,
                                                               modelmatrix = mod,
                                                               coef.dat=coeffs)
        }
        if (isTRUE(simOptions$geneset)) {
          if (verbose) {message(paste0("Fill in with low magnitude Poisson.")) }
          sim.data = .sc.ZINB.FillIn_downsample_counts(sim.options = simOptions,
                                                       phenotype = phenotype,
                                                       modelmatrix = mod,
                                                       coef.dat=coeffs)
        }
      }
    }
  }

  if (simOptions$RNAseq == "bulk") {
    sim.data = .bulk.NB.RNAseq_counts(sim.options = simOptions,
                                      phenotype = phenotype,
                                      modelmatrix = mod,
                                      coef.dat=coeffs)
  }

  sim.counts = sim.data$counts
  sf.value = sim.data$sf
  mus = sim.data$mus
  disps = sim.data$disps
  drops = sim.data$drops

  sim.counts <- apply(sim.counts,2,function(x) {storage.mode(x) <- 'integer'; x})

  # fill in gene pseudonames if missing (only true for in silico)
  if (is.null(rownames(sim.counts))) {
    rownames(sim.counts) <- paste0("G", 1:nrow(sim.counts))
    names(mus) <- names(disps) <- names(drops) <- rownames(sim.counts)
  }
  # add human readable sample names
  if (is.null(colnames(sim.counts))) {
    if(is.null(batch)) {
      ngroup <- phenotype
      colnames(sim.counts) <- paste0("S", 1:ncol(sim.counts), "_", ngroup)
    }
    if(!is.null(batch)) {
      ngroup <- phenotype
      nbatch <- ifelse(batch==-1, "b1", "b2")
      colnames(sim.counts) <- paste0("S", 1:ncol(sim.counts), "_", ngroup, "_", nbatch)
    }
  }

  ## return
  list(counts = sim.counts,
       phenotypes = phenotype,
       mus=mus,
       disps=disps,
       drops=drops,
       batches = batch,
       sf = sf.value,
       simOptions = simOptions)
}


# SINGLE CELL -------------------------------------------------------------

# MEAN-VARIANCE RELATIONSHIP AS IS

#' @importFrom stats rnorm rnbinom
.sc.NB.SampleWReplace_counts <- function(sim.options,
                                         phenotype,
                                         modelmatrix,
                                         coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
      stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
      ncol = nsamples,
      nrow = ngenes,
      dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                      NULL))
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
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

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # effective dispersions
    if(length(disp) == 1) {
      sizevec <- rep(1/disp, ngenes)
    }
    if(length(disp) > 1) {
      sizevec <- disp
    }

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes,
                                   mu = 2 ^ mumat - 1, size = 1/disp),
                    ncol = nsamples, nrow = ngenes)
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
    }

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom rpois
.sc.NB.FillIn_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {

    # # make a low magnitude poisson count table
    # zerocounts = matrix(stats::rpois(ngenes*nsamples, lambda=0.1),
    #                     nrow=ngenes, ncol=nsamples)
    # # estimate nb parameters of this zero count table
    # zero_mus = rowMeans(zerocounts)
    # true.means = ifelse(is.na(zero_mus), 0, zero_mus)

    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    true.means = vector(mode = "numeric", length = ngenes)
    index = sample(1:length(mu), size = length(mu), replace = F)
    true.means[index] = mu[index]
    names(true.means)[index] = names(mu)[index]
    fillnames = sample(x = names(mu), size = ngenes - length(mu), replace = T)
    names(true.means)[-index] = fillnames
    true.means[-index] = min(mu)/2

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the library size factor vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[-index,] = min(mu)/2
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
      stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
      ncol = nsamples,
      nrow = ngenes,
      dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                      NULL))

    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
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

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # effective dispersions
    if(length(disp) == 1) {
      sizevec <- rep(1/disp, ngenes)
    }
    if(length(disp) > 1) {
      sizevec <- disp
    }

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes,
                                   mu = 2 ^ mumat - 1,
                                   size = 1/disp),
                    ncol = nsamples, nrow = ngenes)
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
    }

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom runif approx
.sc.ZINB.SampleWReplace_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))
    muvec = as.vector(mumat)

    meanp0fit.nonamplified = sim.options$meanp0fit.nonamplified
    meanp0fit.amplified = sim.options$meanp0fit.amplified
    meanp0fit = sim.options$meanp0fit
    nonamplified = sim.options$nonamplified

    if(!is.null(meanp0fit.nonamplified) && !is.null(meanp0fit.amplified)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat.nonamplified = abs(lmu)[p0.split==1]
      muvec.dat.amplified = abs(lmu)[p0.split==0]

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

      p0vec = rep(NA, length(lmu))
      p0vec[p0.split==0] = rnorm(n=length(which(p0.split==0)),
                                 mean=predp0.amplified.mean,
                                 sd=predp0.amplified.sd)
      p0vec[p0.split==1] = rnorm(n=length(which(p0.split==1)),
                                 mean=predp0.nonamplified.mean,
                                 sd=predp0.nonamplified.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(any(is.null(meanp0fit.nonamplified), is.null(meanp0fit.amplified)) && !is.null(meanp0fit)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat = abs(lmu)

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

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(is.null(meanp0fit.nonamplified) && is.null(meanp0fit.amplified) && is.null(meanp0fit)) {
      # result count matrix
      counts = matrix(stats::rnbinom(nsamples * ngenes,
                                     mu = 2 ^ mumat - 1,
                                     size = 2 ^ sizevec),
                      ncol = nsamples, nrow = ngenes)
      counts0 = counts == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples
    }
    }

  dimnames(counts) <- list(paste0(rownames(mumat),"_", seq_len(ngenes)), NULL)

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom runif approx rpois
.sc.ZINB.FillIn_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {

    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    true.means = vector(mode = "numeric", length = ngenes)
    index = sample(1:length(mu), size = length(mu), replace = F)
    true.means[index] = mu[index]
    names(true.means)[index] = names(mu)[index]
    fillnames = sample(x = names(mu), size = ngenes - length(mu), replace = T)
    names(true.means)[-index] = fillnames
    true.means[-index] = min(mu)/2

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[-index,] = min(mu)/2
    mumat[mumat < 0] = min(log2(effective.means + 1))
    muvec = as.vector(mumat)

    meanp0fit.nonamplified = sim.options$meanp0fit.nonamplified
    meanp0fit.amplified = sim.options$meanp0fit.amplified
    meanp0fit = sim.options$meanp0fit
    nonamplified = sim.options$nonamplified

    if(!is.null(meanp0fit.nonamplified) && !is.null(meanp0fit.amplified)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat.nonamplified = abs(lmu)[p0.split==1]
      muvec.dat.amplified = abs(lmu)[p0.split==0]

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

      p0vec = rep(NA, length(lmu))
      p0vec[p0.split==0] = rnorm(n=length(which(p0.split==0)),
                                 mean=predp0.amplified.mean,
                                 sd=predp0.amplified.sd)
      p0vec[p0.split==1] = rnorm(n=length(which(p0.split==1)),
                                 mean=predp0.nonamplified.mean,
                                 sd=predp0.nonamplified.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(any(is.null(meanp0fit.nonamplified), is.null(meanp0fit.amplified)) && !is.null(meanp0fit)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat = abs(lmu)

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

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(is.null(meanp0fit.nonamplified) && is.null(meanp0fit.amplified) && is.null(meanp0fit)) {
      # result count matrix
      counts = matrix(stats::rnbinom(nsamples * ngenes,
                                     mu = 2 ^ mumat - 1,
                                     size = 2 ^ sizevec),
                      ncol = nsamples, nrow = ngenes)
      counts0 = counts == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples
    }
    }

  dimnames(counts) <- list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                           NULL)

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

# MEAN-VARIANCE RELATIONSHIP AFTER MEAN DRAW

#' @importFrom stats rnorm rnbinom
.sc.NB.SampleWReplace_downsample_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)


    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # estimate size parameter associated with effective mean values
    lmu = log2(effective.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
      stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
      ncol = nsamples,
      nrow = ngenes,
      dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                      NULL))
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
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

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # effective dispersions
    if(length(disp) == 1) {
      sizevec <- rep(1/disp, ngenes)
    }
    if(length(disp) > 1) {
      sizevec <- disp
    }

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes,
                                   mu = 2 ^ mumat - 1, size = 1/disp),
                    ncol = nsamples, nrow = ngenes)
    }

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom rpois
.sc.NB.FillIn_downsample_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {

    # # make a low magnitude poisson count table
    # zerocounts = matrix(stats::rpois(ngenes*nsamples, lambda=0.1),
    #                     nrow=ngenes, ncol=nsamples)
    # # estimate nb parameters of this zero count table
    # zero_mus = rowMeans(zerocounts)
    # true.means = ifelse(is.na(zero_mus), 0, zero_mus)

    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    true.means = vector(mode = "numeric", length = ngenes)
    index = sample(1:length(mu), size = length(mu), replace = F)
    true.means[index] = mu[index]
    names(true.means)[index] = names(mu)[index]
    fillnames = sample(x = names(mu), size = ngenes - length(mu), replace = T)
    names(true.means)[-index] = fillnames
    true.means[-index] = min(mu)/2

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the library size factor vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # estimate size parameter associated with effective mean values
    lmu = log2(effective.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[-index,] = min(mu)/2
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
      stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
      ncol = nsamples,
      nrow = ngenes,
      dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                      NULL))
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
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

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # effective dispersions
    if(length(disp) == 1) {
      sizevec <- rep(1/disp, ngenes)
    }
    if(length(disp) > 1) {
      sizevec <- disp
    }

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(stats::rnbinom(nsamples * ngenes,
                                   mu = 2 ^ mumat - 1, size = 1/disp),
                    ncol = nsamples, nrow = ngenes)
    counts0 = counts == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples
    }

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom runif approx
.sc.ZINB.SampleWReplace_downsample_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {
    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    index = sample(1:length(mu), size = ngenes, replace = T)
    true.means = mu[index]

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # estimate size parameter associated with effective mean values
    lmu = log2(effective.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))
    muvec = as.vector(mumat)

    meanp0fit.nonamplified = sim.options$meanp0fit.nonamplified
    meanp0fit.amplified = sim.options$meanp0fit.amplified
    meanp0fit = sim.options$meanp0fit
    nonamplified = sim.options$nonamplified

    if(!is.null(meanp0fit.nonamplified) && !is.null(meanp0fit.amplified)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat.nonamplified = abs(lmu)[p0.split==1]
      muvec.dat.amplified = abs(lmu)[p0.split==0]

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

      p0vec = rep(NA, length(lmu))
      p0vec[p0.split==0] = rnorm(n=length(which(p0.split==0)),
                                 mean=predp0.amplified.mean,
                                 sd=predp0.amplified.sd)
      p0vec[p0.split==1] = rnorm(n=length(which(p0.split==1)),
                                 mean=predp0.nonamplified.mean,
                                 sd=predp0.nonamplified.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(any(is.null(meanp0fit.nonamplified), is.null(meanp0fit.amplified)) && !is.null(meanp0fit)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat = abs(lmu)

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

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(is.null(meanp0fit.nonamplified) && is.null(meanp0fit.amplified) && is.null(meanp0fit)) {
      # result count matrix
      counts = matrix(stats::rnbinom(nsamples * ngenes,
                                     mu = 2 ^ mumat - 1,
                                     size = 2 ^ sizevec),
                      ncol = nsamples, nrow = ngenes)

      counts0 = counts == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples
    }
    }

  dimnames(counts) <- list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                           NULL)

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}

#' @importFrom stats rnorm rnbinom runif approx rpois
.sc.ZINB.FillIn_downsample_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

  if (attr(sim.options, 'param.type') == 'estimated') {

    # # make a low magnitude poisson count table
    # zerocounts = matrix(stats::rpois(ngenes*nsamples, lambda=0.1),
    #                     nrow=ngenes, ncol=nsamples)
    # # estimate nb parameters of this zero count table
    # zero_mus = rowMeans(zerocounts)
    # true.means = ifelse(is.na(zero_mus), 0, zero_mus)

    # define NB params
    mu = sim.options$means
    meansizefit = sim.options$meansizefit

    # sample from the mean parameter observed
    true.means = vector(mode = "numeric", length = ngenes)
    index = sample(1:length(mu), size = length(mu), replace = F)
    true.means[index] = mu[index]
    names(true.means)[index] = names(mu)[index]
    fillnames = sample(x = names(mu), size = ngenes - length(mu), replace = T)
    names(true.means)[-index] = fillnames
    true.means[-index] = min(mu)/2



    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # estimate size parameter associated with effective mean values
    lmu = log2(effective.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule=2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule=2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[-index,] = min(mu)/2
    mumat[mumat < 0] = min(log2(effective.means + 1))
    muvec = as.vector(mumat)

    meanp0fit.nonamplified = sim.options$meanp0fit.nonamplified
    meanp0fit.amplified = sim.options$meanp0fit.amplified
    meanp0fit = sim.options$meanp0fit
    nonamplified = sim.options$nonamplified

    if(!is.null(meanp0fit.nonamplified) && !is.null(meanp0fit.amplified)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat.nonamplified = abs(lmu)[p0.split==1]
      muvec.dat.amplified = abs(lmu)[p0.split==0]

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

      p0vec = rep(NA, length(lmu))
      p0vec[p0.split==0] = rnorm(n=length(which(p0.split==0)),
                                 mean=predp0.amplified.mean,
                                 sd=predp0.amplified.sd)
      p0vec[p0.split==1] = rnorm(n=length(which(p0.split==1)),
                                 mean=predp0.nonamplified.mean,
                                 sd=predp0.nonamplified.sd)

      p0vec[p0vec<0] = runif(length(which(p0vec<0)),0,0.05)
      p0vec[p0vec>1] = runif(length(which(p0vec>1)),0.85,1)

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(any(is.null(meanp0fit.nonamplified), is.null(meanp0fit.amplified)) && !is.null(meanp0fit)) {
      # predict the dropout probabilities associated with sampled mean
      p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
      muvec.dat = abs(lmu)

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

      zero.mark = sapply(p0vec, function(x) {
        rbinom(nsamples,  prob = (1 - x), size = 1)
      }, simplify = T)

      p0mat = matrix(t(zero.mark), nrow = ngenes, ncol = nsamples)

      counts0 = p0mat == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples

      # result count matrix
      counts.nb = matrix(stats::rnbinom(nsamples * ngenes,
                                        mu = 2 ^ mumat - 1,
                                        size = 2 ^ sizevec),
                         ncol = nsamples, nrow = ngenes)
      counts = counts.nb * p0mat
    }

    if(is.null(meanp0fit.nonamplified) && is.null(meanp0fit.amplified) && is.null(meanp0fit)) {
      # result count matrix
      counts = matrix(stats::rnbinom(nsamples * ngenes,
                                     mu = 2 ^ mumat - 1,
                                     size = 2 ^ sizevec),
                      ncol = nsamples, nrow = ngenes)
      counts0 = counts == 0
      nn0 = rowSums(!counts0)
      dropout = (nsamples - nn0)/nsamples
    }
    }

  dimnames(counts) <- list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                           NULL)

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}


# BULK --------------------------------------------------------------------

#' @importFrom stats rnorm rbinom rnbinom approx
.bulk.NB.RNAseq_counts <- function(sim.options, phenotype, modelmatrix, coef.dat) {

  set.seed(sim.options$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = sim.options$ngenes
  mod = modelmatrix
  beta = coef.dat

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

    counts0 = p0mat == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    cnts = matrix(
      stats::rnbinom(nsamples * ngenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
      ncol = nsamples, nrow = ngenes,
      dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                      NULL))
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

    counts0 = p0mat == 0
    nn0 = rowSums(!counts0)
    dropout = (nsamples - nn0)/nsamples

    # associate with dispersion parameter
    if (length(sim.options$dispersion) == 1) {
      disp = sim.options$dispersion
    }
    if (is.function(sim.options$dispersion)) {
      disp = sim.options$dispersion(true.means)
    }

    # size factor
    if(is.list(sim.options$size.factors)) {
      all.facs <- c(sim.options$size.factors$n1(length(which(phenotype==-1))),
                    sim.options$size.factors$n2(length(which(phenotype==1))))
    }
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
        stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
      }
      }

    # effective means
    effective.means <- outer(true.means, all.facs, "*")

    # effective dispersions
    if(length(disp) == 1) {
      sizevec <- rep(1/disp, ngenes)
    }
    if(length(disp) > 1) {
      sizevec <- disp
    }

    # make mean expression with beta coefficients added as defined by model matrix
    ind = !apply(mod, 2, function(x) { all(x == 1) })
    mod = cbind(mod[, ind])
    beta = cbind(beta[, ind])
    mumat = log2(effective.means + 1) + beta %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    cnts = matrix(stats::rnbinom(nsamples * ngenes,
                                 mu = 2 ^ mumat - 1,
                                 size = 1/disp),
                  ncol = nsamples, nrow = ngenes)
    counts = p0mat * cnts
    }

  return(list(counts=counts, sf=all.facs, mus=true.means, disps=sizevec, drops=dropout))
}


# Simulate spike-in reads -------------------------------------------------

.simSpike <- function(SpikeOptions, n1, n2, sf = NULL) {
  # sample mean expression values of spike-ins
  predictedMean = rowMeans(SpikeOptions$normCounts) /
    (SpikeOptions$EVGammaThetaEstimates$EGamma * SpikeOptions$EVGammaThetaEstimates$ETheta)

  # create size factor matrix
  if(is.null(sf)) {
    sf =  rep(1, n1+n2)
    sizefactors.mat = .repmat(t(as.matrix(sf)), length(predictedMean), 1)
  }
  if(!is.null(sf)) {
    sizefactors.mat = .repmat(t(as.matrix(sf)), length(predictedMean), 1)
  }


  # simulate spike-ins assuming no biological variance contribution
  spike.cnts <- .simulateCountGenes(Xi=predictedMean,
                                    ETheta = SpikeOptions$EVGammaThetaEstimates$ETheta,
                                    VTheta = SpikeOptions$EVGammaThetaEstimates$VTheta,
                                    EGamma = SpikeOptions$EVGammaThetaEstimates$EGamma,
                                    VGamma =  SpikeOptions$EVGammaThetaEstimates$VGamma,
                                    Aij = sizefactors.mat,
                                    nCell = n1+n2,
                                    BV = 0*(n1+n2))

  return(list(counts=spike.cnts, sf=sf))
}

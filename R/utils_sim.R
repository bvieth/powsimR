# TODO
# change behaviour of minus mean values ?

# Simulate DE between 2 groups (DE of mean) -------------------------------

#' @importFrom stats model.matrix rbinom
.simRNAseq.2grp <- function(simOptions, n1, n2, verbose) {

  set.seed(simOptions$DESetup$sim.seed)

  if(is.null(simOptions$DESetup$bLFC)) {
    ## make group labels for phenotype LFC
    designs <- c(-stats::rbinom(n = n1, size = 1, prob = simOptions$DESetup$p.G),
                   stats::rbinom(n = n2, size = 1, prob = simOptions$DESetup$p.G))
    phenotype <- ifelse(designs == -1, 0, designs)
    batch <- NULL
    ## make model matrix
    modelmatrix = stats::model.matrix(~-1 + phenotype)
    coef.dat = cbind(simOptions$DESetup$pLFC)
  }

  if(!is.null(simOptions$DESetup$bLFC)) {
    if(simOptions$DESetup$bPattern=="uncorrelated") {
      # make group labels for phenotype LFC
      designs <- c(-stats::rbinom(n = n1, size = 1, prob = simOptions$DESetup$p.G),
                   stats::rbinom(n = n2, size = 1, prob = simOptions$DESetup$p.G))
      phenotype <- ifelse(designs == -1, 0, designs)
      # make batch labels for batch LFC
      batch <- rep_len(c(0,1), n1+n2)
      # make model matrix
      modelmatrix = stats::model.matrix(~-1 + phenotype + batch)
      coef.dat = cbind(simOptions$DESetup$pLFC, simOptions$DESetup$bLFC)
    }
    if(simOptions$DESetup$bPattern=="orthogonal") {
      # make group labels for phenotype LFC
      designs <- c(-stats::rbinom(n = n1, size = 1, prob = simOptions$DESetup$p.G),
                   stats::rbinom(n = n2, size = 1, prob = simOptions$DESetup$p.G))
      phenotype <- ifelse(designs == -1, 0, designs)
      # make batch labels for batch LFC
      batch <-  1 - rbinom(n1+n2, size=1, prob=0.5)
      # make model matrix
      modelmatrix = stats::model.matrix(~-1 + phenotype + batch)
      coef.dat = cbind(simOptions$DESetup$pLFC, simOptions$DESetup$bLFC)
    }
    if(simOptions$DESetup$bPattern=="correlated") {
      # make group labels for phenotype LFC
      designs <- c(-stats::rbinom(n = n1, size = 1, prob = simOptions$DESetup$p.G),
                   stats::rbinom(n = n2, size = 1, prob = simOptions$DESetup$p.G))
      phenotype <- ifelse(designs == -1, 0, designs)
      # make batch labels for batch LFC
      flip <- rbinom(n1+n2, size = 1, prob = 0.9)
      batch <- phenotype*flip + -phenotype*(1-flip)
      batch <- ifelse(batch == -1, 0, batch)
      # make model matrix
      modelmatrix = stats::model.matrix(~-1 + phenotype + batch)
      coef.dat = cbind(simOptions$DESetup$pLFC, simOptions$DESetup$bLFC)
    }
  }

  # generate read counts
  # NB distribution
  if (attr(simOptions$estParamRes, 'Distribution') == 'NB'){
    sim.data <- .sc.NB_counts(simOptions = simOptions,
                              phenotype = phenotype,
                              modelmatrix = modelmatrix,
                              coef.dat = coef.dat)
  }

  # ZINB distribution
  if (attr(simOptions$estParamRes, 'Distribution') == 'ZINB'){
    sim.data <- .sc.ZINB_counts(simOptions = simOptions,
                                phenotype = phenotype,
                                modelmatrix = modelmatrix,
                                coef.dat = coef.dat)
  }

  # retrieve output
  sim.counts = sim.data$counts
  sf.value = sim.data$sf
  mus = sim.data$mus
  disps = sim.data$disps
  drops = sim.data$drops

  sim.counts <- apply(sim.counts, 2, function(x) {storage.mode(x) <- 'integer'; x})

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

  invisible(gc())

  ## return
  list(counts = sim.counts,
       mus=mus,
       disps=disps,
       drops=drops,
       designs = designs,
       sf = sf.value,
       simOptions = simOptions)
}

# Simulate NB -------------------------------------------------------------

#' @importFrom stats rnorm rnbinom
#' @importFrom truncnorm rtruncnorm
.sc.NB_counts <- function(simOptions,
                          phenotype,
                          modelmatrix,
                          coef.dat) {

  set.seed(simOptions$DESetup$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = simOptions$DESetup$ngenes

  # define NB params
  means = simOptions$estParamRes$Parameters[[simOptions$DESetup$Draw$MoM]]$means
  means = means[means>0]
  sizes = simOptions$estParamRes$Parameters[[simOptions$DESetup$Draw$MoM]]$size
  meansizefit = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$meansizefit

  # sample from the mean parameter observed
  if(ngenes <= length(means)){
    index = sample(1:length(means), size = ngenes,
                   replace = simOptions$DESetup$SwReplace)
  }
  if(ngenes > length(means)){
    index = sample(1:length(means), size = ngenes, replace = T)
  }
  true.means = means[index]

  # estimate size parameter associated with true mean values
  lmu = log2(true.means)
  predsize.mean = suppressWarnings(approx(meansizefit$x,
                                          meansizefit$y,
                                          xout = lmu, rule=2)$y)
  predsize.sd = suppressWarnings(approx(meansizefit$x,
                                        meansizefit$sd,
                                        xout = lmu, rule=2)$y)
  lsizevec = truncnorm::rtruncnorm(n = length(lmu),
                                  a = min(log2(sizes), na.rm = TRUE),
                                  b = max(log2(sizes), na.rm = TRUE),
                                  mean = predsize.mean,
                                  sd = predsize.sd)
  sizevec = 2^lsizevec

  # size factor
  if (simOptions$SimSetup$LibSize == "equal") {
    all.facs <- rep(1, nsamples)
  }
  if (simOptions$SimSetup$LibSize == "given") {
    if (is.function(simOptions$estParamRes$sf)) {
      all.facs <- simOptions$estParamRes$sf(nsamples)
    }
    if (is.vector(simOptions$estParamRes$sf)) {
      all.facs <- sample(simOptions$estParamRes$sf, nsamples, replace = TRUE)
    }
    if (is.null(simOptions$estParamRes$sf)) {
      stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
    }
  }

  # effective means
  effective.means <- outer(true.means, all.facs, "*")

  # make mean expression with beta coefficients added as defined by model matrix
  ind = !apply(modelmatrix, 2, function(x) { all(x == 1) })
  mod = cbind(modelmatrix[, ind])
  beta = cbind(coef.dat[, ind])
  coef.mat = 2^beta %*% t(mod)
  mumat = effective.means * coef.mat
  # mumat[mumat < 0] = min(log2(effective.means+1))

  # result count matrix
  counts = matrix(
    stats::rnbinom(nsamples * ngenes,
                   mu = mumat,
                   size = sizevec),
    ncol = nsamples,
    nrow = ngenes,
    dimnames = list(paste0(rownames(mumat),"_", seq_len(ngenes)),
                    NULL))

  counts0 = counts == 0
  nn0 = rowSums(!counts0)
  dropout = (nsamples - nn0)/nsamples

  return(list(counts=counts,
              sf=all.facs,
              mus=true.means,
              disps=1/sizevec,
              sizes=sizevec,
              drops=dropout))
}

# Simulate ZINB -----------------------------------------------------------

#' @importFrom stats rnorm rnbinom runif approx
#' @importFrom truncnorm rtruncnorm
.sc.ZINB_counts <- function(simOptions,
                            phenotype,
                            modelmatrix,
                            coef.dat) {

  set.seed(simOptions$DESetup$sim.seed)

  nsamples = nrow(modelmatrix)
  ngenes = simOptions$DESetup$ngenes

  # define ZINB params
  pos.means = simOptions$estParamRes$Parameters[[simOptions$DESetup$Draw$MoM]]$pos.means
  pos.means = pos.means[pos.means>0]
  pos.sizes = simOptions$estParamRes$Parameters[[simOptions$DESetup$Draw$MoM]]$pos.size
  meansizefit = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$meansizefit

  # sample from the positive mean parameter observed
  if(ngenes <= length(pos.means)){
    index = sample(1:length(pos.means), size = ngenes,
                   replace = simOptions$DESetup$SwReplace)
  }
  if(ngenes > length(pos.means)){
    index = sample(1:length(pos.means), size = ngenes, replace = T)
  }
  true.means = pos.means[index]

  # estimate size parameter associated with true mean values
  lmu = log2(true.means)
  predsize.mean = suppressWarnings(approx(meansizefit$x,
                                          meansizefit$y,
                                          xout = lmu, rule=2)$y)
  predsize.sd = suppressWarnings(approx(meansizefit$x,
                                        meansizefit$sd,
                                        xout = lmu, rule=2)$y)
  lsizevec = truncnorm::rtruncnorm(n = length(lmu),
                                  a = min(log2(pos.sizes), na.rm = TRUE),
                                  b = max(log2(pos.sizes), na.rm = TRUE),
                                  mean = predsize.mean,
                                  sd = predsize.sd)
  sizevec = 2^lsizevec

  # size factor
  if (simOptions$SimSetup$LibSize == "equal") {
    all.facs <- rep(1, nsamples)
  }
  if (simOptions$SimSetup$LibSize == "given") {
    if (is.function(simOptions$estParamRes$sf)) {
      all.facs <- simOptions$estParamRes$sf(nsamples)
    }
    if (is.vector(simOptions$estParamRes$sf)) {
      all.facs <- sample(simOptions$estParamRes$sf, nsamples, replace = TRUE)
    }
    if (is.null(simOptions$estParamRes$sf)) {
      stop(message(paste0("You chose to draw from the given size factors,
                            however the sf vector is empty!")))
    }
  }

  # effective means
  effective.means <- outer(true.means, all.facs, "*")

  # make mean expression with beta coefficients added as defined by model matrix
  ind = !apply(modelmatrix, 2, function(x) { all(x == 1) })
  mod = cbind(modelmatrix[, ind])
  beta = cbind(coef.dat[, ind])
  coef.mat = 2^beta %*% t(mod)
  mumat = effective.means * coef.mat
  # mumat[mumat < 0] = min(log2(effective.means+1))

  meang0fit.nonamplified = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$meang0fit.nonamplified
  meang0fit.amplified = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$meang0fit.amplified
  meang0fit = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$meang0fit
  nonamplified = simOptions$estParamRes$Fit[[simOptions$DESetup$Draw$Fit]]$nonamplified

  if(!is.null(meang0fit.nonamplified) && !is.null(meang0fit.amplified)) {
    # predict the dropout probabilities associated with sampled mean
    p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
    muvec.dat.nonamplified = lmu[p0.split==1]
    muvec.dat.amplified = lmu[p0.split==0]

    predp0.amplified.mean = approx(meang0fit.amplified$x,
                                   meang0fit.amplified$y,
                                   xout=muvec.dat.amplified,
                                   rule = 2:1)$y
    predp0.amplified.sd = approx(meang0fit.amplified$x,
                                 meang0fit.amplified$sd,
                                 xout=muvec.dat.amplified,
                                 rule=2:1)$y
    predp0.nonamplified.mean = approx(meang0fit.nonamplified$x,
                                      meang0fit.nonamplified$y,
                                      xout=muvec.dat.nonamplified,
                                      rule=2:1)$y
    predp0.nonamplified.sd = approx(meang0fit.nonamplified$x,
                                    meang0fit.nonamplified$sd,
                                    xout=muvec.dat.nonamplified,
                                    rule=2)$y

    p0vec = rep(NA, length(lmu))
    p0vec[p0.split==0] = truncnorm::rtruncnorm(n=length(which(p0.split==0)),
                                               a = 0, b = 1,
                                               mean=predp0.amplified.mean,
                                               sd=predp0.amplified.sd)
    p0vec[p0.split==1] = truncnorm::rtruncnorm(n=length(which(p0.split==1)),
                                               a = 0, b = 1,
                                               mean=predp0.nonamplified.mean,
                                               sd=predp0.nonamplified.sd)

    p0vec[is.na(p0vec)] = runif(n= length(p0vec[is.na(p0vec)]),
                                min = 0, max = 0.05)

    zero.mark = sapply(p0vec, function(x) {
      rbinom(nsamples,  prob = 1-x, size = 1)
    }, simplify = T)

    p0mat = matrix(t(zero.mark),
                   nrow = ngenes,
                   ncol = nsamples)

  }

  if(any(is.null(meang0fit.nonamplified), is.null(meang0fit.amplified)) &&
     !is.null(meang0fit)) {
    # predict the dropout probabilities associated with sampled mean
    p0.split = rbinom(length(lmu), prob = nonamplified, size = 1)
    muvec.dat = abs(lmu)

    predp0.mean = approx(meang0fit$x,
                         meang0fit$y,
                         xout=muvec.dat,
                         rule = 2:1)$y
    predp0.sd = approx(meang0fit$x,
                       meang0fit$sd,
                       xout=muvec.dat,
                       rule=2:1)$y

    p0vec = truncnorm::rtruncnorm(n=length(muvec.dat),
                                  a = 0, b = 1,
                                  mean=predp0.mean,
                                  sd=predp0.sd)
    p0vec[is.na(p0vec)] = runif(n= length(p0vec[is.na(p0vec)]),
                                min = 0, max = 0.05)

    zero.mark = sapply(p0vec, function(x) {
      rbinom(nsamples,  prob = (1 - x), size = 1)
    }, simplify = T)

    p0mat = matrix(t(zero.mark),
                   nrow = ngenes,
                   ncol = nsamples)
  }

  if(is.null(meang0fit.nonamplified) &&
     is.null(meang0fit.amplified) &&
     is.null(meang0fit)) {
# at the moment this is not implemented, maybe throw an error in estimateParam ?
    p0mat <- NULL

  }

  # make the count matrix
    counts.nb = matrix(
      stats::rnbinom(nsamples * ngenes,
                     mu = mumat,
                     size = sizevec),
      ncol = nsamples,
      nrow = ngenes)

    if(!is.null(p0mat)){

      # check that the overall number of zeroes is correct
      nb0 = counts.nb == 0
      nn0 = rowSums(!nb0)
      dropout.nb = (nsamples - nn0)/nsamples
      p0mat0 = p0mat == 0
      nn0 = rowSums(!p0mat0)
      dropout.p0 = (nsamples - nn0)/nsamples
      change.nb = dropout.nb < dropout.p0

      p0mat.nb = p0mat[change.nb, ]
      p0mat.nb[counts.nb[change.nb, ] == 0 & p0mat[change.nb, ] == 0] = 1
      p0mat[change.nb, ] = p0mat.nb
      counts = counts.nb * p0mat
    }
    if(is.null(p0mat)){
      counts = counts.nb
    }

  dimnames(counts) <- list(paste0(rownames(mumat),"_", seq_len(ngenes)), NULL)

  counts0 = counts == 0
  nn0 = rowSums(!counts0)
  dropout = (nsamples - nn0)/nsamples

  #return object
  return(list(counts=counts,
              sf=all.facs,
              mus=true.means,
              disps=1/sizevec,
              sizes=sizevec,
              drops=dropout))

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
  spike.cnts <- .simulateCountGenes(Xi = predictedMean,
                                    ETheta = SpikeOptions$EVGammaThetaEstimates$ETheta,
                                    VTheta = SpikeOptions$EVGammaThetaEstimates$VTheta,
                                    EGamma = SpikeOptions$EVGammaThetaEstimates$EGamma,
                                    VGamma = SpikeOptions$EVGammaThetaEstimates$VGamma,
                                    Aij = sizefactors.mat,
                                    nCell = n1+n2,
                                    BV = 0*(n1+n2))

  return(list(counts=spike.cnts, sf=sf))
}


# Introduce dropout genes -------------------------------------------------

.dropGene <- function(simOptions, simData){

  set.seed(simOptions$DESetup$sim.seed)

  nsamples = ncol(simData$counts)
  ngenes = simOptions$DESetup$ngenes
  ndrop = floor(ngenes * simOptions$SimSetup$DropRate)
  all.facs = simData$sf

  # define NB params
  means = simOptions$estParamRes$Parameters$DropGene$means
  sizes = simOptions$estParamRes$Parameters$DropGene$size
  meansizefit = simOptions$estParamRes$Fit$DropGene$meansizefit

  # sample from the mean parameter observed
  index = sample(1:length(means), size = ndrop, replace = T)
  true.means = means[index]

  # estimate size parameter associated with true mean values
  lmu = log2(true.means+1)
  predsize.mean = suppressWarnings(approx(meansizefit$x,
                                          meansizefit$y,
                                          xout = lmu, rule=2)$y)
  predsize.sd = suppressWarnings(approx(meansizefit$x,
                                        meansizefit$sd,
                                        xout = lmu, rule=2)$y)
  sizevec = truncnorm::rtruncnorm(n = length(lmu),
                                  a = min(log2(sizes), na.rm = TRUE),
                                  b = max(log2(sizes), na.rm = TRUE),
                                  mean = predsize.mean,
                                  sd = predsize.sd)

  # effective means
  effective.means <- outer(true.means, all.facs, "*")

  d.index = sample(1:ngenes, size = ndrop, replace = F)

  # replace counts with dropouts
  dcounts = simData$counts
  dcounts[d.index,] = suppressWarnings(
    matrix(stats::rnbinom(nsamples * ndrop,
                           mu = 2^effective.means-1,
                           size = 2^sizevec),
            ncol = nsamples,
            nrow = ndrop)
    )
  dcounts[is.na(dcounts)] <- 0
  dcounts <- apply(dcounts, 2, function(x) {storage.mode(x) <- 'integer'; x})

  simData$counts <- dcounts

  return(simData)
}


# Binomial Thinning -------------------------------------------------------

#' @importFrom edgeR thinCounts
.run.thin <- function(countData,
                      Thin,
                      simOptions){

  # convert UMI data to Read data
  if(attr(simOptions$estParamRes, 'Protocol') == "UMI"){
    umireadfit = simOptions$estParamRes$Fit$UmiRead$Fit
    lUMI <- as.vector(log10(countData+1))

    predratio.mean = suppressWarnings(approx(x = umireadfit$x,
                                             y = umireadfit$y,
                                             xout = lUMI,
                                             rule=2)$y)
    predratio.sd = suppressWarnings(approx(x = umireadfit$x,
                                           y = umireadfit$sd,
                                           xout = lUMI,
                                           rule=2)$y)
    predratio.vec = truncnorm::rtruncnorm(n = length(lUMI),
                                          mean = predratio.mean,
                                          sd = predratio.sd,
                                          a = 1)
    predratio.dat <- matrix(predratio.vec,
                            nrow = nrow(countData),
                            ncol = ncol(countData),
                            dimnames = list(rownames(countData), colnames(countData)),
                            byrow = F)

    predread.dat <- log10(countData+1) * predratio.dat
    readData <- floor(10^predread.dat-1)

    # apply binomial thinning
    thinCounts <- edgeR::thinCounts(x = readData, prob = rep(Thin, ncol(readData)))
    DropUMI <- thinCounts - countData < 0
    thinData <- countData
    thinData[DropUMI] <- 0

  }
  if(attr(simOptions$estParamRes, 'Protocol') == "Read"){
    readData <- countData
    thinCounts <- edgeR::thinCounts(x = readData, prob = rep(Thin, ncol(readData)))
    thinData <- thinCounts
  }

  thinData <- apply(thinData, 2, function(x) {storage.mode(x) <- 'integer'; x})

  # return object
  return(thinData)
}





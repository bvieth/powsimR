# set LFC -----------------------------------------------------------------
#' @importFrom methods is

.setFC <- function(input, nDEgenes, k) {
  # two groups
  if(k==2) {
    if (is.vector(input)) { ## vector
      if (length(input) == 1) { ## constant
        lfc = rep(input, nDEgenes)
      } else if (length(input) != nDEgenes) { ## vector
        lfc = sample(input, nDEgenes, replace = TRUE)
      } else {
        lfc = input
      }
    } else if (is.function(input)) { # a function
      lfc = input(nDEgenes)
    }   else {
      stop("Unrecognized form of lfc!\n")
    }
  }

  # multi group
  if(k>2) {
    if(is.list(input)) { # list
      if(!length(input) == k) { stop("The LFC vector has not the same number of entries as n groups specified!")}
      if(all(sapply(input, length)==1)) { # constant lfc
        lfc = sapply(1:length(input), function(i) { rep(i, nDEgenes)})
      }
      if(any(!sapply(input, length)==nDEgenes)) { ## vector of lfcs
        lfc = sapply(1:length(input), function(i) {
          sample(x = input[[i]], size = nDEgenes, replace = TRUE)
        })
      }
      if(all(sapply(input, length)==nDEgenes)) {
        lfc = do.call("cbind", input)
      }
    }
    if(is.function(input)) {
      lfc = input(nDEgenes)
      if(any(c(!ncol(lfc) == k), c(is.vector(lfc)))) {
        stop("The provided LFC function has not the same number of entries as n groups specified!")
        }
    }
    if(all(!is.list(input), !is.function(input))) {
      stop(message("Unrecognized form of lfc!"))
      }
  }
  # return lfc object
  return(lfc)
}

.setMarker <- function(input, nDEgenes, idpool, ngenes, k) {
  lfcs =  matrix(0, nrow = ngenes, ncol = k)
  # generate random id for marker genes and sample group
  Mids <- sample(idpool, nDEgenes, replace = FALSE)
  Kids <- sample(1:k, length(Mids), replace=TRUE)

  if (is.vector(input)) { ## vector
    if (length(input) == 1) { ## constant
      lfc = rep(input, nDEgenes)
    } else if (length(input) != nDEgenes) { ## vector
      lfc = sample(input, nDEgenes, replace = TRUE)
    } else {
      lfc = input
    }
  } else if (is.function(input)) { # a function
    lfc = input(nDEgenes)
  }   else {
    stop("Unrecognized form of lfc!\n")
  }

  for (i in 1:length(Mids)) {
    lfcs[Mids[i], Kids[i]] = lfc[i]
  }

  # return matrix lfc object
  return(list(lfcs=lfcs, id=Mids))
}


# Calculate AIC of model --------------------------------------------------

.myAIC <- function(resids, dat.n, npar, k) {
  sse0 <- rowSums(resids^2)
  aic <- dat.n + dat.n*log(2*pi) + dat.n * log(sse0/dat.n) + k*npar
  return(aic)
}

# ## Calculate goodness of fit chisquare test from fitdistrplus ------------

#' @importFrom fitdistrplus gofstat
.fitdistrplusGOF <- function(fitdistobj) {
  gof <- fitdistrplus::gofstat(f = fitdistobj, discrete = T)
  data.frame(gof.stats=gof$chisq,
             gof.df = gof$chisqdf,
             gof.pval = gof$chisqpvalue,
             stringsAsFactors = F)
}

# Calculate goodness of fit -----------------------------------------------

#' @importFrom stats pchisq
.myGOF <- function(deviances, df.residuals) {
  data.frame(gof.stats = deviances,
             gof.df = df.residuals,
             gof.pvals = pchisq(q=deviances, df=df.residuals, lower.tail = F),
             stringsAsFactors = F)
}

# Poisson Beta Moments Estimation -----------------------------------------

#' @importFrom moments moment
#' @importFrom stats rbeta rpois dpois
.PoissonBetaFit <- function(x.raw, x.norm, nMC = 1e3) {
  # natural estimators of the first three moments
  m1 = moments::moment(x.norm, order=1, central=F)
  m2 = moments::moment(x.norm, order=2, central=T)
  m3 = moments::moment(x.norm, order=3, central=T)

  # parameters of the beta poisson distribution (Hemberg lab)
  bp.gamma = 2*m1 - (m3*(m1^2 + m2))/(- 2*m2^2 + m1*m3)
  bp.alpha = (2*m1*(m1^2*m2 + m3*m1 - m2^2))/(- m3*m1^2 + 4*m1*m2^2 + m3*m2)
  bp.beta = -(2*m2*(m3 + 2*m1*m2)*(m1^2*m2 + m3*m1 - m2^2))/((- 2*m2^2 + m1*m3)*(- m3*m1^2 + 4*m1*m2^2 + m3*m2))

  # parameters of beta-poisson following Marioni estimation
  e1.foo <- function(x) { sum(x)/length(x) }
  e2.foo <- function(x) { sum(x*(x-1))/length(x) }
  e3.foo <- function(x) { sum(x*(x-1)*(x-2))/length(x) }
  # e1, e2, e3 exponentional moments
  e1 = e1.foo(x.norm)
  e2 = e2.foo(x.norm)
  e3 = e3.foo(x.norm)

  r1 = e1
  r2 = e2 / e1
  r3 = e3 / e2

  bp.alpha2 = 2*r1*(r3 - r2) / (r1*r2 - 2*r1*r3 + r2*r3)
  bp.beta2 = 2*(r2 - r1)*(r1 - r3)*(r3 - r2) / ((r1*r2 - 2*r1*r3 + r2*r3) *(r1 - 2*r2 + r3))
  bp.gamma2 = (-r1*r2 + 2*r1*r3 - r2*r3) / (r1 - 2*r2 +r3)

  # calculate the Log Likelihood
  n <- length(x.raw)
  logLike <- rep(0, times=n)
  for (i in 1:n) {
    if (bp.alpha>0 && bp.beta>0 && bp.gamma>0) {
      q <- stats::rbeta(nMC, bp.alpha, bp.beta)
      p <- stats::dpois(x=x.raw[i], lambda=bp.gamma*q)
      logLike[i] <- log(mean(p))
    }
    else {
      logLike[i] <- NA
    }
  }
  LogLikelihood <- sum(logLike)
  estAIC <- -2*LogLikelihood + 2*3
  deviance <- -2*LogLikelihood
  df.residual <- n-2

  # predicted zeroes
  if (bp.alpha>0 && bp.beta>0 && bp.gamma>0) {
    #generate Beta random variables
    y <- stats::rbeta(n = n, bp.alpha, bp.beta)
    #draw counts
    PBcnts <- stats::rpois(n=n, bp.gamma*y)
    predZero <- sum(PBcnts==0)
  }
  else {
    predZero <- NA
  }

  # calculate the Log Likelihood
  logLike2 <- rep(0, times=n)
  for (i in 1:n) {
    if (bp.alpha2>0 && bp.beta2>0 && bp.gamma2>0) {
      q <- stats::rbeta(nMC, bp.alpha2, bp.beta2)
      p <- stats::dpois(x=x.raw[i], lambda=bp.gamma2*q)
      logLike2[i] <- log(mean(p))
    }
    else {
      logLike2[i] <- NA
    }
  }
  LogLikelihood2 <- sum(logLike2)
  estAIC2 <- -2*LogLikelihood2 + 2*3
  deviance2 <- -2*LogLikelihood2
  df.residual2 <- n-2

  # predicted zeroes
  if (bp.alpha2>0 && bp.beta2>0 && bp.gamma2>0) {
    #generate Beta random variables
    y <- stats::rbeta(n = n, bp.alpha2, bp.beta2)
    #draw counts
    PBcnts <- stats::rpois(n=n, bp.gamma2*y)
    predZero2 <- sum(PBcnts==0)
  }
  else {
    predZero2 <- NA
  }

  return(list("bp.alpha.hemberg" = bp.alpha, "bp.beta.hemberg" = bp.beta, "bp.gamma.hemberg" = bp.gamma,
              "bp.alpha.marioni" = bp.alpha2, "bp.beta.marioni" = bp.beta2, "bp.gamma.marioni" = bp.gamma2,
              "LogLikelihood.hemberg" = LogLikelihood, "LogLikelihood.marioni" = LogLikelihood2,
              "Deviance.hemberg" = deviance, "Deviance.marioni" = deviance2,
              "Dd.hemberg" = df.residual, "Dd.marioni" = df.residual2,
              "AIC.hemberg" = estAIC, "AIC.marioni" = estAIC2,
              "PredZero.hemberg"=predZero, "PredZero.marioni"=predZero2))
}


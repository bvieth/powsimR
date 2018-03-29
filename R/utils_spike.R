# Taken from the supplementary code file of
# 'Characterizing noise structure in single-cell RNA-seq distinguishes genuine from technical allelic expression'
# by Jong Kyoung Kim
# published in Nature Communications
# DOI: 10.1038/ncomms9687

# Estimate gamma and theta
#' @importFrom stats lm coefficients nls.control
#' @importFrom minpack.lm nlsLM
.estimateGammaTheta <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  fitData = data.frame(Xi=numberSpikes[,1], Ki=rowMeans(nCountSpikes))
  fit1 <- stats::lm(Ki ~ 0 + Xi, data=fitData)
  gammaTheta = stats::coefficients(fit1)
  PropCell = sapply(1:nrow(nCountSpikes), function(x) {
    sum(nCountSpikes[x,]>0)/ncol(nCountSpikes)
  })
  fitData = data.frame(Xi=numberSpikes[,1], Y=PropCell, Ai=rowMeans(sizeFactorMatrix))
  initialTheta = mean(fitData$Y[fitData$Xi>0.5 & fitData$Xi<5])
  if (initialTheta==0 | is.na(initialTheta)) {
    initialTheta = 0.01
  }
  fit2 <- minpack.lm::nlsLM(Y ~ 1 - (1 - Theta + Theta*exp(-(gammaTheta/Theta)*Ai))^Xi,
                            data=fitData, start=list(Theta=initialTheta),
                            lower=0, upper=1, control=stats::nls.control(warnOnly=TRUE))


  Theta = stats::coefficients(fit2)
  list(gammaTheta=gammaTheta, Theta=Theta)
}

# Initializing V[gamma] and V[theta]
#' @importFrom stats lm coefficients nls.control var
#' @importFrom minpack.lm nlsLM
#' @importFrom MASS fitdistr
#' @importFrom plyr aaply
.estimateVGammaThetaInitial <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  old.o = options()
  options(warn=FALSE)
  gammaThetaSample = plyr::aaply(1:ncol(nCountSpikes), 1, function(x) {
    fitData = data.frame(Xi=numberSpikes[,1], Ki=nCountSpikes[,x])
    fit1 <- stats::lm(Ki ~ 0 + Xi, data=fitData)
    gammaTheta = stats::coefficients(fit1)

    fitData = data.frame(Xi=numberSpikes[,1], Y=nCountSpikes[,x]>0, Ai=sizeFactorMatrix[,x])

    fit2 <- minpack.lm::nlsLM(Y ~ 1 - (1 - Theta + Theta*exp(-(gammaTheta/Theta)*Ai))^Xi,
                              data=fitData, start=list(Theta=0.4),
                              lower=0, upper=1, control=stats::nls.control(warnOnly=TRUE))
    c(coefficients(fit2), gammaTheta/coefficients(fit2))}, .expand=FALSE)
  gammaThetaSample[gammaThetaSample[,1]==1,1] = 0.9999

  abTheta = try(MASS::fitdistr(gammaThetaSample[,1], "beta", list(shape1=0.5, shape2=0.5)), silent=TRUE)
  if (class(abTheta) == "try-error") {
    VTheta = stats::var(gammaThetaSample[,1], na.rm=TRUE)
  } else {
    aTheta = abTheta$estimate[[1]]
    bTheta = abTheta$estimate[[2]]
    VTheta = (aTheta*bTheta) / ( (aTheta+bTheta)^2*(aTheta+bTheta+1))
  }
  abGamma = try(MASS::fitdistr(gammaThetaSample[,2], "gamma"), silent=TRUE)
  if (class(abGamma) == "try-error") {
    VGamma = stats::var(gammaThetaSample[,2], na.rm=TRUE)
  } else {
    aGamma = abGamma$estimate[[1]]
    bGamma = abGamma$estimate[[2]]
    VGamma = aGamma/bGamma^2
  }
  options(old.o)
  c(VGamma, VTheta)
}

# Estimate E[gamma], E[theta], V[Gamma] and V[Theta]
#' @importFrom minpack.lm nlsLM
.estimateEVGammaTheta <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  gammaThetaEstimate = .estimateGammaTheta(nCountSpikes, numberSpikes, sizeFactorMatrix)
  EGamma = gammaThetaEstimate$gammaTheta[[1]] / gammaThetaEstimate$Theta[[1]]
  ETheta = gammaThetaEstimate$Theta[[1]]
  E2Gamma = EGamma^2
  E2Theta = ETheta^2

  varianceGammaThetaEstimate = .estimateVGammaThetaInitial(nCountSpikes, numberSpikes, sizeFactorMatrix)
  VGamma = varianceGammaThetaEstimate[1]
  VTheta = varianceGammaThetaEstimate[2]

  fitData = data.frame(Xi=numberSpikes[,1], Y=apply(nCountSpikes,1,stats::var)/rowMeans(nCountSpikes), Bi=rowMeans(1/sizeFactorMatrix))

  fit <- minpack.lm::nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) +
                 (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
               start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  # second optimization for robust estimates
  fit <- minpack.lm::nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) +
                 (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
               start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  list(EGamma=EGamma, ETheta=ETheta, E2Gamma=E2Gamma, E2Theta=E2Theta, VGamma=VGamma, VTheta=VTheta)
}

.repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


# Simulating spike-in data
#' @importFrom stats rgamma rbeta rbinom lm coefficients
#' @importFrom minpack.lm nlsLM
.simulateCountGenes <- function(Xi, ETheta, VTheta, EGamma, VGamma, Aij, nCell, BV) {

  nGenes = length(Xi)

  if (VGamma > 0) {
    gammaShape = EGamma^2 / VGamma
    gammaScale = VGamma / EGamma
  }

  if (VTheta > 0) {
    thetaA = ETheta*( (ETheta*(1-ETheta))/VTheta - 1)
    thetaB = (1-ETheta)*( (ETheta*(1-ETheta))/VTheta - 1)
  }
  if (VGamma > 0) {
    gammaj = stats::rgamma(nCell*nGenes, shape=gammaShape, scale=gammaScale)
  } else {
    gammaj = rep(EGamma, nCell*nGenes)
  }

  if (VTheta > 0) {
    thetaj = tryCatch({
      stats::rbeta(nCell*nGenes, shape1=thetaA, shape2=thetaB)
    }, warning = function(war) {
      rep(ETheta, nCell*nGenes)
    })
  } else {
    thetaj = rep(ETheta, nCell*nGenes)
  }

  if (sum(BV) > 0) {
    Yi = tryCatch({
      stats::rgamma(nCell*nGenes, shape=Xi^2/BV, scale=BV/Xi)
    }, warning = function(war) {
      Yi = rep(Xi, nCell)
    })
  } else {
    Yi = rep(Xi, nCell)
  }
  Zij = stats::rbinom(nCell*nGenes, round(Yi), thetaj)
  Kij = stats::rpois(nCell*nGenes, gammaj*Aij*Zij)
  Kij = matrix(Kij, nGenes, nCell, byrow=F)
  Kij[Xi==0,] = 0

  # supply spike names and sample dummies
  rownames(Kij) = names(Xi)
  colnames(Kij) = paste0("S", 1:ncol(Kij))

  return(Kij)
}

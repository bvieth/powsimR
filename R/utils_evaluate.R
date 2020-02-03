
# Simulation Parameters ---------------------------------------------------

#' @importFrom MASS rlm
#' @importFrom stats residuals na.exclude
.lfc.evaluate <- function(truth, estimated) {

  # input
  SE <- ((truth - estimated)^2)
  AE <- abs(truth - estimated)
  RMSE <- sqrt(mean(SE, na.rm = T))
  MAE <- mean(AE, na.rm = T)
  AE.uniq <- unique(AE)

  # robust fitting and error of estimation
  fitted <- try(MASS::rlm(AE.uniq ~ 1, na.action = na.exclude), silent = T)
  if(inherits(fitted, 'try-error')){
    resids <- NA
    err.all <- NA
  } else{
    resids <- stats::residuals(fitted)
    err.all <- 2^sqrt(mean(resids^2, na.rm = T)) - 1
  }

  return(c(RMSE.Value=RMSE,
           MAE.Value=MAE,
           RMSE.NA=sum(is.na(SE))/length(SE),
           MAE.NA=sum(is.na(AE))/length(AE),
           ErrorFit=err.all))
}

#' @importFrom MASS rlm
#' @importFrom stats residuals na.exclude
.fiterror.sf <- function(estimated.sf, true.sf){

  # log fold change calculation
  logfold <- log2(estimated.sf) - log2(true.sf)

  # robust fitting and error of estimation
  fitted <- try(MASS::rlm(logfold ~ 1, na.action = na.exclude), silent = T)
  if(inherits(fitted, 'try-error')){
    resids <- NA
    err.all <- NA
  } else{
    resids <- stats::residuals(fitted)
    err.all <- 2^sqrt(mean(resids^2, na.rm = T)) - 1
  }

  return(err.all)
}

#' @importFrom stats aggregate
#' @importFrom tidyr %>%
.ratio.sf <- function(estimated.nsf, true.nsf, group) {
  dat <- data.frame(estimated.nsf, true.nsf, ratio=estimated.nsf/true.nsf, group=as.factor(group))
  dat.proc <- stats::aggregate(ratio ~ group, data=dat, FUN= mean)

  res <- data.frame(group1=as.numeric(dat.proc[1,2]), group2=as.numeric(dat.proc[2,2]))
  return(res)
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
  out <- list(TN=TN,TN.marginal=TN.mar,
              TP=TP, TP.marginal=TP.mar,
              FP=FP, FP.marginal=FP.mar,
              FN=FN, FN.marginal=FN.mar,
              TNR=TNR, TPR=TPR, FPR=FPR, FNR=FNR,FDR=FDR,
              TNR.marginal=TNR.marginal, TPR.marginal=TPR.marginal,
              FPR.marginal=FPR.marginal, FNR.marginal=FNR.marginal,
              FDR.marginal=FDR.marginal)

  return(out)
}

#' @importFrom dplyr mutate
#' @importFrom rlang .data
.roc.calc <- function(df){
  res <- dplyr::mutate(df,
                       FPR = .data$FP / (.data$FP + .data$TN),
                       TNR = .data$TN / (.data$TN + .data$FP),
                       FNR = .data$FN / (.data$FN + .data$TP),
                       PPV = .data$TP / (.data$TP + .data$FP),
                       ACC = (.data$TP + .data$TN) / (.data$TP + .data$TN + .data$FP + .data$FN),
                       F1score = (2*.data$TP) / ((2*.data$TP) + .data$FP + .data$FN),
                       MCC = ( (.data$TP * .data$TN) - (.data$FP * .data$FN) ) /
                         sqrt( (.data$TP + .data$FP)*(.data$TP + .data$FN)*
                                 (.data$TN + .data$FP)*(.data$TN + .data$FN) )
                )

  return(res)
}

#' @importFrom dplyr select group_by summarise_all
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom plotrix std.error
.scores.summary.calc <- function(calc.obj) {

  res <- do.call("rbind", calc.obj) %>%
    tidyr::pivot_longer(-c(.data$Samples, .data$Sim),
                        names_to = "Score") %>%
    dplyr::select(-.data$Sim) %>%
    dplyr::group_by(.data$Samples, .data$Score) %>%
    dplyr::summarise(Mean = mean(.data$value, na.rm = TRUE),
                     SE = plotrix::std.error(.data$value, na.rm = TRUE)) %>%
    data.frame()

  return(res)
}

#' @importFrom dplyr select group_by summarise_all
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom plotrix std.error
.tprvsfdr.summary.calc <- function(calc.obj) {

  res <- do.call("rbind", do.call("rbind", calc.obj)) %>%
    dplyr::select(-c(.data$Sim, .data$TP:.data$FN, .data$NBR)) %>%
    dplyr::group_by(.data$Samples, .data$Threshold) %>%
    dplyr::summarise_all(list(Mean = mean,
                              SE = plotrix::std.error)) %>%
    data.frame()

  return(res)
}

.perf.summary.calc <- function(calc.obj) {

  roc.obj.l <- sapply(names(calc.obj), function(j){
    data.frame(.threshold.avg(perf.obj = calc.obj[[j]], x.name = "FPR", y.name = "TPR"))
  }, USE.NAMES = TRUE, simplify = FALSE)

  roc.obj <- data.table::rbindlist(roc.obj.l, use.names = TRUE, idcol = 'Samples')

  pr.obj.l <- sapply(names(calc.obj), function(j){
    data.frame(.threshold.avg(perf.obj = calc.obj[[j]], x.name = "TPR", y.name = "PPV"))
  }, USE.NAMES = TRUE, simplify = FALSE)

  pr.obj <- data.table::rbindlist(pr.obj.l, use.names = TRUE, idcol = 'Samples')

  res <- list("ROC-Curve" = roc.obj, "PR-Curve" = pr.obj)

  return(res)

}

#' @importFrom stats approxfun
#' @importFrom matrixStats rowSds
.threshold.avg <- function(perf.obj, x.name, y.name) {

  alpha.values <- x.values <- y.values <- NULL

  alpha.tmp <- lapply(1:length(perf.obj), function(i) {
    x <- perf.obj[[i]][, "CUTOFF"]
    x[-1]
  })
  x.tmp <- lapply(1:length(perf.obj), function(i) {
    x <- perf.obj[[i]][, x.name]
    x[-1]
  })
  y.tmp <- lapply(1:length(perf.obj), function(i) {
    x <- perf.obj[[i]][, y.name]
    x[-1]
  })

  alpha.values <- seq(min(unlist(alpha.tmp)),
                      max(unlist(alpha.tmp)),
                      length=max( sapply(alpha.tmp, length)))

  for (i in 1:length(y.tmp)) {
    x.values[[i]] <- stats::approxfun(alpha.tmp[[i]],
                                      x.tmp[[i]],
                                      rule=2, ties=mean)(alpha.values)
    y.values[[i]] <- stats::approxfun(alpha.tmp[[i]],
                                      y.tmp[[i]],
                                      rule=2, ties=mean)(alpha.values)
  }

  res <- matrix(data = cbind(alpha.values,
                             rowMeans(data.frame(x.values)),
                             matrixStats::rowSds(as.matrix(data.frame(x.values))) / sqrt(length(x.values)),
                             rowMeans(data.frame(y.values)) ,
                             matrixStats::rowSds(as.matrix(data.frame(y.values))) / sqrt(length(y.values))),
                nrow = length(alpha.values), ncol = 5,
                dimnames = list(NULL, c("CUTOFF",
                                        paste(rep(x.name, 2), c("Mean", "SE"), sep ="_"),
                                        paste(rep(y.name, 2), c("Mean", "SE"), sep ="_"))
                                )
                )

  return(res)

}

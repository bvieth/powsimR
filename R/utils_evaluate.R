# Simulation Parameters ---------------------------------------------------

#' @importFrom MASS rlm
#' @importFrom stats residuals na.exclude
.lfc.evaluate <- function(truth, estimated) {
  SE <- ((truth - estimated)^2)
  AE <- abs(truth - estimated)
  RMSE <- sqrt(mean(SE, na.rm=T))
  MAE <- mean(AE, na.rm=T)

  folddiff <- abs(truth - estimated)
  # robust fitting
  fitted <- MASS::rlm(folddiff ~ 1, na.action=na.exclude)
  resids <- stats::residuals(fitted)
  # error of estimation
  err.all <- 2^sqrt(mean(resids^2, na.rm = T))-1

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

    # robust fitting
    fitted <- MASS::rlm(logfold ~ 1, na.action=na.exclude)
    resids <- stats::residuals(fitted)

    # error of estimation
    err.all <- 2^sqrt(mean(resids^2,na.rm = T))-1

  return(err.all)
}

#' @importFrom dplyr mutate group_by summarise
#' @importFrom tidyr %>%
.ratio.sf <- function(estimated.nsf, true.nsf, group) {
  dat <- data.frame(estimated.nsf, true.nsf, group)
  dat.proc <- dat %>%
    dplyr::mutate(ratio=estimated.nsf/true.nsf)  %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(avg = mean(ratio))

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
  list(TN=TN,TN.marginal=TN.mar, TP=TP, TP.marginal=TP.mar, FP=FP, FP.marginal=FP.mar, FN=FN, FN.marginal=FN.mar,
       TNR=TNR, TPR=TPR, FPR=FPR, FNR=FNR,FDR=FDR,
       TNR.marginal=TNR.marginal, TPR.marginal=TPR.marginal,
       FPR.marginal=FPR.marginal, FNR.marginal=FNR.marginal,
       FDR.marginal=FDR.marginal)
}

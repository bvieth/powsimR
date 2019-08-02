
# printEvalRes ------------------------------------------------------------

#' @name printEvalDE
#' @aliases printEvalDE
#' @title Summary table of power assessment
#' @description This function takes as input a result object from \code{\link{evaluateDE}}
#' and prints out a table to summarize important error-rates-related quantities.
#' The results are marginalized, meaning that they are averaged quantities
#' over all strata and simulations.
#' This provides a quick view of the marginal results per sample size.
#' @usage printEvalDE(evalRes)
#' @param evalRes The result object from \code{\link{evaluateDE}}.
#' @return A matrix of results per sample size considered (rows). Columns include sample size, specified nomial type I control value (for FDR or p-values), actual error rate, marginal TPR, averaged number of true and false discoveries, and false discovery costs.
#' @author Beate Vieth
#' @seealso \code{\link{simulateDE}}, \code{\link{evaluateDE}}
#' @examples
#' \dontrun{
#' ## for example evaluation result see \code{\link{evaluateDE}}
#' printEvalRes(evalRes=evalres)
#' }
#' @rdname printEvalDE
#' @export
printEvalDE <- function(evalRes) {
  nreps1 <- evalRes$n1
  nreps2 <- evalRes$n2
  alpha.type <- evalRes$alpha.type

  if(alpha.type == "raw") {
    alpha.nam <- "FPR"
    alpha.mar <- rowMeans(evalRes$FPR.marginal, na.rm=TRUE)
  }
  if(alpha.type == "adjusted") {
    alpha.nam <- "FDR"
    alpha.mar <- rowMeans(evalRes$FDR.marginal, na.rm=TRUE)
  }

  res <- matrix(0, nrow=length(nreps1), ncol=5)
  colnames(res) <- c("Sample size group 1", "Sample size group 2", paste(c("Nominal", "Marginal"), alpha.nam),
                     "Marginal TPR")

  res[,1] <- nreps1
  res[,2] <- nreps2
  res[,3] <- evalRes$alpha.nominal
  res[,4] <- alpha.mar
  res[,5] <- rowMeans(evalRes$TPR.marginal, na.rm=TRUE)

  res[, c(4:5)] <- signif(res[,c(4:5)],2)
  print(res)
  return(invisible(res))
}




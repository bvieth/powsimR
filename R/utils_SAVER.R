
# SAVER combine output ----------------------------------------------------

.combine.saver <- function (saver.list)
{
  est <- do.call(rbind, lapply(saver.list, `[[`, 1))
  se <- do.call(rbind, lapply(saver.list, `[[`, 2))
  info <- vector("list", 10)
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff", "lambda.coefs",
                   "total.time")
  info[[1]] <- saver.list[[1]]$info$size.factor
  info.list <- lapply(saver.list, `[[`, 3)
  for (i in 2:9) {
    info[[i]] <- do.call(c, lapply(info.list, `[[`, i))
  }
  info[[10]] <- do.call(sum, (lapply(info.list, `[[`, 10)))
  out <- list(estimate = est, se = se, info = info)
  class(out) <- "saver"
  out
}

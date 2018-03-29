.FindVarGenes <- function (object,
                           mean.function,
                           dispersion.function,
                           set.var.genes = TRUE,
                           x.low.cutoff = 0.1,
                           x.high.cutoff = Inf,
                           y.cutoff = 1,
                           y.high.cutoff = Inf,
                           num.bin = 20,
                           sort.results = TRUE)
{
  data <- object@data
  genes.use <- rownames(x = object@data)

  gene.mean <- rep(x = 0, length(x = genes.use))
  names(x = gene.mean) <- genes.use
  gene.dispersion <- gene.mean
  gene.dispersion.scaled <- gene.mean
  bin.size <- 1000
  max.bin <- floor(x = length(x = genes.use)/bin.size) + 1
  for (i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size *
                                        i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    genes.iter <- genes.use[my.inds]
    data.iter <- data[genes.iter, , drop = F]
    gene.mean[genes.iter] <- apply(X = data.iter,
                                   MARGIN = 1, FUN = mean.function)
    gene.dispersion[genes.iter] <- apply(X = data.iter,
                                         MARGIN = 1, FUN = dispersion.function)
  }

  gene.dispersion[is.na(x = gene.dispersion)] <- 0
  gene.mean[is.na(x = gene.mean)] <- 0
  gene.dispersion[is.infinite(x = gene.dispersion)] <- 0
  gene.mean[is.infinite(x = gene.mean)] <- 0
  data_x_bin <- cut(x = gene.mean, breaks = num.bin)
  names(x = data_x_bin) <- names(x = gene.mean)
  mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                   FUN = mean)
  sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                 FUN = sd)
  gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
  gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
  names(x = gene.dispersion.scaled) <- names(x = gene.mean)
  mv.df <- data.frame(gene.mean, gene.dispersion, gene.dispersion.scaled)
  rownames(x = mv.df) <- rownames(x = data)
  object@hvg.info <- mv.df

  gene.mean <- object@hvg.info[, 1]
  gene.dispersion <- object@hvg.info[, 2]
  gene.dispersion.scaled <- object@hvg.info[, 3]
  names(x = gene.mean) <- names(x = gene.dispersion) <- names(x = gene.dispersion.scaled) <- rownames(x = object@data)
  x.high.cutoff <- median(gene.mean)
  pass.cutoff <- names(x = gene.mean)[which(x = ((gene.mean > x.low.cutoff) &
                                                   (gene.mean < x.high.cutoff)) &
                                              (gene.dispersion.scaled > y.cutoff) &
                                              (gene.dispersion.scaled < y.high.cutoff))]

  if (set.var.genes) {
    object@var.genes <- pass.cutoff
    if (sort.results) {
      object@hvg.info <- object@hvg.info[order(object@hvg.info$gene.dispersion,
                                               decreasing = TRUE), ]
    }
  }
  return(object)
}


# # ORIGINAL ----------------------------------------------------------------
#
# Seurat::FindVariableGenes
# function (object, mean.function = ExpMean, dispersion.function = LogVMR,
#           do.plot = TRUE, set.var.genes = TRUE, x.low.cutoff = 0.1,
#           x.high.cutoff = 8, y.cutoff = 1, y.high.cutoff = Inf, num.bin = 20,
#           do.recalc = TRUE, sort.results = TRUE, do.cpp = TRUE, display.progress = TRUE,
#           ...)
# {
#   parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FindVariableGenes"))]
#   parameters.to.store$mean.function <- as.character(substitute(mean.function))
#   parameters.to.store$dispersion.function <- as.character(substitute(dispersion.function))
#   object <- SetCalcParams(object = object, calculation = "FindVariableGenes",
#                           ... = parameters.to.store)
#   data <- object@data
#   genes.use <- rownames(x = object@data)
#   if (do.recalc) {
#     if (do.cpp) {
#       if (!identical(mean.function, ExpMean)) {
#         warning("No equivalent mean.function implemented in c++ yet, falling back to R version")
#         do.cpp <- FALSE
#       }
#       if (!identical(dispersion.function, LogVMR)) {
#         warning("No equivalent dispersion.function implemented in c++ yet, falling back to R version")
#         do.cpp <- FALSE
#       }
#     }
#     if (do.cpp) {
#       if (class(data) != "dgCMatrix") {
#         data <- as(as.matrix(data), "dgCMatrix")
#       }
#       gene.mean <- FastExpMean(data, display.progress)
#       names(gene.mean) <- genes.use
#       gene.dispersion <- FastLogVMR(data, display.progress)
#       names(gene.dispersion) <- genes.use
#     }
#     if (!do.cpp) {
#       gene.mean <- rep(x = 0, length(x = genes.use))
#       names(x = gene.mean) <- genes.use
#       gene.dispersion <- gene.mean
#       gene.dispersion.scaled <- gene.mean
#       bin.size <- 1000
#       max.bin <- floor(x = length(x = genes.use)/bin.size) +
#         1
#       if (display.progress) {
#         print("Calculating gene dispersion")
#         pb <- txtProgressBar(min = 0, max = max.bin,
#                              style = 3)
#       }
#       for (i in 1:max.bin) {
#         my.inds <- ((bin.size * (i - 1)):(bin.size *
#                                             i - 1)) + 1
#         my.inds <- my.inds[my.inds <= length(x = genes.use)]
#         genes.iter <- genes.use[my.inds]
#         data.iter <- data[genes.iter, , drop = F]
#         gene.mean[genes.iter] <- apply(X = data.iter,
#                                        MARGIN = 1, FUN = mean.function)
#         gene.dispersion[genes.iter] <- apply(X = data.iter,
#                                              MARGIN = 1, FUN = dispersion.function)
#         if (display.progress) {
#           setTxtProgressBar(pb = pb, value = i)
#         }
#       }
#       if (display.progress) {
#         close(con = pb)
#       }
#     }
#     gene.dispersion[is.na(x = gene.dispersion)] <- 0
#     gene.mean[is.na(x = gene.mean)] <- 0
#     data_x_bin <- cut(x = gene.mean, breaks = num.bin)
#     names(x = data_x_bin) <- names(x = gene.mean)
#     mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
#                      FUN = mean)
#     sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
#                    FUN = sd)
#     gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
#     gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
#     names(x = gene.dispersion.scaled) <- names(x = gene.mean)
#     mv.df <- data.frame(gene.mean, gene.dispersion, gene.dispersion.scaled)
#     rownames(x = mv.df) <- rownames(x = data)
#     object@hvg.info <- mv.df
#   }
#   gene.mean <- object@hvg.info[, 1]
#   gene.dispersion <- object@hvg.info[, 2]
#   gene.dispersion.scaled <- object@hvg.info[, 3]
#   names(x = gene.mean) <- names(x = gene.dispersion) <- names(x = gene.dispersion.scaled) <- rownames(x = object@data)
#   pass.cutoff <- names(x = gene.mean)[which(x = ((gene.mean >
#                                                     x.low.cutoff) & (gene.mean < x.high.cutoff)) & (gene.dispersion.scaled >
#                                                                                                       y.cutoff) & (gene.dispersion.scaled < y.high.cutoff))]
#   if (do.plot) {
#     VariableGenePlot(object = object, x.low.cutoff = x.low.cutoff,
#                      x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff,
#                      y.high.cutoff = y.high.cutoff, ...)
#   }
#   if (set.var.genes) {
#     object@var.genes <- pass.cutoff
#     if (sort.results) {
#       object@hvg.info <- object@hvg.info[order(object@hvg.info$gene.dispersion,
#                                                decreasing = TRUE), ]
#     }
#     return(object)
#   }
#   else {
#     return(pass.cutoff)
#   }
# }

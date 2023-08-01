# remove scientific notation ----------------------------------------------

.plain <- function(x,...) {
  format(x, ..., scientific = FALSE, trim = TRUE, digits=1)
}


# ggplot colors -----------------------------------------------------------
#' @importFrom grDevices hcl
.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

# Themes ------------------------------------------------------------------

#' @importFrom ggplot2 theme theme_light theme_classic theme_linedraw element_blank element_text element_rect
#' @importFrom grid unit

.theme_param_qc <- function() {
  ggplot2::theme_light() +
  ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size=10, color='black'),
                 axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                 plot.title = ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white",
                                                          colour = NA))
}

.theme_param_margs <- function() {
  ggplot2::theme_light() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "none",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x = ggplot2::element_text(size=10, color='black'),
                   axis.text.y = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

.theme_param_fit <- function(){
  ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none',
                   axis.text = ggplot2::element_text(size=10),
                   axis.title = ggplot2::element_text(size=12, face="bold"))
}

.theme_eval_de <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "right",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x = ggplot2::element_text(size=10, color='black', angle=45, hjust=1),
                   axis.text.y = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

.theme_eval_sim <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x = ggplot2::element_text(size=10, color='black', angle=45, hjust=1),
                   axis.text.y = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

.theme_eval_time <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1, "lines"),
                   axis.text = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

.theme_eval_dist <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1, "lines"),
                   axis.text = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

.theme_eval_roc <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   legend.key.size = grid::unit(1, "lines"),
                   axis.text = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}

# PlotParam ---------------------------------------------------------------
#' @importFrom ggplot2 ggplot aes_ geom_boxplot geom_point position_jitterdodge scale_fill_manual labs theme element_text element_blank element_rect geom_violin geom_dotplot stat_summary scale_y_continuous scale_y_log10 geom_hline geom_bar facet_grid as_labeller scale_x_discrete stat_density2d scale_fill_gradientn geom_line
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom grid unit
#' @importFrom dplyr bind_rows
#' @importFrom ggpubr ggtexttable ttheme
#' @importFrom grDevices blues9 colorRampPalette
#' @importFrom reshape2 melt
#' @importFrom stats reorder
# sequencing depth with marker for dropout samples
.seqdepth_plot <- function(estParamRes){

  .x = NULL
  lib.size.dat <- data.frame(Seqdepth=estParamRes$Parameters$Raw$seqDepth,
                             Sample=names(estParamRes$Parameters$Raw$seqDepth),
                             Dropout=estParamRes$DropOuts$Sample$totCounts,
                             stringsAsFactors = F)

  if(attr(estParamRes, "RNAseq")=="singlecell" | estParamRes$detectS>12) {

    seqdepth.plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = lib.size.dat,
                          ggplot2::aes_(x = 1, y = quote(Seqdepth), fill=quote(Dropout)),
                          pch = 21,
                          position = ggplot2::position_jitterdodge(dodge.width = 0.5)) +
      ggplot2::geom_boxplot(data = lib.size.dat,
                            ggplot2::aes_(x = 1, y = quote(Seqdepth)),
                            outlier.shape = NA, width = 0.5, alpha=0.5) +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Sequencing Depth") +
      .theme_param_qc()
  }
  if(attr(estParamRes, "RNAseq")=="bulk" | estParamRes$detectS<12) {
    lib.size.dat <- lib.size.dat[order(lib.size.dat$Seqdepth, decreasing = T), ]
    lib.size.dat$Sample <- factor(lib.size.dat$Sample,
                                  levels = lib.size.dat$Sample,
                                  ordered = TRUE)

    seqdepth.plot <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = lib.size.dat,
                        ggplot2::aes_(x = quote(Sample),
                                      y = quote(Seqdepth),
                                      fill = quote(Dropout)),
                        stat="identity", width=.5) +
      ggplot2::geom_hline(yintercept = stats::median(lib.size.dat$Seqdepth),
                          linetype = 2, colour="grey40") +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::scale_x_discrete(limits = rev(levels(lib.size.dat$Sample))) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Sequencing Depth") +
      .theme_param_qc() +
      ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, color='black')) +
      ggplot2::coord_flip()
  }

  # return plot
  return(seqdepth.plot)
}
# library size factor plot
.sf_plot <- function(estParamRes){

  sf.dat <- data.frame(SizeFactor=estParamRes$sf,
                       Sample=names(estParamRes$sf),
                       stringsAsFactors = FALSE)
  sf.max <- max(sf.dat$SizeFactor)*1.05
  if(attr(estParamRes, "RNAseq")=="singlecell" | estParamRes$detectS>12) {
    sf.plot <- ggplot2::ggplot() +
      ggplot2::geom_violin(data = sf.dat,
                           ggplot2::aes_(x = 1, y=quote(SizeFactor)),
                           fill = "grey90", width = 0.8, color = "black") +
      ggplot2::stat_summary(data = sf.dat,
                            ggplot2::aes_(x = 1, y=quote(SizeFactor)),
                            fun = stats::median,
                            fun.min = stats::median,
                            fun.max = stats::median,
                            color = "black",
                            width = 0.5,
                            geom = "crossbar") +
      ggplot2::scale_y_continuous(limits = c(0, sf.max)) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = paste0("Library Size Factors (", estParamRes$normFramework, ")")) +
      .theme_param_qc()
  }
  if(attr(estParamRes, "RNAseq")=="bulk" | estParamRes$detectS<12) {
    sf.plot <- ggplot2::ggplot() +
      ggplot2::geom_dotplot(data = sf.dat,
                            ggplot2::aes_(x = 1, y=quote(SizeFactor)),
                            binaxis='y',
                            stackdir='center',
                            dotsize=1,
                            fill = "grey75") +
      ggplot2::stat_summary(data = sf.dat,
                            ggplot2::aes_(x = 1, y=quote(SizeFactor)),
                            fun = stats::median,
                            fun.min = stats::median,
                            fun.max = stats::median,
                            color = "black",
                            width = 0.5,
                            geom = "crossbar") +
      ggplot2::scale_y_continuous(limits = c(0, sf.max)) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = paste0("Library Size Factors (", estParamRes$normFramework, ")")) +
      .theme_param_qc()
  }

  # return plot
  return(sf.plot)
}
# total gene features
.feat_plot <- function(estParamRes){
  .x = NULL

  totfeatures.dat <- data.frame(TotFeatures=estParamRes$Parameters$Raw$totFeatures,
                                Sample=names(estParamRes$Parameters$Raw$totFeatures),
                                Dropout=estParamRes$DropOuts$Sample$totFeatures,
                                stringsAsFactors = F)

  if(attr(estParamRes, "RNAseq")=="singlecell" | estParamRes$detectS>12) {
    totfeatures.plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = totfeatures.dat,
                          ggplot2::aes_(x = 1, y = quote(TotFeatures), fill=quote(Dropout)),
                          pch = 21,
                          position = ggplot2::position_jitterdodge(dodge.width = 0.5)) +
      ggplot2::geom_boxplot(data = totfeatures.dat,
                            ggplot2::aes_(x = 1, y = quote(TotFeatures)),
                            outlier.shape = NA, width = 0.5, alpha=0.5) +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Detected Genes") +
      .theme_param_qc()
  }
  if(attr(estParamRes, "RNAseq")=="bulk" | estParamRes$detectS<12) {
    totfeatures.dat <- totfeatures.dat[order(totfeatures.dat$TotFeatures, decreasing = T), ]
    totfeatures.dat$Sample <- factor(totfeatures.dat$Sample,
                                     levels = totfeatures.dat$Sample,
                                     ordered = TRUE)

    totfeatures.plot <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = totfeatures.dat,
                        ggplot2::aes_(x = quote(Sample),
                                      y = quote(TotFeatures),
                                      fill = quote(Dropout)),
                        stat="identity", width=.5) +
      ggplot2::geom_hline(yintercept =stats::median(totfeatures.dat$TotFeatures),
                          linetype = 2, colour="grey40") +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::scale_x_discrete(limits = rev(levels(totfeatures.dat$Sample))) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Sequencing Depth") +
      .theme_param_qc() +
      ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, color='black')) +
      ggplot2::coord_flip()
  }

  # return plot
  return(totfeatures.plot)
}
# gene spike ratio
.genespike_ratio_plot <- function(estParamRes){
  genespike.dat <- data.frame(Ratio=estParamRes$DropOuts$Sample$GeneSpikeRatio/100,
                              Sample=names(estParamRes$Parameters$Raw$totFeatures),
                              Dropout=estParamRes$DropOuts$Sample$GeneSpike,
                              stringsAsFactors = F)
  gs.max <- max(genespike.dat$Ratio)*1.05

  if(attr(estParamRes, "RNAseq")=="singlecell" | estParamRes$detectS>12) {

    genespike.plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = genespike.dat,
                          ggplot2::aes_(x=1, y=quote(Ratio), fill=quote(Dropout)),
                          pch = 21,
                          position = ggplot2::position_jitterdodge(dodge.width = 0.5)) +
      ggplot2::geom_boxplot(data = genespike.dat,
                            ggplot2::aes_(x=1, y=quote(Ratio)),
                            outlier.shape = NA, width = 0.5, alpha=0.5) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::scale_y_continuous(limits = c(0, gs.max),
                                  labels = scales::percent_format()) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Gene Spike Count Ratio") +
      .theme_param_qc()
  }

  if(attr(estParamRes, "RNAseq")=="bulk" | estParamRes$detectS<12) {
    genespike.plot <- ggplot2::ggplot() +
      ggplot2::geom_dotplot(data = genespike.dat,
                            ggplot2::aes_(x = 1, y=quote(Ratio), fill=quote(Dropout)),
                            binaxis='y',
                            stackdir='center',
                            dotsize=0.75) +
      ggplot2::stat_summary(data = genespike.dat,
                            ggplot2::aes_(x = 1, y=quote(Ratio)),
                            fun = stats::median,
                            fun.min = stats::median,
                            fun.max = stats::median,
                            color = "black",
                            width = 0.5,
                            geom = "crossbar") +
      ggplot2::scale_y_continuous(limits = c(0, gs.max),
                                  labels = scales::percent_format()) +
      ggplot2::scale_fill_manual(values = c("grey75", "red"),
                                 labels = c("Included", "Outlier")) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = "Gene Spike Count Ratio") +
      .theme_param_qc()
  }

  # return plot
  return(genespike.plot)
}
# marginal distributions for mean, dispersion, dropout
.margs_plot <- function(estParamRes, Distribution){
  param.names <- names(estParamRes$Parameters)[!is.na(estParamRes$Parameters)]
  set.relabels <- c(`Raw` = "Provided",
                    `Full` = "All Genes",
                    `Filtered`="Filtered Genes",
                    `DropGene` = "Dropout Genes",
                    `DropSample` = ifelse(attr(estParamRes, "RNAseq")=="singlecell",
                                          "Cell Outliers", "Sample Outliers"))
  if(Distribution == "NB"){
    param.relabels <- c(`Mean` = "Log Mean",
                        `Dispersion` = "Log Dispersion",
                        `Dropout` = "Gene Dropout Rate")
    margs.L <- sapply(param.names, function(i){
      tmp <- estParamRes$Parameters[[i]]
      data.frame(Mean=log2(tmp$means+1),
                 Dispersion=log2(tmp$dispersion+1),
                 Dropout=tmp$gene.dropout)
    }, simplify = F, USE.NAMES = T)
  }
  if(Distribution == "ZINB"){
    param.relabels <- c(`Mean` = "Log Positive Mean",
                        `Dispersion` = "Log Positive Dispersion",
                        `Dropout` = "Gene Dropout Rate")
    margs.L <- sapply(param.names, function(i){
      tmp <- estParamRes$Parameters[[i]]
      data.frame(Mean=log2(tmp$pos.means+1),
                 Dispersion=log2(tmp$pos.dispersion+1),
                 Dropout=tmp$gene.dropout)
    }, simplify = F, USE.NAMES = T)
  }
  margs.dat <- dplyr::bind_rows(margs.L, .id = "Set")
  margs.dat <- suppressMessages(reshape2::melt(margs.dat))

  margs.plot <- ggplot2::ggplot(margs.dat,
                                ggplot2::aes_(x=quote(Set), y=quote(value))) +
    ggplot2::geom_violin(fill = "#597EB5", alpha = 0.5) +
    ggplot2::stat_summary(fun = stats::median,
                          fun.min = stats::median,
                          fun.max = stats::median,
                          color = "black",
                          width = 0.5,
                          geom = "crossbar") +
    ggplot2::facet_grid(~variable,
                        scales = "free_x",
                        labeller = ggplot2::as_labeller(param.relabels)) +
    ggplot2::scale_x_discrete(labels = set.relabels) +
    ggplot2::labs(x = NULL,
                  y = NULL) +
    .theme_param_margs() +
    ggplot2::coord_flip()

  return(margs.plot)

}
# table with numbers of genes and samples
.estimate_table_print <- function(estParamRes, RNAseq){
  sample.names = ifelse(RNAseq=="singlecell", "Single Cells", "Bulk Samples")
  no.dat <- data.frame(c("Provided", "Detected", "All Genes", "Filtered Genes", "Dropout Genes"),
                       c(estParamRes$totalG,
                         estParamRes$detectG,
                         estParamRes$Parameters$Full$ngenes,
                         estParamRes$Parameters$Filtered$ngenes,
                         ifelse(all(!is.na(estParamRes$Parameters$DropGene)),
                                estParamRes$Parameters$DropGene$ngenes, 0)),
                       c(NA,
                         NA,
                         estParamRes$Fit$Full$estG,
                         estParamRes$Fit$Filtered$estG,
                         ifelse(all(!is.na(estParamRes$Fit$DropGene)),
                                estParamRes$Fit$DropGene$estG, 0)),
                       c(estParamRes$totalS,
                         estParamRes$detectS,
                         estParamRes$Parameters$Full$nsamples,
                         estParamRes$Parameters$Filtered$nsamples,
                         ifelse(all(!is.na(estParamRes$Parameters$DropGene)),
                                estParamRes$Parameters$DropGene$nsamples, NA)),
                       stringsAsFactors = F )
  colnames(no.dat) <- c("Set", "# Genes", "# Genes for Fit", paste0("# ", sample.names))
  no.table <- ggpubr::ggtexttable(no.dat,
                                  rows = NULL,
                                  theme = ggpubr::ttheme("mOrange"))

  return(no.table)
}
# fitting lines
.meanvsdisp_plot <- function(estParamRes, Distribution){
  ..density.. = NULL

  meanvsdisp.dat <- data.frame(Mean=estParamRes$Fit$Filtered$meandispfit$model$x[,"x"],
                               Dispersion=estParamRes$Fit$Filtered$meandispfit$model$y)
  meanvsdisp.fdat <- data.frame(Mean=estParamRes$Fit$Filtered$meandispfit$x,
                                Dispersion=estParamRes$Fit$Filtered$meandispfit$y,
                                Upper=estParamRes$Fit$Filtered$meandispfit$upper,
                                Lower=estParamRes$Fit$Filtered$meandispfit$lower)


  if(Distribution == "NB"){
    mean.name <- "Log Mean"
    disp.name <- "Log Dispersion"
    cdisp <- log2(estParamRes$Parameters$Filtered$common.dispersion+1)
  }
  if(Distribution == "ZINB"){
    mean.name <- "Log Positive Mean"
    disp.name <- "Log Positive Dispersion"
    cdisp <- log2(estParamRes$Parameters$Filtered$pos.common.dispersion+1)
  }

  meanvsdisp.plot <- ggplot2::ggplot(data=meanvsdisp.dat,
                                     ggplot2::aes_(x=quote(Mean), y=quote(Dispersion))) +
    ggplot2::geom_point(size=0.5) +
    ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                                      alpha=ifelse(..density..^0.15<0.4,0,1)),
                            contour=FALSE) +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
    ggplot2::geom_hline(yintercept = cdisp,
                        linetype = 2, colour="grey40") +
    ggplot2::geom_line(data = meanvsdisp.fdat,
                       ggplot2::aes_(x = quote(Mean), y = quote(Dispersion)),
                       linetype=1, size=1.5, colour="orange") +
    ggplot2::geom_line(data = meanvsdisp.fdat,
                       ggplot2::aes_(x = quote(Mean), y = quote(Lower)),
                       linetype=2, size=1, colour="orange") +
    ggplot2::geom_line(data = meanvsdisp.fdat,
                       ggplot2::aes_(x = quote(Mean), y = quote(Upper)),
                       linetype=2, size=1, colour="orange") +
    ggplot2::labs(y=disp.name,
                  x=mean.name) +
    .theme_param_fit()

  return(meanvsdisp.plot)
}
.meanvsdrop_plot <- function(estParamRes, Distribution, RNAseq){
  ..density.. = NULL
  if(Distribution == "NB"){
    meanvsp0.dat <- data.frame(Mean=log2(estParamRes$Parameters$Filtered$means+1),
                               Dropout=estParamRes$Parameters$Filtered$gene.dropout)
    mean.name <- "Log Mean"
  }
  if(Distribution == "ZINB"){
    meanvsp0fit.dat <- data.frame(Mean=estParamRes$Fit$Filtered$meang0fit$x,
                                  Dropout=estParamRes$Fit$Filtered$meang0fit$y)
    meanvsp0.dat <- data.frame(Mean=log2(estParamRes$Parameters$Filtered$pos.means+1),
                               Dropout=estParamRes$Parameters$Filtered$gene.dropout)
    mean.name <- "Log Positive Mean"
  }

  meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat,
                                   ggplot2::aes_(x=quote(Mean), y=quote(Dropout))) +
    ggplot2::geom_point(size=0.5) +
    ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                                      alpha=ifelse(..density..^0.15<0.4,0,1)),
                            contour=FALSE) +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::labs(y="Gene Dropout Rate",
                  x=mean.name) +
    .theme_param_fit()

  if(Distribution == "ZINB"){
    meanvsp0.plot <- meanvsp0.plot +
      ggplot2::geom_vline(xintercept = estParamRes$Fit$Filtered$g0.cut,
                          linetype = 2, colour="grey40") +
      ggplot2::geom_line(data = meanvsp0fit.dat,
                         ggplot2::aes_(x=quote(Mean),
                                       y=quote(Dropout)),
                         linetype=1,
                         size=1.5,
                         colour="orange")
  }

  if(RNAseq == 'bulk'){
    meanvsp0.plot <- meanvsp0.plot +
      ggplot2::geom_vline(xintercept = estParamRes$Fit$Filtered$g0.cut,
                          linetype = 2, colour="grey40")
  }

  return(meanvsp0.plot)

}
# read vs umi
.readumi_plot <- function(estParamRes){
  ..density.. = NULL
  readvsumifit.dat <- data.frame(X=estParamRes$Fit$UmiRead$Fit$x,
                                 Y=estParamRes$Fit$UmiRead$Fit$y,
                                 Lower=estParamRes$Fit$UmiRead$Fit$lower,
                                 Upper=estParamRes$Fit$UmiRead$Fit$upper)
  readvsumi.dat <- data.frame(UMI=estParamRes$Fit$UmiRead$lUMI,
                              Ratio=estParamRes$Fit$UmiRead$lRatio)
  readvsumi.plot <- ggplot2::ggplot(data=readvsumi.dat,
                                    ggplot2::aes_(x=quote(UMI), y=quote(Ratio))) +
    ggplot2::geom_point(size=0.5) +
    ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                                      alpha=ifelse(..density..^0.15<0.4,0,1)),
                            contour=FALSE) +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
    ggplot2::geom_line(data = readvsumifit.dat,
                       ggplot2::aes_(x=quote(X),
                                     y=quote(Y)),
                       linetype=1, size=1.5, colour="orange") +
    ggplot2::geom_line(data = readvsumifit.dat,
                       ggplot2::aes_(x=quote(X),
                                    y=quote(Upper)),
                       linetype=2, size=1, colour="orange") +
    ggplot2::geom_line(data = readvsumifit.dat,
                       ggplot2::aes_(x=quote(X),
                                    y=quote(Lower)),
                       linetype=2, size=1, colour="orange") +
    ggplot2::labs(y=expression(bold(paste("Amplification Rate"))),
                  x=expression(bold(paste(Log[10], " UMI")))) +
    .theme_param_fit()

  return(readvsumi.plot)
}


# PlotEvalROC -------------------------------------------------------------

#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ labs theme scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar theme_minimal expansion
#' @importFrom ggstance geom_linerangeh
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom dplyr group_by summarise ungroup n
#' @importFrom tidyr %>%
# roc curve
.roc_plot <- function(ROCData) {
  roc.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(FPR_Mean),
                                     y = quote(TPR_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         size = 0.75,
                         linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(x = "1 - Specificity (FPR)",
                  y = "Sensitivity (TPR)") +
    .theme_eval_roc()

  return(roc.plot)
}
# pr curve
.pr_plot <- function(ROCData) {
  pr.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(TPR_Mean),
                                     y = quote(PPV_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::scale_x_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(x = "Recall (TPR)", y = "Precision (PPV)") +
    .theme_eval_roc()

  return(pr.plot)
}
# tpr vs fdr curve
.tprvsfdr_plot <- function(ROCData, alpha.nominal) {
  tprvsfdr.nominal <- ROCData %>%
    dplyr::filter(.data$Threshold == alpha.nominal) %>%
    dplyr::mutate(`FDR Control` = ifelse(.data$FDR_Mean <= .data$Threshold, "yes", "no"),
                  FDRLower = .data$FDR_Mean - .data$FDR_SE,
                  FDRUpper = .data$FDR_Mean + .data$FDR_SE,
                  TPRLower = .data$TPR_Mean - .data$TPR_SE,
                  TPRUpper = .data$TPR_Mean + .data$TPR_SE)

  upperlimit <- ifelse(max(tprvsfdr.nominal$FDR_Mean) <= alpha.nominal,
                       max(tprvsfdr.nominal$FDR_Mean)*5,
                       max(tprvsfdr.nominal$FDR_Mean)*2)

  tprvsfdr.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(FDR_Mean),
                                     y = quote(TPR_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::geom_point(data = tprvsfdr.nominal,
                        ggplot2::aes_(x = quote(FDR_Mean),
                                      y = quote(TPR_Mean),
                                      group = quote(Samples),
                                      shape = quote(`FDR Control`)),
                        size = 3) +
    ggplot2::geom_linerange(data = tprvsfdr.nominal,
                            ggplot2::aes_(ymin=quote(TPRLower),
                                          ymax=quote(TPRUpper),
                                          x = quote(FDR_Mean),
                                          group = quote(Samples),
                                          color = quote(Samples)),
    ) +
    ggstance::geom_linerangeh(data = tprvsfdr.nominal,
                              ggplot2::aes_(xmin=quote(FDRLower),
                                            xmax=quote(FDRUpper),
                                            y = quote(TPR_Mean),
                                            group = quote(Samples),
                                            color = quote(Samples))) +
    ggplot2::scale_fill_manual(values = ) +
    ggplot2::scale_x_continuous(limits = c(0,upperlimit),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(y = "TPR",
                  x = "observed FDR") +
    ggplot2::guides() +
    .theme_eval_roc()

  return(tprvsfdr.plot)
}

.summary_table_print <- function(TblData, cutoff) {

  fixed.scores <- c("TPRvsFPR_AUC", "TPRvsPPV_AUC")
  select.scores <- paste(c("ACC", "MCC", "F1score", "TPRvsFDR_pAUC"), cutoff, sep="_")
  scores <- c(select.scores, fixed.scores)

  tbl.dat <- TblData %>%
    dplyr::filter(.data$Score %in% scores) %>%
    dplyr::mutate(`Summary Statistic` = gsub(pattern = paste0("_", cutoff), replacement = "", .data$Score)) %>%
    dplyr::mutate(`Summary Statistic` = gsub(pattern = "_", replacement = " ", .data$`Summary Statistic`)) %>%
    dplyr::mutate(`Summary Statistic` = case_when(.data$`Summary Statistic` == "TPRvsFPR AUC" ~ "ROC AUC",
                                                  .data$`Summary Statistic` == "TPRvsPPV AUC" ~ "PRC AUC",
                                                  TRUE ~ as.character(.data$`Summary Statistic`))) %>%
    tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
    dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2),
                  Mean = round(.data$Mean, digits = 2),
                  SE = round(.data$SE, digits = 2)) %>%
    dplyr::arrange(.data$`Summary Statistic`, .data$SumN) %>%
    tidyr::unite("Value", c("Mean", "SE"), sep = "\u00B1") %>%
    dplyr::select(.data$`Summary Statistic`, .data$Samples, .data$Value) %>%
    tidyr::pivot_wider(id = "Samples", names_from = "Summary Statistic", values_from = "Value")

  summary.tbl <- ggpubr::ggtexttable(tbl.dat,
                                     rows = NULL,
                                     theme = ggpubr::ttheme("lBlackWhite"))

  return(summary.tbl)
}

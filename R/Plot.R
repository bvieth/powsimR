
# plotParam -------------------------------------------------------------
#' @name plotParam
#' @aliases plotParam
#' @title Visualize distributional characteristics of RNAseq experiment
#' @description This function plots the results of the parameter estimation. This includes the absolute and relative sequencing depth (i.e. library size factor) as well as marginal log2(mean+1), log2(dispersion) and dropout. Furthermore, the mean-dispersion relationship with loess fit for simulations is visualized. Lastly, the mean-dropout rate is presented as a smooth scatter plot.
#' @usage plotParam(estParamRes, annot=TRUE)
#' @param estParamRes The output of \code{\link{estimateParam}}.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' plotParam(estParamRes = kolodziejczk_param, annot=TRUE)
#' }
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot aes geom_bar geom_line theme geom_hline labs coord_flip scale_y_continuous geom_point scale_fill_gradientn stat_density2d labs theme_minimal theme_classic element_text ylim
#' @importFrom grDevices blues9 colorRampPalette
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom stats reorder
#' @rdname plotParam
#' @export
plotParam <- function(estParamRes, annot=TRUE) {

  if(estParamRes$RNAseq=="bulk" | estParamRes$estS<15) {
    # library size
    lib.size.dat <- data.frame(Seqdepth=estParamRes$seqDepth,
                               Sample=names(estParamRes$seqDepth))
    libsize.plot <- ggplot2::ggplot(data=lib.size.dat, ggplot2::aes(reorder(Sample, Seqdepth),Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x=NULL, y="Sequencing depth") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estParamRes$sf,
                         Sample=names(estParamRes$sf))
    sf.plot <- ggplot2::ggplot(data=sf.dat, ggplot2::aes(reorder(Sample, SizeFactor),SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x=NULL, y="Library Size Factor") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
    # marginal distributions
    margs.dat <- data.frame(Mean=log2(estParamRes$means+1),
                            Dispersion=log2(estParamRes$dispersion),
                            Dropout=estParamRes$p0)
    margs.dat <- suppressMessages(reshape2::melt(margs.dat))
    margs.plot <- ggplot2::ggplot(margs.dat, ggplot2::aes(value)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=12, face="bold"),
                     strip.text = ggplot2::element_text(size=12, face="bold")) +
      ggplot2::labs(x=NULL, y="Density") +
      ggplot2::facet_wrap(~variable, scales = 'free')
    # mean vs dispersion plot
    meanvsdisp.dat <- data.frame(Means=log2(estParamRes$means+1),
                                 Dispersion=log2(estParamRes$dispersion))
    meanvsdisp.plot <- ggplot2::ggplot(data=meanvsdisp.dat, ggplot2::aes(x=Means, y=Dispersion)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
      ggplot2::geom_hline(yintercept = log2(estParamRes$common.dispersion), linetype = 2, colour="grey40") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$y),
                         linetype=1, size=1.5, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$upper),
                         linetype=2, size=1, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$lower),
                         linetype=2, size=1, colour="orange") +
      ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))),
                    x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold"))
    # mean vs p0 plot
    meanvsp0.dat <- data.frame(Means=log2(estParamRes$means+1),
                               Dropout=estParamRes$p0)
    meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat, ggplot2::aes(x=Means, y=Dropout)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                               alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256))+
      ggplot2::ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold"))
    if(estParamRes$RNAseq=="bulk" && !is.null(estParamRes$p0.cut)) {
      meanvsp0.plot <- meanvsp0.plot + ggplot2::annotate("rect", xmin=0, ymax=1,
                                                         ymin=0, xmax=estParamRes$p0.cut+1,
                                                         fill="red", alpha=0.2)
    }

    top_row <- suppressWarnings(cowplot::plot_grid(libsize.plot,sf.plot,
                                                   labels=c('A', 'B'),
                                                   ncol=2, nrow=1))
    bottom_row <- suppressWarnings(cowplot::plot_grid(meanvsdisp.plot,
                                                      meanvsp0.plot,
                                                      labels=c('D', 'E'),
                                                      ncol=2, nrow=1))
    middle_row <- suppressWarnings(cowplot::plot_grid(margs.plot, labels=c('C'),
                                                      ncol=1, nrow=1))
    p.final <- suppressWarnings(cowplot::plot_grid(top_row, middle_row,
                                                   bottom_row,
                                                   rel_heights = c(1, 1, 1.5),
                                                   ncol=1, nrow=3))
    # annotation under plot
    if (annot) {
      p.final <- cowplot::add_sub(p.final, "A) Sequencing depth per sample with median sequencing depth (grey dashed line).
                                  \nB) Library size normalisation factor per sample with median size factor (grey dashed line).
                                  \nC) Marginal Distribution of mean, dispersion and dropout.
                                  \nD) Local polynomial regression fit between mean and dispersion estimates with variability band per gene (yellow). \nCommon dispersion estimate (grey dashed line).
                                  \nE) Fraction of dropouts versus estimated mean expression per gene.", size=8)
    }
    }

  if(estParamRes$RNAseq=="singlecell" | estParamRes$estS>15) {
    # library size
    lib.size.dat <- data.frame(Seqdepth=estParamRes$seqDepth,
                               Sample=names(estParamRes$seqDepth))
    libsize.plot <- ggplot2::ggplot(lib.size.dat, ggplot2::aes(Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x="Sequencing Depth", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estParamRes$sf, Sample=names(estParamRes$sf))
    sf.plot <- ggplot2::ggplot(sf.dat, ggplot2::aes(SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x="Library Size Factor", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # marginal distributions
    margs.dat <- data.frame(Mean=log2(estParamRes$means+1),
                            Dispersion=log2(estParamRes$dispersion),
                            Dropout=estParamRes$p0)
    margs.dat <- suppressMessages(reshape2::melt(margs.dat))
    margs.plot <- ggplot2::ggplot(margs.dat, ggplot2::aes(value)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=12, face="bold"),
                     strip.text = ggplot2::element_text(size=12, face="bold")) +
      ggplot2::labs(x=NULL, y="Density") +
      ggplot2::facet_wrap(~variable, scales = 'free')
    # mean vs dispersion plot
    meanvsdisp.dat <- data.frame(Means=log2(estParamRes$means+1),
                                 Dispersion=log2(estParamRes$dispersion))
    meanvsdisp.plot <- ggplot2::ggplot(data=meanvsdisp.dat, ggplot2::aes(x=Means, y=Dispersion)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                               alpha=ifelse(..density..^0.15<0.4,0,1)),
                              contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
      ggplot2::geom_hline(yintercept = log2(estParamRes$common.dispersion),
                          linetype = 2, colour="grey40") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$y),
                         linetype=1, size=1.5, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$upper),
                         linetype=2, size=1, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$lower),
                         linetype=2, size=1, colour="orange") +
      ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))),
                    x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold"))
    # mean vs p0 plot
    meanvsp0.dat <- data.frame(Means=log2(estParamRes$means+1),
                               Dropout=estParamRes$p0)
    meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat, ggplot2::aes(x=Means, y=Dropout)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                               alpha=ifelse(..density..^0.15<0.4,0,1)),
                              contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(c("white", grDevices::blues9))(256)) +
      ggplot2::ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none',
                     axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold"))

    top_row <- suppressWarnings(cowplot::plot_grid(libsize.plot,
                                                   sf.plot,
                                                   labels=c('A', 'B'),
                                                   ncol=2, nrow=1))
    bottom_row <- suppressWarnings(cowplot::plot_grid(meanvsdisp.plot,
                                                      meanvsp0.plot,
                                                      labels=c('D', 'E'),
                                                      ncol=2, nrow=1))
    middle_row <- suppressWarnings(cowplot::plot_grid(margs.plot,
                                                      labels=c('C'),
                                                      ncol=1, nrow=1))
    p.final <- suppressWarnings(cowplot::plot_grid(top_row,
                                                   middle_row,
                                                   bottom_row,
                                                   rel_heights = c(1, 1, 1.5),
                                                   ncol=1, nrow=3))
    # annotation under plot
    if (annot) {
      p.final <- cowplot::add_sub(p.final, "A) Sequencing depth per sample with median sequencing depth (grey dashed line).
                                  \nB) Library size normalisation factor per sample with median size factor (grey dashed line).
                                  \nC) Marginal Distribution of mean, dispersion and dropout.
                                  \nD) Local polynomial regression fit between mean and dispersion estimates with variability band per gene (yellow). \nCommon dispersion estimate (grey dashed line).
                                  \nE) Fraction of dropouts versus estimated mean expression per gene.", size=8)
    }
  }

  # draw the plot
  cowplot::ggdraw(p.final)
}

# plotSpike ---------------------------------------------------------------

#' @name plotSpike
#' @aliases plotSpike
#' @title Visualize distributional characteristics of spike-ins
#' @description This function plots the results of the parameter estimation for spike-ins. This includes the absolute and relative sequencing depth (i.e. library size factor), a calibration curve as well as the capture efficiency given as a binomial regression.
#' @usage plotSpike(estSpike, annot=TRUE)
#' @param estSpike The output of \code{\link{estimateParam}}.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' #' ## batch annotation
#' data(scrbseq_spike_cnts)
#' data(scrbseq_spike_info)
#' batch_info <- data.frame(Batch = ifelse(grepl(pattern = "SCRBseqA_",
#' colnames(scrbseq_spike_cnts)), "A", "B"),
#' row.names = colnames(scrbseq_spike_cnts))
#' ## spike information table
#' spike_info <- scrbseq_spike_info[-1,]
#' ## estimation
#' spike_param <- estimateSpike(spikeData = scrbseq_spike_cnts,
#' spikeInfo = spike_info,
#' MeanFragLength = NULL,
#' batchData = batch_info,
#' normalisation = 'depth')
#' ## plotting
#' plotSpike(estSpike = spike_param, annot=TRUE)
#' }
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot aes theme_minimal geom_bar geom_density geom_hline theme labs scale_y_continuous coord_flip geom_pointrange geom_point geom_smooth annotate scale_x_log10 scale_y_log10 annotation_logticks
#' @importFrom dplyr left_join group_by mutate ungroup do summarise
#' @importFrom tidyr "%>%"
#' @importFrom broom glance
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom stats reorder
#' @rdname plotSpike
#' @export
plotSpike <- function(estSpike, annot=TRUE) {

  if(length(estSpike$seqDepth<15)) {
    # library size plot
    lib.size.dat <- data.frame(Seqdepth=estSpike$seqDepth,
                               Sample=names(estSpike$seqDepth))
    libsize.plot <- ggplot2::ggplot(data=lib.size.dat, ggplot2::aes(reorder(Sample, Seqdepth),Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x=NULL, y="Sequencing depth") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estSpike$size.factors,
                         Sample=names(estSpike$size.factors))
    sf.plot <- ggplot2::ggplot(data=sf.dat, ggplot2::aes(reorder(Sample, SizeFactor),SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x=NULL, y="Library Size Factor") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
  }
  if(length(estSpike$seqDepth>=15)) {
    # library size plot
    lib.size.dat <- data.frame(Seqdepth=estSpike$seqDepth,
                               Sample=names(estSpike$seqDepth))
    libsize.plot <- ggplot2::ggplot(lib.size.dat, ggplot2::aes(Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x="Sequencing Depth", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estSpike$size.factors,
                         Sample=names(estSpike$size.factors))
    sf.plot <- ggplot2::ggplot(sf.dat, ggplot2::aes(SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                     axis.title=ggplot2::element_text(size=14, face="bold")) +
      ggplot2::labs(x="Library Size Factor", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain)
  }

  # calibration curve data
  cal.dat <- reshape2::melt(estSpike$normCounts)
  names(cal.dat) <- c("SpikeID", "SampleID", "normCounts")
  cal.info.dat <- cal.dat %>%
    dplyr::left_join(estSpike$Input$spikeInfo, by="SpikeID") %>%
    dplyr::group_by(factor(SpikeInput)) %>%
    dplyr::mutate(Expectation=mean(normCounts),
                  Deviation=sd(normCounts),
                  Error=sd(normCounts)/sqrt(n())) %>%
    dplyr::ungroup()
  limits <- ggplot2::aes(ymax = log10(Expectation) + log10(Error),
                         ymin= log10(Expectation) - log10(Error))
  Calibration = cal.dat %>%
    dplyr::left_join(estSpike$Input$spikeInfo, by="SpikeID") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::do(LmFit = lm(log10(normCounts+1) ~ log10(SpikeInput+1), data = .)) %>%
    broom::glance(LmFit) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(Rsquared=mean(r.squared), RsquaredSE=sd(r.squared))

  # calibration curve plot
  cal.plot <- ggplot2::ggplot(data = cal.info.dat,
                              ggplot2::aes(x=log10(SpikeInput),
                                           y=log10(Expectation))) +
    ggplot2::geom_pointrange(limits) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method='lm',formula=y~x) +
    ggplot2::annotate("text", label = paste0("italic(R) ^ 2 == ",
                                             round(Calibration$Rsquared, digits = 2),
                                             "%+-%",
                                             round(Calibration$RsquaredSE, digits = 2)),
                      parse = T, x = 0.2, y = 4, size = 4) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10(labels=c("0.1","1","10","100","1,000"),breaks=c(0.1,1,10,100,1000)) +
    ggplot2::scale_y_log10(labels=c("0.1","1","10","100","1,000"),breaks=c(0.1,1,10,100,1000)) +
    ggplot2::annotation_logticks(sides = "bl") +
    ggplot2::labs(y=expression(bold(paste(Log[10], " Estimated Expression", sep=""))),
                  x=expression(bold(paste(Log[10], " Spike-In Molecules")))) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                   axis.title=ggplot2::element_text(size=14, face="bold"),
                   axis.line.x = ggplot2::element_line(colour = "black"),
                   axis.line.y = ggplot2::element_line(colour = "black"))

  # capture efficiency data
  capture.dat <- estSpike$CaptureEfficiency %>%
    tibble::rownames_to_column(var = "SpikeID") %>%
    dplyr::select(SpikeID, p_success, hat_p_success, hat_p_success_cilower, hat_p_success_ciupper) %>%
    dplyr::mutate(hat_p_success = ifelse(is.na(hat_p_success), p_success, hat_p_success)) %>%
    dplyr::select( -p_success) %>%
    dplyr::left_join(estSpike$Input$spikeInfo, by="SpikeID")

  capture.plot <- ggplot2::ggplot(data = capture.dat,
                                  ggplot2::aes(x=log10(SpikeInput+1),
                                               y=hat_p_success)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="glm",  method.args = list(family = "binomial"), se=T) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10(labels=c("0.1","1","10","100","1,000"),breaks=c(0.1,1,10,100,1000)) +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::labs(y=expression(bold("Detection Probability")),
                  x=expression(bold(paste(Log[10], " Spike-In Molecules")))) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                   axis.title=ggplot2::element_text(size=14, face="bold"),
                   axis.line.x = ggplot2::element_line(colour = "black"),
                   axis.line.y = ggplot2::element_line(colour = "black"))

  top_row <- suppressWarnings(cowplot::plot_grid(libsize.plot,sf.plot,
                                                 labels=c('A', 'B'),
                                                 ncol=2, nrow=1))
  bottom_row <- suppressWarnings(cowplot::plot_grid(cal.plot,
                                                    capture.plot,
                                                    labels=c('C', 'D'),
                                                    ncol=2, nrow=1))
  p.final <- suppressWarnings(cowplot::plot_grid(top_row,
                                                 bottom_row,
                                                 rel_heights = c(1, 1.5),
                                                 ncol=1, nrow=2))
  # annotation under plot
  if (annot) {
    p.final <- cowplot::add_sub(p.final, "A) Sequencing depth per sample with median sequencing depth (grey dashed line).
                                  \nB) Library size normalisation factor per sample with median size factor (grey dashed line).
                                  \nC) Calibration curve with mean expression estimates and average R squared over all cells.
                                  \nD) Capture efficiency with binomial logistic regression fit over all cells.", size=8)
  }
  # draw the plot
  cowplot::ggdraw(p.final)
}


# plotCounts --------------------------------------------------------------

#' @name plotCounts
#' @aliases plotCounts
#' @title Visualize simulated counts
#' @description This function performs multidimensional scaling of the simulated counts from the pairwise sample distances using variable genes (i.e. variance unequal to zero). Prior to distance calculation, the counts are normalized using the simulated size factors and log2 transformed. In the plot, the samples are annotated by phenotype and batch, if present.
#' @usage plotCounts(simCounts, Distance, Scale, DimReduce, verbose = T)
#' @param simCounts The output of \code{\link{simulateCounts}}.
#' @param Distance The (dis-)similarity measure to be used. This can be "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" for distance measures and "spearman", "pearson" or "kendall" for correlation measures converted into distances, respectively. For more information, see \code{\link[stats]{dist}} and \code{\link[stats]{cor}}.
#' @param Scale A logical vector indicating whether to use scaled log2 transformed counts or not.
#' @param DimReduce The dimension reduction approach to be used. This can be "MDS" \code{\link[stats]{cmdscale}}, "PCA" \code{\link[stats]{prcomp}}, "t-SNE" \code{\link[Rtsne]{Rtsne}}, "ICA" \code{\link[fastICA]{fastICA}} or "LDA" \code{\link[MASS]{lda}}.
#' @param verbose Logical value to indicate whether to print function information. Default is \code{TRUE}.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## not yet
#' }
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal
#' @importFrom stats dist cor as.dist cmdscale prcomp predict
#' @importFrom MASS isoMDS
#' @importFrom Rtsne Rtsne
#' @importFrom fastICA fastICA
#' @importFrom tidyr "%>%" separate
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
#' @rdname plotCounts
#' @export
plotCounts <- function(simCounts, Distance, Scale, DimReduce, verbose = T) {

  # normalize
  norm.counts <- t(t(simCounts$GeneCounts)/simCounts$sf)
  # log2 transform
  lnorm.counts <- log2(norm.counts+1)
  # kick out invariable genes
  drop_genes <- apply(lnorm.counts, 1, function(x) {var(x) < 0.01})
  if(isTRUE(verbose)) {
    message(paste0("Dropping ", sum(drop_genes), " genes out of a total of ",
                   nrow(lnorm.counts), " genes."))
  }
  lnorm.counts <- lnorm.counts[!drop_genes, ]
  # transpose
  vals <- t(lnorm.counts)
  # apply scale
  if(Scale) {
    vals <- scale(vals, scale = TRUE)
  }
  # calculate dissimilarity matrix
  if(Distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") &&
     DimReduce == "MDS") {
    if(isTRUE(verbose)) {
      message(paste0("Calculating distance matrix."))
    }
    dist.mat <- stats::dist(vals, method = Distance)
  }
  if(Distance %in% c("pearson", "kendall", "spearman") &&
     DimReduce == "MDS") {
    if(isTRUE(verbose)) {
      message(paste0("Calculating distance matrix."))
    }
    cor.mat <- stats::cor(vals, method = Distance)
    dist.mat <- stats::as.dist(1-cor.mat)
  }
  # else { stop("Unrecognized form of distance measure!\n") }

  # apply dimension reduction
  if(isTRUE(verbose)) {
    message(paste0("Applying dimension reduction."))
  }
  if(DimReduce == "MDS") {
    mds_out <- stats::cmdscale(d = dist.mat)
  }
  if(DimReduce == "PCA") {
    mds_out <- stats::prcomp(x = vals, center = T, scale = Scale)
    mds_out <- mds_out$x[,c(1:2)]
  }
  if(DimReduce == "t-SNE") {
    # dimension reduction on expression matrix: PCA + t-SNE
    # sample.value <- ncol(lnorm.counts) -2
    # max.perplexity <- sample.value/3
    tsne.res <- Rtsne::Rtsne(X=vals, pca=T, is_distance=F, perplexity=30 )
    mds_out <- tsne.res$Y
    rownames(mds_out) <- colnames(lnorm.counts)
  }
  if(DimReduce == "ICA") {
    ica.res <- fastICA::fastICA(X = vals,
                                n.comp = 2,
                                alg.typ = "deflation",
                                fun = "logcosh",
                                alpha = 1,
                                method = "R",
                                row.norm = FALSE,
                                maxit = 200,
                                tol = 0.0001,
                                verbose = FALSE)
    mds_out <- ica.res$S # ICA components
    rownames(mds_out) <- colnames(lnorm.counts)
  }
  if(DimReduce == "LDA") {
    lda.dat <- vals %>%
      tibble::rownames_to_column(var="Sample") %>%
      tidyr::separate(col = Sample, into = c("SampleID", "Phenotype", "Batch"),
                      extra = "drop", fill = "right", remove = TRUE) %>%
      dplyr::select(-SampleID, -Batch)
    lda.prior <- as.vector(table(lda.dat$Phenotype)/sum(table(lda.dat$Phenotype)))
    lda.res <- MASS::lda(Phenotype ~ ., lda.dat, prior = lda.prior)
    plda.res <- stats::predict(object = lda.res, newdata = lda.dat)
    mds_out <- plda.res$x
    rownames(mds_out) <- colnames(lnorm.counts)
  }
  # else { stop("Unrecognized form of dimension reduction!\n") }

  # collect data to plot
  if(isTRUE(verbose)) {
    message(paste0("Creating plot."))
  }
  colnames(mds_out) <- c("Dimension1", "Dimension2")
  dat.plot <- data.frame(mds_out) %>%
    tibble::rownames_to_column(var="Sample") %>%
    tidyr::separate(col = Sample, into = c("SampleID", "Phenotype", "Batch"),
                    extra = "drop", fill = "right", remove = FALSE)

  # plot
  p1 <- ggplot2::ggplot(data = dat.plot,
                        ggplot2::aes(x = Dimension1, y = Dimension2))
  if(all(!is.na(dat.plot$Batch))) {
    p2 <- p1 + ggplot2::geom_point(ggplot2::aes(shape = Batch, colour = Phenotype),
                                   size = 2, alpha =0.5, data = dat.plot)
  }
  if(all(is.na(dat.plot$Batch))) {
    p2 <- p1 + ggplot2::geom_point(ggplot2::aes(colour = Phenotype),
                                   size = 2, alpha =0.5, data = dat.plot)
  }
  p3 <- p2 + ggplot2::theme_minimal()

  return(p3)

}

# plotEvalDE -------------------------------------------------------------
#' @name plotEvalDE
#' @aliases plotEvalDE
#' @title Visualize power assessment
#' @description This function plots the results of \code{\link{evaluateDE}} for assessing the error rates and sample size requirements.
#' @usage plotEvalDE(evalRes, rate=c('marginal', 'stratified'),
#'                    quick=TRUE, annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateDE}}.
#' @param rate Character vector defining whether the marginal or condtional rates should be plotted. Conditional depends on the choice of stratify.by in \code{\link{evaluateDE}}.
#' @param quick A logical vector. If \code{TRUE}, the TPR and FDR are only plotted. If \code{FALSE}, then all rates are plotted.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## using example data set
#' eval.de <- evaluateDE(simRes = kolodziejczk_simDE)
#' plotEvalDE(evalRes=eval.de, rate ="marginal", quick=T, annot=T)
#' plotEvalDE(evalRes=eval.de, rate ="stratified", quick=T, annot=T)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar theme_minimal
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom tidyr %>%
#' @rdname plotEvalDE
#' @export
plotEvalDE <- function(evalRes, rate=c('marginal', 'stratified'), quick=TRUE, annot=TRUE) {

  rate = match.arg(rate)

  # marginal rates over sample sizes
  if(rate=='marginal') {
    if(quick) {
      dat.marginal <- evalRes[c('TPR.marginal', 'FDR.marginal')]
      names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
      dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalRes[['n1']], " vs ", evalRes[['n2']]))
      dat.marginal.long <- reshape2::melt(dat.marginal)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      dat.marginal.calc <- dat.marginal.long %>% dplyr::group_by(Var1, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      # marginal in one
      grandplot <- ggplot2::ggplot(data = dat.marginal.long, ggplot2::aes(x=Var1, y=value, color=L1)) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=NULL, y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        ggplot2::theme(legend.position='top',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y=ggplot2::element_text(size=10),
                       axis.title=ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"))
      # faceted marginal
      facetplot <-  ggplot2::ggplot(data = dat.marginal.calc, ggplot2::aes(x=Var1, y=Expectation, fill=L1, color=L1)) +
        ggplot2::geom_line(ggplot2::aes(group=L1)) +
        ggplot2::geom_pointrange(limits) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x="Samples", y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::theme(legend.position='none',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y = ggplot2::element_text(size=10),
                       axis.title = ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10, face="bold")) +
        ggplot2::facet_wrap(~L1, scales = 'free', ncol=1) +
        ggplot2::geom_hline(data = refval, ggplot2::aes(yintercept = ref), linetype="dashed", color='grey')
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(grandplot,
                                                     facetplot,
                                                     labels=c('A', 'B'),
                                                     rel_heights = c(1,1.5),
                                                     ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Marginal TPR and FDR per sample size comparison. \nB) Marginal TPR and FDR per sample size comparison with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error).", size=8)
      }

    }
    if(!quick) {
      dat.marginal <- evalRes[grep('*R.marginal', names(evalRes))]
      names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
      dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalRes[['n1']], " vs ", evalRes[['n2']]))
      dat.marginal.long <- reshape2::melt(dat.marginal)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      dat.marginal.calc <- dat.marginal.long %>%
        dplyr::group_by(Var1, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      # marginal in one
      grandplot <- ggplot(data = dat.marginal.long, ggplot2::aes(x=Var1, y=value, color=L1)) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=NULL, y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        ggplot2::theme(legend.position='top',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y=ggplot2::element_text(size=10),
                       axis.title=ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"))
      # faceted marginal
      facetplot <-  ggplot2::ggplot(data = dat.marginal.calc, ggplot2::aes(x=Var1, y=Expectation, fill=L1, color=L1)) +
        ggplot2::geom_line(ggplot2::aes(group=L1)) +
        ggplot2::geom_pointrange(limits) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x='Samples', y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::theme(legend.position='none',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y=ggplot2::element_text(size=10),
                       axis.title=ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10, face="bold")) +
        ggplot2::facet_wrap(~L1, scales = 'free', ncol=2) +
        ggplot2::geom_hline(data = refval, ggplot2::aes(yintercept = ref), linetype="dashed", color='grey')
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(grandplot,
                                                     facetplot,
                                                     labels=c('A', 'B'),
                                                     rel_heights = c(1,2),
                                                     ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Marginal error rates per sample size comparison. \nB) Marginal error rates per sample size comparison with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error).", size=8)
      }
    }
  }

  #stratified rates
  if(rate=='stratified') {
    if(quick){
      dat.stratified <- evalRes[c('TPR', 'FDR')]
      strata <- evalRes$strata.levels
      dat.stratified <- lapply(dat.stratified, "dimnames<-", list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]), NULL))
      dat.stratified.long <- reshape2::melt(dat.stratified)
      dat.stratified.calc <- dat.stratified.long %>%
        dplyr::group_by(Var1, Var2, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      facetplot <-  ggplot2::ggplot(data = dat.stratified.calc, ggplot2::aes(x=Var1, y=Expectation, fill=Var2, color=Var2)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(ggplot2::aes(group=Var2)) +
        ggplot2::geom_pointrange(limits) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=NULL, y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        ggplot2::theme(legend.position='top',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y = ggplot2::element_text(size=10),
                       axis.title = ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10, face="bold")) +
        ggplot2::facet_wrap(~L1, scales = 'free', ncol=2) +
        ggplot2::geom_hline(data = refval, ggplot2::aes(yintercept = ref), linetype="dashed", color='grey')
      # strata genes
      N <- length(evalRes$n1)
      dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],'DEgenes'=evalRes$stratadiffgenes[,N,])
      dat.genes <- lapply(dat.genes, "rownames<-", strata)
      dat.genes.long <- reshape2::melt(dat.genes)
      dat.genes.calc <- dat.genes.long %>% dplyr::group_by(Var1, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      dodge <- ggplot2::position_dodge(width=0.9)
      strataplot <- ggplot2::ggplot(data = dat.genes.calc, ggplot2::aes(x=Var1, y=Expectation, fill=L1)) +
        ggplot2::geom_bar(stat="identity")  +
        ggplot2::theme_minimal() +
        ggplot2::labs(x='Stratum', y="Count") +
        ggplot2::theme(legend.position='right',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y = ggplot2::element_text(size=10),
                       axis.title = ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10)) +
        ggplot2::scale_fill_manual(values=c('grey', 'black'),
                                   breaks = c("DEgenes", "Ngenes"),
                                   labels = c("DE genes", "EE genes"))
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(facetplot,
                                                     strataplot,
                                                     labels=c('A', 'B'),
                                                     rel_heights = c(2,1),
                                                     ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Conditional TPR and FDR per sample size comparison per stratum. \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.", size=8)
      }
    }
    if(!quick) {
      dat.stratified <- evalRes[grep('*R$', names(evalRes))]
      strata <- evalRes$strata.levels
      dat.stratified <- lapply(dat.stratified, "dimnames<-", list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]), NULL))
      dat.stratified.long <- reshape2::melt(dat.stratified)
      dat.stratified.calc <- dat.stratified.long %>%
        dplyr::group_by(Var1, Var2, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      facetplot <-   ggplot2::ggplot(data = dat.stratified.calc, ggplot2::aes(x=Var1, y=Expectation, fill=Var2, color=Var2)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(ggplot2::aes(group=Var2)) +
        ggplot2::geom_pointrange(limits) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=NULL, y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        ggplot2::theme(legend.position = 'top',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y = ggplot2::element_text(size=10),
                       axis.title = ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10, face="bold")) +
        ggplot2::facet_wrap(~L1, scales = 'free_x', ncol=2) +
        ggplot2::geom_hline(data = refval, ggplot2::aes(yintercept = ref), linetype="dashed", color='grey')
      # strata genes
      N <- length(evalRes$n1)
      dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],'DEgenes'=evalRes$stratadiffgenes[,N,])
      dat.genes <- lapply(dat.genes, "rownames<-", strata)
      dat.genes.long <- reshape2::melt(dat.genes)
      dat.genes.calc <- dat.genes.long %>%
        dplyr::group_by(Var1, L1) %>%
        dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>%
        dplyr::ungroup()
      dodge <-  ggplot2::position_dodge(width=0.9)
      strataplot <-  ggplot2::ggplot(data = dat.genes.calc, ggplot2::aes(x=Var1, y=Expectation, fill=L1)) +
        ggplot2::geom_bar(stat="identity")  +
        ggplot2::theme_minimal() +
        ggplot2::labs(x="Stratum", y="Count") +
        ggplot2::theme(legend.position='right',
                       legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size=10, angle=45, hjust=1),
                       axis.text.y=ggplot2::element_text(size=10),
                       axis.title=ggplot2::element_text(size=10, face="bold"),
                       legend.text = ggplot2::element_text(size=10),
                       legend.key.size = grid::unit(1, "cm"),
                       strip.text = ggplot2::element_text(size=10)) +
        ggplot2::scale_fill_manual(values=c('grey', 'black'),
                                   breaks = c("DEgenes", "Ngenes"),
                                   labels = c("DE genes", "EE genes"))
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(facetplot,
                                                     strataplot,
                                                     labels=c('A', 'B'),
                                                     rel_heights = c(3,1),
                                                     ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Conditional error rates over stratum. \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.", size=8)
      }
    }
  }

  # draw final plot
  cowplot::ggdraw(p.final)
}

# plotEvalSim -------------------------------------------------------------

#' @name plotEvalSim
#' @aliases plotEvalSim
#' @title Visualize power assessment
#' @description This function plots the results of \code{\link{evaluateSim}} for assessing the setup performance, i.e. normalisation method performance.
#' @usage plotEvalSim(evalRes, annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateSim}}.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## using example data set
#' eval.sim <- evaluateSim(simRes = kolodziejczk_simDE, timing = T)
#' plotEvalSim(eval.sim)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme element_blank scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar theme_minimal
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom tidyr "%>%"
#' @rdname plotEvalSim
#' @export
plotEvalSim <- function(evalRes, annot=TRUE) {

  if (attr(evalRes, 'Simulation') == 'Flow') {
    # log fold changes
    lfc <- reshape2::melt(evalRes$Log2FoldChange)
    colnames(lfc) <- c("SimNo", "Metric", "Value", "Samples")
    lfc.dat <- lfc %>%
      tidyr::separate(Metric, c("DE-Group", "Metric", "Type"), "_") %>%
      dplyr::filter(!Type=="NAFraction") %>%
      dplyr::group_by(Samples, `DE-Group`, Metric) %>%
      dplyr::summarise(Expectation=mean(Value, na.rm=T),
                       Deviation=sd(Value, na.rm=T),
                       Error=sd(Value, na.rm=T)/sqrt(n())) %>%
      dplyr::ungroup() %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    # to label the x axis from smallest to largest n group!
    lfc.dat$Samples <- factor(lfc.dat$Samples,
                              levels=unique(lfc.dat$Samples[order(lfc.dat$SumN,decreasing = F)]))

    limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
    dodge <- ggplot2::position_dodge(width=0.7)
    lfc.plot <- ggplot2::ggplot(data = lfc.dat, ggplot2::aes(x=Samples,
                                                    y=Expectation,
                                                    fill=`DE-Group`,
                                                    colour=`DE-Group`)) +
      ggplot2::geom_point(position = dodge) +
      # ggplot2::geom_line(ggplot2::aes(group=Metric),position = dodge) +
      ggplot2::geom_pointrange(limits, position = dodge) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::facet_wrap(~Metric, ncol=3, scales="free") +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))

    # size factors
    sf <- reshape2::melt(evalRes$SizeFactors)
    colnames(sf) <- c("SimNo", "Metric", "Value", "Samples")
    sf.stats <- sf %>%
      dplyr::filter(Metric %in% c("MAD", "rRMSE")) %>%
      dplyr::group_by(Samples, Metric) %>%
      dplyr::summarise(Expectation=mean(Value, na.rm=T),
                       Deviation=sd(Value, na.rm=T),
                       Error=sd(Value, na.rm=T)/sqrt(n())) %>%
      dplyr::ungroup() %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    # to label the x axis from smallest to largest n group!
    sf.stats$Samples <- factor(sf.stats$Samples,
                               levels=unique(sf.stats$Samples[order(sf.stats$SumN,decreasing = F)]))
    limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
    dodge <- ggplot2::position_dodge(width=0.7)
    sfstats.plot <- ggplot2::ggplot(data = sf.stats, ggplot2::aes(x=Samples,
                                                         y=Expectation)) +
      ggplot2::geom_point(position = dodge) +
      # ggplot2::geom_line(ggplot2::aes(group=Metric),position = dodge) +
      ggplot2::facet_wrap(~Metric, ncol=2, scales="free") +
      ggplot2::geom_pointrange(limits, position = dodge) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))
    # ratio of size factors per group
    ratio.dat <-  sf %>%
      dplyr::filter(grepl("Group", Metric)) %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    ratio.dat$Samples <- factor(ratio.dat$Samples,
                                levels=unique(ratio.dat$Samples[order(ratio.dat$SumN,decreasing = F)]))

    ratio.plot <- ggplot2::ggplot(data = ratio.dat, ggplot2::aes(x=Samples,
                                                        y=Value,
                                                        color=Metric)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::geom_hline(yintercept=1,linetype="dashed", color='darkgrey') +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))

    # clustering
    bottom_row <- suppressWarnings(cowplot::plot_grid(sfstats.plot,
                                                      ratio.plot,
                                                      labels = c('B', 'C'),
                                                      align = 'hv',
                                                      ncol=2,
                                                      nrow=1,
                                                      rel_widths = c(1.3, 1)))
    p.final <- suppressWarnings(cowplot::plot_grid(lfc.plot,
                                                   bottom_row,
                                                   labels=c('A', ''),
                                                   rel_heights = c(1.3,1),
                                                   ncol=1, nrow=2))
    # annotation under plot
    if(annot) {
      p.final <- cowplot::add_sub(p.final, "A) Mean Absolute Error (MAE), Root Mean Squared Error (RMSE) and robust Root Mean Squared Error (rRMSE) \n for the estimated log2 fold changes of all (ALL), differentially expressed (DE) and equally expressed (EE) genes compared to the true log2 fold changes.
                                  \nB) Median absolute deviation (MAD) and robust Root Mean Squared Error (rRMSE) between estimated and true size factors.
                                  \nC) The average ratio between simulated and true size factors in the two groups of samples.", size=8)
    }

  }

  if (attr(evalRes, 'Simulation') == 'DE') {
    # log fold changes
    lfc <- reshape2::melt(evalRes$Log2FoldChange)
    colnames(lfc) <- c("SimNo", "Metric", "Value", "Samples")
    lfc.dat <- lfc %>%
      tidyr::separate(Metric, c("DE-Group", "Metric", "Type"), "_") %>%
      dplyr::filter(!Type=="NAFraction") %>%
      dplyr::group_by(Samples, `DE-Group`, Metric) %>%
      dplyr::summarise(Expectation=mean(Value, na.rm=T),
                       Deviation=sd(Value, na.rm=T),
                       Error=sd(Value, na.rm=T)/sqrt(n())) %>%
      dplyr::ungroup() %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    # to label the x axis from smallest to largest n group!
    lfc.dat$Samples <- factor(lfc.dat$Samples,
                              levels=unique(lfc.dat$Samples[order(lfc.dat$SumN,decreasing = F)]))

    limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
    dodge <- ggplot2::position_dodge(width=0.7)
    lfc.plot <- ggplot2::ggplot(data = lfc.dat, ggplot2::aes(x=Samples,
                                                     y=Expectation,
                                                     fill=`DE-Group`,
                                                     colour=`DE-Group`)) +
      ggplot2::geom_point(position = dodge) +
      # ggplot2::geom_line(ggplot2::aes(group=Metric),position = dodge) +
      ggplot2::geom_pointrange(limits, position = dodge) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::facet_wrap(~Metric, ncol=3, scales="free") +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))

    # size factors
    sf <- reshape2::melt(evalRes$SizeFactors)
    colnames(sf) <- c("SimNo", "Metric", "Value", "Samples")
    sf.stats <- sf %>%
      dplyr::filter(Metric %in% c("MAD", "rRMSE")) %>%
      dplyr::group_by(Samples, Metric) %>%
      dplyr::summarise(Expectation=mean(Value, na.rm=T),
                       Deviation=sd(Value, na.rm=T),
                       Error=sd(Value, na.rm=T)/sqrt(n())) %>%
      dplyr::ungroup() %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    # to label the x axis from smallest to largest n group!
    sf.stats$Samples <- factor(sf.stats$Samples,
                              levels=unique(sf.stats$Samples[order(sf.stats$SumN,decreasing = F)]))
    limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
    dodge <- ggplot2::position_dodge(width=0.7)
    sfstats.plot <- ggplot2::ggplot(data = sf.stats, ggplot2::aes(x=Samples,
                                                    y=Expectation)) +
      ggplot2::geom_point(position = dodge) +
      # ggplot2::geom_line(ggplot2::aes(group=Metric),position = dodge) +
      ggplot2::facet_wrap(~Metric, ncol=2, scales="free") +
      ggplot2::geom_pointrange(limits, position = dodge) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))
    # ratio of size factors per group
    ratio.dat <-  sf %>%
      dplyr::filter(grepl("Group", Metric)) %>%
      tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
      dplyr::arrange(SumN)
    ratio.dat$Samples <- factor(ratio.dat$Samples,
                               levels=unique(ratio.dat$Samples[order(ratio.dat$SumN,decreasing = F)]))

    ratio.plot <- ggplot2::ggplot(data = ratio.dat, ggplot2::aes(x=Samples,
                                                        y=Value,
                                                        color=Metric)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="Samples", y="Value") +
      ggplot2::geom_hline(yintercept=1,linetype="dashed", color='darkgrey') +
      ggplot2::theme(legend.position='right', legend.title = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=10, angle=45, hjust=1),
                     axis.text.y=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     legend.text = ggplot2::element_text(size=10),
                     legend.key.size = grid::unit(0.5, "cm"),
                     strip.text = ggplot2::element_text(size=10, face="bold"))

    bottom_row <- suppressWarnings(cowplot::plot_grid(sfstats.plot,
                                                      ratio.plot,
                                                      labels = c('B', 'C'),
                                                      align = 'hv',
                                                      ncol=2,
                                                      nrow=1,
                                                      rel_widths = c(1.3, 1)))
    p.final <- suppressWarnings(cowplot::plot_grid(lfc.plot,
                                                   bottom_row,
                                                   labels=c('A', ''),
                                                   rel_heights = c(1.3,1),
                                                   ncol=1, nrow=2))
    # annotation under plot
    if(annot) {
      p.final <- cowplot::add_sub(p.final, "A) Mean Absolute Error (MAE), Root Mean Squared Error (RMSE) and robust Root Mean Squared Error (rRMSE) \n for the estimated log2 fold changes of all (ALL), differentially expressed (DE) and equally expressed (EE) genes compared to the true log2 fold changes.
                                  \nB) Median absolute deviation (MAD) and robust Root Mean Squared Error (rRMSE) between estimated and true size factors.
                                  \nC) The average ratio between simulated and true size factors in the two groups of samples.", size=8)
    }

  }

  # draw final plot
  cowplot::ggdraw(p.final)

}


# plotTime ---------------------------------------------------------------

#' @name plotTime
#' @aliases plotTime
#' @title Visualize computational time
#' @description This function plots the computational running time of simulations.
#' @usage plotTime(simRes, Table=TRUE, annot=TRUE)
#' @param simRes The output of \code{\link{simulateDE}} or \code{\link{simulateFlow}}.
#' @param Table A logical vector. If \code{TRUE}, a table of average computational running time per step and sample size is printed additionally.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## using example data set
#' plotTime(simRes = kolodziejczk_simDE)
#' }
#' @author Beate Vieth
#' @importFrom tidyr "%>%" separate
#' @importFrom dplyr mutate arrange
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes position_dodge geom_point geom_pointrange facet_wrap theme_bw labs theme element_blank element_text guide_legend guides
#' @importFrom grid unit
#' @importFrom cowplot add_sub ggdraw
#' @importFrom matrixStats rowSds
#' @rdname plotTime
#' @export
plotTime <- function(simRes, Table=TRUE, annot=TRUE) {

  # simulation parameters
  Nreps1 = simRes$sim.settings$n1
  Nreps2 = simRes$sim.settings$n2
  time.taken = simRes$time.taken
  nsims = simRes$sim.settings$nsims
  ncores = simRes$sim.settings$NCores
  if(is.null(ncores)) {ncores=1}

  # create output objects
  my.names = paste0(Nreps1, " vs ", Nreps2)
  time.taken.mat <- lapply(1:length(my.names), function(x) {
    data.frame(matrix(NA, nrow = length(rownames(time.taken[,1,]))+1,
                      ncol = 3, dimnames = list(c(rownames(time.taken[,1,]), "Total"),
                                                c("Mean", "SD", "SEM")))
    )
  })
  names(time.taken.mat) <- my.names
  for(j in seq(along=Nreps1)) {
    tmp.time <- time.taken[,j,]
    Total <- colSums(tmp.time, na.rm = T)
    tmp.time <- rbind(tmp.time, Total)
    time.taken.mat[[j]][,"Mean"] <- rowMeans(tmp.time)
    time.taken.mat[[j]][,"SD"] <- matrixStats::rowSds(tmp.time)
    time.taken.mat[[j]][,"SEM"] <- matrixStats::rowSds(tmp.time)/sqrt(nsims)
  }
  time.taken.dat <- do.call('rbind', time.taken.mat)
  # time.taken.dat <- time.taken.dat[!is.na(time.taken.dat$Mean),]
  time.taken.dat <- time.taken.dat %>%
    tibble::rownames_to_column(var="ID") %>%
    tidyr::separate(col = ID, into=c("Samples", "Step"), sep="[.]", remove=TRUE) %>%
    tidyr::separate(Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
    dplyr::mutate(SumN = as.numeric(n1)+as.numeric(n2)) %>%
    dplyr::arrange(SumN)
  # to label the x axis from smallest to largest n group!
  time.taken.dat$Samples <- factor(time.taken.dat$Samples,
                            levels=unique(time.taken.dat$Samples[order(time.taken.dat$SumN,decreasing = F)]))
  time.taken.dat$Step <- factor(time.taken.dat$Step,
                                   levels=c("Preprocess", "Normalisation", "Clustering", "DE", "Moments", "Total"))
  if(isTRUE(Table)) {
    printtime <- time.taken.dat[,c(1,4:7)]
    printtime[,c(3:5)] <- signif(printtime[,c(3:5)],2)
    print(printtime)
  }
  time.taken.dat <- time.taken.dat[!is.na(time.taken.dat$Mean),]
  # plot
  limits <- ggplot2::aes(ymax = Mean + SEM, ymin= Mean - SEM)
  dodge <- ggplot2::position_dodge(width=0.7)
  p.final <- ggplot2::ggplot(data = time.taken.dat, ggplot2::aes(x=Step,
                                                          y=Mean,
                                                          colour=Samples)) +
    ggplot2::geom_point(position = dodge) +
    # ggplot2::geom_line(ggplot2::aes(group=Metric),position = dodge) +
    ggplot2::geom_pointrange(limits, position = dodge) +
    ggplot2::facet_wrap(~Step, nrow=1, scales="free") +
    ggplot2::theme_bw() +
    ggplot2::labs(x=NULL, y="Time in minutes") +
    ggplot2::theme(legend.position='bottom', legend.title = ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_text(size=10),
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_text(size=10, face="bold"),
                   legend.text = ggplot2::element_text(size=10),
                   legend.key.size = grid::unit(0.5, "cm"),
                   strip.text = ggplot2::element_text(size=10, face="bold")) +
    ggplot2::guides(colour=ggplot2::guide_legend(nrow=1))

  if(isTRUE(annot)) {
    subtitle <- paste0("The average time in minutes per simulation step and sample size.\nThe number of cores was set to ", ncores, ".")
    p.final <- cowplot::add_sub(p.final, subtitle, size=8)
  }

  # draw final plot
  cowplot::ggdraw(p.final)

}

# plotEvalDist ------------------------------------------------------------

#' @name plotEvalDist
#' @aliases plotEvalDist
#' @title Visualize distribution assessment
#' @description This function plots the results of \code{\link{evaluateDist}} to assess goodness-of-fit testing.
#' @usage plotEvalDist(evalDist, annot=TRUE)
#' @param evalDist The output of \code{\link{evaluateDist}}.
#' @param annot A logical vector. If \code{TRUE}, a short description of the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## using example data set
#' data(kolodziejczk_cnts)
#' evaldist <- evaluateDist(countData = kolodziejczk_cnts,
#' RNAseq = "singlecell", normalisation="scran",
#' frac.genes=1, min.meancount = 0.1,
#' max.dropout=0.7, min.libsize=1000,
#' verbose = TRUE)
#' plotEvalDist(evaldist, annot = TRUE)
#' }
#' @author Beate Vieth, Ines Hellmann
#' @importFrom ggplot2 ggplot aes ylab xlab theme scale_y_continuous geom_boxplot position_dodge geom_bar theme_minimal scale_x_discrete labs ggtitle coord_flip
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr filter group_by summarise mutate bind_rows slice select ungroup  left_join
#' @importFrom tidyr %>% separate gather spread
#' @importFrom utils stack
#' @rdname plotEvalDist
#' @export
plotEvalDist <- function(evalDist, annot=TRUE){
  # define naming labels for plotting
  #                    "Multiple"='Multiple',
  dist_labels <- c("None"='None',
                   "ZIP"='Zero-Inflated \n Poisson',
                   "Poisson"='Poisson',
                   "ZINB"="Zero-Inflated \n Negative Binomial",
                   "NB"="Negative \n Binomial"
  )
  label_names<-c("PoiBeta" ="Beta-Poisson",
                 'zifpois' = "Zero-Inflated \n Poisso",
                 'pois' = "Poisson",
                 "zifnbinom" = "Zero-Inflated \n Negative Binomial",
                 "nbinom" = "Negative \n Binomial"
  )
  combi.ind <- c("1_1_1_1"='Multiple',
                 "0_0_0_0"='None',
                 "0_1_0_0"='Poisson',
                 "1_0_0_0"="NB",
                 "0_0_1_0"="ZINB",
                 "0_0_0_1"='ZIP',
                 "1_1_0_0"='Multiple',
                 "1_0_0_1"='Multiple',
                 "1_0_1_0"="Multiple",
                 "0_0_1_1"="Multiple",
                 "1_0_1_1"='Multiple',
                 "1_1_1_0"='Multiple',
                 "0_1_1_1"='Multiple'
  )
  combi.ind.df <- data.frame(combi.ind, stringsAsFactors=F) %>%
    tibble::rownames_to_column(var = "nbinom_pois_zifnbinom_zifpois") %>%
    dplyr::rename(Distribution=combi.ind)

  # extract the GOF results from object
  gofres <- evalDist$GOF_res
  # reshape the table
  gofres <- cbind(rn = rownames(gofres), utils::stack(gofres))
  gofres <- gofres %>%
    tidyr::separate(ind, c("distribution", "framework", 'type'), "_", remove = F)
  # extract the observed zero table
  obszero <- evalDist$ObservedZeros
  obszero <- cbind(rn = rownames(obszero), obszero)
  obszero <- obszero[which(obszero$rn %in% gofres$rn), ]

  # GOF p-value based on chisquare test
  gof.pval <- gofres %>%
    dplyr::filter(type=='gofpval', framework=='standard') %>%
    dplyr::mutate(values, TestRes=ifelse(values>0.05, 1, 0)) %>%
    dplyr::select(rn, TestRes, distribution) %>%
    tidyr::spread( distribution, TestRes) %>%
    tidyr::unite(nbinom_pois_zifnbinom_zifpois, nbinom, pois, zifnbinom,zifpois, remove = F) %>%
    dplyr::full_join(combi.ind.df, by='nbinom_pois_zifnbinom_zifpois') %>%
    dplyr::mutate(Distribution=ifelse(is.na(Distribution), 'Multiple', Distribution)) %>%
    dplyr::group_by(Distribution) %>%
    dplyr::summarise(Total=length(nbinom_pois_zifnbinom_zifpois)) %>%
    dplyr::mutate(Percentage=Total/sum(Total)) %>%
    dplyr::filter(Distribution != 'Multiple')
  p.chisquare <- ggplot2::ggplot(gof.pval, aes(x = Distribution, y=Percentage)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = 'Distribution', y = "Percentage")  +
    geom_bar(stat='identity', width=0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::scale_x_discrete(labels= dist_labels, limits=names(dist_labels)) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle('Goodness-of-fit statistic')

  # AIC (lowest value AND intersect with GOF p-value)
  AIC.calc <- gofres %>%
    dplyr::filter(framework %in%c('Marioni', 'standard'), type %in% c("gofpval","aic")) %>%
    dplyr::select(-ind ) %>%
    tidyr::spread(type,values) %>%
    dplyr::group_by(rn) %>%
    dplyr::mutate(minAIC= (aic==min(aic)), GOF = gofpval>0.05) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(distribution) %>%
    dplyr::summarise(`Total Lowest\nAIC`= sum(minAIC,na.rm = T),
                     `Total Lowest\nAIC + \nGOF \np >0.05` = sum((minAIC & GOF), na.rm=T)) %>%
    dplyr::mutate(Total.Lowest=sum(`Total Lowest\nAIC`),Total.Sub=sum(`Total Lowest\nAIC + \nGOF \np >0.05`)) %>%
    dplyr::mutate(`Percentage Lowest\nAIC`=`Total Lowest\nAIC`/Total.Lowest, `Percentage Lowest\nAIC + \nGOF \np >0.05`=`Total Lowest\nAIC + \nGOF \np >0.05`/Total.Sub) %>%
    dplyr::select(-Total.Lowest, -Total.Sub) %>%
    gather(variable, value, `Total Lowest\nAIC`:`Percentage Lowest\nAIC + \nGOF \np >0.05`, factor_key=FALSE) %>%
    mutate(Type=ifelse(grepl(pattern='Percentage', variable), 'Percentage', 'Total'), Set=sub(".+? ", "", variable)) %>%
    dplyr::filter(Type=='Percentage')
  AIC.calc$distribution <- factor(AIC.calc$distribution, levels = c('PoiBeta', 'zifpois', 'pois', 'zifnbinom', 'nbinom'))
  p.aic <- ggplot2::ggplot(data=AIC.calc, aes(x=Set, y=value, fill=distribution)) +
    ggplot2::geom_bar(stat='identity', width=0.7, colour="black") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::scale_fill_brewer(labels= label_names, limits=names(label_names), palette="Set1") +
    ggplot2::xlab(NULL) +
    ggplot2::ylab('Percentage') +
    ggplot2::labs(fill=NULL) +
    ggplot2::guides(fill=guide_legend(reverse=TRUE)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle("Akaike Information Criterion") +
    ggplot2::coord_flip()

  # predicted zeroes vs observed zeros
  predzero.dat <- gofres %>%
    dplyr::filter(type %in% c("gofpval","predzero"),
                  framework %in% c("Marioni", 'standard')) %>%
    dplyr::select(-ind ) %>%
    tidyr::spread(type,values) %>%
    dplyr::left_join(obszero, by='rn') %>%
    dplyr::mutate(`Obs. vs. Pred.\n Zeros`=ObsZero-predzero,
                  `Obs. vs. Pred.\n Zeros, GOF p>0.05`=ifelse(gofpval>0.05,ObsZero-predzero,NA)) %>%
    tidyr::gather(comp,divergentzero,7:8)
  dd = 1:5
  zpdat <- dplyr::filter(predzero.dat, comp == "Obs. vs. Pred.\n Zeros")
  p.zero <- ggplot2::ggplot(zpdat, aes(x = distribution, y = divergentzero)) +
    geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Distribution", y = "Observed - Predicted Zeros") +
    ggplot2::scale_x_discrete(labels = label_names[dd], limits = names(label_names[dd])) +
    ggplot2::ggtitle("Dropouts") +
    ggplot2::theme(legend.position = "none", legend.title = element_blank()) +
    coord_flip()

  # Best fit by LRT / Vuong
  model.pval <- gofres %>%
    dplyr::filter(type=='gofpval', framework=='standard') %>%
    dplyr::mutate(values, TestRes=ifelse(values>0.05, 1, 0)) %>%
    dplyr::select(rn, TestRes, distribution) %>%
    tidyr::spread( distribution, TestRes) %>%
    na.omit() %>%
    tidyr::unite(nbinom_pois_zifnbinom_zifpois, nbinom, pois, zifnbinom,zifpois, remove = F) %>%
    mutate_if(is.factor, as.character)

  model.calc <- gofres %>%
    dplyr::filter(distribution %in% c('LRT', 'Vuong'), values<0.05) %>%
    dplyr::select(rn, type) %>%
    mutate_if(is.factor, as.character)
  out <- split(model.calc, f = model.calc$type)
  out2 <- lapply(out, function(x) {
    left_join(x, model.pval, by='rn') %>% na.omit()
  })
  out2[['NBPoisson']] <- out2[['NBPoisson']] %>% dplyr::filter(nbinom==1, pois==1)  %>% select(rn, type)
  out2[['ZNB']] <- out2[['ZNB']] %>% dplyr::filter(nbinom==1, zifnbinom==1)  %>% select(rn, type)
  out2[['ZNBZPoisson']] <- out2[['ZNBZPoisson']] %>% dplyr::filter(zifnbinom==1, zifpois==1)  %>% select(rn, type)
  out2[['ZPoisson']] <- out2[['ZPoisson']] %>% dplyr::filter(pois==1, zifpois==1)  %>% select(rn, type)
  out3 <- do.call('rbind', out2)
  ModelTest <- out3 %>%
    dplyr::group_by(type) %>%
    count() %>% mutate(Percentage=n/sum(n))
  p.modeltest <- ggplot2::ggplot(data = ModelTest, aes(x = type, y = Percentage)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(x = "Test", y = "Percentage") +
    ggplot2::scale_x_discrete(labels = c(NBPoisson = "Negative Binomial \n> Poisson",
                                         ZNB = "Zero-inflated Negative Binomial \n> Negative Binomial",
                                         ZPoisson = "Zero-inflated Poisson \n> Poisson",
                                         ZNBZPoisson = "Zero-inflated Negative Binomial \n> Zero-inflated Poisson"),
                              limits = c("ZNBZPoisson", "ZPoisson", "NBPoisson", "ZNB")) +
    ggplot2::ggtitle("Model Comparisons") +
    ggplot2::coord_flip()

  p.final <- cowplot::ggdraw() +
    cowplot::draw_plot(p.chisquare, x= 0 , y= 0.5,  width = 0.5, height = 0.5) +
    cowplot::draw_plot(p.aic, x=0.5, y=0.5, width = 0.5, height=0.5) +
    cowplot::draw_plot(p.zero, x=0, y= 0,  width = 0.5, height=0.5) +
    cowplot::draw_plot(p.modeltest, x=0.5, y= 0,  width = 0.5, height=0.5) +
    cowplot::draw_plot_label( c("A","B","C", "D"), x=c(0, 0.5, 0, 0.5), y=c(1,1,0.5, 0.5), size=15, vjust=1)

  if (annot) {
    p.final <- cowplot::add_sub(p.final, "A) Goodness-of-fit of the model assessed with a chi-square test based on residual deviance and degrees of freedom.
                                \nB) Akaike Information Criterion per gene: Model with the lowest AIC. Model with the lowest AIC and passed goodness-of-fit statistic test.
                                \nC) Observed versus predicted dropouts per model and gene plotted without outliers.
                                \nD) Model Assessment based on LRT for nested models and Vuong test for nonnested models.",
                                size = 8)
  }
  cowplot::ggdraw(p.final)
}

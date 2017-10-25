
# plotParam -------------------------------------------------------------
#' @name plotParam
#' @aliases plotParam
#' @title Visualize
#' @description This function plots the results of the parameter estimation. This includes the absolute and relative sequencing depth (i.e. library size factor) and marginal mean, dispersion and dropout. Also the mean-dispersion relationship with loess fit for simulations is visualized. Lastly the mean-dropout rate is presented as a smooth scatter plot.
#' @usage plotParam(estParamRes, annot=TRUE)
#' @param estParamRes The output of \code{\link{estimateParam}} or a string specifying the name of precalculated estimates (see details).
#' @param annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## for example see \code{\link{insilicoNBParam}}
#' }
#' @details Precalculated negative binomial parameters for simulations are available via \url{}.
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot aes geom_bar geom_line theme geom_hline labs coord_flip scale_y_continuous geom_point scale_fill_gradientn stat_density2d labs theme_minimal theme_classic element_text ylim
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom ggExtra ggMarginal
#' @importFrom ggthemes theme_base
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
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=10, face="bold")) +
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
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=10, face="bold")) +
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
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=8, face="bold"),
                     strip.text = ggplot2::element_text(size=8, face="bold")) +
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
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
      ggplot2::geom_hline(yintercept = log2(estParamRes$common.dispersion), linetype = 2, colour="grey40") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$y),
                         linetype=1, size=1.5, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$upper),
                         linetype=2, size=1, colour="orange") +
      ggplot2::geom_line(ggplot2::aes(x=estParamRes$meandispfit$x, y=estParamRes$meandispfit$lower),
                         linetype=2, size=1, colour="orange") +
      ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))),
                    x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none',
                     axis.text=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=12, face="bold"))
    # mean vs p0 plot
    meanvsp0.dat <- data.frame(Means=log2(estParamRes$means+1),
                               Dropout=estParamRes$p0)
    meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat, ggplot2::aes(x=Means, y=Dropout)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25,
                                               alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+
      ggplot2::ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none',
                     axis.text=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=12, face="bold"))
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
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=10, face="bold")) +
      ggplot2::labs(x="Sequencing Depth", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estParamRes$sf, Sample=names(estParamRes$sf))
    sf.plot <- ggplot2::ggplot(sf.dat, ggplot2::aes(SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=10, face="bold")) +
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
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=10, face="bold"),
                     strip.text = ggplot2::element_text(size=10, face="bold")) +
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
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
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
      ggplot2::theme(legend.position='none',
                     axis.text=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=12, face="bold"))
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
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
      ggplot2::ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none',
                     axis.text=ggplot2::element_text(size=10),
                     axis.title=ggplot2::element_text(size=12, face="bold"))

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
#' ## for example simres object see \code{\link{simulateDE}}
#' evalres <- evaluateDE(simRes=simres,
#' alpha.type="adjusted",
#' MTC="BH", alpha.nominal=0.1,
#' stratify.by="mean", filter.by="none",
#' target.by="lfc", delta=0)
#' plotEvalRes(evalRes=evalres, rate ="marginal", quick=T, annot=T)
#' plotEvalRes(evalRes=evalres, rate ="stratified", quick=T, annot=T)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar
#' @importFrom ggthemes theme_base
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom tidyr %>%
#' @rdname plotEvalDE
#' @export
plotEvalDE <- function(evalRes, rate=c('marginal', 'stratified'), quick=TRUE, annot=TRUE) {

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
      grandplot <- ggplot2::ggplot(data = dat.marginal.long, ggplot2::aes(x=Var1, y=value, fill=L1)) +
        ggplot2::geom_boxplot() +
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
      grandplot <- ggplot(data = dat.marginal.long, ggplot2::aes(x=Var1, y=value, fill=L1)) +
        ggplot2::geom_boxplot() +
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
        ggthemes::theme_base() +
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
#' ## not yet
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
#' ## not yet
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



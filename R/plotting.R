
# plotParam -------------------------------------------------------------
#' @name plotParam
#' @aliases plotParam
#' @title Visualize
#' @description This function plots the results of the parameter estimation. This includes the absolute and relative sequencing depth (i.e. library size factor) and marginal mean, dispersion and dropout. Also the mean-dispersion relationship with loess fit for simulations is visualized. Lastly the mean-dropout rate is presented as a smooth scatter plot.
#' @usage plotParam(estParam.out, annot=TRUE)
#' @param estParam.out The output of \code{\link{estimateParam}} or a string specifying the name of precalculated estimates (see details).
#' @param annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## for example see \code{\link{insilicoNBParam}}
#' }
#' @details Precalculated negative binomial parameters for simulations are available via \url{}.
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot geom_bar geom_line theme geom_hline xlab ylab coord_flip scale_y_continuous geom_point scale_fill_gradientn stat_density2d labs theme_minimal theme_classic
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom ggExtra ggMarginal
#' @importFrom ggthemes theme_base
#' @importFrom stats reorder
#' @rdname plotParam
#' @export
plotParam <- function(estParam.out, annot=TRUE) {

  if(estParam.out$RNAseq=="bulk" | estParam.out$estS<15) {
    # library size
    lib.size.dat <- data.frame(Seqdepth=estParam.out$seqDepth,
                               Sample=names(estParam.out$seqDepth))
    libsize.plot <- ggplot2::ggplot(data=lib.size.dat, aes(reorder(Sample, Seqdepth),Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold")) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Sequencing depth") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estParam.out$sf,
                         Sample=names(estParam.out$sf))
    sf.plot <- ggplot2::ggplot(data=sf.dat, aes(reorder(Sample, SizeFactor),SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold")) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Library Size Factor") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip()
    # marginal distributions
    margs.dat <- data.frame(Mean=log2(estParam.out$means+1),
                            Dispersion=log2(estParam.out$dispersion),
                            Dropout=estParam.out$p0)
    margs.dat <- suppressMessages(reshape2::melt(margs.dat))
    margs.plot <- ggplot2::ggplot(margs.dat, aes(value)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"), strip.text = element_text(size=8, face="bold")) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Density") +
      ggplot2::facet_wrap(~variable, scales = 'free')
    # mean vs dispersion plot
    meanvsdisp.dat <- data.frame(Means=log2(estParam.out$means+1),
                                 Dispersion=log2(estParam.out$dispersion))
    meanvsdisp.plot <- ggplot2::ggplot(data=meanvsdisp.dat, aes(x=Means, y=Dispersion)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
      ggplot2::geom_hline(yintercept = log2(estParam.out$common.dispersion), linetype = 2, colour="grey40") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$y), linetype=1, size=1.5, colour="orange") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$upper), linetype=2, size=1, colour="orange") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$lower), linetype=2, size=1, colour="orange") +
      ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))), x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"))
    # mean vs p0 plot
    meanvsp0.dat <- data.frame(Means=log2(estParam.out$means+1),
                               Dropout=estParam.out$p0)
    meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat, aes(x=Means, y=Dropout)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+
      ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"))
    if(estParam.out$RNAseq=="bulk" && !is.null(estParam.out$p0.cut)) {
      meanvsp0.plot <- meanvsp0.plot + ggplot2::annotate("rect", xmin=0, ymax=1, ymin=0, xmax=estParam.out$p0.cut+1, fill="red", alpha=0.2)
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

  if(estParam.out$RNAseq=="singlecell" | estParam.out$estS>15) {
    # library size
    lib.size.dat <- data.frame(Seqdepth=estParam.out$seqDepth,
                               Sample=names(estParam.out$seqDepth))
    libsize.plot <- ggplot2::ggplot(lib.size.dat, aes(Seqdepth)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(lib.size.dat$Seqdepth), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold")) +
      ggplot2::xlab("Sequencing Depth") +
      ggplot2::ylab("Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # size factor plot
    sf.dat <- data.frame(SizeFactor=estParam.out$sf, Sample=names(estParam.out$sf))
    sf.plot <- ggplot2::ggplot(sf.dat, aes(SizeFactor)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold")) +
      ggplot2::xlab("Library Size Factor") +
      ggplot2::ylab("Density") +
      ggplot2::scale_x_continuous(labels=.plain)
    # marginal distributions
    margs.dat <- data.frame(Mean=log2(estParam.out$means+1),
                            Dispersion=log2(estParam.out$dispersion),
                            Dropout=estParam.out$p0)
    margs.dat <- suppressMessages(reshape2::melt(margs.dat))
    margs.plot <- ggplot2::ggplot(margs.dat, aes(value)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_density() +
      ggplot2::theme(axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"), strip.text = element_text(size=8, face="bold")) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Density") +
      ggplot2::facet_wrap(~variable, scales = 'free')
    # mean vs dispersion plot
    meanvsdisp.dat <- data.frame(Means=log2(estParam.out$means+1),
                                 Dispersion=log2(estParam.out$dispersion))
    meanvsdisp.plot <- ggplot2::ggplot(data=meanvsdisp.dat, aes(x=Means, y=Dispersion)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
      ggplot2::geom_hline(yintercept = log2(estParam.out$common.dispersion), linetype = 2, colour="grey40") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$y), linetype=1, size=1.5, colour="orange") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$upper), linetype=2, size=1, colour="orange") +
      ggplot2::geom_line(aes(x=estParam.out$meandispfit$x, y=estParam.out$meandispfit$lower), linetype=2, size=1, colour="orange") +
      ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))), x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"))
    # mean vs p0 plot
    meanvsp0.dat <- data.frame(Means=log2(estParam.out$means+1),
                                Dropout=estParam.out$p0)
    meanvsp0.plot <- ggplot2::ggplot(data=meanvsp0.dat, aes(x=Means, y=Dropout)) +
      ggplot2::theme_classic() +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + ggplot2::geom_point(size=0.5) +
      ggplot2::stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=ifelse(..density..^0.15<0.4,0,1)), contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+ ylim(c(0,1)) +
      ggplot2::labs(y="Dropout Fraction", x=expression(bold(paste(Log[2], " (Mean)")))) +
      ggplot2::theme(legend.position='none', axis.text=element_text(size=6), axis.title=element_text(size=8, face="bold"))

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


# plotEvalRes -------------------------------------------------------------
#' @name plotEvalRes
#' @aliases plotEvalRes
#' @title Visualize power assessment
#' @description This function plots the results of \code{\link{evaluateSim}} for assessing the error rates and sample size requirements.
#' @usage plotEvalRes(evalRes, rate=c('marginal', 'stratified'),
#'                    quick=TRUE, annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateSim}}.
#' @param rate Character vector defining whether the marginal or condtional rates should be plotted. Conditional depends on the choice of stratify.by in \code{\link{evaluateSim}}.
#' @param quick A logical vector. If \code{TRUE}, the TPR and FDR are only plotted. If \code{FALSE}, then all rates are plotted.
#' @param annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## for example simres object see \code{\link{simulateDE}}
#' evalres <- evaluateSim(simRes=simres,
#' alpha.type="adjusted",
#' MTC="BH", alpha.nominal=0.1,
#' stratify.by="mean", filter.by="none",
#' target.by="lfc", delta=0)
#' plotEvalRes(evalRes=evalres, rate ="marginal", quick=T, annot=T)
#' plotEvalRes(evalRes=evalres, rate ="stratified", quick=T, annot=T)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes ylab xlab theme scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar
#' @importFrom ggthemes theme_base
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom tidyr %>%
#' @rdname plotEvalRes
#' @export
plotEvalRes <- function(evalRes, rate=c('marginal', 'stratified'), quick=TRUE, annot=TRUE) {

  # marginal rates over sample sizes
  if(rate=='marginal') {
    if(quick) {
      dat.marginal <- evalRes[c('TPR.marginal', 'FDR.marginal')]
      names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
      dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalRes[['n1']], " vs ", evalRes[['n2']]))
      dat.marginal.long <- reshape2::melt(dat.marginal)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      dat.marginal.calc <- dat.marginal.long %>% dplyr::group_by(Var1, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      # marginal in one
      grandplot <- ggplot2::ggplot(data = dat.marginal.long, aes(x=Var1, y=value, fill=L1)) + ggplot2::geom_boxplot() + ggthemes::theme_base() + ggplot2::xlab(NULL) + ggplot2::ylab('Rate') + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) + ggplot2::theme(legend.position='top', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"))
      # faceted marginal
      facetplot <-  ggplot2::ggplot(data = dat.marginal.calc, aes(x=Var1, y=Expectation, fill=L1, color=L1)) + ggplot2::geom_line(aes(group=L1)) + ggplot2::geom_pointrange(limits) + ggthemes::theme_base() + ggplot2::xlab('Samples') + ggplot2::ylab('Rate') + ggplot2::scale_y_continuous(labels = scales::percent) + ggplot2::theme(legend.position='none', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10, face="bold")) + ggplot2::facet_wrap(~L1, scales = 'free', ncol=1) + ggplot2::geom_hline(data = refval, aes(yintercept = ref), linetype="dashed", color='grey')
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(grandplot, facetplot, labels=c('A', 'B'), rel_heights = c(1,1.5), ncol=1, nrow=2))
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
      dat.marginal.calc <- dat.marginal.long %>% dplyr::group_by(Var1, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      # marginal in one
      grandplot <- ggplot(data = dat.marginal.long, aes(x=Var1, y=value, fill=L1)) + ggplot2::geom_boxplot() + ggthemes::theme_base() + ggplot2::xlab(NULL) + ggplot2::ylab('Rate') + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) + ggplot2::theme(legend.position='top', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"))
      # faceted marginal
      facetplot <-  ggplot2::ggplot(data = dat.marginal.calc, aes(x=Var1, y=Expectation, fill=L1, color=L1)) + ggplot2::geom_line(aes(group=L1)) + ggplot2::geom_pointrange(limits) + ggthemes::theme_base() + ggplot2::xlab('Samples') + ggplot2::ylab('Rate') + ggplot2::scale_y_continuous(labels = scales::percent) + ggplot2::theme(legend.position='none', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10, face="bold")) + ggplot2::facet_wrap(~L1, scales = 'free', ncol=2) + ggplot2::geom_hline(data = refval, aes(yintercept = ref), linetype="dashed", color='grey')
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(grandplot, facetplot, labels=c('A', 'B'), rel_heights = c(1,2), ncol=1, nrow=2))
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
      dat.stratified.calc <- dat.stratified.long %>% dplyr::group_by(Var1, Var2, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      limits <- aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      facetplot <-  ggplot2::ggplot(data = dat.stratified.calc, aes(x=Var1, y=Expectation, fill=Var2, color=Var2)) + ggplot2::geom_point() + ggplot2::geom_line(aes(group=Var2)) +  ggplot2::geom_pointrange(limits) + ggthemes::theme_base() + ggplot2::xlab(NULL) + ggplot2::ylab('Rate') + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) + ggplot2::theme(legend.position='top', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10, face="bold")) + ggplot2::facet_wrap(~L1, scales = 'free', ncol=2) + ggplot2::geom_hline(data = refval, aes(yintercept = ref), linetype="dashed", color='grey')
      # strata genes
      N <- length(evalRes$n1)
      dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],'DEgenes'=evalRes$stratadiffgenes[,N,])
      dat.genes <- lapply(dat.genes, "rownames<-", strata)
      dat.genes.long <- reshape2::melt(dat.genes)
      dat.genes.calc <- dat.genes.long %>% dplyr::group_by(Var1, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      dodge <- ggplot2::position_dodge(width=0.9)
      strataplot <- ggplot2::ggplot(data = dat.genes.calc, aes(x=Var1, y=Expectation, fill=L1)) + ggplot2::geom_bar(stat="identity")  + ggthemes::theme_base() + ggplot2::xlab('Stratum') + ggplot2::ylab('Count') + ggplot2::theme(legend.position='right', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10)) + ggplot2::scale_fill_manual(values=c('grey', 'black'),breaks = c("DEgenes", "Ngenes"), labels = c("DE genes", "EE genes"))
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(facetplot, strataplot, labels=c('A', 'B'), rel_heights = c(2,1), ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Conditional TPR and FDR per sample size comparison per stratum. \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.", size=8)
      }
    }
    if(!quick) {
      dat.stratified <- evalRes[grep('*R$', names(evalRes))]
      strata <- evalRes$strata.levels
      dat.stratified <- lapply(dat.stratified, "dimnames<-", list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]), NULL))
      dat.stratified.long <- reshape2::melt(dat.stratified)
      dat.stratified.calc <- dat.stratified.long %>% dplyr::group_by(Var1, Var2, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      limits <- ggplot2::aes(ymax = Expectation + Deviation, ymin= Expectation - Deviation)
      refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))
      facetplot <-   ggplot2::ggplot(data = dat.stratified.calc, aes(x=Var1, y=Expectation, fill=Var2, color=Var2)) +  ggplot2::geom_point() + ggplot2::geom_line(aes(group=Var2)) +   ggplot2::geom_pointrange(limits) + ggthemes::theme_base() +  ggplot2::xlab(NULL) +  ggplot2::ylab('Rate') +  ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  ggplot2::theme(legend.position='top', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10, face="bold")) +  ggplot2::facet_wrap(~L1, scales = 'free_x', ncol=2) +  ggplot2::geom_hline(data = refval, aes(yintercept = ref), linetype="dashed", color='grey')
      # strata genes
      N <- length(evalRes$n1)
      dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],'DEgenes'=evalRes$stratadiffgenes[,N,])
      dat.genes <- lapply(dat.genes, "rownames<-", strata)
      dat.genes.long <- reshape2::melt(dat.genes)
      dat.genes.calc <- dat.genes.long %>% dplyr::group_by(Var1, L1) %>% dplyr::summarise(Expectation=mean(value), Deviation=sd(value), Error=sd(value)/sqrt(n())) %>% dplyr::ungroup()
      dodge <-  ggplot2::position_dodge(width=0.9)
      strataplot <-  ggplot2::ggplot(data = dat.genes.calc, aes(x=Var1, y=Expectation, fill=L1)) +  ggplot2::geom_bar(stat="identity")  + ggthemes::theme_base() +  ggplot2::xlab('Stratum') +  ggplot2::ylab('Count') +  ggplot2::theme(legend.position='right', legend.title = element_blank(), axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10), axis.title=element_text(size=10, face="bold"), legend.text = element_text(size=10), legend.key.size = unit(1, "cm"), strip.text = element_text(size=10)) +  ggplot2::scale_fill_manual(values=c('grey', 'black'),breaks = c("DEgenes", "Ngenes"), labels = c("DE genes", "EE genes"))
      # annotation under plot
      p.final <- suppressWarnings(cowplot::plot_grid(facetplot, strataplot, labels=c('A', 'B'), rel_heights = c(3,1), ncol=1, nrow=2))
      if(annot) {
        p.final <- cowplot::add_sub(p.final, "A) Conditional error rates over stratum. \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.", size=8)
      }
    }
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
#' ## for example see \code{\link{evaluateDist}}
#' @author Beate Vieth, Ines Hellmann
#' @importFrom ggplot2 ggplot aes ylab xlab theme scale_y_continuous geom_boxplot position_dodge geom_bar theme_minimal scale_x_discrete labs ggtitle coord_flip
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr filter group_by summarise mutate bind_rows slice select ungroup  left_join
#' @importFrom tidyr %>% separate gather spread
#' @rdname plotEvalDist
#' @export
plotEvalDist <- function(evalDist, annot=TRUE){
  # TODO use NBGOF results for bulk data (ie small number of samples) and nsims bigger than 1
  # extract the GOF results from object
  label_names<-c('pois' = "Poisson",
                 "nbinom" = "Negative Binomial",
                 'zifpois' = "Zero-Inflated \n Poisson",
                 "zifnbinom" = "Zero-Inflated \n Negative Binomial",
                 "PoiBeta" ="Beta-Poisson")
  gofres <- evalDist$GOF_res
  # reshape the table
  gofres <- cbind(rn = rownames(gofres), stack(gofres))
  gofres <- gofres %>% tidyr::separate(ind, c("distribution", "framework", 'type'), "_", remove = F)
  # extract the observed zero table
  obszero <- evalDist$ObservedZeros
  obszero <- cbind(rn = rownames(obszero), obszero)
  obszero <- obszero[which(obszero$rn %in% gofres$rn), ]
  # GOF p-value based on chisquare test
  gof.pval <- gofres %>% dplyr::filter(type=='gofpval', framework=='standard') %>%
    dplyr::group_by(distribution) %>%
    dplyr::summarise(percgenes=sum(values>0.05,na.rm = T)/
                       sum(is.na(values)==F) )
  # barplot of % acceptable fit (GOF chisquare test above 5%)
  gofpval.plot <- ggplot2::ggplot(data=gof.pval, aes(x=distribution,y=percgenes)) +
    ggplot2::geom_bar(stat="identity", width = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::labs(x = "Distribution", y = "Percentage") +
    ggplot2::scale_x_discrete(labels=label_names, limits=names(label_names)) +
    ggplot2::ggtitle("Goodness-of-fit statistic")+
    ggplot2::theme(axis.title.y = element_blank())

  # AIC (lowest value AND intersect with GOF p-value)
  AIC.calc <- gofres %>%
    dplyr::filter(framework %in%c('Marioni', 'standard'), type %in% c("gofpval","aic")) %>%
    dplyr::select(-ind ) %>%
    tidyr::spread(type,values) %>%
    dplyr::group_by(rn) %>%
    dplyr::mutate(minAIC= (aic==min(aic)), GOF = gofpval>0.05) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(distribution) %>%
    dplyr::summarise(`Lowest\nAIC`= sum(minAIC,na.rm = T),
                     `Lowest\nAIC + \nGOF \np >0.05` = sum((minAIC & GOF), na.rm=T)) %>%
    tidyr::gather(Set,totalgenes,2:3) %>%
    dplyr::group_by(Set) %>%
    dplyr::mutate(Percentage=totalgenes/sum(totalgenes))
  # barplot of lowest AIC with GOF test passed
  AIC.plot <- ggplot2::ggplot(data=AIC.calc, aes(x=distribution,y=percgenes, fill=Set)) +
    ggplot2::geom_bar(stat="identity", width = 0.7,position="dodge") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::labs(x = "Distribution", y = "Percentage")  +
    ggplot2::scale_x_discrete(labels= label_names, limits=names(label_names)) +
    ggplot2::theme(legend.position = 'left', legend.title = element_blank(),axis.title.y = element_blank()) +
    ggplot2::ggtitle('Akaike Information Criterion')

  # predicted zeroes vs observed zeros
  predzero.dat <- gofres %>% dplyr::filter(type %in% c("gofpval","predzero"), framework %in% c("Marioni", 'standard')) %>%
    dplyr::select(-ind ) %>%
    tidyr::spread(type,values) %>%
    dplyr::left_join(obszero, by='rn') %>%
    dplyr::mutate(`Obs. vs. Pred.\n Zeros`=ObsZero-predzero,
                  `Obs. vs. Pred.\n Zeros, GOF p>0.05`=ifelse(gofpval>0.05,ObsZero-predzero,NA)) %>%
    tidyr::gather(comp,divergentzero,7:8)

  dd=1:5
  zpdat<- dplyr::filter(predzero.dat,comp=="Obs. vs. Pred.\n Zeros")
  zero.plot <- ggplot2::ggplot( zpdat,aes(x=distribution, y=divergentzero)) +
    geom_boxplot( width = 0.7, position=position_dodge(width=0.8), outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Distribution", y = "Observed - Predicted Zeros")  +
    ggplot2::scale_x_discrete(labels=label_names[dd], limits=names(label_names[dd]))  +
    ggplot2::ggtitle('Dropouts')+
    ggplot2::theme(legend.position = 'none', legend.title = element_blank(),axis.title.y = element_blank())
  # Best fit by LRT / Vuong
  ModelTest <- gofres %>% dplyr::filter(distribution %in% c('LRT', 'Vuong')) %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(percgenes=sum(values<0.05,na.rm = T)/
                       sum(is.na(values)==F) )

  ModelTest.plot <- ggplot2::ggplot(data=ModelTest, aes(x=type,y=percgenes)) +
    ggplot2::geom_bar(stat="identity", width = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::labs(x = "Test", y = "Percentage") +
    ggplot2::scale_x_discrete(labels=c('NBPoisson' = "Negative-Binomial \n> Poisson",
                                       "ZNB" = "Zero-inflated Negative-Binomial \n> Negative-Binomial",
                                       'ZPoisson' = "Zero-Inflated Poisson \n> Poisson",
                                       "ZNBZPoisson" = "Zero-Inflated Negative Binomial \n> Zero-Inflated Poisson"),
                              limits=c('NBPoisson', 'ZNB','ZPoisson', 'ZNBZPoisson')) +
    ggplot2::ggtitle('Model Comparisons')+
    ggplot2::theme(axis.title.y = element_blank())
  # combine into one plot
  ### flip Plots
  ModelTest.plot <- ModelTest.plot + ggplot2::coord_flip()
  gofpval.plot <- gofpval.plot + ggplot2::coord_flip()
  AIC.plot    <- AIC.plot + ggplot2::coord_flip()
  zero.plot   <- zero.plot + ggplot2::coord_flip()

  p.final <- suppressWarnings(cowplot::plot_grid(gofpval.plot, AIC.plot, zero.plot, ModelTest.plot, labels=c('A', 'B', 'C', 'D'), ncol=2, nrow=2))

  # annotation under plot
  if(annot) {
    p.final <- cowplot::add_sub(p.final, "A) Goodness-of-fit of the model assessed with a chi-square test based on residual deviance and degrees of freedom. \n
                                          B) Akaike Information Criterion per gene: Model with the lowest AIC. Model with the lowest AIC and passed goodness-of-fit statistic test. \n
                                         C) Observed versus predicted dropouts per model and gene plotted without outliers. \n
                                         D) Model Assessment based on LRT for nested models and Vuong test for nonnested models.",
                                size=8)
  }
  # draw final plot
  cowplot::ggdraw(p.final)
}

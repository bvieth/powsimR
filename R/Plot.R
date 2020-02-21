
# plotParam -------------------------------------------------------------
#' @name plotParam
#' @aliases plotParam
#' @title Visualize distributional characteristics of RNA-seq experiment
#' @description This function plots the results of the parameter estimation. This includes the absolute and relative sequencing depth (i.e. library size factor) as well as marginal log mean, log dispersion and dropout. Furthermore, the mean-dispersion relationship with loess fit for simulations is visualized. Lastly, the mean-dropout rate is presented as a smooth scatter plot.
#' @usage plotParam(estParamRes, Annot=TRUE)
#' @param estParamRes The output of \code{\link{estimateParam}}.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # using example data set
#' data("CELseq2_Gene_UMI_Counts")
#' data("CELseq2_Gene_Read_Counts")
#' Batches <- data.frame(Batch = sapply(strsplit(colnames(CELseq2_Gene_UMI_Counts), "_"), "[[", 1),
#'                   stringsAsFactors = FALSE, row.names = colnames(CELseq2_Gene_UMI_Counts))
#' data("GeneLengths_mm10")
#' data("CELseq2_SpikeIns_UMI_Counts")
#' data("CELseq2_SpikeInfo")
#' # estimation
#' estparam <-  estimateParam(countData = CELseq2_Gene_UMI_Counts,
#' readData = CELseq2_Gene_Read_Counts,
#' batchData = Batches,
#' spikeData = NULL,
#' spikeInfo = NULL,
#' Lengths = GeneLengths_mm10,
#' MeanFragLengths = NULL,
#' Distribution = 'NB',
#' RNAseq = 'singlecell',
#' Protocol = 'UMI',
#' Normalisation = 'scran',
#' GeneFilter = 0.1,
#' SampleFilter = 3,
#' sigma = 1.96,
#' NCores = NULL,
#' verbose = TRUE)
#' # plotting
#' plotParam(estParamRes = estparam, Annot=TRUE)
#' }
#' @author Beate Vieth
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom stats reorder
#' @rdname plotParam
#' @export
plotParam <- function(estParamRes, Annot=TRUE) {

  # QC plots
  seqdepth.plot <- .seqdepth_plot(estParamRes=estParamRes)
  libsize.plot <- .sf_plot(estParamRes=estParamRes)
  genefeatures.plot <- .feat_plot(estParamRes=estParamRes)
  if(all(!is.na(estParamRes$DropOuts$Sample$GeneSpikeRatio))){
    gs.plot <- .genespike_ratio_plot(estParamRes)
  } else {
    gs.plot <- NULL
  }

  # Marginal Distributions
  margs.plot <- .margs_plot(estParamRes = estParamRes,
                            Distribution = attr(estParamRes, "Distribution"))

  # Table with numbers of genes / samples
  no.table <- .estimate_table_print(estParamRes = estParamRes,
                                    RNAseq = attr(estParamRes, "RNAseq"))

 # Fitting Lines
 # mean-disp and mean-p0
 meanvsdisp.plot <- .meanvsdisp_plot(estParamRes = estParamRes,
                                     Distribution = attr(estParamRes, "Distribution"))
 meanvsdrop.plot <- .meanvsdrop_plot(estParamRes = estParamRes,
                                     Distribution = attr(estParamRes, "Distribution"),
                                     RNAseq = attr(estParamRes, 'RNAseq'))

 # read-umi
 if(all(!is.na(estParamRes$Fit$UmiRead))) {
   readvsumi.plot <- .readumi_plot(estParamRes = estParamRes)
 } else {
   readvsumi.plot <- NULL
 }

 # combine the plots into one output
 if(!is.null(gs.plot)){
   top_row <- suppressWarnings(cowplot::plot_grid(seqdepth.plot,
                                                  libsize.plot,
                                                  genefeatures.plot,
                                                  gs.plot,
                                                  align = 'hv',
                                                  ncol=4, nrow=1))
 }
 if(is.null(gs.plot)){
   top_row <- suppressWarnings(cowplot::plot_grid(seqdepth.plot,
                                                  libsize.plot,
                                                  genefeatures.plot,
                                                  align = 'hv',
                                                  ncol=3, nrow=1))
 }

  middle_row <- suppressWarnings(cowplot::plot_grid(margs.plot,
                                                    no.table,
                                                    labels=c('B', 'C'),
                                                    rel_widths = c(0.7,0.3),
                                                    ncol=2,
                                                    nrow=1))
  if(!is.null(readvsumi.plot)){
    bottom_row <- suppressWarnings(cowplot::plot_grid(meanvsdisp.plot,
                                                      meanvsdrop.plot,
                                                      readvsumi.plot,
                                                      labels=c('D', 'E', 'F'),
                                                      ncol=3, nrow=1))
  }
  if(is.null(readvsumi.plot)){
    bottom_row <- suppressWarnings(cowplot::plot_grid(meanvsdisp.plot,
                                                      meanvsdrop.plot,
                                                      labels=c('D', 'E'),
                                                      ncol=2, nrow=1))
  }

  p.final <- suppressWarnings(cowplot::plot_grid(top_row,
                                                 middle_row,
                                                 bottom_row,
                                                 labels=c('A', NULL, NULL),
                                                 ncol=1, nrow=3))
  # Annotation under plot
  if (Annot) {
    annottext.a <-  c("A) Quality Control Metrics: Sequencing depth; Library size factors with median (black line) for the filtered data set; Detected genes; Ratio of gene to spike-in counts (if spike-ins were provided). Outliers are marked in red.")
    annottext.b <- c("\nB) Marginal Distribution of gene mean, dispersion and dropout rate per estimation set.")
    annottext.c <- c("\nC) Number of genes and samples per estimation set. Provided by the user; Detected = number of genes and samples with at least one count; All = number of genes for which mean, dispersion and dropout could be estimated using non-outlying samples. \nFiltered = number of genes above filter threshold for which mean, dispersion and dropout could be estimated using non-outlying samples. Dropout Genes = number of genes filtered out due to gene dropout rate.")
    annottext.d <- c("\nD) Local polynomial regression fit between mean and dispersion estimates with variability band per gene (yellow). Common dispersion estimate (grey dashed line).")
    annottext.e <- c("\nE) Fraction of dropouts versus estimated mean expression per gene.")
    if(is.null(readvsumi.plot)) {
      annottext.f <- NULL
    }
    if(!is.null(readvsumi.plot)) {
      annottext.f <-     annottext.f <- c("\nF) Local polynomial regression fit between UMI and read counts with variability band per gene (yellow).")
    }

    annottext <- paste0(annottext.a, annottext.b, annottext.c, annottext.d, annottext.e, annottext.f)
    p.final <- cowplot::add_sub(plot = p.final,
                                label = annottext,
                                x = 0, hjust = 0, size=8)
  }

  # draw the plot
  cowplot::ggdraw(p.final)
}

# plotSpike ---------------------------------------------------------------

#' @name plotSpike
#' @aliases plotSpike
#' @title Visualize distributional characteristics of spike-ins
#' @description This function plots the results of the parameter estimation for spike-ins. This includes the absolute and relative sequencing depth (i.e. library size factor), a calibration curve as well as the capture efficiency given as a binomial regression.
#' @usage plotSpike(estSpike, Annot = TRUE)
#' @param estSpike The output of \code{\link{estimateSpike}}.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # using example data set
#' data("SCRBseq_SpikeIns_Read_Counts")
#' data("SCRBseq_SpikeInfo")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SCRBseq_SpikeIns_Read_Counts), "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SCRBseq_SpikeIns_Read_Counts))
#' # estimation
#' spikeparam <- estimateSpike(spikeData = SCRBseq_SpikeIns_Read_Counts,
#'                              spikeInfo = SCRBseq_SpikeInfo,
#'                              MeanFragLength = NULL,
#'                              batchData = Batches,
#'                              Normalisation = 'depth')
#' # plotting
#' plotSpike(estSpike = spikeparam, Annot = TRUE)
#' }
#' @author Beate Vieth
#' @importFrom ggplot2 ggplot aes_ theme_minimal geom_bar geom_density geom_hline theme labs scale_y_continuous coord_flip geom_pointrange geom_point geom_smooth annotate scale_x_log10 scale_y_log10 annotation_logticks
#' @importFrom dplyr left_join group_by mutate ungroup do summarise n
#' @importFrom tidyr "%>%"
#' @importFrom broom glance
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom stats reorder
#' @rdname plotSpike
#' @export
plotSpike <- function(estSpike, Annot = TRUE) {

  # library size plot
  lib.size.dat <- data.frame(Seqdepth=estSpike$seqDepth,
                             Sample=names(estSpike$seqDepth))
  lib.size.dat <- lib.size.dat[order(lib.size.dat$Seqdepth),]
  lib.size.dat$Sample <- factor(lib.size.dat$Sample, levels = lib.size.dat$Sample)

  # size factor plot
  sf.dat <- data.frame(SizeFactor=estSpike$size.factors,
                       Sample=names(estSpike$size.factors))
  sf.dat <- sf.dat[order(sf.dat$SizeFactor),]
  sf.dat$Sample <- factor(sf.dat$Sample, levels = sf.dat$Sample)

  if(length(estSpike$seqDepth)<15) {
    # library size plot
    libsize.plot <- ggplot2::ggplot(data=lib.size.dat,
                                    ggplot2::aes_(x = quote(Sample),
                                                  y = quote(Seqdepth))) +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(lib.size.dat$Seqdepth),
                          linetype = 2, colour="grey40") +
      ggplot2::labs(x=NULL, y="Sequencing Depth") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip() +
      .theme_eval_time()
    # size factor plot
    sf.plot <- ggplot2::ggplot(data=sf.dat,
                               ggplot2::aes_(x = quote(Sample),
                                             y = quote(SizeFactor))) +
      ggplot2::geom_bar(stat="identity",width=.5) +
      ggplot2::geom_hline(yintercept = median(sf.dat$SizeFactor),
                          linetype = 2, colour="grey40") +
      ggplot2::labs(x=NULL, y="Library Size Factor") +
      ggplot2::scale_y_continuous(labels=.plain) +
      ggplot2::coord_flip() +
      .theme_eval_time()
  }
  if(length(estSpike$seqDepth)>=15) {
    # library size plot
    libsize.plot <- ggplot2::ggplot(data = lib.size.dat,
                                    ggplot2::aes_(quote(Seqdepth))) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(lib.size.dat$Seqdepth),
                          linetype = 2, colour="grey40") +
      ggplot2::labs(x="Sequencing Depth", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain) +
      .theme_eval_time()
    # size factor plot
    sf.plot <- ggplot2::ggplot(data = sf.dat,
                               ggplot2::aes_(quote(SizeFactor))) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = median(sf.dat$SizeFactor), linetype = 2, colour="grey40") +
      ggplot2::labs(x="Library Size Factor", y="Density") +
      ggplot2::scale_x_continuous(labels=.plain) +
      .theme_eval_time()
  }

  # calibration curve data
  cal.dat <- reshape2::melt(estSpike$normCounts)
  names(cal.dat) <- c("SpikeID", "SampleID", "normCounts")
  cal.info.dat <- cal.dat %>%
    dplyr::left_join(estSpike$FilteredInput$spikeInfo, by="SpikeID") %>%
    dplyr::mutate(FSpike = factor(.data$SpikeInput)) %>%
    dplyr::group_by(.data$FSpike) %>%
    dplyr::mutate(Expectation=mean(.data$normCounts),
                  LExpectation=mean(log10(.data$normCounts)),
                  Deviation=sd(.data$normCounts),
                  Error=sd(.data$normCounts)/sqrt(dplyr::n()),
                  LError=sd(log10(.data$normCounts))/sqrt(dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LSpikeInput = log10(.data$SpikeInput),
                  LExpectation = log10(.data$Expectation))

  LmFit <- NULL
  Calibration <- cal.dat %>%
    dplyr::left_join(estSpike$FilteredInput$spikeInfo, by="SpikeID") %>%
    dplyr::group_by(.data$SampleID) %>%
    dplyr::do(LmFit = lm(log10(normCounts+1) ~ log10(SpikeInput+1), data = .data)) %>%
    broom::glance(LmFit) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(Rsquared=mean(.data$r.squared),
                     RsquaredSE=sd(.data$r.squared))

  # calibration curve plot
  cal.plot <- ggplot2::ggplot(data = cal.info.dat,
                              ggplot2::aes_(x=quote(LSpikeInput),
                                           y=quote(LExpectation))) +
    ggplot2::geom_pointrange(data = cal.info.dat,
                             ggplot2::aes_(ymax = ~ LExpectation + LError,
                                           ymin = ~ LExpectation - LError)) +
    ggplot2::geom_point(data = cal.info.dat,
                        ggplot2::aes_(x=quote(LSpikeInput),
                                      y=quote(LExpectation))) +
    ggplot2::geom_smooth(method='lm',formula=y~x) +
    ggplot2::annotate("text", label = paste0("italic(R) ^ 2 == ",
                                             round(Calibration$Rsquared, digits = 2),
                                             "%+-%",
                                             round(Calibration$RsquaredSE, digits = 2)),
                      parse = T, x = 0.2, y = 4, size = 4) +
    ggplot2::annotation_logticks(sides = "bl") +
    ggplot2::labs(y=expression(bold(paste(Log[10], " Estimated Expression", sep=""))),
                  x=expression(bold(paste(Log[10], " Spike-In Molecules")))) +
    .theme_eval_time()

  # capture efficiency data
  capture.dat <- estSpike$CaptureEfficiency$`Spike-In` %>%
    tibble::rownames_to_column(var = "SpikeID") %>%
    dplyr::select(.data$SpikeID, .data$p_success,
                  .data$hat_p_success_cilower, .data$hat_p_success_ciupper) %>%
    dplyr::left_join(estSpike$FilteredInput$spikeInfo, by="SpikeID") %>%
    dplyr::mutate(LSpikeInput = log10(.data$SpikeInput))

  capture.plot <- ggplot2::ggplot(data = capture.dat,
                                  ggplot2::aes_(x=quote(LSpikeInput),
                                               y=quote(p_success))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="glm",
                         method.args = list(family = "binomial"), se=T) +
    ggplot2::scale_x_log10(labels=c("0.1","1","10","100","1,000"),breaks=c(0.1,1,10,100,1000)) +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::labs(y=expression(bold("Detection Probability")),
                  x=expression(bold(paste(Log[10], " Spike-In Molecules")))) +
    .theme_eval_time()

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
  if (Annot) {
    annot.text <- c("A) Sequencing depth per sample with median sequencing depth (grey dashed line). \nB) Library size normalisation factor per sample with median size factor (grey dashed line). \nC) Calibration curve with mean expression estimates and average R squared over all cells. \nD) Capture efficiency with binomial logistic regression fit over all cells.")
    p.final <- cowplot::add_sub(plot = p.final,
                                label = annot.text,
                                x = 0, hjust = 0, size=8)

  }
  # draw the plot
  cowplot::ggdraw(p.final)
}

# plotEvalDE -------------------------------------------------------------

#' @name plotEvalDE
#' @aliases plotEvalDE
#' @title Visualize power assessment
#' @description This function plots the results of \code{\link{evaluateDE}} for assessing the error rates and sample size requirements.
#' @usage plotEvalDE(evalRes, rate=c('marginal', 'conditional'),
#'                    quick=TRUE, Annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateDE}}.
#' @param rate Character vector defining whether the \code{"marginal"} or \code{"conditional"} rates should be plotted. Conditional depends on the choice of stratify.by in \code{\link{evaluateDE}}.
#' @param quick A logical vector. If \code{TRUE}, the TPR and FDR are only plotted. If \code{FALSE}, then all rates are plotted.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # estimate gene parameters
#' data("Bulk_Read_Counts")
#' data("GeneLengths_hg19")
#' estparam_gene <- estimateParam(countData = Bulk_Read_Counts,
#'                                readData = NULL,
#'                                batchData = NULL,
#'                                spikeData = NULL, spikeInfo = NULL,
#'                                Lengths = GeneLengths_hg19, MeanFragLengths = NULL,
#'                                RNAseq = 'bulk', Protocol = 'Read',
#'                                Distribution = 'NB', Normalisation = "MR",
#'                                GeneFilter = 0.25, SampleFilter = 3,
#'                                sigma = 1.96, NCores = NULL, verbose = TRUE)
#' # define log fold change
#' p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 2, rate = 2)
#' # set up simulations
#' setupres <- Setup(ngenes = 10000, nsims = 10,
#'                   p.DE = 0.1, pLFC = p.lfc,
#'                   n1 = c(3,6,12), n2 = c(3,6,12),
#'                   Thinning = c(1,0.9,0.8), LibSize = 'given',
#'                   estParamRes = estparam_gene,
#'                   estSpikeRes = NULL,
#'                   DropGenes = FALSE,
#'                   sim.seed = 4379, verbose = TRUE)
#' # run simulation
#' simres <- simulateDE(SetupRes = setupres,
#'                      Prefilter = NULL, Imputation = NULL,
#'                      Normalisation = 'MR', Label = 'none',
#'                      DEmethod = "limma-trend", DEFilter = FALSE,
#'                      NCores = NULL, verbose = TRUE)
#' # DE evaluation
#' evalderes <- evaluateDE(simRes = simres, alpha.type="adjusted",
#'                         MTC='BH', alpha.nominal=0.05,
#'                         stratify.by = "mean", filter.by = "none")
#' plotEvalDE(evalderes, rate = "marginal", quick = FALSE, Annot = TRUE)
#' plotEvalDE(evalderes, rate = "conditional", quick = FALSE, Annot = TRUE)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ labs theme scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar theme_minimal
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup n
#' @importFrom tidyr %>%
#' @rdname plotEvalDE
#' @export
plotEvalDE <- function(evalRes, rate=c('marginal', 'conditional'), quick=TRUE, Annot=TRUE) {

  rate = match.arg(rate)

  pal = structure(c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"),
                  names = c("TPR", "FNR", "TNR", "FPR", "FDR"),
                  levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))

  # marginal rates over sample sizes
  if(rate=='marginal') {
    if(quick) {
      dat.marginal <- evalRes[c('TPR.marginal', 'FDR.marginal')]
      annot.text <- c("A) Marginal TPR and FDR per sample size setup (mean +/- standard error). \nB) Marginal TPR and FDR per sample size setup with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error).")
    }

    if(!quick) {
      dat.marginal <- evalRes[grep('*R.marginal', names(evalRes))]
      annot.text <- c("A) Marginal error rates per sample size setup (mean +/- standard error). \nB) Marginal error rates per sample size comparison with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error).")
    }
    names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
    dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalRes[['n1']], " vs ", evalRes[['n2']]))
    dat.marginal.long <- reshape2::melt(dat.marginal)
    dat.marginal.long$L1 <- factor(dat.marginal.long$L1,
                                   levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))
    refval <- data.frame(L1 = c("FDR", "TPR"), ref = c(evalRes$alpha.nominal, 0.8))

      # marginal in one
      grandplot <- ggplot2::ggplot(data = dat.marginal.long,
                                   ggplot2::aes_(x = quote(Var1),
                                                y = quote(value),
                                                color = quote(L1))) +
        stat_summary(fun.data = "mean_se",
                     size = 0.5,
                     position = position_dodge(width = 0.5)) +
        ggplot2::labs(x = "Sample Size Setup", y = "Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                    limits = c(0,1)) +
        ggplot2::scale_color_manual(values = pal) +
        .theme_eval_de()

      # faceted marginal
      facetplot <-  ggplot2::ggplot() +
        ggplot2::stat_summary(data = dat.marginal.long,
                              ggplot2::aes_(x = quote(Var1),
                                           y = quote(value),
                                           color = quote(L1)),
                              fun.data = "mean_se", size = 0.5) +
        ggplot2::stat_summary(data = dat.marginal.long,
                              ggplot2::aes_(x = quote(Var1),
                                           y = quote(value),
                                           group = quote(L1),
                                           color = quote(L1)),
                              fun.y = mean, geom="line") +
        ggplot2::geom_hline(data = refval,
                            ggplot2::aes_(yintercept = quote(ref)),
                            linetype="dashed",
                            color='black') +
        ggplot2::labs(x="Sample Size Setup", y="Rate") +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                    limits = c(0,1)) +
        ggplot2::scale_color_manual(values = pal) +
        ggplot2::facet_wrap(~L1, scales = 'fixed', ncol = 2) +
        .theme_eval_de() +
        ggplot2::theme(legend.position = "none")


      p.final <- suppressWarnings(cowplot::plot_grid(grandplot,
                                                     facetplot,
                                                     labels=c('A', 'B'),
                                                     label_size = 12,
                                                     rel_heights = c(1,1.75),
                                                     ncol = 1, nrow = 2))
  }

  #stratified rates
  if(rate=='conditional') {
    # stratum name
    stratum.name <- dplyr::case_when(evalRes$stratify.by == "mean" ~ "(Log Mean Expression)",
                                     evalRes$stratify.by == "dispersion" ~ "(Log Dispersion)",
                                     evalRes$stratify.by == "lfc" ~ "(Log Fold Change)",
                                     evalRes$stratify.by == "dropout" ~ "(Gene Dropout Rate)")
    # strata genes
    strata <- evalRes$strata.levels
    N <- length(evalRes$n1)
    dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],
                      'DEgenes'=evalRes$stratadiffgenes[,N,])
    dat.genes <- lapply(dat.genes, "rownames<-", strata)
    dat.genes.long <- reshape2::melt(dat.genes)
    dat.genes.calc <- dat.genes.long %>%
      dplyr::group_by(.data$Var1, .data$L1) %>%
      dplyr::summarise(Expectation=mean(.data$value),
                       Deviation=sd(.data$value),
                       Error=sd(.data$value)/sqrt(dplyr::n())) %>%
      dplyr::ungroup()

    refval <- data.frame(L1 = c("FDR", "TPR"),
                         ref = c(evalRes$alpha.nominal, 0.8))
    if(quick){
      dat.stratified <- evalRes[c('TPR', 'FDR')]
      annot.text <- c("A) Conditional TPR and FDR per sample size comparison per stratum (mean +/- standard error). \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.")
    }
    if(!quick){
      dat.stratified <- evalRes[grep('*R$', names(evalRes))]
      annot.text <- c("A) Conditional error rates per sample size comparison per stratum (mean +/- standard error). \nB) Number of equally (EE) and differentially expressed (DE) genes per stratum.")
    }
    dat.stratified <- lapply(dat.stratified, "dimnames<-",
                             list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]),
                                  NULL))
    dat.stratified.long <- reshape2::melt(dat.stratified)
    dat.stratified.long$L1 <- factor(dat.stratified.long$L1,
                                   levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))

    # rates
    facetplot <-  ggplot2::ggplot() +
      ggplot2::stat_summary(data = dat.stratified.long,
                            ggplot2::aes_(x = quote(Var1),
                                         y = quote(value),
                                         color = quote(Var2)),
                            fun.data = "mean_se", size = 0.5) +
      ggplot2::stat_summary(data = dat.stratified.long,
                            ggplot2::aes_(x = quote(Var1),
                                         y = quote(value),
                                         group = quote(Var2),
                                         color = quote(Var2)),
                            fun.y = mean, geom="line") +
      ggplot2::geom_hline(data = refval,
                          ggplot2::aes_(yintercept = quote(ref)),
                          linetype="dashed",
                          color='black') +
      ggplot2::labs(x=paste0("Stratum ", stratum.name),
                    y="Rate") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                  limits = c(0,1)) +
      ggplot2::facet_wrap(~ L1, scales = 'fixed', nrow = 3) +
      .theme_eval_de()

    # genes / stratum
    strataplot <- ggplot2::ggplot(data = dat.genes.calc,
                                  ggplot2::aes_(x=quote(Var1),
                                               y=quote(Expectation),
                                               fill=quote(L1))) +
      ggplot2::geom_bar(stat="identity")  +
      ggplot2::labs(x=paste0("Stratum ", stratum.name),
                    y="Count") +
      ggplot2::scale_fill_manual(values=c('grey', 'black'),
                                 breaks = c("DEgenes", "Ngenes"),
                                 labels = c("DE genes", "EE genes")) +
      .theme_eval_de()

    p.final <- suppressWarnings(cowplot::plot_grid(facetplot,
                                                   strataplot,
                                                   labels=c('A', 'B'),
                                                   rel_heights = c(2,1),
                                                   ncol=1, nrow=2))
  }

  # annotation under plot
  if(Annot) {
    p.final <- cowplot::add_sub(plot = p.final,
                                label = annot.text,
                                x = 0, hjust = 0, size=8)
  }

  # draw final plot
  cowplot::ggdraw(p.final)
}

# plotEvalSim -------------------------------------------------------------

#' @name plotEvalSim
#' @aliases plotEvalSim
#' @title Visualize power assessment
#' @description This function plots the results of \code{\link{evaluateSim}} for assessing the setup performance, i.e. normalisation method performance.
#' @usage plotEvalSim(evalRes, Annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateSim}}.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # estimate gene parameters
#' data("SmartSeq2_Gene_Read_Counts")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SmartSeq2_Gene_Read_Counts))
#' data("GeneLengths_mm10")
#' estparam_gene <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
#'                                readData = NULL,
#'                                batchData = Batches,
#'                                spikeData = NULL, spikeInfo = NULL,
#'                                Lengths = GeneLengths_mm10, MeanFragLengths = NULL,
#'                                RNAseq = 'singlecell', Protocol = 'Read',
#'                                Distribution = 'ZINB', Normalisation = "scran",
#'                                GeneFilter = 0.1, SampleFilter = 3,
#'                                sigma = 1.96, NCores = NULL, verbose = TRUE)
#' # define log fold change
#' p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#' # set up simulations
#' setupres <- Setup(ngenes = 10000, nsims = 10,
#'                   p.DE = 0.1, pLFC = p.lfc,
#'                   n1 = c(20,50,100), n2 = c(30,60,120),
#'                   Thinning = c(1,0.9,0.8), LibSize = 'given',
#'                   estParamRes = estparam_gene,
#'                   estSpikeRes = NULL,
#'                   DropGenes = FALSE,
#'                   sim.seed = 66437, verbose = TRUE)
#' # run simulation
#' simres <- simulateDE(SetupRes = setupres,
#'                      Prefilter = "FreqFilter",
#'                      Imputation = NULL,
#'                      Normalisation = 'scran', Label = 'none',
#'                      DEmethod = "limma-trend", DEFilter = FALSE,
#'                      NCores = NULL, verbose = TRUE)
#' # evaluation
#' evalsimres <- evaluateSim(simRes = simres)
#' plotEvalSim(evalRes = evalsimres, Annot = TRUE)
#' }
#' @author Beate Vieth
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ labs theme element_blank scale_y_continuous geom_line geom_hline geom_pointrange facet_wrap geom_boxplot position_dodge scale_fill_manual geom_bar theme_minimal
#' @importFrom grid unit
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr group_by summarise ungroup n
#' @importFrom tidyr "%>%"
#' @rdname plotEvalSim
#' @export
plotEvalSim <- function(evalRes, Annot=TRUE) {

    # log fold changes
    lfc <- reshape2::melt(evalRes$LogFoldChange)
    colnames(lfc) <- c("SimNo", "Metric", "Value", "Samples")
    lfc.dat <- lfc %>%
      tidyr::separate(.data$Metric, c("DE-Group", "Metric", "Type"), "_") %>%
      dplyr::filter(! .data$Type=="NAFraction") %>%
      tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2)) %>%
      dplyr::arrange(.data$SumN)
    # to label the x axis from smallest to largest n group!
    lfc.dat$Samples <- factor(lfc.dat$Samples,
                              levels=unique(lfc.dat$Samples[order(lfc.dat$SumN,decreasing = F)]))

    lfc.plot <- ggplot2::ggplot(data = lfc.dat,
                                ggplot2::aes_(x=quote(Samples),
                                             y=quote(Value),
                                             fill=quote(`DE-Group`),
                                             colour=quote(`DE-Group`))) +
      stat_summary(fun.data = "mean_se",
                   size = 0.5,
                   position = position_dodge(width = 0.5)) +
      ggplot2::labs(x="Sample Size Setup", y="Value") +
      ggplot2::facet_wrap(~Metric, ncol=1, scales="free") +
      ggplot2::expand_limits(y = 0) +
      .theme_eval_sim()

    # size factors
    sf <- reshape2::melt(evalRes$SizeFactors)
    colnames(sf) <- c("SimNo", "Metric", "Value", "Samples")
    sf.stats <- sf %>%
      dplyr::filter(.data$Metric %in% c("MAD", "rRMSE")) %>%
      tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2)) %>%
      dplyr::arrange(.data$SumN)
    # to label the x axis from smallest to largest n group!
    sf.stats$Samples <- factor(sf.stats$Samples,
                              levels=unique(sf.stats$Samples[order(sf.stats$SumN,decreasing = F)]))
    sfstats.plot <- ggplot2::ggplot(data = sf.stats,
                                    ggplot2::aes_(x=quote(Samples),
                                                 y=quote(Value))) +
      stat_summary(fun.data = "mean_se",
                   size = 0.5,
                   position = position_dodge(width = 0.5)) +
      ggplot2::facet_wrap(~Metric, ncol=1, scales="free") +
      ggplot2::labs(x="Sample Size Setup", y="Value") +
      ggplot2::expand_limits(y = 0) +
      .theme_eval_sim()

    # ratio of size factors per group
    ratio.dat <-  sf %>%
      dplyr::filter(grepl("Group", .data$Metric)) %>%
      tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
      dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2)) %>%
      dplyr::arrange(.data$SumN)
    ratio.dat$Samples <- factor(ratio.dat$Samples,
                               levels=unique(ratio.dat$Samples[order(ratio.dat$SumN,decreasing = F)]))

    ratio.plot <- ggplot2::ggplot(data = ratio.dat,
                                  ggplot2::aes_(x=quote(Samples),
                                                y=quote(Value),
                                                color=quote(Metric))) +
      stat_summary(fun.data = "mean_se",
                   size = 0.5,
                   position = position_dodge(width = 0.5)) +
      ggplot2::geom_hline(yintercept=1,linetype="dashed", color='darkgrey') +
      ggplot2::labs(x="Sample Size Setup", y="Value") +
      .theme_eval_sim()

    right_col <- suppressWarnings(cowplot::plot_grid(sfstats.plot,
                                                      ratio.plot,
                                                      labels = c('B', 'C'),
                                                      align = 'hv',
                                                      ncol=1,
                                                      nrow=2,
                                                     rel_heights = c(1.5, 1)))
    p.final <- suppressWarnings(cowplot::plot_grid(lfc.plot,
                                                   right_col,
                                                   labels=c('A', ''),
                                                   rel_widths = c(1.25, 1),
                                                   ncol=2, nrow=1))
    # annotation under plot
    if(Annot) {
      annot.text <- c("All values are mean +/- standard error. A) Mean Absolute Error (MAE), Root Mean Squared Error (RMSE) and robust Root Mean Squared Error (rRMSE) for the estimated log fold changes of all (ALL), differentially expressed (DE) and equally expressed (EE) genes compared to the true log fold changes. \nB) Median absolute deviation (MAD) and robust Root Mean Squared Error (rRMSE) between estimated and simulated size factors. \nC) The average ratio between simulated and estimated size factors in the two groups per sample size setup.")
      p.final <- cowplot::add_sub(plot = p.final,
                                  label = annot.text,
                                  x = 0, hjust = 0, size=8)
    }

  # draw final plot
  cowplot::ggdraw(p.final)

}


# plotTime ---------------------------------------------------------------

#' @name plotTime
#' @aliases plotTime
#' @title Visualize computational time
#' @description This function plots the computational running time of the simulations.
#' @usage plotTime(evalRes, Table=TRUE, Annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateSim}}.
#' @param Table A logical vector. If \code{TRUE}, a table of average running time in seconds per sample size setup and pipeline step is printed.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # estimate gene parameters
#' data("SmartSeq2_Gene_Read_Counts")
#' Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
#'                      stringsAsFactors = F,
#'                      row.names = colnames(SmartSeq2_Gene_Read_Counts))
#' data("GeneLengths_mm10")
#' estparam_gene <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
#'                                readData = NULL,
#'                                batchData = Batches,
#'                                spikeData = NULL, spikeInfo = NULL,
#'                                Lengths = GeneLengths_mm10, MeanFragLengths = NULL,
#'                                RNAseq = 'singlecell', Protocol = 'Read',
#'                                Distribution = 'ZINB', Normalisation = "scran",
#'                                GeneFilter = 0.1, SampleFilter = 3,
#'                                sigma = 1.96, NCores = NULL, verbose = TRUE)
#' # define log fold change
#' p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#' # set up simulations
#' setupres <- Setup(ngenes = 10000, nsims = 10,
#'                   p.DE = 0.1, pLFC = p.lfc,
#'                   n1 = c(20,50,100), n2 = c(30,60,120),
#'                   Thinning = c(1,0.9,0.8), LibSize = 'given',
#'                   estParamRes = estparam_gene,
#'                   estSpikeRes = NULL,
#'                   DropGenes = FALSE,
#'                   sim.seed = 66437, verbose = TRUE)
#' # run simulation
#' simres <- simulateDE(SetupRes = setupres,
#'                      Prefilter = "FreqFilter",
#'                      Imputation = NULL,
#'                      Normalisation = 'scran', Label = 'none',
#'                      DEmethod = "limma-trend", DEFilter = FALSE,
#'                      NCores = NULL, verbose = TRUE)
#'
#' # evaluation
#' evalsimres <- evaluateSim(simRes = simres)
#' plotEvalSim(evalRes = evalsimres, Annot = TRUE)
#' plotTime(evalRes = evalsimres, Annot = TRUE)
#' }
#' @author Beate Vieth
#' @importFrom tidyr "%>%" separate
#' @importFrom dplyr case_when bind_rows filter mutate arrange
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes_ position_dodge geom_point geom_pointrange facet_wrap theme_bw labs theme element_blank element_text guide_legend guides scale_fill_manual expand_scale alpha
#' @importFrom grid unit
#' @importFrom cowplot add_sub ggdraw
#' @importFrom matrixStats rowSds
#' @rdname plotTime
#' @export
plotTime <- function(evalRes, Table=TRUE, Annot=TRUE) {

  # proc time object summarised
  time.taken = evalRes$Timing
  # pipeline choices
  pipeline = evalRes$Pipeline
  prefilter = ifelse(is.null(pipeline$Prefilter), "no", pipeline$Prefilter)
  impute = ifelse(is.null(pipeline$Imputation), "no", pipeline$Imputation)
  normalisation = pipeline$Normalisation
  label = dplyr::case_when(pipeline$Label == "none" ~ "no labels",
                           pipeline$Label == "known" ~ "known group labels",
                           pipeline$Label == "clustering" ~ "labels derived by clustering")
  demethod = pipeline$DEmethod
  defilter = ifelse(pipeline$DEFilter == FALSE, "raw", "preprocessed")
  ncores = ifelse(is.null(pipeline$NCores), "1 core", paste0(pipeline$NCores, " cores"))
  # de settings
  desetup = evalRes$DESetup
  pde = desetup$p.DE
  ngenes = desetup$ngenes
  nsims = desetup$nsims

  # make plotting data object
  time.taken <- sapply(names(time.taken), function(i){
    time.taken[[i]] %>%
      tibble::rownames_to_column(var = "Step")
  }, USE.NAMES = TRUE, simplify = F)

  time.dat <- dplyr::bind_rows(time.taken, .id = "Samples") %>%
    tidyr::separate(col = .data$Step, into = c("Step", "Type")) %>%
    dplyr::filter(.data$Type == "Elapsed") %>%
    tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
    dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2)) %>%
    dplyr::arrange(.data$SumN) %>%
    dplyr::mutate(Lower = .data$Mean - .data$SEM,
                  Upper = .data$Mean + .data$SEM)

  time.dat$Samples <- factor(time.dat$Samples,
                                   levels=unique(time.dat$Samples[order(time.dat$SumN,decreasing = F)]))
  time.dat$Step <- factor(time.dat$Step,
                                levels=c("Simulation", "Preprocessing", "Normalisation", "DE", "Moments", "Total"))

  cols = .gg_color_hue(length(unique(time.dat$Samples)))

  p.final <- ggplot2::ggplot() +
    ggplot2::geom_bar(data = time.dat,
                      ggplot2::aes_(x = quote(Step), y = quote(Mean),
                                    color = quote(Samples),
                                    fill = quote(Samples)),
             stat="identity",
             position=position_dodge(width = 0.75), width=0.5) +
    ggplot2::geom_errorbar(data = time.dat,
                           ggplot2::aes_(ymin=quote(Lower),
                                         ymax=quote(Upper),
                                         x = quote(Step),
                                         group = quote(Samples),
                                         color = quote(Samples)),
                  position = position_dodge(width = 0.75), width = 0.5) +
    ggplot2::scale_fill_manual(values = ggplot2::alpha(cols, 0.5)) +
    ggplot2::scale_color_manual(values = ggplot2::alpha(cols, 1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(0, 0.1))) +
    ggplot2::labs(x = NULL, y = "Elapsed Time (seconds)") +
    ggplot2::coord_flip() +
    .theme_eval_time()

  if(isTRUE(Table)) {
    printtime <- time.dat[,c(1, 4, 6)]
    printtime[,c(3)] <- round(printtime[,c(3)], digits = 2)
    print(printtime %>%
            tidyr::pivot_wider(names_from = "Samples", values_from = "Mean") %>%
            data.frame(check.names = FALSE))
  }

  if(isTRUE(Annot)) {
    annot.text <- paste0("Time in seconds (mean +/- standard error) per pipeline step and sample size setup.\nExpression of ", ngenes, " genes was simulated ", nsims, " times. \nThe pipeline choices: ", prefilter, " filtering; ", impute, " imputation method; ", normalisation, " as normalisation method; ",  demethod, " for DE-testing on ", defilter, " counts. ", ncores, " utilized. ")
    p.final <- cowplot::add_sub(plot = p.final,
                                label = annot.text,
                                x = 0, hjust = 0, size=8)
  }

  # draw final plot
  cowplot::ggdraw(p.final)

}

# plotEvalROC -------------------------------------------------------------

#' @name plotEvalROC
#' @aliases plotEvalROC
#' @title Visualize error rate curves and associated summary statistics
#' @description This function plots the results of \code{\link{evaluateROC}} for assessing relative operating characteristic curves and summary statistics.
#' @usage plotEvalROC(evalRes,
#' cutoff=c('liberal', 'conservative'),
#' Annot=TRUE)
#' @param evalRes The output of \code{\link{evaluateROC}}.
#' @param cutoff Character vector defining whether the \code{"liberal"} or \code{"conservative"} FDR control is considered.
#' @param Annot A logical vector. If \code{TRUE}, a short figure legend under the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # estimate gene parameters
#' data("CELseq2_Gene_UMI_Counts")
#' estparam_gene <- estimateParam(countData = CELseq2_Gene_UMI_Counts,
#'                                readData = NULL,
#'                                batchData = NULL,
#'                                spikeData = NULL, spikeInfo = NULL,
#'                                Lengths = NULL, MeanFragLengths = NULL,
#'                                RNAseq = 'singlecell', Protocol = 'UMI',
#'                                Distribution = 'NB', Normalisation = "scran",
#'                                GeneFilter = 0.1, SampleFilter = 3,
#'                                sigma = 1.96, NCores = NULL, verbose = TRUE)
#' # define log2 fold change
#' p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#' # set up simulations
#' setupres <- Setup(ngenes = 10000, nsims = 10,
#'                   p.DE = 0.1, pLFC = p.lfc,
#'                   n1 = c(20,50,100), n2 = c(30,60,120),
#'                   Thinning = NULL, LibSize = 'equal',
#'                   estParamRes = estparam_gene,
#'                   estSpikeRes = NULL,
#'                   DropGenes = FALSE,
#'                   sim.seed = 34269, verbose = TRUE)
#' # run simulation
#' simres <- simulateDE(SetupRes = setupres,
#'                      Prefilter = "FreqFilter", Imputation = NULL,
#'                      Normalisation = 'scran', Label = 'none',
#'                      DEmethod = "limma-trend", DEFilter = FALSE,
#'                      NCores = NULL, verbose = TRUE)
#' # evaluation
#' evalrocres <- evaluateROC(simRes = simres,
#'                           alpha.type = "adjusted",
#'                           MTC = 'BH', alpha.nominal = 0.05,
#'                           raw = FALSE)
#' # plot evaluation
#' plotEvalROC(evalRes = evalrocres, cutoff = "conservative", Annot = TRUE)
#' }
#' @author Beate Vieth
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @rdname plotEvalROC
#' @export
plotEvalROC <- function(evalRes, cutoff=c('liberal', 'conservative'), Annot=TRUE) {

  cutoff = match.arg(cutoff)
  cutoff = ifelse(cutoff == "liberal", "lib", "conv")

  # ROC curve
  roc.plot <- .roc_plot(ROCData = evalRes$Performances$`ROC-Curve`)

  # PR curve
  pr.plot <- .pr_plot(ROCData = evalRes$Performances$`PR-Curve`)

  # TPR vs FDR curve
  alpha.nominal <- as.numeric(evalRes$Settings["alpha.nominal"])
  tprvsfdr.plot <-  .tprvsfdr_plot(ROCData = evalRes$TPRvsFDR,
                                   alpha.nominal = alpha.nominal)

  # table with summary statistics
  summary.tbl <- .summary_table_print(TblData = evalRes$Scores,
                                      cutoff = cutoff)


  curve.legend <- cowplot::get_legend(tprvsfdr.plot)
  # combine the plots
  top_row <- suppressWarnings(cowplot::plot_grid(roc.plot + ggplot2::theme(legend.position = "none"),
                                                 pr.plot + ggplot2::theme(legend.position = "none"),
                                                 labels= LETTERS[1:2],
                                                 ncol=2, nrow=1))
  bottom_row <- suppressWarnings(cowplot::plot_grid(tprvsfdr.plot + ggplot2::theme(legend.position = "none"),
                                                    summary.tbl,
                                                    labels= LETTERS[3:4],
                                                    ncol=2, nrow=1,
                                                    rel_widths = c(0.4, 0.6)))

  p.combined <- suppressWarnings(cowplot::plot_grid(top_row,
                                                 bottom_row,
                                                 ncol=1, nrow=2))

  p.final <- suppressWarnings(cowplot::plot_grid(p.combined,
                                                 curve.legend, rel_heights = c(1,0.1),
                                                 ncol=1, nrow=2))

  # annotation under plot
  if(Annot) {
    delta.text <- ifelse(evalRes$Settings$delta != 0,
                         paste0("considering genes with at least ", evalRes$Settings$delta, " as biologically meaningful DE genes.\n"),
                         " considering all DE genes.\n")
    settings.text <- paste0(evalRes$Settings$alpha.type, " p-values with nominal level equal to ", evalRes$Settings$alpha.nominal, delta.text)
    annot.text <- c("A) Receiver-Operator-Characteristics (ROC) Curve per sample size setup. \nB) Precision-Recall (PR) Curve per sample size setup. \nC) TPR versus observed FDR per sample size setup. The filling of the point indicates whether FDR is controlled at the chosen nominal level. \nD) Summary Statistics per sample size setup rounded to two digits.")
    p.final <- cowplot::add_sub(plot = p.final,
                                label = paste0(settings.text, annot.text),
                                x = 0, hjust = 0, size=8)
  }

  # draw final plot
  cowplot::ggdraw(p.final)
}

# plotEvalDist ------------------------------------------------------------

#' @name plotEvalDist
#' @aliases plotEvalDist
#' @title Visualize distribution assessment
#' @description This function plots the results of \code{\link{evaluateDist}} to assess goodness-of-fit testing.
#' @usage plotEvalDist(evalDistRes, Annot=TRUE)
#' @param evalDistRes The output of \code{\link{evaluateDist}}.
#' @param Annot A logical vector. If \code{TRUE}, a short description of the plot is included.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ## using example data set, but run it for fraction of genes
#' data("CELseq2_Gene_UMI_Counts")
#' evalDistRes <- evaluateDist(countData = CELseq2_Gene_UMI_Counts, batchData = NULL,
#'                             spikeData = NULL, spikeInfo = NULL,
#'                             Lengths = NULL, MeanFragLengths = NULL,
#'                             RNAseq = "singlecell", Protocol = "UMI",
#'                             Normalisation = "scran",
#'                             GeneFilter = 0.1, SampleFilter = 3,
#'                             FracGenes = 0.1,
#'                             verbose = TRUE)
#' plotEvalDist(evalDistRes)
#' }
#' @author Beate Vieth, Ines Hellmann
#' @importFrom ggplot2 ggplot aes_ theme scale_y_continuous geom_boxplot position_dodge geom_bar theme_minimal scale_x_discrete labs coord_flip
#' @importFrom scales percent
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom dplyr filter group_by summarise mutate mutate_if bind_rows slice select ungroup  left_join count
#' @importFrom tidyr %>% separate pivot_longer pivot_wider
#' @importFrom utils stack
#' @rdname plotEvalDist
#' @export
plotEvalDist <- function(evalDistRes, Annot=TRUE){

  dist_labels <- c("None"='None',
                   "ZIP"='Zero-Inflated \n Poisson',
                   "Poisson"='Poisson',
                   "ZINB"="Zero-Inflated \n Negative Binomial",
                   "NB"="Negative \n Binomial")
  label_names<-c("PoiBeta" ="Beta-Poisson",
                 'zifpois' = "Zero-Inflated \n Poisson",
                 'pois' = "Poisson",
                 "zifnbinom" = "Zero-Inflated \n Negative Binomial",
                 "nbinom" = "Negative \n Binomial")
  set_names <- c("plowestaic" = "Percentage Lowest\nAIC",
                 "plowestgoodaic" = "Percentage Lowest\nAIC + \nGOF \np > 0.05")
  zero_names <- c("diffzero" = "Obs. vs. Pred.\n Zeros",
                  "diffzerogood" = "Obs. vs. Pred.\n Zeros, GOF p>0.05")
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
                 "0_1_1_1"='Multiple')
  combi.ind.df <- data.frame(Distribution = combi.ind,
                             nbinom_pois_zifnbinom_zifpois = names(combi.ind),
                             row.names = NULL,
                             stringsAsFactors=F)

  # extract the GOF results from object
  gofres <- evalDistRes$GOF
  # reshape the table
  gofres <- cbind(rn = rownames(gofres), utils::stack(gofres))
  gofres <- gofres %>%
    tidyr::separate(.data$ind, c("distribution", "framework", 'type'),
                    "_", remove = F)
  # extract the observed zero table
  obszero <- evalDistRes$ObservedZeros
  obszero <- cbind(rn = rownames(obszero), obszero)
  obszero <- obszero[which(obszero$rn %in% gofres$rn), ]

  # GOF p-value based on chisquare test
  gof.pval <- gofres %>%
    dplyr::filter(.data$type=='gofpval', .data$framework=='standard') %>%
    dplyr::mutate(TestRes=ifelse(.data$values>0.05, 1, 0)) %>%
    dplyr::select(.data$rn, .data$TestRes, .data$distribution) %>%
    tidyr::spread(.data$distribution, .data$TestRes) %>%
    tidyr::unite(col = "nbinom_pois_zifnbinom_zifpois",
                 .data$nbinom, .data$pois, .data$zifnbinom, .data$zifpois, remove = F) %>%
    dplyr::full_join(combi.ind.df, by='nbinom_pois_zifnbinom_zifpois') %>%
    dplyr::mutate(Distribution=ifelse(is.na(.data$Distribution), 'Multiple', .data$Distribution)) %>%
    dplyr::group_by(.data$Distribution) %>%
    dplyr::summarise(Total=length(.data$nbinom_pois_zifnbinom_zifpois)) %>%
    dplyr::mutate(Percentage=.data$Total/sum(.data$Total)) %>%
    dplyr::filter(.data$Distribution != 'Multiple')
  p.chisquare <- ggplot2::ggplot(gof.pval,
                                 ggplot2::aes_(x = quote(Distribution),
                                      y=quote(Percentage))) +
    ggplot2::geom_bar(stat='identity', width=0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent,
                                limits = c(0,1)) +
    ggplot2::scale_x_discrete(labels= dist_labels,
                              limits=names(dist_labels)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = 'Distribution', y = "Percentage",
                  title = 'Goodness-of-fit statistic') +
    .theme_eval_dist()


  # AIC (lowest value AND intersect with GOF p-value)
  AIC.calc <- gofres %>%
    dplyr::filter(.data$framework %in%c('Marioni', 'standard'),
                  .data$type %in% c("gofpval","aic")) %>%
    dplyr::select(-.data$ind) %>%
    tidyr::spread(.data$type, .data$values) %>%
    dplyr::group_by(.data$rn) %>%
    dplyr::mutate(minAIC=(.data$aic==min(.data$aic)),
                  GOF = .data$gofpval>0.05) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$distribution) %>%
    dplyr::summarise(nlowestaic = sum(.data$minAIC,na.rm = T),
                     nlowestgoodaic = sum((.data$minAIC & .data$GOF), na.rm=T)) %>%
    dplyr::mutate(TotalLowest=sum(.data$nlowestaic),
                  TotalSub=sum(.data$nlowestgoodaic)) %>%
    dplyr::mutate(plowestaic=.data$nlowestaic/.data$TotalLowest,
                  plowestgoodaic=.data$nlowestgoodaic/.data$TotalSub) %>%
    dplyr::select(-.data$TotalLowest, -.data$TotalSub) %>%
    tidyr::pivot_longer(cols = .data$nlowestaic:.data$plowestgoodaic,
                        names_to = 'variable', values_to = 'value') %>%
    dplyr::mutate(Type=ifelse(grepl(pattern='p', .data$variable), 'Percentage', 'Total')) %>%
    dplyr::filter(.data$Type=='Percentage')
  AIC.calc$distribution <- factor(AIC.calc$distribution,
                                  levels = c('PoiBeta', 'zifpois', 'pois', 'zifnbinom', 'nbinom'))

  p.aic <- ggplot2::ggplot(data=AIC.calc,
                           ggplot2::aes_(x=quote(variable),
                                y=quote(value),
                                fill=quote(distribution))) +
    ggplot2::geom_bar(stat='identity', width=0.7, colour="black") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::scale_x_discrete(labels = set_names) +
    ggplot2::scale_fill_brewer(labels= label_names,
                               limits=names(label_names), palette="Set1") +
    ggplot2::labs(x = NULL, y = 'Percentage',
                  fill = NULL, title = "Akaike Information Criterion") +
    ggplot2::guides(fill=guide_legend(reverse=TRUE)) +
    ggplot2::coord_flip() +
    .theme_eval_dist()

  # predicted zeroes vs observed zeros
  predzero.dat <- gofres %>%
    dplyr::filter(.data$type %in% c("gofpval","predzero"),
                  .data$framework %in% c("Marioni", 'standard')) %>%
    dplyr::select(-.data$ind) %>%
    tidyr::pivot_wider(names_from = .data$type, values_from = .data$values) %>%
    dplyr::left_join(obszero, by='rn') %>%
    dplyr::mutate(diffzero=.data$ObsZero-.data$predzero,
                  diffzerogood=ifelse(.data$gofpval>0.05,.data$ObsZero-.data$predzero,NA)) %>%
    tidyr::pivot_longer(c(.data$diffzero, .data$diffzerogood),
                        names_to = "comp", values_to = "divergentzero")
  zpdat <- dplyr::filter(predzero.dat,
                         .data$comp == "diffzero")
  p.zero <- ggplot2::ggplot(zpdat,
                            ggplot2::aes_(x = quote(distribution),
                                          y = quote(divergentzero))) +
    geom_boxplot(width = 0.7, position = position_dodge(width = 0.8),
                 outlier.shape = NA) +
    ggplot2::labs(x = "Distribution",
                  y = "Observed - Predicted Zeros",
                  title = "Dropouts") +
    ggplot2::scale_x_discrete(labels = label_names,
                              limits = names(label_names)) +
    ggplot2::coord_flip() +
    .theme_eval_dist()

  # Best fit by LRT / Vuong
  model.pval <- gofres %>%
    dplyr::filter(.data$type=='gofpval', .data$framework=='standard') %>%
    dplyr::mutate(TestRes=ifelse(.data$values>0.05, 1, 0)) %>%
    dplyr::select(.data$rn, .data$TestRes, .data$distribution) %>%
    tidyr::pivot_wider(names_from = .data$distribution,
                       values_from = .data$TestRes) %>%
    na.omit() %>%
    tidyr::unite("nbinom_pois_zifnbinom_zifpois",
                 .data$nbinom, .data$pois,
                 .data$zifnbinom, .data$zifpois, remove = F) %>%
    dplyr::mutate_if(is.factor, as.character)

  model.calc <- gofres %>%
    dplyr::filter(.data$distribution %in% c('LRT', 'Vuong'),
                  .data$values<0.05) %>%
    dplyr::select(.data$rn, .data$type) %>%
    dplyr::mutate_if(is.factor, as.character)
  out <- split(model.calc, f = model.calc$type)
  out2 <- lapply(out, function(x) {
    dplyr::left_join(x, model.pval, by='rn') %>%
      na.omit()
  })
  out2[['NBPoisson']] <- out2[['NBPoisson']] %>%
    dplyr::filter(.data$nbinom==1, .data$pois==1)  %>%
    dplyr::select(.data$rn, .data$type)
  out2[['ZNB']] <- out2[['ZNB']] %>%
    dplyr::filter(.data$nbinom==1, .data$zifnbinom==1)  %>%
    dplyr::select(.data$rn, .data$type)
  out2[['ZNBZPoisson']] <- out2[['ZNBZPoisson']] %>%
    dplyr::filter(.data$zifnbinom==1, .data$zifpois==1)  %>%
    dplyr::select(.data$rn, .data$type)
  out2[['ZPoisson']] <- out2[['ZPoisson']] %>%
    dplyr::filter(.data$pois==1, .data$zifpois==1)  %>%
    dplyr::select(.data$rn, .data$type)
  out3 <- do.call('rbind', out2)
  ModelTest <- out3 %>%
    dplyr::group_by(.data$type) %>%
    dplyr::count() %>%
    mutate(Percentage=.data$n/sum(.data$n))
  p.modeltest <- ggplot2::ggplot(data = ModelTest,
                                 ggplot2::aes_(x = quote(type),
                                               y = quote(Percentage))) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(x = "Test", y = "Percentage", title="Model Comparisons") +
    ggplot2::scale_x_discrete(labels = c(NBPoisson = "Negative Binomial \n> Poisson",
                                         ZNB = "Zero-inflated Negative Binomial \n> Negative Binomial",
                                         ZPoisson = "Zero-inflated Poisson \n> Poisson",
                                         ZNBZPoisson = "Zero-inflated Negative Binomial \n> Zero-inflated Poisson"),
                              limits = c("ZNBZPoisson", "ZPoisson", "NBPoisson", "ZNB")) +
    ggplot2::coord_flip() +
    .theme_eval_dist()

  p.final <- cowplot::ggdraw() +
    cowplot::draw_plot(p.chisquare, x= 0 , y= 0.5,  width = 0.5, height = 0.5) +
    cowplot::draw_plot(p.aic, x=0.5, y=0.5, width = 0.5, height=0.5) +
    cowplot::draw_plot(p.zero, x=0, y= 0,  width = 0.5, height=0.5) +
    cowplot::draw_plot(p.modeltest, x=0.5, y= 0,  width = 0.5, height=0.5) +
    cowplot::draw_plot_label( c("A","B","C", "D"),
                              x=c(0, 0.5, 0, 0.5),
                              y=c(1,1,0.5, 0.5), size=15, vjust=1)

  if (Annot) {
    p.final <- cowplot::add_sub(p.final, "A) Goodness-of-fit of the model assessed with a chi-square test based on residual deviance and degrees of freedom.
                                \nB) Akaike Information Criterion per gene: Model with the lowest AIC. Model with the lowest AIC and passed goodness-of-fit statistic test.
                                \nC) Observed versus predicted dropouts per model and gene plotted without outliers.
                                \nD) Model Assessment based on LRT for nested models and Vuong test for nonnested models.",
                                size = 8)
  }
  cowplot::ggdraw(p.final)
}

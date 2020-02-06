## ----knitr_init, echo=FALSE, results="asis", cache=FALSE----------------------
library(knitr)
library(rmdformats)
## Global options
options(max.print = "75")
opts_chunk$set(echo = FALSE,
	             cache = FALSE,
               prompt = FALSE,
               tidy = FALSE,
               comment = NA,
               message = FALSE,
               warning = FALSE)
opts_knit$set(width = 75)

## ----DevTools, echo = TRUE, eval = FALSE--------------------------------------
#  install.packages('devtools')
#  library(devtools)

## ----dep, echo = TRUE, eval = FALSE-------------------------------------------
#  ipak <- function(pkg, repository=c('CRAN', 'Bioconductor', 'github')){
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    # new.pkg <- pkg
#    if (length(new.pkg)) {
#      if(repository=='CRAN') {
#        install.packages(new.pkg, dependencies = TRUE)
#      }
#      if(repository=='Bioconductor') {
#        if(strsplit(version[['version.string']], ' ')[[1]][3] > "3.5.0"){
#          if (!requireNamespace("BiocManager")){
#              install.packages("BiocManager")
#          }
#              BiocManager::install(new.pkg, dependencies=TRUE, ask=FALSE)
#        }
#        if(strsplit(version[['version.string']], ' ')[[1]][3] < "3.5.0"){
#                source("https://bioconductor.org/biocLite.R")
#        biocLite(new.pkg, dependencies=TRUE, ask=FALSE)
#        }
#      }
#      if(repository=='github') {
#        devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE, dependencies=TRUE)
#      }
#    }
#  }
#  
#  # CRAN PACKAGES
#  cranpackages <- c("broom", "cobs", "cowplot",
#                    "data.table", "doParallel", "dplyr", "DrImpute",
#                    "fastICA", "fitdistrplus", "foreach", "future",
#                    "gamlss.dist", "ggplot2", "ggpubr", "grDevices",
#                    "grid", "Hmisc", "kernlab", "MASS", "magrittr", "MBESS", "Matrix",
#                    "matrixStats", "mclust", "methods", "minpack.lm", "moments", "msir",
#                    "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl",
#                    "reshape2", "Rmagic", "rsvd", "Rtsne", "scales", "Seurat", "snow", "sctransform",
#                    "stats", "tibble", "tidyr", "truncnorm", "VGAM", "ZIM", "zoo")
#  ipak(cranpackages, repository='CRAN')
#  
#  # BIOCONDUCTOR
#  biocpackages <- c("bayNorm", "baySeq", "BiocGenerics", "BiocParallel",
#                    "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA",
#                    "limma", "Linnorm", "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq",
#                    "S4Vectors", "scater", "scDD", "scde", "scone", "scran", "SCnorm",
#                    "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
#  ipak(biocpackages, repository='Bioconductor')
#  
#  # GITHUB
#  githubpackages <- c('nghiavtr/BPSC', 'cz-ye/DECENT',
#                      'mohuangx/SAVER', 'statOmics/zingeR')
#  ipak(githubpackages, repository = 'github')

## ----depcheck, echo = TRUE, eval = FALSE--------------------------------------
#  powsimRdeps <- data.frame(Package = c(cranpackages,
#                                        biocpackages,
#                                        sapply(strsplit(githubpackages, "/"), "[[", 2)),
#                            stringsAsFactors = F)
#  
#  ip <- as.data.frame(installed.packages()[,c(1,3:4)], stringsAsFactors = F)
#  
#  ip.check <- cbind(powsimRdeps,
#                    Version = ip[match(powsimRdeps$Package, rownames(ip)),"Version"])
#  
#  table(is.na(ip.check$Version))  # all should be FALSE

## ----install1, echo=T, eval=F, tidy=T-----------------------------------------
#  devtools::install_github('bvieth/powsimR',
#                           build_vignettes = TRUE,
#                           dependencies=FALSE)
#  library("powsimR")

## ----install2, echo=T, eval=F, tidy=T-----------------------------------------
#  devtools::install_github("bvieth/powsimR")

## ----schematic, fig.cap="PowsimR schematic overview. (A) Estimation (B) Simulation (C) Evaluation.", echo=F, eval=T, include=T----
knitr::include_graphics("powsimRvignetteschematic.png")

## ----geneparams, echo=T, eval=F, include=T------------------------------------
#  data("CELseq2_Gene_UMI_Counts")
#  Batches <- data.frame(Batch = sapply(strsplit(colnames(CELseq2_Gene_UMI_Counts), "_"), "[[", 1),
#                        stringsAsFactors = FALSE, row.names = colnames(CELseq2_Gene_UMI_Counts))
#  data("GeneLengths_mm10")
#  
#  # estimation
#  estparam_gene <- estimateParam(countData = CELseq2_Gene_UMI_Counts,
#                            readData = NULL,
#                            batchData = Batches,
#                            spikeData = NULL,
#                            spikeInfo = NULL,
#                            Lengths = GeneLengths, MeanFragLengths = NULL,
#                            RNAseq = 'singlecell', Protocol = 'UMI',
#                            Distribution = 'NB', Normalisation = "scran",
#                            GeneFilter = 0.1, SampleFilter = 3,
#                            sigma = 1.96, NCores = NULL, verbose = TRUE)
#  
#  # plotting
#  plotParam(estParamRes = estparam_gene, Annot = T)

## ----geneparamsplot, echo = FALSE, eval = TRUE, include = TRUE, fig.wide = TRUE, fig.cap="Estimated parameters for Ziegenhain data set. A) Quality Control Metrics: Sequencing depth; Library size factors with median (black line) for the filtered data set; Detected genes; Ratio of gene to spike-in counts (if spike-ins were provided). Outliers are marked in red. B) Marginal Distribution of gene mean, dispersion and dropout rate per estimation set. C) Number of genes and samples per estimation set. Provided by the user; Detected = number of genes and samples with at least one count; All = number of genes for which mean, dispersion and dropout could be estimated using non-outlying samples. \nFiltered = number of genes above filter threshold for which mean, dispersion and dropout could be estimated using non-outlying samples. Dropout Genes = number of genes filtered out due to dropout rate. D) Local polynomial regression fit between mean and dispersion estimates with variability band per gene (yellow). Common dispersion estimate (grey dashed line). E) Fraction of dropouts versus estimated mean expression per gene."----
knitr::include_graphics("estparam_gene_celseq2.png")

## ----spikeparams, echo=T, eval=F, include=T, warning=F------------------------
#  data("SmartSeq2_SpikeIns_Read_Counts")
#  data("SmartSeq2_SpikeInfo")
#  Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_SpikeIns_Read_Counts), "_"), "[[", 1),
#                         stringsAsFactors = F,
#                         row.names = colnames(SmartSeq2_SpikeIns_Read_Counts))
#  # estimation
#  estparam_spike <- estimateSpike(spikeData = SmartSeq2_SpikeIns_Read_Counts,
#  spikeInfo = SmartSeq2_SpikeInfo,
#  MeanFragLength = NULL,
#  batchData = Batches,
#  Normalisation = 'depth')
#  
#  # plotting
#  plotSpike(estparam_spike)
#  

## ----spikeplot, echo=F, warning=F, eval=T, include=T, fig.height = 7, fig.width=10, fig.align='centered', fig.cap="Estimated parameters for the spike-ins added to Smartseq2 libraries in Ziegenhain dataset. (A) Sequencing depth per sample with median sequencing depth (grey dashed line). (B) Library size normalisation factor per sample with median size factor (grey dashed line). (C) Calibration curve with mean expression estimates and average R squared over all cells. (D) Capture efficiency with binomial logistic regression fit over all cells."----
knitr::include_graphics("estparam_spike_celseq2.png")

## ----lfcs, echo=F, eval=T, include=T, fig.cap="Examples of Log Fold Changes following a gamma, uniform and normal distribution."----
knitr::include_graphics("lfcdist.jpeg")

## ----simsetup, echo = TRUE, eval = FALSE--------------------------------------
#  # define log fold change
#  p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#  # set up simulations
#  setupres <- Setup(ngenes = 10000, nsims = 25,
#                    p.DE = 0.05, pLFC = p.lfc,
#                    n1 = c(48, 96, 384, 800), n2 = c(48, 96, 384, 800),
#                    Thinning = NULL, LibSize = 'equal',
#                    estParamRes = estparam_gene,
#                    estSpikeRes = estparam_spike,
#                    DropGenes = TRUE,
#                    sim.seed = 5299, verbose = TRUE)
#  

## ----simrun, eval=F, echo=T---------------------------------------------------
#  simres <- simulateDE(SetupRes = setupres,
#                       Prefilter = NULL, Imputation = NULL,
#                       Normalisation = 'scran',
#                       DEmethod = "limma-trend", DEFilter = FALSE,
#                       NCores = NULL, verbose = TRUE)

## ----evalderes, echo = T, eval=F----------------------------------------------
#  evalderes = evaluateDE(simRes = simres,
#                       alpha.type = 'adjusted',
#                       MTC = 'BH',
#                       alpha.nominal = 0.1,
#                       stratify.by = 'mean',
#                       filter.by = 'none',
#                       strata.filtered = 1,
#                       target.by = 'lfc',
#                       delta = 0)

## ----evaldeplot1, echo=F, eval=T, fig.cap="Marginal Error Rates. (A) Marginal TPR and FDR per sample size comparison. (B) Marginal TPR and FDR per sample size comparison with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error)."----
knitr::include_graphics("evalderes_marginal_celseq2.png")

## ----evaldeplot2, echo=F, eval=T, fig.cap="Stratified Error Rates. (A) Conditional TPR and FDR per sample size comparison per stratum. (B) Number of equally (EE) and differentially expressed (DE) genes per stratum."----
knitr::include_graphics("evalderes_conditional_celseq2.png")

## ----evalrocres, echo = T, eval=F---------------------------------------------
#  evalrocres = evaluateROC(simRes = simres,
#                           alpha.type="adjusted",
#                           MTC='BH',
#                           alpha.nominal = 0.1)
#  
#  plotEvalROC(evalRes = evalrocres, cutoff = "liberal")

## ----evalrocresplot, echo=F, eval=T, fig.cap="A) Receiver-Operator-Characteristics (ROC) Curve per sample size setup. \nB) Precision-Recall (PR) Curve per sample size setup. \nC) TPR versus observed FDR per sample size setup. The filling of the point indicates whether FDR is controlled at the chosen nominal level. \nD) Summary Statistics per sample size setup rounded to two digits."----
knitr::include_graphics("evalrocres_liberal_celseq2.png")

## ----evalsimres, echo = T, eval=F---------------------------------------------
#  evalsimres = evaluateSim(simRes = simres)
#  
#  plotEvalSim(evalRes = evalsimres)
#  
#  plotTime(evalRes = evalsimres)
#  

## ----evalsimresplot, echo=F, eval=T, fig.cap="Pipeline Evaluation. A) Mean Absolute Error (MAE), Root Mean Squared Error (RMSE) and robust Root Mean Squared Error (rRMSE) for the estimated log fold changes of all (ALL), differentially expressed (DE) and equally expressed (EE) genes compared to the true log fold changes. \nB) Median absolute deviation (MAD) and robust Root Mean Squared Error (rRMSE) between estimated and simulated size factors. \nC) The average ratio between simulated and estimated size factors in the two groups per sample size setup. All values are mean +/- standard error."----
knitr::include_graphics("evalsimres_celseq2.png")

## ----evaltimeplot, echo=F, eval=T, fig.cap="Computational Run Time in seconds per simulateDE() pipeline step."----
knitr::include_graphics("evaltime_celseq2.png")

## ----online_repos, echo = T, eval = F-----------------------------------------
#  # Install and load the R package
#  BiocManager::install("recount")
#  library('recount')
#  
#  # Download the data set
#  url <- download_study('SRP060416')
#  
#  # Load the data
#  load(file.path('SRP060416', 'rse_gene.Rdata'))
#  
#  # count table
#  cnts <- assay(rse_gene)
#  # sample annotation
#  sample.info <- data.frame(colData(rse_gene)@listData,
#                            stringsAsFactors=F)
#  # gene annotation
#  gene.info <- data.frame(GeneID=rowData(rse_gene)@listData$gene_id,
#                          GeneLength=rowData(rse_gene)@listData$bp_length,
#                          stringsAsFactors=F)

## ----evaldist, eval=F, echo=T-------------------------------------------------
#  data("SmartSeq2_Gene_Read_Counts")
#  evalDistRes <- evaluateDist(countData = SmartSeq2_Gene_Read_Counts,
#                              batchData = NULL,
#                              spikeData = NULL, spikeInfo = NULL,
#                              Lengths = NULL, MeanFragLengths = NULL,
#                              RNAseq = "singlecell", Protocol = "UMI",
#                              Normalisation = "scran",
#                              GeneFilter = 0.1, SampleFilter = 3,
#                              FracGenes = 0.1,
#                              verbose = TRUE)
#  plotEvalDist(evalDistRes)

## ----evaldistplot, echo=F, eval=T, include=T, fig.wide = TRUE, fig.cap="Distribution Evaluation. A) Goodness-of-fit of the model assessed with a Chi-Square Test based on residual deviance and degrees of freedom. B) Akaike Information Criterion per gene: Model with the lowest AIC. Model with the lowest AIC and passed goodness-of-fit statistic test.  C) Observed versus predicted dropouts per model and gene plotted without outliers. D) Model Assessment based on LRT for nested models and Vuong test for nonnested models. "----
knitr::include_graphics("evaldist_smartseq2.png")

## ----thinning, echo = TRUE, eval = FALSE--------------------------------------
#  data("CELseq2_Gene_UMI_Counts")
#  data("CELseq2_Gene_Read_Counts")
#  Batches <- data.frame(Batch = sapply(strsplit(colnames(CELseq2_Gene_UMI_Counts), "_"), "[[", 1),
#                        stringsAsFactors = FALSE, row.names = colnames(CELseq2_Gene_UMI_Counts))
#  data("GeneLengths_mm10")
#  
#  # estimation
#  estparam_gene <- estimateParam(countData = CELseq2_Gene_UMI_Counts,
#                                 readData = CELseq2_Gene_Read_Counts,
#                                 batchData = Batches,
#                                 spikeData = NULL,
#                                 spikeInfo = NULL,
#                                 Lengths = GeneLengths, MeanFragLengths = NULL,
#                                 RNAseq = 'singlecell', Protocol = 'UMI',
#                                 Distribution = 'NB', Normalisation = "scran",
#                                 GeneFilter = 0.1, SampleFilter = 3,
#                                 sigma = 1.96, NCores = NULL, verbose = TRUE)
#  
#  plotParam(estParamRes = estparam_gene)
#  
#  # define log fold change
#  p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
#  
#  # set up simulations
#  setupres <- Setup(ngenes = 10000, nsims = 25,
#                    p.DE = 0.1, pLFC = p.lfc,
#                    n1 = c(500, 1000), n2 = c(500, 1000),
#                    Thinning = c(1, 0.5), LibSize = 'given',
#                    estParamRes = estparam_gene,
#                    estSpikeRes = NULL,
#                    DropGenes = TRUE,
#                    sim.seed = 5299, verbose = TRUE)
#  
#  # run simulations
#  simres <- simulateDE(SetupRes = setupres,
#                       Prefilter = "FreqFilter", Imputation = NULL,
#                       Normalisation = 'scran',
#                       DEmethod = "limma-trend", DEFilter = FALSE,
#                       NCores = NULL, verbose = TRUE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()


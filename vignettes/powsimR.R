## ----knitset, eval=TRUE, include=FALSE, cache=FALSE------------------------
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=50), 
                      fig.align = 'center', 
                      message=FALSE, error=FALSE, warning=FALSE)

## ----dep, echo = TRUE, eval = FALSE, tidy = FALSE--------------------------
#  ipak <- function(pkg, repository=c('CRAN', 'Bioconductor', 'github')){
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    if (length(new.pkg)) {
#      if(repository=='CRAN') {
#        install.packages(new.pkg, dependencies = TRUE)
#      }
#      if(repository=='Bioconductor') {
#        source("https://bioconductor.org/biocLite.R")
#        biocLite(new.pkg, dependencies=TRUE, ask=FALSE)
#      }
#      if(repository=='github') {
#        devtools::install_github(new.pkg, build_vignettes = FALSE, dependencies=TRUE)
#      }
#    }
#  }
#  
#  # CRAN PACKAGES
#  cranpackages <- c("bbmle", "broom", "cluster", "cobs", "cowplot",
#                    "data.table", "doParallel", "dplyr", "drc", "DrImpute", "fastICA", "fitdistrplus",
#                    "foreach", "gamlss.dist", "ggExtra", "ggplot2", "ggthemes", "grDevices",
#                    "glmnet", "grid", "gtools", "Hmisc", "kernlab", "MASS",
#                    "matrixStats", "mclust", "methods", "minpack.lm", "moments", "msir",
#                    "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl",
#                    "reshape2", "ROCR", "Rtsne", "scales", "Seurat", "snow",
#                    "stats", "tibble", "tidyr", "VGAM", "ZIM")
#  ipak(cranpackages, repository='CRAN')
#  
#  # BIOCONDUCTOR
#  biocpackages <- c("AnnotationDbi", "baySeq", "Biobase", "BiocGenerics",
#                    "BiocParallel", "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "limma",
#                    "Linnorm", "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq",
#                    "S4Vectors", "scater", "scDD", "scde", "scone", "scran", "SCnorm",
#                    "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
#  ipak(biocpackages, repository='Bioconductor')
#  
#  # GITHUB
#  githubpackages <- c('nghiavtr/BPSC', 'VCCRI/cidr', 'cz-ye/DECENT',
#                      'mohuangx/SAVER', 'statOmics/zingeR')
#  ipak(githubpackages, repository = 'github')

## ----install, echo=TRUE, eval=FALSE, tidy=FALSE----------------------------
#  devtools::install_github('bvieth/powsimR', build_vignettes = TRUE,
#                           dependencies = FALSE)

## ----load, echo=TRUE, eval=T, tidy=T---------------------------------------
library(powsimR)

## ----schematic, fig.cap="PowsimR schematic overview. (A) Estimation (B) Simulation (C) Evaluation.", echo=F, eval=T, include=T, fig.wide = T----
knitr::include_graphics("powsimR-vignette-schematic.png")

## ----params, echo=T, eval=T, include=T-------------------------------------
data("kolodziejczk_cnts")
kolodziejczk_cnts <- kolodziejczk_cnts[, grep('standard',
                                              colnames(kolodziejczk_cnts))]
TwoiLIF.params <- estimateParam(countData=kolodziejczk_cnts,
                                batchData = NULL,
                                spikeData = NULL,
                                spikeInfo = NULL,
                                Lengths = NULL,
                                MeanFragLengths = NULL,
                                Distribution = 'ZINB',
                                RNAseq = 'singlecell',
                                normalisation = 'scran',
                                sigma = 1.96,
                                NCores = NULL)

## ----paramsplot, echo=T, eval=T, include=T, fig.height = 7, fig.width=10, fig.cap="Estimated parameters for Kolodziejczyk data set. (A) Sequencing depth per sample with median sequencing depth (grey dashed line). (B) Library size normalisation factor per sample with median size factor (grey dashed line). (C) Marginal Distribution of log2(mean), log2(dispersion) and dropout. (D) Local polynomial regression fit between log2(mean) and log2(dispersion) estimates with variability band per gene (yellow). Common dispersion estimate (grey dashed line). E) Fraction of dropouts versus estimated mean expression per gene."----
plotParam(TwoiLIF.params, annot = F)

## ----simeval, fig.cap="Comparison of estimated and simulated read counts. (A) Dispersion versus Mean. (B) Dropout versus Mean.", echo=F, eval=T, include=T, out.width = "95%"----
knitr::include_graphics("simeval.jpeg")

## ----spikeparams, echo=T, eval=T, include=T--------------------------------
data("scrbseq_spike_cnts")
data("scrbseq_spike_info")
batch_info <- data.frame(Batch = ifelse(grepl(pattern = "SCRBseqA_", colnames(scrbseq_spike_cnts)), "A", "B"), row.names = colnames(scrbseq_spike_cnts))
## spike information table
spike_info <- scrbseq_spike_info[-1,]
## estimation
spike.param <- estimateSpike(spikeData = scrbseq_spike_cnts,
                             spikeInfo = spike_info,
                             MeanFragLength = NULL,
                             batchData = batch_info,
                             normalisation = 'depth')

## ----spikeplot, fig.cap="Estimated parameters for the spike-ins added to SCRBseq libraries in Ziegenhain dataset. (A) Sequencing depth per sample with median sequencing depth (grey dashed line). (B) Library size normalisation factor per sample with median size factor (grey dashed line). (C) Calibration curve with mean expression estimates and average R squared over all cells. (D) Capture efficiency with binomial logistic regression fit over all cells.", echo=T, eval=T, include=T, fig.height = 7, fig.width=10, fig.align='centered'----
plotSpike(spike.param, annot = F)

## ----lfcs, fig.cap="Log2 fold change examples for gamma, uniform and normal distribution.", echo=F, eval=T, include=T, out.width = "95%"----
knitr::include_graphics("lfcdist.jpeg")

## ----simsetup, echo = TRUE, eval = FALSE-----------------------------------
#  lfc.gamma = function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, 3, 3)
#  de.opts = DESetup(ngenes=10000, nsims=25,
#                    p.DE=0.2, pLFC=lfc.gamma,
#                    sim.seed = 58673)
#  sim.opts = SimSetup(desetup = de.opts,
#                      params = TwoiLIF.params,
#                      spike=NULL,
#                      size.factors='equal',
#                      downsample=FALSE, geneset = FALSE)

## ----simrun, eval=F, echo=T------------------------------------------------
#  simDE = simulateDE(n1 = c(24,48,96,192,384,800),
#                     n2 = c(24,48,96,192,384,800),
#                     sim.settings = sim.opts,
#                     DEmethod = "limma-trend",
#                     normalisation = "scran",
#                     Preclust = FALSE,
#                     Prefilter = NULL,
#                     Impute = NULL,
#                     spikeIns = FALSE,
#                     NCores = NULL,
#                     verbose = TRUE)

## ---- echo = T, eval=T-----------------------------------------------------
data("kolodziejczk_simDE")
simDE = kolodziejczk_simDE
evalDE = evaluateDE(simRes = simDE,
                     alpha.type = 'adjusted',
                     MTC = 'BH',
                     alpha.nominal = 0.1,
                     stratify.by = 'mean',
                     filter.by = 'none',
                     strata.filtered = 1,
                     target.by = 'lfc',
                     delta = 0)

## ----evalplot1, echo=T, eval=T, fig.cap="Marginal Error Rates. (A) Marginal TPR and FDR per sample size comparison. (B) Marginal TPR and FDR per sample size comparison with dashed line indicating nominal alpha level (type I error) and nominal 1-beta level, i.e. 80% power (type II error)."----
plotEvalDE(evalRes = evalDE,
            rate='marginal',
            quick=TRUE, annot=FALSE)

## ----evalplot2, echo=T, eval=T, fig.cap="Stratified Error Rates. (A) Conditional TPR and FDR per sample size comparison per stratum. (B) Number of equally (EE) and differentially expressed (DE) genes per stratum."----
plotEvalDE(evalRes = evalDE,
            rate='stratified',
            quick=TRUE, annot=FALSE)

## ----twogroup, echo = TRUE, eval = TRUE------------------------------------
plfc.foo = function(x) sample(c(-1,1), size=x, prob = c(0.25,0.75),replace=T)*
  rgamma(x, 2, 4)
blfc.foo = function(x) rnorm(x, sd = 0.25)
simcounts.2grp <- simulateCounts(n=c(120, 100),
                                 ngenes=10000,
                                 p.DE=0.1, pLFC=plfc.foo,
                                 p.B=0.1, bLFC=blfc.foo, 
                                 bPattern="uncorrelated",
                                 p.M=NULL, mLFC=NULL,
                                 params=kolodziejczk_param,
                                 spike=NULL,
                                 spikeIns=FALSE,
                                 size.factors="given",
                                 downsample=F,
                                 geneset=F,
                                 sim.seed=NULL,
                                 verbose=TRUE)

## ----plot2grp, echo=T, eval=T----------------------------------------------
plotCounts(simCounts = simcounts.2grp, Distance = "euclidean", Scale = T, DimReduce = "PCA", verbose = T)

## ----multigroup, echo = TRUE, eval = TRUE----------------------------------
if(length(grep("MBESS",installed.packages()))==0){
   install.packages("MBESS", dependencies = TRUE)
}
if(length(grep("mvtnorm",installed.packages()))==0){
   install.packages("mvtnorm", dependencies = TRUE)
 }
cor.lfc <- matrix(c(1,0.5,0.7,0.5,1,0.95,0.7,0.95,1), nrow=3, ncol=3)
v.lfc <- c(4,1,1)
cov.lfc <- MBESS::cor2cov(cor.lfc,v.lfc)

plfc.foo = function(x) {
  mu.tmp = stats::rnorm(n = ncol(cov.lfc), mean = 0, sd = 0.5)
  mvtnorm::rmvnorm(x, mean = mu.tmp, sigma = cov.lfc)
}

simcounts.3grp <- simulateCounts(n=c(100, 50, 50),
                                 ngenes=10000,
                                 p.DE=0.1, pLFC=plfc.foo,
                                 p.B=NULL, bLFC=NULL, 
                                 bPattern="uncorrelated",
                                 p.M=NULL, mLFC=NULL,
                                 params=kolodziejczk_param,
                                 spike=NULL,
                                 spikeIns=FALSE,
                                 size.factors="given",
                                 downsample=F,
                                 geneset=F,
                                 sim.seed=NULL,
                                 verbose=TRUE)

## ----plot3grp, echo=T, eval=T----------------------------------------------
plotCounts(simCounts = simcounts.3grp, Distance = "euclidean", Scale = F, DimReduce = "PCA", verbose = T)

## ----online_repos, echo = T, eval = F--------------------------------------
#  # Install and load the R package
#  source('http://bioconductor.org/biocLite.R')
#  biocLite('recount')
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
#  sample.info <- data.frame(colData(rse_gene)@listData, stringsAsFactors=F)
#  # gene annotation
#  gene.info <- data.frame(GeneID=rowData(rse_gene)@listData$gene_id, GeneLength=rowData(rse_gene)@listData$bp_length, stringsAsFactors=F)

## ----evaldist, eval=F, echo=T----------------------------------------------
#  library("powsimRDev")
#  data("kolodziejczk_cnts")
#  kolodziejczk_cnts <- kolodziejczk_cnts[, grep('standard',
#                                                colnames(kolodziejczk_cnts))]
#  TwoiLIF.dist <- evaluateDist(countData = kolodziejczk_cnts,
#                               batchData = NULL,
#                               spikeData = NULL,
#                               spikeInfo = NULL,
#                               Lengths = NULL,
#                               MeanFragLengths = NULL,
#                               RNAseq = 'singlecell',
#                               normalisation = 'scran',
#                               frac.genes = 0.2,
#                               min.meancount = 0.1,
#                               max.dropout = 0.8,
#                               min.libsize = 1000,
#                               verbose = TRUE)

## ----evaldistplot, echo=F, eval=T, include=T, fig.wide = TRUE, fig.cap="Distribution Evaluation. A) Goodness-of-fit of the model assessed with a Chi-Square Test based on residual deviance and degrees of freedom. B) Akaike Information Criterion per gene: Model with the lowest AIC. Model with the lowest AIC and passed goodness-of-fit statistic test.  C) Observed versus predicted dropouts per model and gene plotted without outliers. D) Model Assessment based on LRT for nested models and Vuong test for nonnested models. "----
knitr::include_graphics("evaldist.png")

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()


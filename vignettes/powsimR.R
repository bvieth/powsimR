## ----knitset, eval=TRUE, include=FALSE, cache=FALSE---------------------------
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=50), 
                      fig.align = 'center', 
                      message=FALSE, error=FALSE, warning=FALSE)

## ----dep, echo = TRUE, eval = FALSE, tidy = FALSE-----------------------------
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

## ----install, echo=TRUE, eval=FALSE, tidy=FALSE-------------------------------
#  devtools::install_github('bvieth/powsimR',
#                           build_vignettes = TRUE,
#                           dependencies = FALSE)
#  library("powsimR")


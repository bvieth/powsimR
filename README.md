
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

NEWS
----

Version 1.1.0 released with the following changes / additions (2018-03-29):

-   simulation of batch effects (see options `p.B`, `bLFC` and `bPattern` in `DESetup` and `simulateCounts`)
-   simulation of spike-in expression (see `estimateSpike` , `plotSpike` and option `spikeIns` in `simulateDE` and `simulateCounts`)
-   simulation of multiple sample groups (e.g. single cell populations) with `simulateCounts`
-   imputation and prefiltering options prior to normalisation in DE power simulations added (scImpute, scone, Seurat, DrImpute, SAVER)
-   additional normalisation options and DE tools (esp. for single cells) included in `simulateDE`
-   evaluation of simulation setup using estimated versus true value comparisons of library size factors and log2 fold changes in `evaluateSim` and `plotEvalSim`
-   increased flexibility in preprocessing for distribution evaluation in `evaluateDist`

Installation Guide
------------------

You can install powsimR from GitHub using devtools:

``` r
library(devtools)
install_github("bvieth/powsimR")
```

If the default installation does not work, you can install powsimR by first installing its dependencies:

``` r
ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) {
        if (repository == "CRAN") {
            install.packages(new.pkg, dependencies = TRUE)
        }
        if (repository == "Bioconductor") {
            source("https://bioconductor.org/biocLite.R")
            biocLite(new.pkg, dependencies = TRUE, ask = FALSE)
        }
        if (repository == "github") {
            devtools::install_github(new.pkg, build_vignettes = FALSE, dependencies = TRUE)
        }
    }
}

# CRAN PACKAGES
cranpackages <- c("bbmle", "broom", "cluster", "cobs", "cowplot", "data.table", 
    "doParallel", "dplyr", "drc", "DrImpute", "fastICA", "fitdistrplus", "foreach", 
    "gamlss.dist", "ggExtra", "ggplot2", "ggthemes", "grDevices", "glmnet", 
    "grid", "gtools", "Hmisc", "kernlab", "MASS", "matrixStats", "mclust", "methods", 
    "minpack.lm", "moments", "msir", "NBPSeq", "nonnest2", "parallel", "penalized", 
    "plyr", "pscl", "reshape2", "ROCR", "Rtsne", "scales", "Seurat", "snow", 
    "stats", "tibble", "tidyr", "VGAM", "ZIM")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("AnnotationDbi", "BASiCS", "baySeq", "Biobase", "BiocGenerics", 
    "BiocParallel", "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "limma", "Linnorm", 
    "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", 
    "scDD", "scde", "scone", "scran", "SingleCellExperiment", "SummarizedExperiment", 
    "zinbwave")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("nghiavtr/BPSC", "VCCRI/cidr", "cz-ye/DECENT", "mohuangx/SAVER", 
    "rhondabacher/SCnorm", "statOmics/zingeR")
ipak(githubpackages, repository = "github")
```

After installing the dependencies, powsimR can be installed by using devtools as well.

``` r
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
library("powsimR")
```

Some users have experienced issues installing powsimR due to Tex compilation errors. If that is the case, you can leave out building the vignette.

User Guide
----------

For examples and tips on using the package, please see the vignette [here](https://github.com/bvieth/powsimR/tree/master/vignettes/powsimR.Rmd) or open it in R by typing

``` r
browseVignettes("powsimR")
```

Citation
--------

Please use the following entry for citing powsimR.

``` r
citation("powsimR")
```

powsimR is published in [Bioinformatics](https://doi.org/10.1101/117150). A preprint paper is also on [bioRxiv](https://doi.org/10.1101/117150).

Notes
-----

Please send bug reports and feature requests by opening a new issue on [this page](https://github.com/bvieth/powsimR/issues).

Note that the error "maximal number of DLLs reached..." might occur due to the loading of many shared objects by Bioconductor packages. Restarting the R session after installing dependencies / powsimR will help.

Starting with R version 3.4.0, one can set the environmental variable 'R\_MAX\_NUM\_DLLS' to a higher number. See `?Startup()` for more information.

R Session Info
--------------

``` r
sessionInfo()
#> R version 3.4.4 (2018-03-15)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 16.04.4 LTS
#> 
#> Matrix products: default
#> BLAS: /usr/lib/openblas-base/libblas.so.3
#> LAPACK: /usr/lib/libopenblasp-r0.2.18.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_3.4.4  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
#>  [5] formatR_1.5     tools_3.4.4     htmltools_0.3.6 yaml_2.1.18    
#>  [9] Rcpp_0.12.16    stringi_1.1.7   rmarkdown_1.9   knitr_1.20     
#> [13] stringr_1.3.0   digest_0.6.15   evaluate_0.10.1
```

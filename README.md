
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

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
cranpackages <- c("methods", "stats", "matrixStats", "Rtsne", "moments", "minpack.lm", 
    "glmnet", "cluster", "mclust", "MASS", "gtools", "doParallel", "parallel", 
    "snow", "reshape2", "plyr", "dplyr", "tidyr", "tibble", "data.table", "ggplot2", 
    "grid", "ggthemes", "ggExtra", "cowplot", "scales", "cobs", "msir", "drc", 
    "DrImpute", "VGAM", "NBPSeq")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("S4Vectors", "DEDS", "AnnotationDbi", "Biobase", "BiocGenerics", 
    "SummarizedExperiment", "BiocParallel", "RUVSeq", "scran", "scater", "SingleCellExperiment", 
    "Linnorm", "edgeR", "limma", "DESeq2", "baySeq", "NOISeq", "EBSeq", "MAST", 
    "scde", "scDD", "ROTS", "monocle", "IHW", "qvalue")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("nghiavtr/BPSC", "rhondabacher/SCnorm", "catavallejos/BASiCS")
ipak(githubpackages, repository = "github")
```

After installing the dependencies, powsimR can be installed by using devtools as well.

``` r
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
```

Some users have experienced issues installing powsimR due to Tex compilation errors. If that is the case, you can leave out building the vignette.

User Guide
----------

For examples and tips on using the package, please see the vignette PDF [here](https://github.com/bvieth/powsimR/tree/master/vignettes/powsimR.pdf) or open it in R by typing

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

Starting with R version 3.4.0, one can set the environmental variable 'R\_MAX\_NUM\_DLLS' to a higher number. See `?Startup()` or the [vignette](https://github.com/bvieth/powsimR/tree/master/vignettes/powsimR.pdf) for more information.

R Session Info
--------------

``` r
sessionInfo()
#> R version 3.4.2 (2017-09-28)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 14.04.5 LTS
#> 
#> Matrix products: default
#> BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
#> LAPACK: /usr/lib/lapack/liblapack.so.3.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_GB.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_3.4.2  backports_1.1.1 magrittr_1.5    rprojroot_1.2  
#>  [5] formatR_1.5     tools_3.4.2     htmltools_0.3.6 yaml_2.1.14    
#>  [9] Rcpp_0.12.13    stringi_1.1.5   rmarkdown_1.6   knitr_1.17     
#> [13] stringr_1.2.0   digest_0.6.12   evaluate_0.10.1
```

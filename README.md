
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `powsimR` <br/> Power analysis for bulk and <br/> single cell RNA-seq experiments <img src="vignettes/powsimR.png" align="right" width="200" />

Please also consult my Github Page of
[powsimR](https://bvieth.github.io/powsimR/) made with
[pkgdown](http://pkgdown.r-lib.org/index.html)!

## :arrow_double_down: Installation Guide

For the installation, `devtools` is needed.

``` r
install.packages("devtools")
library(devtools)
```

I have compiled a lite version of powsimR package on the lite branch
which only includes our recommended tools for normalisation, imputation
and DE-testing. You can install this version using devtools:

``` r
install_github("bvieth/powsimR", ref = "lite", dependencies = TRUE, upgrade = "ask")
```

I recommend to install first the dependencies manually and then powsimR.
If you plan to use MAGIC for imputation, then please follow their
[instruction](https://github.com/KrishnaswamyLab/MAGIC) to install the
python implementation before installing powsimR-master.

``` r
ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    # new.pkg <- pkg
    if (length(new.pkg)) {
        if (repository == "CRAN") {
            install.packages(new.pkg, dependencies = TRUE)
        }
        if (repository == "Bioconductor") {
            if (strsplit(version[["version.string"]], " ")[[1]][3] > "4.0.0") {
                if (!requireNamespace("BiocManager")) {
                  install.packages("BiocManager")
                }
                BiocManager::install(new.pkg, dependencies = TRUE, ask = FALSE)
            }
            if (strsplit(version[["version.string"]], " ")[[1]][3] < "3.6.0") {
                stop(message("powsimR depends on packages and functions that are only available in R 4.0.0 and higher."))
            }
        }
        if (repository == "github") {
            remotes::install_github(new.pkg, build_vignettes = FALSE, force = FALSE,
                dependencies = TRUE)
        }
    }
}

# CRAN PACKAGES
cranpackages <- c("broom", "cobs", "cowplot", "data.table", "doParallel", "dplyr",
    "DrImpute", "fastICA", "fitdistrplus", "foreach", "future", "gamlss.dist", "ggplot2",
    "ggpubr", "ggstance", "grDevices", "grid", "Hmisc", "kernlab", "MASS", "magrittr",
    "MBESS", "Matrix", "matrixStats", "mclust", "methods", "minpack.lm", "moments",
    "msir", "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl", "reshape2",
    "Rmagic", "rsvd", "Rtsne", "scales", "Seurat", "snow", "sctransform", "stats",
    "tibble", "tidyr", "truncnorm", "VGAM", "ZIM", "zoo")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("bayNorm", "baySeq", "BiocGenerics", "BiocParallel", "DESeq2",
    "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", "Linnorm", "MAST", "monocle", "NOISeq",
    "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", "scDD", "scde", "scone", "scran",
    "SCnorm", "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("cz-ye/DECENT", "nghiavtr/BPSC", "mohuangx/SAVER", "statOmics/zingeR",
    "Vivianstats/scImpute")
ipak(githubpackages, repository = "github")
```

To check whether all dependencies are installed, you can run the
following lines:

``` r

powsimRdeps <- data.frame(Package = c(cranpackages, 
                                      biocpackages, 
                                      sapply(strsplit(githubpackages, "/"), "[[", 2)), 
                          stringsAsFactors = F)

ip <- as.data.frame(installed.packages()[,c(1,3:4)], stringsAsFactors = F)

ip.check <- cbind(powsimRdeps, 
                  Version = ip[match(powsimRdeps$Package, rownames(ip)),"Version"]) 

table(is.na(ip.check$Version))  # all should be FALSE
```

After installing the dependencies, powsimR can be installed by using
devtools as well.

``` r
remotes::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
library("powsimR")
```

Alternative, you can try to install powsimR and its dependencies
directly using devtools:

``` r
remotes::install_github("bvieth/powsimR")
```

## :book: User Guide

For examples and tips on using the package, please consult the vignette
after successful installation by

``` r
browseVignettes("powsimR")
```

Some users have experienced issues installing powsimR due to vignette
compilation errors or because they are missing the necessary R packages
to build the vignette, i.e. knitr and rmdformats. If that is the case,
you can either install these dependencies or leave out building the
vignette (by setting build_vignettes to FALSE) and read it on my Github
Page of
[powsimR](https://bvieth.github.io/powsimR/articles/powsimR.html) or
download it as a html file
[here](https://github.com/bvieth/powsimR/blob/master/vignettes/powsimR.html).

### DLLs and ulimit

Note that the error “maximal number of DLLs reached…” might occur due to
the loading of many shared objects by Bioconductor packages. Restarting
the R session after installing dependencies / powsimR will help.
Starting with R version 3.4.0, one can set the environmental variable
‘R_MAX_NUM_DLLS’ to a higher number. See `?Startup()` for more
information. I recommend to increase the maximum number of DLLs that can
be loaded to 500. The environmental variable R_MAX_NUM_DLLS can be set
in R_HOME/etc/Renviron prior to starting R. For that locate the Renviron
file and add the following line: R_MAX_NUM_DLLS=xy where xy is the
number of DLLs. On my Ubuntu machine, the Renviron file is in
/usr/lib/R/etc/ and I can set it to 500.

In addition, the user limits for open files (unix: ulimit) might have to
be set to a higher number to accomodate the increase in DLLs. Please
check out the help pages for
[MACs](https://gist.github.com/tombigel/d503800a282fcadbee14b537735d202c)
and
[Linux](https://glassonionblog.wordpress.com/2013/01/27/increase-ulimit-and-file-descriptors-limit/)
for guidance.

## :scroll: Citation

Please use the following entry for citing powsimR.

``` r
citation("powsimR")
```

powsimR is published in
[Bioinformatics](https://doi.org/10.1093/bioinformatics/btx435). A
preprint paper is also on [bioRxiv](https://doi.org/10.1101/117150).

## :incoming_envelope: Notes

Please send bug reports and feature requests by opening a new issue on
[this page](https://github.com/bvieth/powsimR/issues). I try to keep up
to date with new developments / changes of methods implemented in
powsimR, but if you encounter run errors while using a certain tool
(e.g. for imputation), then I appreciate if you can post this as an
[issue](https://github.com/bvieth/powsimR/issues).

## `R` Session Info

``` r
library(powsimR)
#> Loading required package: gamlss.dist
#> Loading required package: MASS
sessionInfo()
#> R version 4.3.1 (2023-06-16)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.6 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so;  LAPACK version 3.7.1
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] powsimR_0.0.900   gamlss.dist_6.0-5 MASS_7.3-60      
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.3.1               later_1.3.1                
#>   [3] bitops_1.0-7                tibble_3.2.1               
#>   [5] lpsymphony_1.28.1           minpack.lm_1.2-3           
#>   [7] rpart_4.1.19                lifecycle_1.0.3            
#>   [9] rstatix_0.7.2               msir_1.3.3                 
#>  [11] edgeR_3.42.4                doParallel_1.0.17          
#>  [13] lattice_0.21-8              backports_1.4.1            
#>  [15] magrittr_2.0.3              limma_3.56.2               
#>  [17] Hmisc_5.1-0                 rmarkdown_2.23             
#>  [19] yaml_2.3.7                  plotrix_3.8-2              
#>  [21] metapod_1.8.0               shinyBS_0.61.1             
#>  [23] httpuv_1.6.11               cowplot_1.1.1              
#>  [25] abind_1.4-5                 zlibbioc_1.46.0            
#>  [27] quadprog_1.5-8              GenomicRanges_1.52.0       
#>  [29] purrr_1.0.1                 BiocGenerics_0.46.0        
#>  [31] RCurl_1.98-1.12             nnet_7.3-19                
#>  [33] sandwich_3.0-2              GenomeInfoDbData_1.2.10    
#>  [35] IRanges_2.34.1              S4Vectors_0.38.1           
#>  [37] irlba_2.3.5.1               moments_0.14.1             
#>  [39] MatrixModels_0.5-2          dqrng_0.3.0                
#>  [41] fitdistrplus_1.1-11         DelayedMatrixStats_1.22.1  
#>  [43] codetools_0.2-19            DelayedArray_0.26.7        
#>  [45] scuttle_1.10.1              DT_0.28                    
#>  [47] tidyselect_1.2.0            ScaledMatrix_1.8.1         
#>  [49] matrixStats_1.0.0           stats4_4.3.1               
#>  [51] base64enc_0.1-3             BiocNeighbors_1.18.0       
#>  [53] ellipsis_0.3.2              Formula_1.2-5              
#>  [55] survival_3.5-5              iterators_1.0.14           
#>  [57] foreach_1.5.2               tools_4.3.1                
#>  [59] Rcpp_1.0.11                 glue_1.6.2                 
#>  [61] mnormt_2.1.1                gridExtra_2.3              
#>  [63] xfun_0.39                   DESeq2_1.40.2              
#>  [65] SAVER_1.1.2                 qvalue_2.32.0              
#>  [67] MatrixGenerics_1.12.2       GenomeInfoDb_1.36.1        
#>  [69] dplyr_1.1.2                 shinydashboard_0.7.2       
#>  [71] formatR_1.14                fastmap_1.1.1              
#>  [73] bluster_1.10.0              fansi_1.0.4                
#>  [75] SparseM_1.81                truncnorm_1.0-9            
#>  [77] digest_0.6.33               rsvd_1.0.5                 
#>  [79] R6_2.5.1                    mime_0.12                  
#>  [81] colorspace_2.1-0            UpSetR_1.4.0               
#>  [83] utf8_1.2.3                  tidyr_1.3.0                
#>  [85] generics_0.1.3              data.table_1.14.8          
#>  [87] htmlwidgets_1.6.2           S4Arrays_1.0.5             
#>  [89] pkgconfig_2.0.3             gtable_0.3.3               
#>  [91] SingleCellExperiment_1.22.0 XVector_0.40.0             
#>  [93] cobs_1.3-5                  htmltools_0.5.5            
#>  [95] lavaan_0.6-16               carData_3.0-5              
#>  [97] scales_1.2.1                Biobase_2.60.0             
#>  [99] scran_1.28.2                knitr_1.43                 
#> [101] rstudioapi_0.15.0           reshape2_1.4.4             
#> [103] checkmate_2.2.0             zoo_1.8-12                 
#> [105] stringr_1.5.0               parallel_4.3.1             
#> [107] nonnest2_0.5-5              foreign_0.8-84             
#> [109] pillar_1.9.0                grid_4.3.1                 
#> [111] vctrs_0.6.3                 slam_0.1-50                
#> [113] pscl_1.5.5.1                promises_1.2.0.1           
#> [115] ggpubr_0.6.0                car_3.1-2                  
#> [117] BiocSingular_1.16.0         IHW_1.28.0                 
#> [119] beachmat_2.16.0             xtable_1.8-4               
#> [121] cluster_2.1.4               htmlTable_2.4.1            
#> [123] evaluate_0.21               pbivnorm_0.6.0             
#> [125] mvtnorm_1.2-2               cli_3.6.1                  
#> [127] locfit_1.5-9.8              compiler_4.3.1             
#> [129] rlang_1.1.1                 crayon_1.5.2               
#> [131] ggsignif_0.6.4              fdrtool_1.2.17             
#> [133] mclust_6.0.0                plyr_1.8.8                 
#> [135] stringi_1.7.12              BiocParallel_1.34.2        
#> [137] munsell_0.5.0               CompQuadForm_1.4.3         
#> [139] quantreg_5.96               Matrix_1.6-0               
#> [141] sparseMatrixStats_1.12.2    ggplot2_3.4.2              
#> [143] statmod_1.5.0               shiny_1.7.4.1              
#> [145] SummarizedExperiment_1.30.2 ROCR_1.0-11                
#> [147] igraph_1.5.0.1              broom_1.0.5                
#> [149] ggstance_0.3.6              iCOBRA_1.28.1
```

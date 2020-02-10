
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powsimR: Power analysis for bulk and single cell RNA-seq experiments

Please also consult my Github Page of
[powsimR](https://bvieth.github.io/powsimR/) made with
[pkgdown](http://pkgdown.r-lib.org/index.html)\!

## Installation Guide

For the installation, the R package `devtools` is needed.

``` r
install.packages("devtools")
library(devtools)
```

I recommend to install first the dependencies manually and then powsimR.
If you plan to use MAGIC for imputation, then please follow their
[instruction](https://github.com/KrishnaswamyLab/MAGIC) to install the
python implementation before installing
powsimR.

``` r
ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    # new.pkg <- pkg
    if (length(new.pkg)) {
        if (repository == "CRAN") {
            install.packages(new.pkg, dependencies = TRUE)
        }
        if (repository == "Bioconductor") {
            if (strsplit(version[["version.string"]], " ")[[1]][3] > "3.6.0") {
                if (!requireNamespace("BiocManager")) {
                  install.packages("BiocManager")
                }
                BiocManager::install(new.pkg, dependencies = TRUE, ask = FALSE)
            }
            if (strsplit(version[["version.string"]], " ")[[1]][3] < "3.6.0") {
                stop(message("powsimR depends on packages that are only available in R 3.6.0 and higher."))
            }
        }
        if (repository == "github") {
            devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE, 
                dependencies = TRUE)
        }
    }
}

# CRAN PACKAGES
cranpackages <- c("broom", "cobs", "cowplot", "data.table", "doParallel", "dplyr", 
    "DrImpute", "fastICA", "fitdistrplus", "foreach", "future", "gamlss.dist", "ggplot2", 
    "ggpubr", "grDevices", "grid", "Hmisc", "kernlab", "MASS", "magrittr", "MBESS", 
    "Matrix", "matrixStats", "mclust", "methods", "minpack.lm", "moments", "msir", 
    "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl", "reshape2", "Rmagic", 
    "rsvd", "Rtsne", "scales", "Seurat", "snow", "sctransform", "stats", "tibble", 
    "tidyr", "truncnorm", "VGAM", "ZIM", "zoo")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("bayNorm", "baySeq", "BiocGenerics", "BiocParallel", "DEDS", "DESeq2", 
    "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", "Linnorm", "MAST", "monocle", "NOISeq", 
    "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", "scDD", "scde", "scone", "scran", 
    "SCnorm", "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("nghiavtr/BPSC", "cz-ye/DECENT", "mohuangx/SAVER", "statOmics/zingeR")
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
devtools as
well.

``` r
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
library("powsimR")
```

Alternative, you can try to install powsimR and its dependencies
directly using devtools:

``` r
devtools::install_github("bvieth/powsimR")
```

## User Guide

For examples and tips on using the package, please consult the vignette
after successful installation by

``` r
browseVignettes("powsimR")
```

Some users have experienced issues installing powsimR due to vignette
compilation errors. If that is the case, you can leave out building the
vignette (by setting build\_vignettes to FALSE) and read it on my Github
Page of
[powsimR](https://bvieth.github.io/powsimR/articles/powsimR.html) or
download it as a html file
[here](https://github.com/bvieth/powsimR/tree/master/inst/doc/powsimR.html).

### DLLs and ulimit

Note that the error “maximal number of DLLs reached…” might occur due to
the loading of many shared objects by Bioconductor packages. Restarting
the R session after installing dependencies / powsimR will help.
Starting with R version 3.4.0, one can set the environmental variable
‘R\_MAX\_NUM\_DLLS’ to a higher number. See `?Startup()` for more
information. I recommend to increase the maximum number of DLLs that can
be loaded to 500. The environmental variable R\_MAX\_NUM\_DLLS can be
set in R\_HOME/etc/Renviron prior to starting R. For that locate the
Renviron file and add the following line: R\_MAX\_NUM\_DLLS=xy where xy
is the number of DLLs. On my Ubuntu machine, the Renviron file is in
/usr/lib/R/etc/ and I can set it to 500.

In addition, the user limits for open files (unix: ulimit) might have to
be set to a higher number to accomodate the increase in DLLs. Please
check out the help pages for
[MACs](https://gist.github.com/tombigel/d503800a282fcadbee14b537735d202c)
and
[Linux](https://glassonionblog.wordpress.com/2013/01/27/increase-ulimit-and-file-descriptors-limit/)
for guidance.

## Citation

Please use the following entry for citing powsimR.

``` r
citation("powsimR")
```

powsimR is published in
[Bioinformatics](https://doi.org/10.1093/bioinformatics/btx435). A
preprint paper is also on [bioRxiv](https://doi.org/10.1101/117150).

## Notes

Please send bug reports and feature requests by opening a new issue on
[this page](https://github.com/bvieth/powsimR/issues).

## R Session Info

``` r
library(powsimR)
#> Loading required package: gamlss.dist
#> Loading required package: MASS
#> Warning: replacing previous import 'DECENT::lrTest' by 'MAST::lrTest' when
#> loading 'powsimR'
#> Warning: replacing previous import 'penalized::predict' by 'stats::predict' when
#> loading 'powsimR'
#> Warning: replacing previous import 'zinbwave::glmWeightedF' by
#> 'zingeR::glmWeightedF' when loading 'powsimR'
sessionInfo()
#> R version 3.6.2 (2019-12-12)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
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
#> other attached packages:
#> [1] powsimR_1.2.0     gamlss.dist_5.1-6 MASS_7.3-51.5    
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.2.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            lattice_0.20-38            
#>   [5] vctrs_0.2.2                 fastICA_1.2-2              
#>   [7] mgcv_1.8-31                 penalized_0.9-51           
#>   [9] blob_1.2.1                  survival_3.1-8             
#>  [11] Rmagic_2.0.3                later_1.0.0                
#>  [13] nloptr_1.2.1                DBI_1.1.0                  
#>  [15] R.utils_2.9.2               SingleCellExperiment_1.8.0 
#>  [17] rappdirs_0.3.1              Linnorm_2.10.0             
#>  [19] dqrng_0.2.1                 jpeg_0.1-8.1               
#>  [21] zlibbioc_1.32.0             MatrixModels_0.4-1         
#>  [23] htmlwidgets_1.5.1           mvtnorm_1.0-12             
#>  [25] future_1.16.0               UpSetR_1.4.0               
#>  [27] parallel_3.6.2              scater_1.14.6              
#>  [29] irlba_2.3.3                 DEoptimR_1.0-8             
#>  [31] Rcpp_1.0.3                  KernSmooth_2.23-16         
#>  [33] DT_0.12                     promises_1.1.0             
#>  [35] gdata_2.18.0                DDRTree_0.1.5              
#>  [37] DelayedArray_0.12.2         limma_3.42.2               
#>  [39] vegan_2.5-6                 Hmisc_4.3-1                
#>  [41] ShortRead_1.44.3            apcluster_1.4.8            
#>  [43] RSpectra_0.16-0             msir_1.3.2                 
#>  [45] mnormt_1.5-6                digest_0.6.23              
#>  [47] png_0.1-7                   qlcMatrix_0.9.7            
#>  [49] sctransform_0.2.1           cowplot_1.0.0              
#>  [51] pkgconfig_2.0.3             docopt_0.6.1               
#>  [53] DelayedMatrixStats_1.8.0    ggbeeswarm_0.6.0           
#>  [55] iterators_1.0.12            minqa_1.2.4                
#>  [57] lavaan_0.6-5                reticulate_1.14            
#>  [59] SummarizedExperiment_1.16.1 spam_2.5-1                 
#>  [61] beeswarm_0.2.3              modeltools_0.2-22          
#>  [63] xfun_0.12                   zoo_1.8-7                  
#>  [65] tidyselect_1.0.0            ZIM_1.1.0                  
#>  [67] reshape2_1.4.3              purrr_0.3.3                
#>  [69] kernlab_0.9-29              EDASeq_2.20.0              
#>  [71] viridisLite_0.3.0           snow_0.4-3                 
#>  [73] rtracklayer_1.46.0          rlang_0.4.4                
#>  [75] hexbin_1.28.1               glue_1.3.1                 
#>  [77] RColorBrewer_1.1-2          fpc_2.2-5                  
#>  [79] matrixStats_0.55.0          stringr_1.4.0              
#>  [81] fields_10.3                 ggsignif_0.6.0             
#>  [83] DESeq2_1.26.0               SparseM_1.78               
#>  [85] httpuv_1.5.2                class_7.3-15               
#>  [87] BPSC_0.99.2                 BiocNeighbors_1.4.1        
#>  [89] annotate_1.64.0             jsonlite_1.6.1             
#>  [91] XVector_0.26.0              bit_1.1-15.1               
#>  [93] mime_0.9                    gridExtra_2.3              
#>  [95] gplots_3.0.1.2              Rsamtools_2.2.1            
#>  [97] zingeR_0.1.0                stringi_1.4.5              
#>  [99] gmodels_2.18.1              bitops_1.0-6               
#> [101] maps_3.3.0                  RSQLite_2.2.0              
#> [103] tidyr_1.0.2                 pheatmap_1.0.12            
#> [105] data.table_1.12.8           DEDS_1.60.0                
#> [107] rstudioapi_0.11             GenomicAlignments_1.22.1   
#> [109] nlme_3.1-144                qvalue_2.18.0              
#> [111] scran_1.14.6                fastcluster_1.1.25         
#> [113] locfit_1.5-9.1              scone_1.10.0               
#> [115] listenv_0.8.0               cobs_1.3-4                 
#> [117] R.oo_1.23.0                 prabclus_2.3-2             
#> [119] dbplyr_1.4.2                segmented_1.1-0            
#> [121] BiocGenerics_0.32.0         lifecycle_0.1.0            
#> [123] ROTS_1.14.0                 munsell_0.5.0              
#> [125] R.methodsS3_1.7.1           moments_0.14               
#> [127] hwriter_1.3.2               caTools_1.18.0             
#> [129] codetools_0.2-16            coda_0.19-3                
#> [131] Biobase_2.46.0              GenomeInfoDb_1.22.0        
#> [133] vipor_0.4.5                 htmlTable_1.13.3           
#> [135] bayNorm_1.4.14              lsei_1.2-0                 
#> [137] rARPACK_0.11-0              xtable_1.8-4               
#> [139] SAVER_1.1.2                 ROCR_1.0-7                 
#> [141] diptest_0.75-7              formatR_1.7                
#> [143] lpsymphony_1.14.0           abind_1.4-5                
#> [145] FNN_1.1.3                   RANN_2.6.1                 
#> [147] askpass_1.1                 sparsesvd_0.2              
#> [149] CompQuadForm_1.4.3          GenomicRanges_1.38.0       
#> [151] tibble_2.1.3                ggdendro_0.1-20            
#> [153] cluster_2.1.0               future.apply_1.4.0         
#> [155] Matrix_1.2-18               prettyunits_1.1.1          
#> [157] shinyBS_0.61                NOISeq_2.30.0              
#> [159] shinydashboard_0.7.1        mclust_5.4.5               
#> [161] igraph_1.2.4.2              ggstance_0.3.3             
#> [163] slam_0.1-47                 testthat_2.3.1             
#> [165] doSNOW_1.0.18               htmltools_0.4.0            
#> [167] BiocFileCache_1.10.2        yaml_2.2.1                 
#> [169] GenomicFeatures_1.38.1      XML_3.99-0.3               
#> [171] ggpubr_0.2.4                DrImpute_1.0               
#> [173] foreign_0.8-75              fitdistrplus_1.0-14        
#> [175] BiocParallel_1.20.1         aroma.light_3.16.0         
#> [177] bit64_0.9-7                 rngtools_1.5               
#> [179] doRNG_1.8.2                 foreach_1.4.8              
#> [181] robustbase_0.93-5           outliers_0.14              
#> [183] Biostrings_2.54.0           combinat_0.0-8             
#> [185] rsvd_1.0.2                  iCOBRA_1.14.0              
#> [187] memoise_1.1.0               evaluate_0.14              
#> [189] VGAM_1.1-2                  nonnest2_0.5-2             
#> [191] geneplotter_1.64.0          permute_0.9-5              
#> [193] curl_4.3                    fdrtool_1.2.15             
#> [195] acepack_1.4.1               edgeR_3.28.0               
#> [197] checkmate_2.0.0             npsurv_0.4-0               
#> [199] truncnorm_1.0-8             DECENT_1.1.0               
#> [201] tensorA_0.36.1              ellipse_0.4.1              
#> [203] ggplot2_3.2.1               ggrepel_0.8.1              
#> [205] scDD_1.10.0                 tools_3.6.2                
#> [207] sandwich_2.5-1              magrittr_1.5               
#> [209] RCurl_1.98-1.1              pbivnorm_0.6.0             
#> [211] bayesm_3.1-4                EBSeq_1.26.0               
#> [213] httr_1.4.1                  assertthat_0.2.1           
#> [215] rmarkdown_2.1               boot_1.3-24                
#> [217] globals_0.12.5              R6_2.4.1                   
#> [219] Rhdf5lib_1.8.0              nnet_7.3-12                
#> [221] progress_1.2.2              genefilter_1.68.0          
#> [223] gtools_3.8.1                statmod_1.4.33             
#> [225] BiocSingular_1.2.1          rhdf5_2.30.1               
#> [227] splines_3.6.2               colorspace_1.4-1           
#> [229] amap_0.8-18                 generics_0.0.2             
#> [231] stats4_3.6.2                NBPSeq_0.3.0               
#> [233] base64enc_0.1-3             compositions_1.40-3        
#> [235] baySeq_2.20.0               pillar_1.4.3               
#> [237] HSMMSingleCell_1.6.0        GenomeInfoDbData_1.2.2     
#> [239] plyr_1.8.5                  dotCall64_1.0-0            
#> [241] gtable_0.3.0                SCnorm_1.8.2               
#> [243] monocle_2.14.0              knitr_1.28                 
#> [245] RcppArmadillo_0.9.850.1.0   latticeExtra_0.6-29        
#> [247] biomaRt_2.42.0              IRanges_2.20.2             
#> [249] fastmap_1.0.1               doParallel_1.0.15          
#> [251] pscl_1.5.2                  flexmix_2.3-15             
#> [253] quantreg_5.54               AnnotationDbi_1.48.0       
#> [255] broom_0.5.4                 openssl_1.4.1              
#> [257] scales_1.1.0                arm_1.10-1                 
#> [259] backports_1.1.5             plotrix_3.7-7              
#> [261] IHW_1.14.0                  S4Vectors_0.24.3           
#> [263] densityClust_0.3            lme4_1.1-21                
#> [265] hms_0.5.3                   DESeq_1.38.0               
#> [267] Rtsne_0.15                  dplyr_0.8.4                
#> [269] shiny_1.4.0                 grid_3.6.2                 
#> [271] lazyeval_0.2.2              Formula_1.2-3              
#> [273] blockmodeling_0.3.4         crayon_1.3.4               
#> [275] MAST_1.12.0                 RUVSeq_1.20.0              
#> [277] viridis_0.5.1               rpart_4.1-15               
#> [279] zinbwave_1.8.0              compiler_3.6.2
```

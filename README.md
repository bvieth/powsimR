
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

Please also consult my Github Page of [powsimR](https://bvieth.github.io/powsimR/) made with [pkgdown](http://pkgdown.r-lib.org/index.html)!

Installation Guide
------------------

For the installation, the R package `devtools` is needed.

``` r
install.packages("devtools")
library(devtools)
```

I recommend to install first the dependencies manually and then powsimR:

``` r
ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    # new.pkg <- pkg
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
cranpackages <- c("bbmle", "broom", "cobs", "cowplot", "data.table", "devtools", 
    "doParallel", "dplyr", "drc", "DrImpute", "fastICA", "fitdistrplus", "foreach", 
    "gamlss.dist", "ggExtra", "ggplot2", "ggthemes", "grDevices", "glmnet", 
    "grid", "gtools", "Hmisc", "kernlab", "MASS", "matrixStats", "mclust", "methods", 
    "minpack.lm", "moments", "msir", "NBPSeq", "nonnest2", "parallel", "penalized", 
    "plyr", "pscl", "reshape2", "ROCR", "Rtsne", "scales", "Seurat", "snow", 
    "stats", "tibble", "tidyr", "VGAM", "ZIM")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("AnnotationDbi", "baySeq", "Biobase", "BiocGenerics", "BiocParallel", 
    "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", "Linnorm", 
    "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", 
    "scDD", "scde", "scone", "scran", "SCnorm", "SingleCellExperiment", "SummarizedExperiment", 
    "zinbwave")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("nghiavtr/BPSC", "cz-ye/DECENT", "mohuangx/SAVER", "statOmics/zingeR")
ipak(githubpackages, repository = "github")
```

After installing the dependencies, powsimR can be installed by using devtools as well.

``` r
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
library("powsimR")
```

Alternative, you can try to install powsimR and its dependencies directly using devtools:

``` r
devtools::install_github("bvieth/powsimR")
```

User Guide
----------

For examples and tips on using the package, please consult the vignette after successful installation by

``` r
browseVignettes("powsimR")
```

Some users have experienced issues installing powsimR due to vignette compilation errors. If that is the case, you can leave out building the vignette and read it on my Github Page of [powsimR](https://bvieth.github.io/powsimR/articles/powsimR.html) or download it as a html file [here](https://github.com/bvieth/powsimR/tree/master/inst/doc/powsimR.html).

### DLLs and ulimit

Note that the error "maximal number of DLLs reached..." might occur due to the loading of many shared objects by Bioconductor packages. Restarting the R session after installing dependencies / powsimR will help. Starting with R version 3.4.0, one can set the environmental variable 'R\_MAX\_NUM\_DLLS' to a higher number. See `?Startup()` for more information. I recommend to increase the maximum number of DLLs that can be loaded to at least 500. The environmental variable R\_MAX\_NUM\_DLLS can be set in R\_HOME/etc/Renviron prior to starting R. For that locate the Renviron file and add the following line: R\_MAX\_NUM\_DLLS=xy where xy is the number of DLLs. On my Ubuntu machine, the Renviron file is in /usr/lib/R/etc/ and I can set it to 500.

In addition, the user limits for open files (unix: ulimit) might have to be set to a higher number to accomodate the increase in DLLs. Please check out the help pages for [MACs](https://gist.github.com/tombigel/d503800a282fcadbee14b537735d202c) and [Linux](https://glassonionblog.wordpress.com/2013/01/27/increase-ulimit-and-file-descriptors-limit/) for guidance.

Citation
--------

Please use the following entry for citing powsimR.

``` r
citation("powsimR")
```

powsimR is published in [Bioinformatics](https://doi.org/10.1093/bioinformatics/btx435). A preprint paper is also on [bioRxiv](https://doi.org/10.1101/117150).

Notes
-----

Please send bug reports and feature requests by opening a new issue on [this page](https://github.com/bvieth/powsimR/issues).

R Session Info
--------------

``` r
library(powsimR)
#> Loading required package: gamlss.dist
#> Loading required package: MASS
#> Warning: replacing previous import 'DECENT::lrTest' by 'MAST::lrTest' when
#> loading 'powsimR'
#> Warning: replacing previous import 'parallel::makeCluster' by
#> 'snow::makeCluster' when loading 'powsimR'
#> Warning: replacing previous import 'parallel::stopCluster' by
#> 'snow::stopCluster' when loading 'powsimR'
#> Warning: replacing previous import 'penalized::predict' by 'stats::predict'
#> when loading 'powsimR'
#> Warning: replacing previous import 'zinbwave::glmWeightedF' by
#> 'zingeR::glmWeightedF' when loading 'powsimR'
sessionInfo()
#> R version 3.5.1 (2018-07-02)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.1 LTS
#> 
#> Matrix products: default
#> BLAS: /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
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
#> [1] powsimR_1.1.3     gamlss.dist_5.0-6 MASS_7.3-50      
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            pbapply_1.3-4              
#>   [5] haven_1.1.2                 lattice_0.20-35            
#>   [7] fastICA_1.2-1               mgcv_1.8-24                
#>   [9] penalized_0.9-51            blob_1.1.1                 
#>  [11] survival_2.42-6             later_0.7.4                
#>  [13] nloptr_1.0.4                DBI_1.0.0                  
#>  [15] R.utils_2.7.0               SingleCellExperiment_1.2.0 
#>  [17] Linnorm_2.4.0               bindr_0.1.1                
#>  [19] zlibbioc_1.26.0             MatrixModels_0.4-1         
#>  [21] pspline_1.0-18              pcaMethods_1.72.0          
#>  [23] SDMTools_1.1-221            htmlwidgets_1.2            
#>  [25] mvtnorm_1.0-8               hdf5r_1.0.0                
#>  [27] UpSetR_1.3.3                parallel_3.5.1             
#>  [29] scater_1.8.4                irlba_2.3.2                
#>  [31] DEoptimR_1.0-8              lars_1.2                   
#>  [33] Rcpp_0.12.18                KernSmooth_2.23-15         
#>  [35] DT_0.4                      promises_1.0.1             
#>  [37] gdata_2.18.0                DDRTree_0.1.5              
#>  [39] DelayedArray_0.6.5          limma_3.36.3               
#>  [41] vegan_2.5-2                 Hmisc_4.1-1                
#>  [43] ShortRead_1.38.0            apcluster_1.4.7            
#>  [45] RSpectra_0.13-1             msir_1.3.1                 
#>  [47] mnormt_1.5-5                digest_0.6.16              
#>  [49] png_0.1-7                   qlcMatrix_0.9.7            
#>  [51] cowplot_0.9.3               glmnet_2.0-16              
#>  [53] pkgconfig_2.0.2             docopt_0.6                 
#>  [55] DelayedMatrixStats_1.2.0    ggbeeswarm_0.6.0           
#>  [57] iterators_1.0.10            minqa_1.2.4                
#>  [59] lavaan_0.6-2                reticulate_1.10            
#>  [61] SummarizedExperiment_1.10.1 spam_2.2-0                 
#>  [63] beeswarm_0.2.3              modeltools_0.2-22          
#>  [65] zoo_1.8-3                   tidyselect_0.2.4           
#>  [67] ZIM_1.1.0                   reshape2_1.4.3             
#>  [69] purrr_0.2.5                 kernlab_0.9-27             
#>  [71] ica_1.0-2                   pcaPP_1.9-73               
#>  [73] EDASeq_2.14.1               viridisLite_0.3.0          
#>  [75] snow_0.4-2                  rtracklayer_1.40.6         
#>  [77] rlang_0.2.2                 hexbin_1.27.2              
#>  [79] manipulateWidget_0.10.0     glue_1.3.0                 
#>  [81] metap_1.0                   RColorBrewer_1.1-2         
#>  [83] registry_0.5                fpc_2.1-11.1               
#>  [85] matrixStats_0.54.0          stringr_1.3.1              
#>  [87] pkgmaker_0.27               fields_9.6                 
#>  [89] DESeq2_1.20.0               SparseM_1.77               
#>  [91] gbRd_0.4-11                 httpuv_1.4.5               
#>  [93] class_7.3-14                BPSC_0.99.2                
#>  [95] RMTstat_0.3                 annotate_1.58.0            
#>  [97] webshot_0.5.0               jsonlite_1.5               
#>  [99] XVector_0.20.0              bit_1.1-14                 
#> [101] mime_0.5                    gridExtra_2.3              
#> [103] gplots_3.0.1                Rsamtools_1.32.3           
#> [105] zingeR_0.1.0                stringi_1.2.4              
#> [107] gmodels_2.18.1              gsl_1.9-10.3               
#> [109] bitops_1.0-6                Rdpack_0.9-0               
#> [111] maps_3.3.0                  RSQLite_2.1.1              
#> [113] tidyr_0.8.1                 pheatmap_1.0.10            
#> [115] data.table_1.11.4           DEDS_1.54.0                
#> [117] energy_1.7-5                rstudioapi_0.7             
#> [119] GenomicAlignments_1.16.0    nlme_3.1-137               
#> [121] qvalue_2.12.0               scran_1.8.4                
#> [123] fastcluster_1.1.25          scone_1.4.0                
#> [125] locfit_1.5-9.1              miniUI_0.1.1.1             
#> [127] cobs_1.3-3                  R.oo_1.22.0                
#> [129] prabclus_2.2-6              segmented_0.5-3.0          
#> [131] BiocGenerics_0.26.0         readxl_1.1.0               
#> [133] ROTS_1.8.0                  cellranger_1.1.0           
#> [135] munsell_0.5.0               R.methodsS3_1.7.1          
#> [137] moments_0.14                hwriter_1.3.2              
#> [139] caTools_1.17.1.1            codetools_0.2-15           
#> [141] coda_0.19-1                 Biobase_2.40.0             
#> [143] GenomeInfoDb_1.16.0         vipor_0.4.5                
#> [145] lmtest_0.9-36               htmlTable_1.12             
#> [147] rARPACK_0.11-0              lsei_1.2-0                 
#> [149] xtable_1.8-3                SAVER_1.1.0                
#> [151] ROCR_1.0-7                  diptest_0.75-7             
#> [153] formatR_1.5                 lpsymphony_1.8.0           
#> [155] abind_1.4-5                 FNN_1.1.2.1                
#> [157] RANN_2.6                    sparsesvd_0.1-4            
#> [159] CompQuadForm_1.4.3          GenomicRanges_1.32.6       
#> [161] bibtex_0.4.2                rgl_0.99.16                
#> [163] tibble_1.4.2                ggdendro_0.1-20            
#> [165] cluster_2.0.7-1             Seurat_2.3.4               
#> [167] Matrix_1.2-14               prettyunits_1.0.2          
#> [169] shinyBS_0.61                ggridges_0.5.0             
#> [171] NOISeq_2.24.0               shinydashboard_0.7.0       
#> [173] mclust_5.4.1                igraph_1.2.2               
#> [175] slam_0.1-43                 testthat_2.0.0             
#> [177] doSNOW_1.0.16               htmltools_0.3.6            
#> [179] yaml_2.2.0                  GenomicFeatures_1.32.2     
#> [181] XML_3.98-1.16               DrImpute_1.0               
#> [183] foreign_0.8-71              withr_2.1.2                
#> [185] fitdistrplus_1.0-11         BiocParallel_1.14.2        
#> [187] aroma.light_3.10.0          bit64_0.9-7                
#> [189] rngtools_1.3.1              doRNG_1.7.1                
#> [191] foreach_1.4.4               robustbase_0.93-2          
#> [193] outliers_0.14               scde_2.8.0                 
#> [195] Biostrings_2.48.0           combinat_0.0-8             
#> [197] rsvd_0.9                    iCOBRA_1.8.0               
#> [199] memoise_1.1.0               evaluate_0.11              
#> [201] VGAM_1.0-6                  nonnest2_0.5-1             
#> [203] forcats_0.3.0               rio_0.5.10                 
#> [205] geneplotter_1.58.0          permute_0.9-4              
#> [207] curl_3.2                    fdrtool_1.2.15             
#> [209] trimcluster_0.1-2.1         acepack_1.4.1              
#> [211] edgeR_3.22.3                checkmate_1.8.5            
#> [213] npsurv_0.4-0                tensorA_0.36.1             
#> [215] DECENT_0.99.1               ellipse_0.4.1              
#> [217] ggplot2_3.0.0               rjson_0.2.20               
#> [219] openxlsx_4.1.0              ggrepel_0.8.0              
#> [221] distillery_1.0-4            dtw_1.20-1                 
#> [223] scDD_1.4.0                  rprojroot_1.3-2            
#> [225] stabledist_0.7-1            Lmoments_1.2-3             
#> [227] tools_3.5.1                 sandwich_2.5-0             
#> [229] magrittr_1.5                RCurl_1.95-4.11            
#> [231] proxy_0.4-22                car_3.0-2                  
#> [233] pbivnorm_0.6.0              ape_5.1                    
#> [235] bayesm_3.1-0.1              EBSeq_1.20.0               
#> [237] httr_1.3.1                  assertthat_0.2.0           
#> [239] rmarkdown_1.10              boot_1.3-20                
#> [241] R6_2.2.2                    Rhdf5lib_1.2.1             
#> [243] nnet_7.3-12                 progress_1.2.0             
#> [245] tximport_1.8.0              genefilter_1.62.0          
#> [247] gtools_3.8.1                statmod_1.4.30             
#> [249] Rook_1.1-1                  rhdf5_2.24.0               
#> [251] splines_3.5.1               carData_3.0-1              
#> [253] colorspace_1.3-2            amap_0.8-16                
#> [255] stats4_3.5.1                NBPSeq_0.3.0               
#> [257] compositions_1.40-2         base64enc_0.1-3            
#> [259] baySeq_2.14.0               pillar_1.3.0               
#> [261] HSMMSingleCell_0.114.0      bindrcpp_0.2.2             
#> [263] GenomeInfoDbData_1.1.0      plyr_1.8.4                 
#> [265] extRemes_2.0-9              dotCall64_1.0-0            
#> [267] gtable_0.2.0                zip_1.0.0                  
#> [269] SCnorm_1.2.1                monocle_2.8.0              
#> [271] knitr_1.20                  RcppArmadillo_0.9.100.5.0  
#> [273] latticeExtra_0.6-28         biomaRt_2.36.1             
#> [275] IRanges_2.14.11             ADGofTest_0.3              
#> [277] copula_0.999-18             crosstalk_1.0.0            
#> [279] Cairo_1.5-9                 doParallel_1.0.11          
#> [281] pscl_1.5.2                  flexmix_2.3-14             
#> [283] quantreg_5.36               AnnotationDbi_1.42.1       
#> [285] broom_0.5.0                 scales_1.0.0               
#> [287] arm_1.10-1                  backports_1.1.2            
#> [289] IHW_1.8.0                   S4Vectors_0.18.3           
#> [291] densityClust_0.3            lme4_1.1-18-1              
#> [293] brew_1.0-6                  hms_0.4.2                  
#> [295] DESeq_1.32.0                Rtsne_0.13                 
#> [297] dplyr_0.7.6                 shiny_1.1.0                
#> [299] grid_3.5.1                  numDeriv_2016.8-1          
#> [301] bbmle_1.0.20                lazyeval_0.2.1             
#> [303] dynamicTreeCut_1.63-1       Formula_1.2-3              
#> [305] tsne_0.1-3                  blockmodeling_0.3.1        
#> [307] crayon_1.3.4                MAST_1.6.1                 
#> [309] RUVSeq_1.14.0               viridis_0.5.1              
#> [311] rpart_4.1-13                zinbwave_1.2.0             
#> [313] compiler_3.5.1
```

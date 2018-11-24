
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
#> Warning: replacing previous import 'SingleCellExperiment::weights' by
#> 'stats::weights' when loading 'SCnorm'
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
#> [1] powsimR_1.1.4     gamlss.dist_5.1-0 MASS_7.3-51.1    
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            pbapply_1.3-4              
#>   [5] haven_2.0.0                 lattice_0.20-38            
#>   [7] fastICA_1.2-1               mgcv_1.8-26                
#>   [9] penalized_0.9-51            blob_1.1.1                 
#>  [11] survival_2.43-1             later_0.7.5                
#>  [13] nloptr_1.2.1                DBI_1.0.0                  
#>  [15] R.utils_2.7.0               SingleCellExperiment_1.4.0 
#>  [17] Linnorm_2.6.0               bindr_0.1.1                
#>  [19] zlibbioc_1.28.0             MatrixModels_0.4-1         
#>  [21] pspline_1.0-18              pcaMethods_1.74.0          
#>  [23] SDMTools_1.1-221            htmlwidgets_1.3            
#>  [25] mvtnorm_1.0-8               hdf5r_1.0.1                
#>  [27] UpSetR_1.3.3                parallel_3.5.1             
#>  [29] scater_1.10.0               irlba_2.3.2                
#>  [31] DEoptimR_1.0-8              lars_1.2                   
#>  [33] Rcpp_1.0.0                  KernSmooth_2.23-15         
#>  [35] DT_0.5                      promises_1.0.1             
#>  [37] gdata_2.18.0                DDRTree_0.1.5              
#>  [39] DelayedArray_0.8.0          limma_3.38.2               
#>  [41] vegan_2.5-3                 Hmisc_4.1-1                
#>  [43] ShortRead_1.40.0            apcluster_1.4.7            
#>  [45] RSpectra_0.13-1             msir_1.3.1                 
#>  [47] mnormt_1.5-5                digest_0.6.18              
#>  [49] png_0.1-7                   qlcMatrix_0.9.7            
#>  [51] cowplot_0.9.3               glmnet_2.0-16              
#>  [53] pkgconfig_2.0.2             docopt_0.6.1               
#>  [55] DelayedMatrixStats_1.4.0    ggbeeswarm_0.6.0           
#>  [57] iterators_1.0.10            minqa_1.2.4                
#>  [59] lavaan_0.6-3                reticulate_1.10            
#>  [61] SummarizedExperiment_1.12.0 spam_2.2-0                 
#>  [63] beeswarm_0.2.3              modeltools_0.2-22          
#>  [65] zoo_1.8-4                   tidyselect_0.2.5           
#>  [67] ZIM_1.1.0                   reshape2_1.4.3             
#>  [69] purrr_0.2.5                 kernlab_0.9-27             
#>  [71] ica_1.0-2                   pcaPP_1.9-73               
#>  [73] EDASeq_2.16.0               viridisLite_0.3.0          
#>  [75] snow_0.4-3                  rtracklayer_1.42.0         
#>  [77] rlang_0.3.0.1               hexbin_1.27.2              
#>  [79] manipulateWidget_0.10.0     glue_1.3.0                 
#>  [81] metap_1.0                   RColorBrewer_1.1-2         
#>  [83] registry_0.5                fpc_2.1-11.1               
#>  [85] matrixStats_0.54.0          stringr_1.3.1              
#>  [87] pkgmaker_0.27               fields_9.6                 
#>  [89] DESeq2_1.22.1               SparseM_1.77               
#>  [91] gbRd_0.4-11                 httpuv_1.4.5               
#>  [93] class_7.3-14                BPSC_0.99.2                
#>  [95] BiocNeighbors_1.0.0         RMTstat_0.3                
#>  [97] annotate_1.60.0             webshot_0.5.1              
#>  [99] jsonlite_1.5                XVector_0.22.0             
#> [101] bit_1.1-14                  mime_0.6                   
#> [103] gridExtra_2.3               gplots_3.0.1               
#> [105] Rsamtools_1.34.0            zingeR_0.1.0               
#> [107] stringi_1.2.4               gmodels_2.18.1             
#> [109] gsl_1.9-10.3                bitops_1.0-6               
#> [111] Rdpack_0.10-1               maps_3.3.0                 
#> [113] RSQLite_2.1.1               tidyr_0.8.2                
#> [115] pheatmap_1.0.10             data.table_1.11.8          
#> [117] DEDS_1.56.0                 energy_1.7-5               
#> [119] rstudioapi_0.8              GenomicAlignments_1.18.0   
#> [121] nlme_3.1-137                qvalue_2.14.0              
#> [123] scran_1.10.1                fastcluster_1.1.25         
#> [125] scone_1.6.0                 locfit_1.5-9.1             
#> [127] miniUI_0.1.1.1              cobs_1.3-3                 
#> [129] R.oo_1.22.0                 prabclus_2.2-6             
#> [131] segmented_0.5-3.0           BiocGenerics_0.28.0        
#> [133] readxl_1.1.0                ROTS_1.10.0                
#> [135] cellranger_1.1.0            munsell_0.5.0              
#> [137] R.methodsS3_1.7.1           moments_0.14               
#> [139] hwriter_1.3.2               caTools_1.17.1.1           
#> [141] codetools_0.2-15            coda_0.19-2                
#> [143] Biobase_2.42.0              GenomeInfoDb_1.18.1        
#> [145] vipor_0.4.5                 lmtest_0.9-36              
#> [147] htmlTable_1.12              rARPACK_0.11-0             
#> [149] lsei_1.2-0                  xtable_1.8-3               
#> [151] SAVER_1.1.1                 ROCR_1.0-7                 
#> [153] diptest_0.75-7              formatR_1.5                
#> [155] lpsymphony_1.10.0           abind_1.4-5                
#> [157] FNN_1.1.2.1                 RANN_2.6                   
#> [159] sparsesvd_0.1-4             CompQuadForm_1.4.3         
#> [161] GenomicRanges_1.34.0        bibtex_0.4.2               
#> [163] rgl_0.99.16                 tibble_1.4.2               
#> [165] ggdendro_0.1-20             cluster_2.0.7-1            
#> [167] Seurat_2.3.4                Matrix_1.2-15              
#> [169] prettyunits_1.0.2           shinyBS_0.61               
#> [171] ggridges_0.5.1              NOISeq_2.26.0              
#> [173] shinydashboard_0.7.1        mclust_5.4.2               
#> [175] igraph_1.2.2                slam_0.1-43                
#> [177] testthat_2.0.1              doSNOW_1.0.16              
#> [179] htmltools_0.3.6             yaml_2.2.0                 
#> [181] GenomicFeatures_1.34.1      XML_3.98-1.16              
#> [183] DrImpute_1.0                foreign_0.8-71             
#> [185] withr_2.1.2                 fitdistrplus_1.0-11        
#> [187] BiocParallel_1.16.0         aroma.light_3.12.0         
#> [189] bit64_0.9-7                 rngtools_1.3.1             
#> [191] doRNG_1.7.1                 foreach_1.4.4              
#> [193] robustbase_0.93-3           outliers_0.14              
#> [195] scde_2.10.0                 Biostrings_2.50.1          
#> [197] combinat_0.0-8              rsvd_1.0.0                 
#> [199] iCOBRA_1.10.0               memoise_1.1.0              
#> [201] evaluate_0.12               VGAM_1.0-6                 
#> [203] nonnest2_0.5-2              forcats_0.3.0              
#> [205] rio_0.5.10                  geneplotter_1.60.0         
#> [207] permute_0.9-4               curl_3.2                   
#> [209] fdrtool_1.2.15              trimcluster_0.1-2.1        
#> [211] acepack_1.4.1               edgeR_3.24.0               
#> [213] checkmate_1.8.5             npsurv_0.4-0               
#> [215] tensorA_0.36.1              DECENT_0.99.1              
#> [217] ellipse_0.4.1               ggplot2_3.1.0              
#> [219] rjson_0.2.20                openxlsx_4.1.0             
#> [221] ggrepel_0.8.0               distillery_1.0-4           
#> [223] dtw_1.20-1                  scDD_1.6.0                 
#> [225] rprojroot_1.3-2             stabledist_0.7-1           
#> [227] Lmoments_1.2-3              tools_3.5.1                
#> [229] sandwich_2.5-0              magrittr_1.5               
#> [231] RCurl_1.95-4.11             proxy_0.4-22               
#> [233] car_3.0-2                   pbivnorm_0.6.0             
#> [235] ape_5.2                     bayesm_3.1-0.1             
#> [237] EBSeq_1.22.0                httr_1.3.1                 
#> [239] assertthat_0.2.0            rmarkdown_1.10             
#> [241] boot_1.3-20                 R6_2.3.0                   
#> [243] Rhdf5lib_1.4.0              nnet_7.3-12                
#> [245] progress_1.2.0              genefilter_1.64.0          
#> [247] gtools_3.8.1                statmod_1.4.30             
#> [249] Rook_1.1-1                  HDF5Array_1.10.0           
#> [251] rhdf5_2.26.0                splines_3.5.1              
#> [253] carData_3.0-2               colorspace_1.3-2           
#> [255] amap_0.8-16                 stats4_3.5.1               
#> [257] NBPSeq_0.3.0                compositions_1.40-2        
#> [259] base64enc_0.1-3             baySeq_2.16.0              
#> [261] pillar_1.3.0                HSMMSingleCell_1.2.0       
#> [263] bindrcpp_0.2.2              GenomeInfoDbData_1.2.0     
#> [265] plyr_1.8.4                  extRemes_2.0-9             
#> [267] dotCall64_1.0-0             gtable_0.2.0               
#> [269] zip_1.0.0                   SCnorm_1.4.0               
#> [271] monocle_2.10.0              knitr_1.20                 
#> [273] RcppArmadillo_0.9.200.4.0   latticeExtra_0.6-28        
#> [275] biomaRt_2.38.0              IRanges_2.16.0             
#> [277] ADGofTest_0.3               copula_0.999-18            
#> [279] crosstalk_1.0.0             Cairo_1.5-9                
#> [281] doParallel_1.0.14           pscl_1.5.2                 
#> [283] flexmix_2.3-14              quantreg_5.36              
#> [285] AnnotationDbi_1.44.0        broom_0.5.0                
#> [287] scales_1.0.0                arm_1.10-1                 
#> [289] backports_1.1.2             IHW_1.10.0                 
#> [291] S4Vectors_0.20.1            densityClust_0.3           
#> [293] lme4_1.1-19                 brew_1.0-6                 
#> [295] hms_0.4.2                   DESeq_1.34.0               
#> [297] Rtsne_0.15                  dplyr_0.7.8                
#> [299] shiny_1.2.0                 grid_3.5.1                 
#> [301] numDeriv_2016.8-1           bbmle_1.0.20               
#> [303] lazyeval_0.2.1              dynamicTreeCut_1.63-1      
#> [305] Formula_1.2-3               tsne_0.1-3                 
#> [307] blockmodeling_0.3.1         crayon_1.3.4               
#> [309] MAST_1.8.0                  RUVSeq_1.16.0              
#> [311] viridis_0.5.1               rpart_4.1-13               
#> [313] zinbwave_1.4.0              compiler_3.5.1
```

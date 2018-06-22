
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

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

Some users have experienced issues installing powsimR due to vignette compilation errors. If that is the case, you can leave out building the vignette and read it [here](https://github.com/bvieth/powsimR/tree/master/inst/doc/powsimR.html) instead.

### DLLs and ulimit

Note that the error "maximal number of DLLs reached..." might occur due to the loading of many shared objects by Bioconductor packages. Restarting the R session after installing dependencies / powsimR will help. Starting with R version 3.4.0, one can set the environmental variable 'R\_MAX\_NUM\_DLLS' to a higher number. See `?Startup()` for more information. I recommend to increase the maximum number of DLLs that can be loaded to at least 500. The environmental variable R\_MAX\_NUM\_DLLS can be set in R\_HOME/etc/Renviron prior to starting R. For that locate the Renviron file and add the following line: R\_MAX\_NUM\_DLLS=xy where xy is the number of DLLs. On my Ubuntu machine, the Renviron file is in /usr/lib/R/etc/ and I can set it to 500.

In addition, the user limits for open files (unix: ulimit) might have to be set to a higher number to accomodate the increase in DLLs. Please check out the help pages for [MACs](https://gist.github.com/tombigel/d503800a282fcadbee14b537735d202c) and [Linux](https://glassonionblog.wordpress.com/2013/01/27/increase-ulimit-and-file-descriptors-limit/) for guidance.

User Guide
----------

For examples and tips on using the package, please see the vignette [here](https://github.com/bvieth/powsimR/tree/master/inst/doc/powsimR.html).

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
#> other attached packages:
#> [1] powsimR_1.1.2     gamlss.dist_5.0-6 MASS_7.3-50      
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0             softImpute_1.4            
#>   [3] minpack.lm_1.2-1           pbapply_1.3-4             
#>   [5] haven_1.1.1                lattice_0.20-35           
#>   [7] fastICA_1.2-1              mgcv_1.8-24               
#>   [9] penalized_0.9-50           blob_1.1.1                
#>  [11] survival_2.42-3            prodlim_2018.04.18        
#>  [13] later_0.7.3                nloptr_1.0.4              
#>  [15] DBI_1.0.0                  R.utils_2.6.0             
#>  [17] SingleCellExperiment_1.0.0 Linnorm_2.2.0             
#>  [19] bindr_0.1.1                zlibbioc_1.24.0           
#>  [21] MatrixModels_0.4-1         pspline_1.0-18            
#>  [23] pcaMethods_1.70.0          pls_2.6-0                 
#>  [25] SDMTools_1.1-221           htmlwidgets_1.2           
#>  [27] mvtnorm_1.0-8              hdf5r_1.0.0               
#>  [29] UpSetR_1.3.3               tclust_1.4-1              
#>  [31] parallel_3.4.4             scater_1.6.3              
#>  [33] irlba_2.3.2                DEoptimR_1.0-8            
#>  [35] lars_1.2                   Rcpp_0.12.17              
#>  [37] KernSmooth_2.23-15         DT_0.4                    
#>  [39] promises_1.0.1             gdata_2.18.0              
#>  [41] DDRTree_0.1.5              DelayedArray_0.4.1        
#>  [43] limma_3.34.9               vegan_2.5-2               
#>  [45] CVST_0.2-2                 Hmisc_4.1-1               
#>  [47] ShortRead_1.36.1           apcluster_1.4.7           
#>  [49] RSpectra_0.13-1            msir_1.3.1                
#>  [51] mnormt_1.5-5               ranger_0.10.1             
#>  [53] digest_0.6.15              RMySQL_0.10.15            
#>  [55] png_0.1-7                  qlcMatrix_0.9.7           
#>  [57] cowplot_0.9.2              glmnet_2.0-16             
#>  [59] pkgconfig_2.0.1            docopt_0.4.5              
#>  [61] gower_0.1.2                ggbeeswarm_0.6.0          
#>  [63] iterators_1.0.9            minqa_1.2.4               
#>  [65] lavaan_0.6-1               reticulate_1.8            
#>  [67] SummarizedExperiment_1.8.1 spam_2.2-0                
#>  [69] beeswarm_0.2.3             modeltools_0.2-21         
#>  [71] RcppNumerical_0.3-2        zoo_1.8-2                 
#>  [73] tidyselect_0.2.4           ZIM_1.0.3                 
#>  [75] reshape2_1.4.3             purrr_0.2.5               
#>  [77] kernlab_0.9-26             ica_1.0-2                 
#>  [79] pcaPP_1.9-73               EDASeq_2.12.0             
#>  [81] viridisLite_0.3.0          snow_0.4-2                
#>  [83] rtracklayer_1.38.3         rlang_0.2.1               
#>  [85] hexbin_1.27.2              manipulateWidget_0.10.0   
#>  [87] glue_1.2.0                 metap_0.9                 
#>  [89] RColorBrewer_1.1-2         registry_0.5              
#>  [91] fpc_2.1-11                 matrixStats_0.53.1        
#>  [93] stringr_1.3.1              pkgmaker_0.27             
#>  [95] lava_1.6.1                 fields_9.6                
#>  [97] DESeq2_1.18.1              recipes_0.1.3             
#>  [99] SparseM_1.77               httpuv_1.4.4.1            
#> [101] class_7.3-14               BPSC_0.99.1               
#> [103] RMTstat_0.3                annotate_1.56.2           
#> [105] webshot_0.5.0              jsonlite_1.5              
#> [107] XVector_0.18.0             bit_1.1-14                
#> [109] mime_0.5                   gridExtra_2.3             
#> [111] gplots_3.0.1               Rsamtools_1.30.0          
#> [113] zingeR_0.1.0               stringi_1.2.3             
#> [115] gmodels_2.16.2             RcppRoll_0.3.0            
#> [117] gsl_1.9-10.3               bitops_1.0-6              
#> [119] maps_3.3.0                 RSQLite_2.1.1             
#> [121] tidyr_0.8.1                pheatmap_1.0.10           
#> [123] data.table_1.11.4          DEDS_1.52.0               
#> [125] energy_1.7-4               rstudioapi_0.7            
#> [127] GenomicAlignments_1.14.2   sfsmisc_1.1-2             
#> [129] nlme_3.1-137               qvalue_2.10.0             
#> [131] scran_1.6.9                fastcluster_1.1.25        
#> [133] scone_1.2.0                locfit_1.5-9.1            
#> [135] miniUI_0.1.1.1             cobs_1.3-3                
#> [137] R.oo_1.22.0                prabclus_2.2-6            
#> [139] segmented_0.5-3.0          BiocGenerics_0.24.0       
#> [141] readxl_1.1.0               dimRed_0.1.0              
#> [143] timeDate_3043.102          ROTS_1.6.0                
#> [145] cellranger_1.1.0           munsell_0.5.0             
#> [147] R.methodsS3_1.7.1          moments_0.14              
#> [149] hwriter_1.3.2              caTools_1.17.1            
#> [151] codetools_0.2-15           coda_0.19-1               
#> [153] Biobase_2.38.0             magic_1.5-8               
#> [155] GenomeInfoDb_1.14.0        diffusionMap_1.1-0        
#> [157] vipor_0.4.5                lmtest_0.9-36             
#> [159] htmlTable_1.12             rARPACK_0.11-0            
#> [161] xtable_1.8-2               SAVER_0.4.0               
#> [163] ROCR_1.0-7                 diptest_0.75-7            
#> [165] formatR_1.5                scatterplot3d_0.3-41      
#> [167] lpsymphony_1.7.1           abind_1.4-5               
#> [169] FNN_1.1                    RANN_2.5.1                
#> [171] sparsesvd_0.1-4            CompQuadForm_1.4.3        
#> [173] GenomicRanges_1.30.3       bibtex_0.4.2              
#> [175] rgl_0.99.16                tibble_1.4.2              
#> [177] ggdendro_0.1-20            cluster_2.0.7-1           
#> [179] Seurat_2.3.2               Matrix_1.2-14             
#> [181] prettyunits_1.0.2          shinyBS_0.61              
#> [183] lubridate_1.7.4            ggridges_0.5.0            
#> [185] NOISeq_2.22.1              shinydashboard_0.7.0      
#> [187] mclust_5.4                 igraph_1.2.1              
#> [189] slam_0.1-43                testthat_2.0.0            
#> [191] doSNOW_1.0.16              geometry_0.3-6            
#> [193] htmltools_0.3.6            yaml_2.1.19               
#> [195] GenomicFeatures_1.30.3     XML_3.98-1.11             
#> [197] ModelMetrics_1.1.0         DrImpute_1.0              
#> [199] foreign_0.8-70             withr_2.1.2               
#> [201] fitdistrplus_1.0-9         BiocParallel_1.12.0       
#> [203] aroma.light_3.8.0          bit64_0.9-7               
#> [205] rngtools_1.3.1             doRNG_1.6.6               
#> [207] foreach_1.4.4              robustbase_0.93-0         
#> [209] outliers_0.14              scde_2.6.0                
#> [211] Biostrings_2.46.0          combinat_0.0-8            
#> [213] iCOBRA_1.6.0               memoise_1.1.0             
#> [215] evaluate_0.10.1            VGAM_1.0-5                
#> [217] nonnest2_0.5-1             forcats_0.3.0             
#> [219] rio_0.5.10                 geneplotter_1.56.0        
#> [221] permute_0.9-4              caret_6.0-80              
#> [223] curl_3.2                   fdrtool_1.2.15            
#> [225] trimcluster_0.1-2          acepack_1.4.1             
#> [227] edgeR_3.20.9               checkmate_1.8.5           
#> [229] tensorA_0.36               DECENT_0.2.0              
#> [231] ellipse_0.4.1              ggplot2_2.2.1             
#> [233] rjson_0.2.20               openxlsx_4.1.0            
#> [235] ggrepel_0.8.0              distillery_1.0-4          
#> [237] dtw_1.20-1                 scDD_1.2.0                
#> [239] rprojroot_1.3-2            stabledist_0.7-1          
#> [241] Lmoments_1.2-3             tools_3.4.4               
#> [243] sandwich_2.4-0             magrittr_1.5              
#> [245] RCurl_1.95-4.10            proxy_0.4-22              
#> [247] car_3.0-0                  pbivnorm_0.6.0            
#> [249] ape_5.1                    bayesm_3.1-0.1            
#> [251] EBSeq_1.18.0               httr_1.3.1                
#> [253] assertthat_0.2.0           rmarkdown_1.10            
#> [255] boot_1.3-20                R6_2.2.2                  
#> [257] nnet_7.3-12                progress_1.2.0            
#> [259] tximport_1.6.0             genefilter_1.60.0         
#> [261] gtools_3.5.0               statmod_1.4.30            
#> [263] Rook_1.1-1                 rhdf5_2.22.0              
#> [265] splines_3.4.4              carData_3.0-1             
#> [267] colorspace_1.3-2           amap_0.8-16               
#> [269] stats4_3.4.4               NBPSeq_0.3.0              
#> [271] compositions_1.40-2        base64enc_0.1-3           
#> [273] baySeq_2.12.0              pillar_1.2.3              
#> [275] HSMMSingleCell_0.112.0     bindrcpp_0.2.2            
#> [277] GenomeInfoDbData_1.0.0     plyr_1.8.4                
#> [279] extRemes_2.0-9             dotCall64_0.9-5.2         
#> [281] gtable_0.2.0               zip_1.0.0                 
#> [283] SCnorm_1.0.0               monocle_2.6.4             
#> [285] psych_1.8.4                knitr_1.20                
#> [287] RcppArmadillo_0.8.500.0    latticeExtra_0.6-28       
#> [289] biomaRt_2.34.2             IRanges_2.12.0            
#> [291] ADGofTest_0.3              copula_0.999-18           
#> [293] Cairo_1.5-9                crosstalk_1.0.0           
#> [295] doParallel_1.0.11          pscl_1.5.2                
#> [297] flexmix_2.3-14             quantreg_5.36             
#> [299] AnnotationDbi_1.40.0       broom_0.4.4               
#> [301] scales_0.5.0               arm_1.10-1                
#> [303] backports_1.1.2            IHW_1.6.0                 
#> [305] S4Vectors_0.16.0           densityClust_0.3          
#> [307] ipred_0.9-6                lme4_1.1-17               
#> [309] brew_1.0-6                 hms_0.4.2                 
#> [311] DESeq_1.30.0               Rtsne_0.13                
#> [313] dplyr_0.7.5                shiny_1.1.0               
#> [315] ddalpha_1.3.3              grid_3.4.4                
#> [317] numDeriv_2016.8-1          bbmle_1.0.20              
#> [319] lazyeval_0.2.1             dynamicTreeCut_1.63-1     
#> [321] Formula_1.2-3              tsne_0.1-3                
#> [323] blockmodeling_0.3.1        crayon_1.3.4              
#> [325] DRR_0.0.3                  MAST_1.4.1                
#> [327] RUVSeq_1.12.0              viridis_0.5.1             
#> [329] rpart_4.1-13               zinbwave_1.0.0            
#> [331] compiler_3.4.4
```

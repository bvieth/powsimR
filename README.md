
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

NEWS
----

**Version 1.1.2** released with the following fixes (2018-06-20):

-   `estimateParam` error fixed concerning expression cleanup.
-   precompiled vignette in inst/doc/.

**Version 1.1.1** released with the following changes / additions (2018-04-19):

-   `simulateDE` now with the option to perform DE testing on filtered/imputed counts (option `DEFilter`)

**Version 1.1.0** released with the following changes / additions (2018-03-29):

-   simulation of batch effects (see options `p.B`, `bLFC` and `bPattern` in `DESetup` and `simulateCounts`)
-   simulation of spike-in expression (see `estimateSpike` , `plotSpike` and option `spikeIns` in `simulateDE` and `simulateCounts`)
-   simulation of multiple sample groups (e.g. single cell populations) with `simulateCounts`
-   imputation and prefiltering options prior to normalisation in DE power simulations added (scImpute, scone, Seurat, DrImpute, SAVER)
-   additional normalisation options and DE tools (esp. for single cells) included in `simulateDE`
-   evaluation of simulation setup using estimated versus true value comparisons of library size factors and log2 fold changes in `evaluateSim` and `plotEvalSim`
-   increased flexibility in preprocessing for distribution evaluation in `evaluateDist`

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
cranpackages <- c("bbmle", "broom", "cluster", "cobs", "cowplot", "data.table", 
    "devtools", "doParallel", "dplyr", "drc", "DrImpute", "fastICA", "fitdistrplus", 
    "foreach", "gamlss.dist", "ggExtra", "ggplot2", "ggthemes", "grDevices", 
    "glmnet", "grid", "gtools", "Hmisc", "kernlab", "MASS", "matrixStats", "mclust", 
    "methods", "minpack.lm", "moments", "msir", "NBPSeq", "nonnest2", "parallel", 
    "penalized", "plyr", "pscl", "reshape2", "ROCR", "Rtsne", "scales", "Seurat", 
    "snow", "stats", "tibble", "tidyr", "VGAM", "ZIM")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("AnnotationDbi", "baySeq", "Biobase", "BiocGenerics", "BiocParallel", 
    "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", "Linnorm", 
    "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", 
    "scDD", "scde", "scone", "scran", "SCnorm", "SingleCellExperiment", "SummarizedExperiment", 
    "zinbwave")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB
githubpackages <- c("nghiavtr/BPSC", "VCCRI/cidr", "cz-ye/DECENT", "mohuangx/SAVER", 
    "statOmics/zingeR")
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

For examples and tips on using the package, please see the vignette [here](https://github.com/bvieth/powsimR/tree/master/vignettes/powsimR.html)

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
#>  [45] CVST_0.2-2                 RcppParallel_4.4.0        
#>  [47] Hmisc_4.1-1                ShortRead_1.36.1          
#>  [49] apcluster_1.4.7            RSpectra_0.13-1           
#>  [51] msir_1.3.1                 mnormt_1.5-5              
#>  [53] ranger_0.10.1              digest_0.6.15             
#>  [55] RMySQL_0.10.15             png_0.1-7                 
#>  [57] qlcMatrix_0.9.7            cidr_0.1.5                
#>  [59] cowplot_0.9.2              glmnet_2.0-16             
#>  [61] pkgconfig_2.0.1            docopt_0.4.5              
#>  [63] gower_0.1.2                ggbeeswarm_0.6.0          
#>  [65] iterators_1.0.9            minqa_1.2.4               
#>  [67] lavaan_0.6-1               reticulate_1.8            
#>  [69] SummarizedExperiment_1.8.1 spam_2.2-0                
#>  [71] beeswarm_0.2.3             modeltools_0.2-21         
#>  [73] RcppNumerical_0.3-2        zoo_1.8-2                 
#>  [75] tidyselect_0.2.4           clusterCrit_1.2.7         
#>  [77] ZIM_1.0.3                  reshape2_1.4.3            
#>  [79] purrr_0.2.5                kernlab_0.9-26            
#>  [81] ica_1.0-2                  pcaPP_1.9-73              
#>  [83] EDASeq_2.12.0              viridisLite_0.3.0         
#>  [85] snow_0.4-2                 rtracklayer_1.38.3        
#>  [87] rlang_0.2.1                hexbin_1.27.2             
#>  [89] manipulateWidget_0.10.0    NbClust_3.0               
#>  [91] glue_1.2.0                 metap_0.9                 
#>  [93] RColorBrewer_1.1-2         registry_0.5              
#>  [95] fpc_2.1-11                 matrixStats_0.53.1        
#>  [97] stringr_1.3.1              pkgmaker_0.27             
#>  [99] lava_1.6.1                 fields_9.6                
#> [101] DESeq2_1.18.1              recipes_0.1.3             
#> [103] SparseM_1.77               httpuv_1.4.4.1            
#> [105] class_7.3-14               BPSC_0.99.1               
#> [107] RMTstat_0.3                annotate_1.56.2           
#> [109] webshot_0.5.0              jsonlite_1.5              
#> [111] XVector_0.18.0             bit_1.1-14                
#> [113] mime_0.5                   gridExtra_2.3             
#> [115] gplots_3.0.1               Rsamtools_1.30.0          
#> [117] zingeR_0.1.0               stringi_1.2.3             
#> [119] gmodels_2.16.2             RcppRoll_0.3.0            
#> [121] gsl_1.9-10.3               bitops_1.0-6              
#> [123] maps_3.3.0                 RSQLite_2.1.1             
#> [125] tidyr_0.8.1                pheatmap_1.0.10           
#> [127] data.table_1.11.4          DEDS_1.52.0               
#> [129] energy_1.7-4               rstudioapi_0.7            
#> [131] GenomicAlignments_1.14.2   sfsmisc_1.1-2             
#> [133] nlme_3.1-137               qvalue_2.10.0             
#> [135] scran_1.6.9                fastcluster_1.1.25        
#> [137] scone_1.2.0                locfit_1.5-9.1            
#> [139] miniUI_0.1.1.1             cobs_1.3-3                
#> [141] R.oo_1.22.0                prabclus_2.2-6            
#> [143] segmented_0.5-3.0          BiocGenerics_0.24.0       
#> [145] readxl_1.1.0               dimRed_0.1.0              
#> [147] timeDate_3043.102          ROTS_1.6.0                
#> [149] cellranger_1.1.0           munsell_0.5.0             
#> [151] R.methodsS3_1.7.1          moments_0.14              
#> [153] hwriter_1.3.2              caTools_1.17.1            
#> [155] codetools_0.2-15           coda_0.19-1               
#> [157] Biobase_2.38.0             magic_1.5-8               
#> [159] GenomeInfoDb_1.14.0        diffusionMap_1.1-0        
#> [161] vipor_0.4.5                lmtest_0.9-36             
#> [163] htmlTable_1.12             rARPACK_0.11-0            
#> [165] xtable_1.8-2               SAVER_0.4.0               
#> [167] ROCR_1.0-7                 diptest_0.75-7            
#> [169] formatR_1.5                scatterplot3d_0.3-41      
#> [171] lpsymphony_1.7.1           abind_1.4-5               
#> [173] FNN_1.1                    RANN_2.5.1                
#> [175] sparsesvd_0.1-4            CompQuadForm_1.4.3        
#> [177] GenomicRanges_1.30.3       bibtex_0.4.2              
#> [179] rgl_0.99.16                tibble_1.4.2              
#> [181] ggdendro_0.1-20            cluster_2.0.7-1           
#> [183] Seurat_2.3.2               Matrix_1.2-14             
#> [185] prettyunits_1.0.2          shinyBS_0.61              
#> [187] lubridate_1.7.4            ggridges_0.5.0            
#> [189] NOISeq_2.22.1              shinydashboard_0.7.0      
#> [191] mclust_5.4                 igraph_1.2.1              
#> [193] RcppEigen_0.3.3.4.0        slam_0.1-43               
#> [195] testthat_2.0.0             doSNOW_1.0.16             
#> [197] geometry_0.3-6             htmltools_0.3.6           
#> [199] yaml_2.1.19                GenomicFeatures_1.30.3    
#> [201] XML_3.98-1.11              ModelMetrics_1.1.0        
#> [203] DrImpute_1.0               foreign_0.8-70            
#> [205] withr_2.1.2                fitdistrplus_1.0-9        
#> [207] BiocParallel_1.12.0        aroma.light_3.8.0         
#> [209] bit64_0.9-7                rngtools_1.3.1            
#> [211] doRNG_1.6.6                foreach_1.4.4             
#> [213] robustbase_0.93-0          outliers_0.14             
#> [215] scde_2.6.0                 Biostrings_2.46.0         
#> [217] combinat_0.0-8             iCOBRA_1.6.0              
#> [219] memoise_1.1.0              evaluate_0.10.1           
#> [221] VGAM_1.0-5                 nonnest2_0.5-1            
#> [223] forcats_0.3.0              rio_0.5.10                
#> [225] geneplotter_1.56.0         permute_0.9-4             
#> [227] caret_6.0-80               curl_3.2                  
#> [229] fdrtool_1.2.15             trimcluster_0.1-2         
#> [231] acepack_1.4.1              edgeR_3.20.9              
#> [233] checkmate_1.8.5            tensorA_0.36              
#> [235] DECENT_0.2.0               ellipse_0.4.1             
#> [237] rjson_0.2.20               ggplot2_2.2.1             
#> [239] openxlsx_4.1.0             ggrepel_0.8.0             
#> [241] distillery_1.0-4           ade4_1.7-11               
#> [243] dtw_1.20-1                 scDD_1.2.0                
#> [245] rprojroot_1.3-2            stabledist_0.7-1          
#> [247] Lmoments_1.2-3             tools_3.4.4               
#> [249] sandwich_2.4-0             magrittr_1.5              
#> [251] RCurl_1.95-4.10            proxy_0.4-22              
#> [253] car_3.0-0                  pbivnorm_0.6.0            
#> [255] ape_5.1                    bayesm_3.1-0.1            
#> [257] EBSeq_1.18.0               httr_1.3.1                
#> [259] assertthat_0.2.0           rmarkdown_1.10            
#> [261] boot_1.3-20                R6_2.2.2                  
#> [263] nnet_7.3-12                tximport_1.6.0            
#> [265] progress_1.2.0             genefilter_1.60.0         
#> [267] gtools_3.5.0               statmod_1.4.30            
#> [269] Rook_1.1-1                 rhdf5_2.22.0              
#> [271] splines_3.4.4              carData_3.0-1             
#> [273] colorspace_1.3-2           amap_0.8-16               
#> [275] stats4_3.4.4               NBPSeq_0.3.0              
#> [277] compositions_1.40-2        base64enc_0.1-3           
#> [279] baySeq_2.12.0              pillar_1.2.3              
#> [281] HSMMSingleCell_0.112.0     bindrcpp_0.2.2            
#> [283] GenomeInfoDbData_1.0.0     plyr_1.8.4                
#> [285] extRemes_2.0-9             dotCall64_0.9-5.2         
#> [287] gtable_0.2.0               zip_1.0.0                 
#> [289] SCnorm_1.0.0               monocle_2.6.4             
#> [291] psych_1.8.4                knitr_1.20                
#> [293] RcppArmadillo_0.8.500.0    latticeExtra_0.6-28       
#> [295] biomaRt_2.34.2             IRanges_2.12.0            
#> [297] ADGofTest_0.3              copula_0.999-18           
#> [299] Cairo_1.5-9                crosstalk_1.0.0           
#> [301] doParallel_1.0.11          pscl_1.5.2                
#> [303] flexmix_2.3-14             quantreg_5.36             
#> [305] AnnotationDbi_1.40.0       broom_0.4.4               
#> [307] scales_0.5.0               arm_1.10-1                
#> [309] backports_1.1.2            IHW_1.6.0                 
#> [311] S4Vectors_0.16.0           densityClust_0.3          
#> [313] ipred_0.9-6                lme4_1.1-17               
#> [315] brew_1.0-6                 hms_0.4.2                 
#> [317] DESeq_1.30.0               Rtsne_0.13                
#> [319] dplyr_0.7.5                shiny_1.1.0               
#> [321] ddalpha_1.3.3              grid_3.4.4                
#> [323] numDeriv_2016.8-1          bbmle_1.0.20              
#> [325] lazyeval_0.2.1             dynamicTreeCut_1.63-1     
#> [327] Formula_1.2-3              tsne_0.1-3                
#> [329] blockmodeling_0.3.1        crayon_1.3.4              
#> [331] DRR_0.0.3                  MAST_1.4.1                
#> [333] RUVSeq_1.12.0              viridis_0.5.1             
#> [335] rpart_4.1-13               zinbwave_1.0.0            
#> [337] compiler_3.4.4
```

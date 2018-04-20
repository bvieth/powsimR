
<!-- README.md is generated from README.Rmd. Please edit that file -->
powsimR: Power analysis for bulk and single cell RNA-seq experiments
====================================================================

NEWS
----

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

Some users have experienced issues installing powsimR due to vignette compilation errors. If that is the case, you can leave out building the vignette.

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
#> [1] powsimR_1.1.1     gamlss.dist_5.0-4 MASS_7.3-49      
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0             softImpute_1.4            
#>   [3] minpack.lm_1.2-1           pbapply_1.3-4             
#>   [5] haven_1.1.1                lattice_0.20-35           
#>   [7] fastICA_1.2-1              mgcv_1.8-23               
#>   [9] penalized_0.9-50           blob_1.1.1                
#>  [11] survival_2.42-3            prodlim_2018.04.18        
#>  [13] nloptr_1.0.4               DBI_0.8                   
#>  [15] R.utils_2.6.0              SingleCellExperiment_1.0.0
#>  [17] Linnorm_2.2.0              bindr_0.1.1               
#>  [19] zlibbioc_1.24.0            MatrixModels_0.4-1        
#>  [21] pspline_1.0-18             pcaMethods_1.70.0         
#>  [23] SDMTools_1.1-221           htmlwidgets_1.2           
#>  [25] mvtnorm_1.0-7              UpSetR_1.3.3              
#>  [27] tclust_1.3-1               parallel_3.4.4            
#>  [29] scater_1.6.3               irlba_2.3.2               
#>  [31] DEoptimR_1.0-8             lars_1.2                  
#>  [33] Rcpp_0.12.16               KernSmooth_2.23-15        
#>  [35] DT_0.4                     gdata_2.18.0              
#>  [37] DDRTree_0.1.5              DelayedArray_0.4.1        
#>  [39] limma_3.34.9               vegan_2.5-1               
#>  [41] CVST_0.2-1                 RcppParallel_4.4.0        
#>  [43] Hmisc_4.1-1                ShortRead_1.36.1          
#>  [45] apcluster_1.4.5            RSpectra_0.12-0           
#>  [47] msir_1.3.1                 mnormt_1.5-5              
#>  [49] ranger_0.9.0               digest_0.6.15             
#>  [51] RMySQL_0.10.14             png_0.1-7                 
#>  [53] qlcMatrix_0.9.5            cidr_0.1.5                
#>  [55] cowplot_0.9.2              glmnet_2.0-16             
#>  [57] pkgconfig_2.0.1            gower_0.1.2               
#>  [59] ggbeeswarm_0.6.0           iterators_1.0.9           
#>  [61] minqa_1.2.4                lavaan_0.5-23.1097        
#>  [63] SummarizedExperiment_1.8.1 spam_2.1-4                
#>  [65] beeswarm_0.2.3             modeltools_0.2-21         
#>  [67] RcppNumerical_0.3-2        zoo_1.8-1                 
#>  [69] tidyselect_0.2.4           clusterCrit_1.2.7         
#>  [71] ZIM_1.0.3                  reshape2_1.4.3            
#>  [73] purrr_0.2.4                kernlab_0.9-25            
#>  [75] ica_1.0-1                  pcaPP_1.9-73              
#>  [77] EDASeq_2.12.0              viridisLite_0.3.0         
#>  [79] snow_0.4-2                 rtracklayer_1.38.3        
#>  [81] rlang_0.2.0                hexbin_1.27.2             
#>  [83] manipulateWidget_0.9.0     NbClust_3.0               
#>  [85] glue_1.2.0                 metap_0.8                 
#>  [87] RColorBrewer_1.1-2         registry_0.5              
#>  [89] fpc_2.1-11                 matrixStats_0.53.1        
#>  [91] stringr_1.3.0              pkgmaker_0.22             
#>  [93] lava_1.6.1                 fields_9.6                
#>  [95] DESeq2_1.18.1              recipes_0.1.2             
#>  [97] SparseM_1.77               httpuv_1.3.6.2            
#>  [99] class_7.3-14               BPSC_0.99.1               
#> [101] RMTstat_0.3                annotate_1.56.2           
#> [103] jsonlite_1.5               XVector_0.18.0            
#> [105] bit_1.1-12                 mime_0.5                  
#> [107] gridExtra_2.3              gplots_3.0.1              
#> [109] Rsamtools_1.30.0           zingeR_0.1.0              
#> [111] stringi_1.1.7              gmodels_2.16.2            
#> [113] RcppRoll_0.2.2             gsl_1.9-10.3              
#> [115] quadprog_1.5-5             bitops_1.0-6              
#> [117] maps_3.3.0                 RSQLite_2.1.0             
#> [119] tidyr_0.8.0                pheatmap_1.0.8            
#> [121] data.table_1.10.4-3        DEDS_1.52.0               
#> [123] energy_1.7-2               rstudioapi_0.7            
#> [125] GenomicAlignments_1.14.2   sfsmisc_1.1-2             
#> [127] nlme_3.1-137               qvalue_2.10.0             
#> [129] scran_1.6.9                fastcluster_1.1.24        
#> [131] scone_1.2.0                locfit_1.5-9.1            
#> [133] miniUI_0.1.1               cobs_1.3-3                
#> [135] R.oo_1.21.0                prabclus_2.2-6            
#> [137] segmented_0.5-3.0          BiocGenerics_0.24.0       
#> [139] readxl_1.0.0               dimRed_0.1.0              
#> [141] timeDate_3043.102          ROTS_1.6.0                
#> [143] cellranger_1.1.0           munsell_0.4.3             
#> [145] R.methodsS3_1.7.1          moments_0.14              
#> [147] hwriter_1.3.2              caTools_1.17.1            
#> [149] codetools_0.2-15           coda_0.19-1               
#> [151] Biobase_2.38.0             magic_1.5-8               
#> [153] GenomeInfoDb_1.14.0        diffusionMap_1.1-0        
#> [155] vipor_0.4.5                lmtest_0.9-36             
#> [157] htmlTable_1.11.2           rARPACK_0.11-0            
#> [159] xtable_1.8-2               SAVER_0.4.0               
#> [161] ROCR_1.0-7                 diptest_0.75-7            
#> [163] formatR_1.5                scatterplot3d_0.3-41      
#> [165] lpsymphony_1.7.1           abind_1.4-5               
#> [167] FNN_1.1                    RANN_2.5.1                
#> [169] CompQuadForm_1.4.3         GenomicRanges_1.30.3      
#> [171] rgl_0.99.16                tibble_1.4.2              
#> [173] ggdendro_0.1-20            cluster_2.0.7-1           
#> [175] Seurat_2.3.0               Matrix_1.2-14             
#> [177] prettyunits_1.0.2          shinyBS_0.61              
#> [179] lubridate_1.7.4            ggridges_0.5.0            
#> [181] NOISeq_2.22.1              shinydashboard_0.7.0      
#> [183] mclust_5.4                 igraph_1.2.1              
#> [185] RcppEigen_0.3.3.4.0        slam_0.1-42               
#> [187] testthat_2.0.0             doSNOW_1.0.16             
#> [189] geometry_0.3-6             htmltools_0.3.6           
#> [191] yaml_2.1.18                GenomicFeatures_1.30.3    
#> [193] XML_3.98-1.11              ModelMetrics_1.1.0        
#> [195] DrImpute_1.0               foreign_0.8-69            
#> [197] withr_2.1.2                fitdistrplus_1.0-9        
#> [199] BiocParallel_1.12.0        aroma.light_3.8.0         
#> [201] bit64_0.9-7                rngtools_1.2.4            
#> [203] doRNG_1.6.6                foreach_1.4.4             
#> [205] robustbase_0.92-8          outliers_0.14             
#> [207] scde_2.6.0                 Biostrings_2.46.0         
#> [209] combinat_0.0-8             iCOBRA_1.6.0              
#> [211] memoise_1.1.0              evaluate_0.10.1           
#> [213] VGAM_1.0-5                 nonnest2_0.5-1            
#> [215] forcats_0.3.0              rio_0.5.10                
#> [217] geneplotter_1.56.0         permute_0.9-4             
#> [219] caret_6.0-79               curl_3.2                  
#> [221] fdrtool_1.2.15             trimcluster_0.1-2         
#> [223] acepack_1.4.1              edgeR_3.20.9              
#> [225] checkmate_1.8.5            tensorA_0.36              
#> [227] DECENT_0.2.0               ellipse_0.4.1             
#> [229] ggplot2_2.2.1              rjson_0.2.15              
#> [231] openxlsx_4.0.17            ggrepel_0.7.0             
#> [233] distillery_1.0-4           ade4_1.7-11               
#> [235] dtw_1.18-1                 scDD_1.2.0                
#> [237] rprojroot_1.3-2            stabledist_0.7-1          
#> [239] Lmoments_1.2-3             tools_3.4.4               
#> [241] sandwich_2.4-0             magrittr_1.5              
#> [243] RCurl_1.95-4.10            proxy_0.4-22              
#> [245] car_3.0-0                  pbivnorm_0.6.0            
#> [247] ape_5.1                    bayesm_3.1-0.1            
#> [249] EBSeq_1.18.0               httr_1.3.1                
#> [251] assertthat_0.2.0           rmarkdown_1.9             
#> [253] boot_1.3-20                R6_2.2.2                  
#> [255] nnet_7.3-12                progress_1.1.2            
#> [257] tximport_1.6.0             genefilter_1.60.0         
#> [259] gtools_3.5.0               statmod_1.4.30            
#> [261] Rook_1.1-1                 rhdf5_2.22.0              
#> [263] splines_3.4.4              carData_3.0-1             
#> [265] colorspace_1.3-2           amap_0.8-14               
#> [267] stats4_3.4.4               NBPSeq_0.3.0              
#> [269] compositions_1.40-1        base64enc_0.1-3           
#> [271] baySeq_2.12.0              pillar_1.2.1              
#> [273] sn_1.5-1                   HSMMSingleCell_0.112.0    
#> [275] bindrcpp_0.2.2             GenomeInfoDbData_1.0.0    
#> [277] plyr_1.8.4                 extRemes_2.0-8            
#> [279] dotCall64_0.9-5.2          gtable_0.2.0              
#> [281] SCnorm_1.0.0               monocle_2.6.4             
#> [283] psych_1.8.3.3              knitr_1.20                
#> [285] RcppArmadillo_0.8.400.0.0  latticeExtra_0.6-28       
#> [287] biomaRt_2.34.2             IRanges_2.12.0            
#> [289] ADGofTest_0.3              copula_0.999-18           
#> [291] crosstalk_1.0.0            Cairo_1.5-9               
#> [293] doParallel_1.0.11          pscl_1.5.2                
#> [295] flexmix_2.3-14             quantreg_5.35             
#> [297] AnnotationDbi_1.40.0       broom_0.4.4               
#> [299] scales_0.5.0               arm_1.10-1                
#> [301] backports_1.1.2            IHW_1.6.0                 
#> [303] S4Vectors_0.16.0           densityClust_0.3          
#> [305] ipred_0.9-6                lme4_1.1-17               
#> [307] brew_1.0-6                 DESeq_1.30.0              
#> [309] Rtsne_0.13                 dplyr_0.7.4               
#> [311] shiny_1.0.5                ddalpha_1.3.2             
#> [313] grid_3.4.4                 numDeriv_2016.8-1         
#> [315] bbmle_1.0.20               lazyeval_0.2.1            
#> [317] dynamicTreeCut_1.63-1      Formula_1.2-2             
#> [319] tsne_0.1-3                 blockmodeling_0.3.0       
#> [321] DRR_0.0.3                  MAST_1.4.1                
#> [323] RUVSeq_1.12.0              viridis_0.5.1             
#> [325] rpart_4.1-13               zinbwave_1.0.0            
#> [327] compiler_3.4.4
```

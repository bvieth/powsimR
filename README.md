
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
python implementation beforre installing
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
            if (strsplit(version[["version.string"]], " ")[[1]][3] > "3.5.0") {
                if (!requireNamespace("BiocManager")) {
                  install.packages("BiocManager")
                }
                BiocManager::install(new.pkg, dependencies = TRUE, ask = FALSE)
            }
            if (strsplit(version[["version.string"]], " ")[[1]][3] < "3.5.0") {
                source("https://bioconductor.org/biocLite.R")
                biocLite(new.pkg, dependencies = TRUE, ask = FALSE)
            }
        }
        if (repository == "github") {
            devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE, 
                dependencies = TRUE)
        }
    }
}

# CRAN PACKAGES
cranpackages <- c("broom", "cobs", "cowplot", "data.table", "devtools", "doParallel", 
    "dplyr", "drc", "DrImpute", "fastICA", "fitdistrplus", "foreach", "gamlss.dist", 
    "ggExtra", "ggplot2", "ggthemes", "grDevices", "glmnet", "grid", "gtools", 
    "Hmisc", "kernlab", "MASS", "MBESS", "matrixStats", "mclust", "methods", 
    "minpack.lm", "moments", "msir", "NBPSeq", "nonnest2", "parallel", "penalized", 
    "plyr", "pscl", "reshape2", "Rmagic", "rsvd", "Rtsne", "scales", "Seurat", 
    "snow", "stats", "tibble", "tidyr", "VGAM", "ZIM")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR
biocpackages <- c("AnnotationDbi", "bayNorm", "baySeq", "Biobase", "BiocGenerics", 
    "BiocParallel", "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", 
    "Linnorm", "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq", "S4Vectors", 
    "scater", "scDD", "scde", "scone", "scran", "SCnorm", "SingleCellExperiment", 
    "SummarizedExperiment", "zinbwave")
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
vignette and read it on my Github Page of
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
#> Registered S3 method overwritten by 'R.oo':
#>   method        from       
#>   throw.default R.methodsS3
#> Warning: replacing previous import 'parallel::makeCluster' by
#> 'snow::makeCluster' when loading 'powsimR'
#> Warning: replacing previous import 'parallel::stopCluster' by
#> 'snow::stopCluster' when loading 'powsimR'
#> Warning: replacing previous import 'penalized::predict' by 'stats::predict'
#> when loading 'powsimR'
#> Warning: replacing previous import 'zinbwave::glmWeightedF' by
#> 'zingeR::glmWeightedF' when loading 'powsimR'
sessionInfo()
#> R version 3.6.0 (2019-04-26)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
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
#> [1] powsimR_1.1.5     gamlss.dist_5.1-4 MASS_7.3-51.4    
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            lattice_0.20-38            
#>   [5] vctrs_0.2.0                 fastICA_1.2-2              
#>   [7] mgcv_1.8-28                 penalized_0.9-51           
#>   [9] blob_1.2.0                  survival_2.44-1.1          
#>  [11] Rmagic_1.4.0                later_0.8.0                
#>  [13] nloptr_1.2.1                DBI_1.0.0                  
#>  [15] R.utils_2.9.0               SingleCellExperiment_1.6.0 
#>  [17] Linnorm_2.8.0               dqrng_0.2.1                
#>  [19] zlibbioc_1.30.0             MatrixModels_0.4-1         
#>  [21] pspline_1.0-18              htmlwidgets_1.3            
#>  [23] mvtnorm_1.0-11              future_1.14.0              
#>  [25] UpSetR_1.4.0                parallel_3.6.0             
#>  [27] scater_1.12.2               irlba_2.3.3                
#>  [29] DEoptimR_1.0-8              lars_1.2                   
#>  [31] Rcpp_1.0.2                  KernSmooth_2.23-15         
#>  [33] DT_0.7                      promises_1.0.1             
#>  [35] gdata_2.18.0                DDRTree_0.1.5              
#>  [37] DelayedArray_0.10.0         limma_3.40.6               
#>  [39] vegan_2.5-5                 Hmisc_4.2-0                
#>  [41] ShortRead_1.42.0            apcluster_1.4.7            
#>  [43] RSpectra_0.15-0             msir_1.3.2                 
#>  [45] mnormt_1.5-5                digest_0.6.20              
#>  [47] qlcMatrix_0.9.7             sctransform_0.2.0          
#>  [49] cowplot_1.0.0               glmnet_2.0-18              
#>  [51] pkgconfig_2.0.2             docopt_0.6.1               
#>  [53] DelayedMatrixStats_1.6.0    ggbeeswarm_0.6.0           
#>  [55] iterators_1.0.12            minqa_1.2.4                
#>  [57] lavaan_0.6-4                reticulate_1.13            
#>  [59] SummarizedExperiment_1.14.0 spam_2.2-2                 
#>  [61] beeswarm_0.2.3              modeltools_0.2-22          
#>  [63] xfun_0.8                    zoo_1.8-6                  
#>  [65] tidyselect_0.2.5            ZIM_1.1.0                  
#>  [67] reshape2_1.4.3              purrr_0.3.2                
#>  [69] kernlab_0.9-27              pcaPP_1.9-73               
#>  [71] EDASeq_2.18.0               viridisLite_0.3.0          
#>  [73] snow_0.4-3                  rtracklayer_1.44.2         
#>  [75] rlang_0.4.0                 hexbin_1.27.3              
#>  [77] glue_1.3.1                  RColorBrewer_1.1-2         
#>  [79] registry_0.5-1              fpc_2.2-3                  
#>  [81] matrixStats_0.54.0          stringr_1.4.0              
#>  [83] pkgmaker_0.27               fields_9.8-3               
#>  [85] ggsignif_0.5.0              DESeq2_1.24.0              
#>  [87] SparseM_1.77                httpuv_1.5.1               
#>  [89] class_7.3-15                BPSC_0.99.2                
#>  [91] BiocNeighbors_1.2.0         annotate_1.62.0            
#>  [93] jsonlite_1.6                XVector_0.24.0             
#>  [95] bit_1.1-14                  mime_0.7                   
#>  [97] gridExtra_2.3               gplots_3.0.1.1             
#>  [99] Rsamtools_2.0.0             zingeR_0.1.0               
#> [101] stringi_1.4.3               gmodels_2.18.1             
#> [103] gsl_2.1-6                   bitops_1.0-6               
#> [105] maps_3.3.0                  RSQLite_2.1.2              
#> [107] tidyr_0.8.3                 pheatmap_1.0.12            
#> [109] data.table_1.12.2           DEDS_1.58.0                
#> [111] energy_1.7-6                rstudioapi_0.10            
#> [113] GenomicAlignments_1.20.1    nlme_3.1-140               
#> [115] qvalue_2.16.0               scran_1.12.1               
#> [117] fastcluster_1.1.25          locfit_1.5-9.1             
#> [119] scone_1.8.0                 listenv_0.7.0              
#> [121] cobs_1.3-3                  R.oo_1.22.0                
#> [123] prabclus_2.3-1              segmented_1.0-0            
#> [125] BiocGenerics_0.30.0         ROTS_1.12.0                
#> [127] munsell_0.5.0               R.methodsS3_1.7.1          
#> [129] moments_0.14                hwriter_1.3.2              
#> [131] caTools_1.17.1.2            codetools_0.2-16           
#> [133] coda_0.19-3                 Biobase_2.44.0             
#> [135] GenomeInfoDb_1.20.0         vipor_0.4.5                
#> [137] htmlTable_1.13.1            bayNorm_1.2.0              
#> [139] lsei_1.2-0                  rARPACK_0.11-0             
#> [141] xtable_1.8-4                SAVER_1.1.1                
#> [143] ROCR_1.0-7                  diptest_0.75-7             
#> [145] formatR_1.7                 lpsymphony_1.12.0          
#> [147] abind_1.4-5                 FNN_1.1.3                  
#> [149] RANN_2.6.1                  sparsesvd_0.2              
#> [151] CompQuadForm_1.4.3          GenomicRanges_1.36.0       
#> [153] bibtex_0.4.2                tibble_2.1.3               
#> [155] ggdendro_0.1-20             cluster_2.1.0              
#> [157] future.apply_1.3.0          zeallot_0.1.0              
#> [159] Matrix_1.2-17               prettyunits_1.0.2          
#> [161] shinyBS_0.61                NOISeq_2.28.0              
#> [163] shinydashboard_0.7.1        mclust_5.4.5               
#> [165] igraph_1.2.4.1              slam_0.1-45                
#> [167] testthat_2.2.1              doSNOW_1.0.18              
#> [169] htmltools_0.3.6             yaml_2.2.0                 
#> [171] GenomicFeatures_1.36.4      XML_3.98-1.20              
#> [173] ggpubr_0.2.1                DrImpute_1.0               
#> [175] foreign_0.8-71              withr_2.1.2                
#> [177] fitdistrplus_1.0-14         BiocParallel_1.18.0        
#> [179] aroma.light_3.14.0          bit64_0.9-7                
#> [181] rngtools_1.4                doRNG_1.7.1                
#> [183] foreach_1.4.7               robustbase_0.93-5          
#> [185] outliers_0.14               Biostrings_2.52.0          
#> [187] combinat_0.0-8              rsvd_1.0.2                 
#> [189] iCOBRA_1.12.1               memoise_1.1.0              
#> [191] evaluate_0.14               VGAM_1.1-1                 
#> [193] nonnest2_0.5-2              geneplotter_1.62.0         
#> [195] permute_0.9-5               fdrtool_1.2.15             
#> [197] acepack_1.4.1               edgeR_3.26.5               
#> [199] checkmate_1.9.4             npsurv_0.4-0               
#> [201] truncnorm_1.0-8             DECENT_1.1.0               
#> [203] tensorA_0.36.1              ellipse_0.4.1              
#> [205] ggplot2_3.2.0               ggrepel_0.8.1              
#> [207] scDD_1.8.0                  tools_3.6.0                
#> [209] stabledist_0.7-1            sandwich_2.5-1             
#> [211] magrittr_1.5                RCurl_1.95-4.12            
#> [213] pbivnorm_0.6.0              bayesm_3.1-2               
#> [215] EBSeq_1.24.0                httr_1.4.0                 
#> [217] assertthat_0.2.1            rmarkdown_1.14             
#> [219] boot_1.3-23                 globals_0.12.4             
#> [221] R6_2.4.0                    Rhdf5lib_1.6.0             
#> [223] nnet_7.3-12                 progress_1.2.2             
#> [225] genefilter_1.66.0           gtools_3.8.1               
#> [227] statmod_1.4.32              BiocSingular_1.0.0         
#> [229] rhdf5_2.28.0                splines_3.6.0              
#> [231] colorspace_1.4-1            amap_0.8-17                
#> [233] generics_0.0.2              stats4_3.6.0               
#> [235] NBPSeq_0.3.0                base64enc_0.1-3            
#> [237] compositions_1.40-2         baySeq_2.18.0              
#> [239] pillar_1.4.2                HSMMSingleCell_1.4.0       
#> [241] GenomeInfoDbData_1.2.1      plyr_1.8.4                 
#> [243] dotCall64_1.0-0             gtable_0.3.0               
#> [245] SCnorm_1.6.0                monocle_2.12.0             
#> [247] knitr_1.23                  RcppArmadillo_0.9.600.4.0  
#> [249] latticeExtra_0.6-28         biomaRt_2.40.3             
#> [251] IRanges_2.18.1              ADGofTest_0.3              
#> [253] copula_0.999-19.1           doParallel_1.0.14          
#> [255] pscl_1.5.2                  flexmix_2.3-15             
#> [257] quantreg_5.42.1             AnnotationDbi_1.46.0       
#> [259] broom_0.5.2                 scales_1.0.0               
#> [261] arm_1.10-1                  backports_1.1.4            
#> [263] IHW_1.12.0                  S4Vectors_0.22.0           
#> [265] densityClust_0.3            lme4_1.1-21                
#> [267] blme_1.0-4                  hms_0.5.0                  
#> [269] DESeq_1.36.0                Rtsne_0.15                 
#> [271] dplyr_0.8.3                 shiny_1.3.2                
#> [273] grid_3.6.0                  numDeriv_2016.8-1.1        
#> [275] bbmle_1.0.20                lazyeval_0.2.2             
#> [277] dynamicTreeCut_1.63-1       Formula_1.2-3              
#> [279] blockmodeling_0.3.4         crayon_1.3.4               
#> [281] MAST_1.10.0                 RUVSeq_1.18.0              
#> [283] viridis_0.5.1               rpart_4.1-15               
#> [285] compiler_3.6.0              zinbwave_1.6.0
```

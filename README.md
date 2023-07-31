
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `powsimR` <br/> Power analysis for bulk and <br/> single cell RNA-seq experiments <img src="vignettes/powsimR.png" align="right" width="200" />

Please also consult my Github Page of
[powsimR](https://bvieth.github.io/powsimR/) made with
[pkgdown](http://pkgdown.r-lib.org/index.html)!

## :arrow_double_down: Installation Guide

For the installation, the R package `devtools` is needed.

``` r
install.packages("devtools")
library(devtools)
```

I recommend to install first the dependencies manually and then powsimR.
If you plan to use MAGIC for imputation, then please follow their
[instruction](https://github.com/KrishnaswamyLab/MAGIC) to install the
python implementation before installing powsimR.

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
            devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE,
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
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
library("powsimR")
```

Alternative, you can try to install powsimR and its dependencies
directly using devtools:

``` r
devtools::install_github("bvieth/powsimR")
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
#> Registered S3 method overwritten by 'gdata':
#>   method         from  
#>   reorder.factor gplots
#> Warning: replacing previous import 'DECENT::lrTest' by 'MAST::lrTest' when
#> loading 'powsimR'
#> Warning: replacing previous import 'penalized::predict' by 'stats::predict' when
#> loading 'powsimR'
#> Warning: replacing previous import 'zinbwave::glmWeightedF' by
#> 'zingeR::glmWeightedF' when loading 'powsimR'
sessionInfo()
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.6 LTS
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
#> [1] powsimR_1.2.3     gamlss.dist_6.0-1 MASS_7.3-54      
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.2.0              softImpute_1.4-1           
#>   [3] minpack.lm_1.2-1            lattice_0.20-45            
#>   [5] vctrs_0.3.8                 fastICA_1.2-3              
#>   [7] mgcv_1.8-38                 penalized_0.9-51           
#>   [9] blob_1.2.2                  survival_3.2-13            
#>  [11] prodlim_2019.11.13          Rmagic_2.0.3               
#>  [13] later_1.3.0                 nloptr_1.2.2.3             
#>  [15] DBI_1.1.1                   R.utils_2.11.0             
#>  [17] rappdirs_0.3.3              SingleCellExperiment_1.16.0
#>  [19] Linnorm_2.18.0              dqrng_0.3.0                
#>  [21] jpeg_0.1-9                  zlibbioc_1.40.0            
#>  [23] MatrixModels_0.5-0          htmlwidgets_1.5.4          
#>  [25] mvtnorm_1.1-3               future_1.23.0              
#>  [27] UpSetR_1.4.0                parallel_4.1.2             
#>  [29] scater_1.22.0               irlba_2.3.3                
#>  [31] DEoptimR_1.0-9              Rcpp_1.0.7                 
#>  [33] KernSmooth_2.23-20          DT_0.20                    
#>  [35] promises_1.2.0.1            gdata_2.18.0               
#>  [37] DDRTree_0.1.5               DelayedArray_0.20.0        
#>  [39] limma_3.50.0                vegan_2.5-7                
#>  [41] Hmisc_4.6-0                 ShortRead_1.52.0           
#>  [43] apcluster_1.4.8             RSpectra_0.16-0            
#>  [45] msir_1.3.3                  mnormt_2.0.2               
#>  [47] digest_0.6.28               png_0.1-7                  
#>  [49] bluster_1.4.0               qlcMatrix_0.9.7            
#>  [51] sctransform_0.3.2           cowplot_1.1.1              
#>  [53] pkgconfig_2.0.3             docopt_0.7.1               
#>  [55] DelayedMatrixStats_1.16.0   gower_0.2.2                
#>  [57] ggbeeswarm_0.6.0            iterators_1.0.13           
#>  [59] minqa_1.2.4                 lavaan_0.6-9               
#>  [61] reticulate_1.22             SummarizedExperiment_1.24.0
#>  [63] spam_2.7-0                  beeswarm_0.4.0             
#>  [65] modeltools_0.2-23           xfun_0.28                  
#>  [67] zoo_1.8-9                   tidyselect_1.1.1           
#>  [69] ZIM_1.1.0                   reshape2_1.4.4             
#>  [71] purrr_0.3.4                 kernlab_0.9-29             
#>  [73] EDASeq_2.28.0               viridisLite_0.4.0          
#>  [75] snow_0.4-4                  rtracklayer_1.54.0         
#>  [77] rlang_0.4.12                hexbin_1.28.2              
#>  [79] glue_1.5.0                  RColorBrewer_1.1-2         
#>  [81] fpc_2.2-9                   matrixStats_0.61.0         
#>  [83] MatrixGenerics_1.6.0        stringr_1.4.0              
#>  [85] lava_1.6.10                 fields_13.3                
#>  [87] ggsignif_0.6.3              DESeq2_1.34.0              
#>  [89] recipes_0.1.17              SparseM_1.81               
#>  [91] httpuv_1.6.3                class_7.3-19               
#>  [93] BPSC_0.99.2                 BiocNeighbors_1.12.0       
#>  [95] annotate_1.72.0             jsonlite_1.7.2             
#>  [97] XVector_0.34.0              tmvnsim_1.0-2              
#>  [99] bit_4.0.4                   mime_0.12                  
#> [101] gridExtra_2.3               gplots_3.1.1               
#> [103] Rsamtools_2.10.0            zingeR_0.1.0               
#> [105] stringi_1.7.5               gmodels_2.18.1             
#> [107] rhdf5filters_1.6.0          bitops_1.0-7               
#> [109] maps_3.4.0                  RSQLite_2.2.8              
#> [111] tidyr_1.1.4                 pheatmap_1.0.12            
#> [113] data.table_1.14.2           rstudioapi_0.13            
#> [115] GenomicAlignments_1.30.0    nlme_3.1-153               
#> [117] qvalue_2.26.0               scran_1.22.1               
#> [119] fastcluster_1.2.3           locfit_1.5-9.4             
#> [121] scone_1.18.0                listenv_0.8.0              
#> [123] cobs_1.3-4                  R.oo_1.24.0                
#> [125] prabclus_2.3-2              segmented_1.3-4            
#> [127] dbplyr_2.1.1                BiocGenerics_0.40.0        
#> [129] lifecycle_1.0.1             timeDate_3043.102          
#> [131] ROTS_1.22.0                 munsell_0.5.0              
#> [133] hwriter_1.3.2               R.methodsS3_1.8.1          
#> [135] moments_0.14                caTools_1.18.2             
#> [137] codetools_0.2-18            coda_0.19-4                
#> [139] Biobase_2.54.0              GenomeInfoDb_1.30.0        
#> [141] vipor_0.4.5                 htmlTable_2.3.0            
#> [143] bayNorm_1.12.0              rARPACK_0.11-0             
#> [145] xtable_1.8-4                SAVER_1.1.2                
#> [147] ROCR_1.0-11                 diptest_0.76-0             
#> [149] formatR_1.11                lpsymphony_1.22.0          
#> [151] abind_1.4-5                 FNN_1.1.3                  
#> [153] parallelly_1.29.0           RANN_2.6.1                 
#> [155] sparsesvd_0.2               CompQuadForm_1.4.3         
#> [157] BiocIO_1.4.0                GenomicRanges_1.46.1       
#> [159] tibble_3.1.6                ggdendro_0.1.22            
#> [161] cluster_2.1.2               future.apply_1.8.1         
#> [163] Matrix_1.3-4                ellipsis_0.3.2             
#> [165] prettyunits_1.1.1           shinyBS_0.61               
#> [167] lubridate_1.8.0             NOISeq_2.38.0              
#> [169] shinydashboard_0.7.2        mclust_5.4.8               
#> [171] igraph_1.2.9                ggstance_0.3.5             
#> [173] slam_0.1-49                 testthat_3.1.0             
#> [175] doSNOW_1.0.19               htmltools_0.5.2            
#> [177] BiocFileCache_2.2.0         GenomicFeatures_1.46.1     
#> [179] yaml_2.2.1                  utf8_1.2.2                 
#> [181] XML_3.99-0.8                ModelMetrics_1.2.2.2       
#> [183] ggpubr_0.4.0                DrImpute_1.0               
#> [185] foreign_0.8-81              withr_2.4.2                
#> [187] scuttle_1.4.0               fitdistrplus_1.1-6         
#> [189] BiocParallel_1.28.2         aroma.light_3.24.0         
#> [191] bit64_4.0.5                 foreach_1.5.1              
#> [193] robustbase_0.93-9           outliers_0.14              
#> [195] Biostrings_2.62.0           combinat_0.0-8             
#> [197] rsvd_1.0.5                  ScaledMatrix_1.2.0         
#> [199] iCOBRA_1.22.1               memoise_2.0.1              
#> [201] evaluate_0.14               VGAM_1.1-5                 
#> [203] nonnest2_0.5-5              geneplotter_1.72.0         
#> [205] permute_0.9-5               caret_6.0-90               
#> [207] curl_4.3.2                  fdrtool_1.2.17             
#> [209] fansi_0.5.0                 conquer_1.2.1              
#> [211] edgeR_3.36.0                checkmate_2.0.0            
#> [213] cachem_1.0.6                truncnorm_1.0-8            
#> [215] tensorA_0.36.2              DECENT_1.1.0               
#> [217] ellipse_0.4.2               rjson_0.2.20               
#> [219] metapod_1.2.0               ggplot2_3.3.5              
#> [221] rstatix_0.7.0               ggrepel_0.9.1              
#> [223] scDD_1.18.0                 tools_4.1.2                
#> [225] sandwich_3.0-1              magrittr_2.0.1             
#> [227] RCurl_1.98-1.5              car_3.0-12                 
#> [229] pbivnorm_0.6.0              bayesm_3.1-4               
#> [231] xml2_1.3.2                  EBSeq_1.34.0               
#> [233] httr_1.4.2                  assertthat_0.2.1           
#> [235] rmarkdown_2.11              Rhdf5lib_1.16.0            
#> [237] boot_1.3-28                 globals_0.14.0             
#> [239] R6_2.5.1                    nnet_7.3-16                
#> [241] progress_1.2.2              genefilter_1.76.0          
#> [243] KEGGREST_1.34.0             gtools_3.9.2               
#> [245] statmod_1.4.36              beachmat_2.10.0            
#> [247] BiocSingular_1.10.0         rhdf5_2.38.0               
#> [249] splines_4.1.2               carData_3.0-4              
#> [251] colorspace_2.0-2            amap_0.8-18                
#> [253] generics_0.1.1              stats4_4.1.2               
#> [255] NBPSeq_0.3.0                compositions_2.0-2         
#> [257] base64enc_0.1-3             baySeq_2.28.0              
#> [259] pillar_1.6.4                HSMMSingleCell_1.14.0      
#> [261] GenomeInfoDbData_1.2.7      plyr_1.8.6                 
#> [263] dotCall64_1.0-1             gtable_0.3.0               
#> [265] SCnorm_1.16.0               monocle_2.22.0             
#> [267] restfulr_0.0.13             knitr_1.36                 
#> [269] RcppArmadillo_0.10.7.3.0    latticeExtra_0.6-29        
#> [271] biomaRt_2.50.1              IRanges_2.28.0             
#> [273] fastmap_1.1.0               doParallel_1.0.16          
#> [275] pscl_1.5.5                  flexmix_2.3-17             
#> [277] quantreg_5.86               AnnotationDbi_1.56.2       
#> [279] broom_0.7.10                filelock_1.0.2             
#> [281] scales_1.1.1                arm_1.12-2                 
#> [283] backports_1.4.0             plotrix_3.8-2              
#> [285] IHW_1.22.0                  S4Vectors_0.32.3           
#> [287] densityClust_0.3            ipred_0.9-12               
#> [289] lme4_1.1-27.1               hms_1.1.1                  
#> [291] Rtsne_0.15                  dplyr_1.0.7                
#> [293] shiny_1.7.1                 grid_4.1.2                 
#> [295] Formula_1.2-4               blockmodeling_1.0.5        
#> [297] crayon_1.4.2                MAST_1.20.0                
#> [299] RUVSeq_1.28.0               pROC_1.18.0                
#> [301] sparseMatrixStats_1.6.0     viridis_0.6.2              
#> [303] rpart_4.1-15                zinbwave_1.16.0            
#> [305] compiler_4.1.2
```

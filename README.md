
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `powsimR` <br/> Power analysis for bulk and <br/> single cell RNA-seq experiments <img src="vignettes/powsimR.png" align="right" width="200" />

Please also consult my Github Page of
[powsimR](https://bvieth.github.io/powsimR/) made with
[pkgdown](http://pkgdown.r-lib.org/index.html)\!

## :arrow\_double\_down: Installation Guide

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

## :book: User Guide

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

## :scroll: Citation

Please use the following entry for citing powsimR.

``` r
citation("powsimR")
```

powsimR is published in
[Bioinformatics](https://doi.org/10.1093/bioinformatics/btx435). A
preprint paper is also on [bioRxiv](https://doi.org/10.1101/117150).

## :incoming\_envelope: Notes

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
#> Warning: replacing previous import 'DECENT::lrTest' by 'MAST::lrTest' when
#> loading 'powsimR'
#> Registered S3 methods overwritten by 'lme4':
#>   method                          from
#>   cooks.distance.influence.merMod car 
#>   influence.merMod                car 
#>   dfbeta.influence.merMod         car 
#>   dfbetas.influence.merMod        car
#> Warning: replacing previous import 'penalized::predict' by 'stats::predict' when
#> loading 'powsimR'
#> Warning: replacing previous import 'zinbwave::glmWeightedF' by
#> 'zingeR::glmWeightedF' when loading 'powsimR'
sessionInfo()
#> R version 4.0.1 (2020-06-06)
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
#> [1] powsimR_1.2.2     gamlss.dist_5.1-6 MASS_7.3-51.6    
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.2.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            lattice_0.20-41            
#>   [5] haven_2.3.1                 vctrs_0.3.1                
#>   [7] fastICA_1.2-2               mgcv_1.8-31                
#>   [9] penalized_0.9-51            blob_1.2.1                 
#>  [11] survival_3.2-3              Rmagic_2.0.3               
#>  [13] later_1.1.0.1               nloptr_1.2.2.1             
#>  [15] DBI_1.1.0                   R.utils_2.9.2              
#>  [17] SingleCellExperiment_1.10.1 rappdirs_0.3.1             
#>  [19] Linnorm_2.12.0              dqrng_0.2.1                
#>  [21] jpeg_0.1-8.1                zlibbioc_1.34.0            
#>  [23] MatrixModels_0.4-1          htmlwidgets_1.5.1          
#>  [25] mvtnorm_1.1-1               future_1.17.0              
#>  [27] UpSetR_1.4.0                parallel_4.0.1             
#>  [29] scater_1.16.1               irlba_2.3.3                
#>  [31] DEoptimR_1.0-8              Rcpp_1.0.4.6               
#>  [33] KernSmooth_2.23-17          DT_0.13                    
#>  [35] promises_1.1.1              gdata_2.18.0               
#>  [37] DDRTree_0.1.5               DelayedArray_0.14.0        
#>  [39] limma_3.44.3                vegan_2.5-6                
#>  [41] Hmisc_4.4-0                 ShortRead_1.46.0           
#>  [43] apcluster_1.4.8             RSpectra_0.16-0            
#>  [45] msir_1.3.2                  mnormt_2.0.0               
#>  [47] digest_0.6.25               png_0.1-7                  
#>  [49] qlcMatrix_0.9.7             sctransform_0.2.1          
#>  [51] cowplot_1.0.0               pkgconfig_2.0.3            
#>  [53] docopt_0.6.1                DelayedMatrixStats_1.10.0  
#>  [55] ggbeeswarm_0.6.0            iterators_1.0.12           
#>  [57] minqa_1.2.4                 lavaan_0.6-6               
#>  [59] reticulate_1.16             SummarizedExperiment_1.18.1
#>  [61] spam_2.5-1                  beeswarm_0.2.3             
#>  [63] modeltools_0.2-23           xfun_0.14                  
#>  [65] zoo_1.8-8                   tidyselect_1.1.0           
#>  [67] ZIM_1.1.0                   reshape2_1.4.4             
#>  [69] purrr_0.3.4                 kernlab_0.9-29             
#>  [71] EDASeq_2.22.0               viridisLite_0.3.0          
#>  [73] snow_0.4-3                  rtracklayer_1.48.0         
#>  [75] rlang_0.4.6                 hexbin_1.28.1              
#>  [77] glue_1.4.1                  RColorBrewer_1.1-2         
#>  [79] fpc_2.2-5                   matrixStats_0.56.0         
#>  [81] stringr_1.4.0               fields_10.3                
#>  [83] ggsignif_0.6.0              DESeq2_1.28.1              
#>  [85] SparseM_1.78                httpuv_1.5.4               
#>  [87] class_7.3-17                BPSC_0.99.2                
#>  [89] BiocNeighbors_1.6.0         annotate_1.66.0            
#>  [91] jsonlite_1.6.1              XVector_0.28.0             
#>  [93] tmvnsim_1.0-2               bit_1.1-15.2               
#>  [95] mime_0.9                    gridExtra_2.3              
#>  [97] gplots_3.0.3                Rsamtools_2.4.0            
#>  [99] zingeR_0.1.0                stringi_1.4.6              
#> [101] gmodels_2.18.1              bitops_1.0-6               
#> [103] maps_3.3.0                  RSQLite_2.2.0              
#> [105] tidyr_1.1.0                 pheatmap_1.0.12            
#> [107] data.table_1.12.8           rstudioapi_0.11            
#> [109] GenomicAlignments_1.24.0    nlme_3.1-148               
#> [111] qvalue_2.20.0               scran_1.16.0               
#> [113] fastcluster_1.1.25          locfit_1.5-9.4             
#> [115] scone_1.12.0                listenv_0.8.0              
#> [117] cobs_1.3-4                  R.oo_1.23.0                
#> [119] prabclus_2.3-2              dbplyr_1.4.4               
#> [121] segmented_1.1-0             BiocGenerics_0.34.0        
#> [123] readxl_1.3.1                lifecycle_0.2.0            
#> [125] ROTS_1.16.0                 munsell_0.5.0              
#> [127] cellranger_1.1.0            R.methodsS3_1.8.0          
#> [129] moments_0.14                hwriter_1.3.2              
#> [131] caTools_1.18.0              codetools_0.2-16           
#> [133] coda_0.19-3                 Biobase_2.48.0             
#> [135] GenomeInfoDb_1.24.2         vipor_0.4.5                
#> [137] htmlTable_1.13.3            bayNorm_1.6.0              
#> [139] rARPACK_0.11-0              xtable_1.8-4               
#> [141] SAVER_1.1.2                 ROCR_1.0-11                
#> [143] diptest_0.75-7              formatR_1.7                
#> [145] lpsymphony_1.16.0           abind_1.4-5                
#> [147] FNN_1.1.3                   RANN_2.6.1                 
#> [149] askpass_1.1                 sparsesvd_0.2              
#> [151] CompQuadForm_1.4.3          GenomicRanges_1.40.0       
#> [153] tibble_3.0.1                ggdendro_0.1-20            
#> [155] cluster_2.1.0               future.apply_1.5.0         
#> [157] Matrix_1.2-18               ellipsis_0.3.1             
#> [159] prettyunits_1.1.1           shinyBS_0.61               
#> [161] NOISeq_2.31.0               shinydashboard_0.7.1       
#> [163] mclust_5.4.6                igraph_1.2.5               
#> [165] ggstance_0.3.4              slam_0.1-47                
#> [167] testthat_2.3.2              doSNOW_1.0.18              
#> [169] htmltools_0.5.0             BiocFileCache_1.12.0       
#> [171] yaml_2.2.1                  GenomicFeatures_1.40.0     
#> [173] XML_3.99-0.3                ggpubr_0.3.0               
#> [175] DrImpute_1.0                foreign_0.8-80             
#> [177] fitdistrplus_1.1-1          BiocParallel_1.22.0        
#> [179] aroma.light_3.18.0          bit64_0.9-7                
#> [181] foreach_1.5.0               robustbase_0.93-6          
#> [183] outliers_0.14               Biostrings_2.56.0          
#> [185] combinat_0.0-8              rsvd_1.0.3                 
#> [187] iCOBRA_1.16.0               memoise_1.1.0              
#> [189] evaluate_0.14               VGAM_1.1-3                 
#> [191] nonnest2_0.5-3              forcats_0.5.0              
#> [193] rio_0.5.16                  geneplotter_1.66.0         
#> [195] permute_0.9-5               curl_4.3                   
#> [197] fdrtool_1.2.15              acepack_1.4.1              
#> [199] edgeR_3.30.3                checkmate_2.0.0            
#> [201] truncnorm_1.0-8             DECENT_1.1.0               
#> [203] tensorA_0.36.1              ellipse_0.4.2              
#> [205] ggplot2_3.3.1               openxlsx_4.1.5             
#> [207] rstatix_0.5.0               ggrepel_0.8.2              
#> [209] scDD_1.12.0                 tools_4.0.1                
#> [211] sandwich_2.5-1              magrittr_1.5               
#> [213] RCurl_1.98-1.2              car_3.0-8                  
#> [215] pbivnorm_0.6.0              bayesm_3.1-4               
#> [217] EBSeq_1.28.0                httr_1.4.1                 
#> [219] assertthat_0.2.1            rmarkdown_2.2              
#> [221] boot_1.3-25                 globals_0.12.5             
#> [223] R6_2.4.1                    Rhdf5lib_1.10.0            
#> [225] nnet_7.3-14                 progress_1.2.2             
#> [227] genefilter_1.70.0           gtools_3.8.2               
#> [229] statmod_1.4.34              BiocSingular_1.4.0         
#> [231] rhdf5_2.32.0                splines_4.0.1              
#> [233] carData_3.0-4               colorspace_1.4-1           
#> [235] amap_0.8-18                 generics_0.0.2             
#> [237] stats4_4.0.1                NBPSeq_0.3.0               
#> [239] base64enc_0.1-3             compositions_1.40-5        
#> [241] baySeq_2.22.0               pillar_1.4.4               
#> [243] HSMMSingleCell_1.8.0        GenomeInfoDbData_1.2.3     
#> [245] plyr_1.8.6                  dotCall64_1.0-0            
#> [247] gtable_0.3.0                zip_2.0.4                  
#> [249] SCnorm_1.10.0               monocle_2.16.0             
#> [251] knitr_1.28                  RcppArmadillo_0.9.900.1.0  
#> [253] latticeExtra_0.6-29         biomaRt_2.44.0             
#> [255] IRanges_2.22.2              fastmap_1.0.1              
#> [257] doParallel_1.0.15           pscl_1.5.5                 
#> [259] flexmix_2.3-15              quantreg_5.55              
#> [261] AnnotationDbi_1.50.0        broom_0.5.6                
#> [263] openssl_1.4.1               scales_1.1.1               
#> [265] arm_1.11-1                  backports_1.1.8            
#> [267] plotrix_3.7-8               IHW_1.16.0                 
#> [269] S4Vectors_0.26.1            densityClust_0.3           
#> [271] lme4_1.1-23                 hms_0.5.3                  
#> [273] DESeq_1.39.0                Rtsne_0.15                 
#> [275] dplyr_1.0.0                 shiny_1.4.0.2              
#> [277] grid_4.0.1                  Formula_1.2-3              
#> [279] blockmodeling_0.3.6         crayon_1.3.4               
#> [281] MAST_1.14.0                 RUVSeq_1.22.0              
#> [283] viridis_0.5.1               rpart_4.1-15               
#> [285] zinbwave_1.10.0             compiler_4.0.1
```

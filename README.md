
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
#> Warning: replacing previous import 'parallel::makeCluster' by
#> 'snow::makeCluster' when loading 'powsimR'
#> Warning: replacing previous import 'parallel::stopCluster' by
#> 'snow::stopCluster' when loading 'powsimR'
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
#> [1] powsimR_1.1.4     gamlss.dist_5.1-5 MASS_7.3-51.5    
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            pbapply_1.4-2              
#>   [5] lattice_0.20-38             vctrs_0.2.2                
#>   [7] fastICA_1.2-2               mgcv_1.8-31                
#>   [9] penalized_0.9-51            blob_1.2.1                 
#>  [11] survival_3.1-8              later_1.0.0                
#>  [13] nloptr_1.2.1                DBI_1.1.0                  
#>  [15] R.utils_2.9.2               SingleCellExperiment_1.8.0 
#>  [17] rappdirs_0.3.1              uwot_0.1.5                 
#>  [19] Linnorm_2.10.0              dqrng_0.2.1                
#>  [21] jpeg_0.1-8.1                zlibbioc_1.32.0            
#>  [23] MatrixModels_0.4-1          pcaMethods_1.78.0          
#>  [25] SDMTools_1.1-221.2          htmlwidgets_1.5.1          
#>  [27] mvtnorm_1.0-12              future_1.16.0              
#>  [29] UpSetR_1.4.0                leiden_0.3.2               
#>  [31] parallel_3.6.2              scater_1.14.6              
#>  [33] irlba_2.3.3                 DEoptimR_1.0-8             
#>  [35] lars_1.2                    Rcpp_1.0.3                 
#>  [37] KernSmooth_2.23-16          DT_0.11                    
#>  [39] promises_1.1.0              gdata_2.18.0               
#>  [41] DDRTree_0.1.5               DelayedArray_0.12.2        
#>  [43] limma_3.42.0                vegan_2.5-6                
#>  [45] RcppParallel_4.4.4          Hmisc_4.3-0                
#>  [47] ShortRead_1.44.1            apcluster_1.4.8            
#>  [49] RSpectra_0.16-0             msir_1.3.2                 
#>  [51] mnormt_1.5-5                digest_0.6.23              
#>  [53] png_0.1-7                   qlcMatrix_0.9.7            
#>  [55] sctransform_0.2.1           cowplot_1.0.0              
#>  [57] pkgconfig_2.0.3             docopt_0.6.1               
#>  [59] DelayedMatrixStats_1.8.0    ggbeeswarm_0.6.0           
#>  [61] iterators_1.0.12            minqa_1.2.4                
#>  [63] lavaan_0.6-5                reticulate_1.14            
#>  [65] SummarizedExperiment_1.16.1 beeswarm_0.2.3             
#>  [67] spam_2.5-1                  modeltools_0.2-22          
#>  [69] xfun_0.12                   zoo_1.8-7                  
#>  [71] tidyselect_1.0.0            ZIM_1.1.0                  
#>  [73] reshape2_1.4.3              purrr_0.3.3                
#>  [75] kernlab_0.9-29              ica_1.0-2                  
#>  [77] EDASeq_2.20.0               snow_0.4-3                 
#>  [79] viridisLite_0.3.0           rtracklayer_1.46.0         
#>  [81] rlang_0.4.4                 hexbin_1.28.0              
#>  [83] glue_1.3.1                  metap_1.3                  
#>  [85] RColorBrewer_1.1-2          fpc_2.2-4                  
#>  [87] matrixStats_0.55.0          stringr_1.4.0              
#>  [89] fields_10.2                 DESeq2_1.26.0              
#>  [91] SparseM_1.78                gbRd_0.4-11                
#>  [93] mutoss_0.1-12               httpuv_1.5.2               
#>  [95] class_7.3-15                BPSC_0.99.2                
#>  [97] RMTstat_0.3                 BiocNeighbors_1.4.1        
#>  [99] TH.data_1.0-10              annotate_1.64.0            
#> [101] jsonlite_1.6.1              XVector_0.26.0             
#> [103] bit_1.1-15.1                mime_0.8                   
#> [105] gridExtra_2.3               gplots_3.0.1.2             
#> [107] Rsamtools_2.2.1             zingeR_0.1.0               
#> [109] stringi_1.4.5               gmodels_2.18.1             
#> [111] bitops_1.0-6                Rdpack_0.11-1              
#> [113] maps_3.3.0                  RSQLite_2.2.0              
#> [115] tidyr_1.0.2                 pheatmap_1.0.12            
#> [117] data.table_1.12.8           DEDS_1.60.0                
#> [119] rstudioapi_0.10             GenomicAlignments_1.22.1   
#> [121] nlme_3.1-143                qvalue_2.18.0              
#> [123] scran_1.14.5                fastcluster_1.1.25         
#> [125] scone_1.10.0                locfit_1.5-9.1             
#> [127] listenv_0.8.0               cobs_1.3-4                 
#> [129] R.oo_1.23.0                 prabclus_2.3-2             
#> [131] segmented_1.1-0             dbplyr_1.4.2               
#> [133] BiocGenerics_0.32.0         lifecycle_0.1.0            
#> [135] ROTS_1.14.0                 munsell_0.5.0              
#> [137] R.methodsS3_1.7.1           moments_0.14               
#> [139] hwriter_1.3.2               caTools_1.18.0             
#> [141] codetools_0.2-16            coda_0.19-3                
#> [143] Biobase_2.46.0              vipor_0.4.5                
#> [145] GenomeInfoDb_1.22.0         lmtest_0.9-37              
#> [147] htmlTable_1.13.3            rARPACK_0.11-0             
#> [149] lsei_1.2-0                  xtable_1.8-4               
#> [151] SAVER_1.1.2                 ROCR_1.0-7                 
#> [153] diptest_0.75-7              formatR_1.7                
#> [155] lpsymphony_1.14.0           abind_1.4-5                
#> [157] FNN_1.1.3                   RANN_2.6.1                 
#> [159] askpass_1.1                 sparsesvd_0.2              
#> [161] CompQuadForm_1.4.3          GenomicRanges_1.38.0       
#> [163] bibtex_0.4.2.2              RcppAnnoy_0.0.14           
#> [165] tibble_2.1.3                ggdendro_0.1-20            
#> [167] cluster_2.1.0               future.apply_1.4.0         
#> [169] Seurat_3.1.2                Matrix_1.2-18              
#> [171] prettyunits_1.1.1           shinyBS_0.61               
#> [173] ggridges_0.5.2              NOISeq_2.30.0              
#> [175] shinydashboard_0.7.1        mclust_5.4.5               
#> [177] igraph_1.2.4.2              multtest_2.42.0            
#> [179] slam_0.1-47                 TFisher_0.2.0              
#> [181] testthat_2.3.1              htmltools_0.4.0            
#> [183] BiocFileCache_1.10.2        yaml_2.2.1                 
#> [185] GenomicFeatures_1.38.1      plotly_4.9.1               
#> [187] XML_3.99-0.3                DrImpute_1.0               
#> [189] foreign_0.8-75              fitdistrplus_1.0-14        
#> [191] BiocParallel_1.20.1         aroma.light_3.16.0         
#> [193] bit64_0.9-7                 rngtools_1.5               
#> [195] doRNG_1.8.2                 multcomp_1.4-12            
#> [197] foreach_1.4.7               robustbase_0.93-5          
#> [199] outliers_0.14               scde_2.14.0                
#> [201] Biostrings_2.54.0           combinat_0.0-8             
#> [203] rsvd_1.0.2                  iCOBRA_1.14.0              
#> [205] memoise_1.1.0               evaluate_0.14              
#> [207] VGAM_1.1-2                  nonnest2_0.5-2             
#> [209] geneplotter_1.64.0          permute_0.9-5              
#> [211] curl_4.3                    fdrtool_1.2.15             
#> [213] acepack_1.4.1               edgeR_3.28.0               
#> [215] checkmate_1.9.4             npsurv_0.4-0               
#> [217] tensorA_0.36.1              DECENT_1.1.0               
#> [219] ellipse_0.4.1               rjson_0.2.20               
#> [221] ggplot2_3.2.1               ggrepel_0.8.1              
#> [223] distillery_1.0-6            scDD_1.10.0                
#> [225] Lmoments_1.3-1              tools_3.6.2                
#> [227] sandwich_2.5-1              magrittr_1.5               
#> [229] RCurl_1.98-1.1              pbivnorm_0.6.0             
#> [231] ape_5.3                     bayesm_3.1-4               
#> [233] EBSeq_1.26.0                httr_1.4.1                 
#> [235] assertthat_0.2.1            rmarkdown_2.1              
#> [237] Rhdf5lib_1.8.0              boot_1.3-24                
#> [239] globals_0.12.5              R6_2.4.1                   
#> [241] nnet_7.3-12                 progress_1.2.2             
#> [243] genefilter_1.68.0           gtools_3.8.1               
#> [245] statmod_1.4.33              Rook_1.1-1                 
#> [247] BiocSingular_1.2.1          rhdf5_2.30.1               
#> [249] splines_3.6.2               colorspace_1.4-1           
#> [251] amap_0.8-18                 generics_0.0.2             
#> [253] stats4_3.6.2                NBPSeq_0.3.0               
#> [255] compositions_1.40-3         base64enc_0.1-3            
#> [257] baySeq_2.20.0               pillar_1.4.3               
#> [259] sn_1.5-5                    HSMMSingleCell_1.6.0       
#> [261] GenomeInfoDbData_1.2.2      plyr_1.8.5                 
#> [263] extRemes_2.0-11             dotCall64_1.0-0            
#> [265] gtable_0.3.0                bdsmatrix_1.3-4            
#> [267] SCnorm_1.8.2                monocle_2.14.0             
#> [269] knitr_1.27                  RcppArmadillo_0.9.800.4.0  
#> [271] latticeExtra_0.6-29         biomaRt_2.42.0             
#> [273] IRanges_2.20.2              fastmap_1.0.1              
#> [275] Cairo_1.5-10                doParallel_1.0.15          
#> [277] pscl_1.5.2                  flexmix_2.3-15             
#> [279] quantreg_5.54               AnnotationDbi_1.48.0       
#> [281] broom_0.5.4                 openssl_1.4.1              
#> [283] scales_1.1.0                arm_1.10-1                 
#> [285] backports_1.1.5             plotrix_3.7-7              
#> [287] IHW_1.14.0                  S4Vectors_0.24.3           
#> [289] densityClust_0.3            lme4_1.1-21                
#> [291] brew_1.0-6                  hms_0.5.3                  
#> [293] DESeq_1.38.0                Rtsne_0.15                 
#> [295] dplyr_0.8.4                 shiny_1.4.0                
#> [297] grid_3.6.2                  numDeriv_2016.8-1.1        
#> [299] bbmle_1.0.22                lazyeval_0.2.2             
#> [301] Formula_1.2-3               tsne_0.1-3                 
#> [303] blockmodeling_0.3.4         crayon_1.3.4               
#> [305] MAST_1.12.0                 RUVSeq_1.20.0              
#> [307] viridis_0.5.1               rpart_4.1-15               
#> [309] zinbwave_1.8.0              compiler_3.6.2
```

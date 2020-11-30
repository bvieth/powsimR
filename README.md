
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `powsimR` <br/> Power analysis for bulk and <br/> single cell RNA-seq experiments <img src="vignettes/powsimR.png" align="right" width="200" />

Please also consult my Github Page of
[powsimR](https://bvieth.github.io/powsimR/) made with
[pkgdown](http://pkgdown.r-lib.org/index.html)!

## :arrow\_double\_down: Installation Guide

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
vignette (by setting build\_vignettes to FALSE) and read it on my Github
Page of
[powsimR](https://bvieth.github.io/powsimR/articles/powsimR.html) or
download it as a html file
[here](https://github.com/bvieth/powsimR/blob/master/vignettes/powsimR.html).

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
be set to a higher number to accommodate the increase in DLLs. Please
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
#> Registered S3 method overwritten by 'gdata':
#>   method         from  
#>   reorder.factor gplots
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
#> R version 4.0.3 (2020-10-10)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.5 LTS
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
#> [1] powsimR_1.2.3     gamlss.dist_5.1-7 MASS_7.3-53      
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.2.0              softImpute_1.4             
#>   [3] minpack.lm_1.2-1            lattice_0.20-41            
#>   [5] haven_2.3.1                 vctrs_0.3.5                
#>   [7] fastICA_1.2-2               mgcv_1.8-33                
#>   [9] penalized_0.9-51            blob_1.2.1                 
#>  [11] survival_3.2-7              Rmagic_2.0.3               
#>  [13] later_1.1.0.1               nloptr_1.2.2.2             
#>  [15] DBI_1.1.0                   R.utils_2.10.1             
#>  [17] SingleCellExperiment_1.12.0 rappdirs_0.3.1             
#>  [19] Linnorm_2.14.0              dqrng_0.2.1                
#>  [21] jpeg_0.1-8.1                zlibbioc_1.36.0            
#>  [23] MatrixModels_0.4-1          htmlwidgets_1.5.2          
#>  [25] mvtnorm_1.1-1               future_1.20.1              
#>  [27] UpSetR_1.4.0                parallel_4.0.3             
#>  [29] scater_1.18.3               irlba_2.3.3                
#>  [31] DEoptimR_1.0-8              Rcpp_1.0.5                 
#>  [33] KernSmooth_2.23-18          DT_0.16                    
#>  [35] promises_1.1.1              gdata_2.18.0               
#>  [37] DDRTree_0.1.5               DelayedArray_0.16.0        
#>  [39] limma_3.46.0                vegan_2.5-7                
#>  [41] Hmisc_4.4-2                 ShortRead_1.48.0           
#>  [43] apcluster_1.4.8             RSpectra_0.16-0            
#>  [45] msir_1.3.2                  mnormt_2.0.2               
#>  [47] digest_0.6.27               png_0.1-7                  
#>  [49] bluster_1.0.0               qlcMatrix_0.9.7            
#>  [51] sctransform_0.3.1           cowplot_1.1.0              
#>  [53] pkgconfig_2.0.3             docopt_0.7.1               
#>  [55] DelayedMatrixStats_1.12.0   ggbeeswarm_0.6.0           
#>  [57] iterators_1.0.13            minqa_1.2.4                
#>  [59] lavaan_0.6-7                reticulate_1.18            
#>  [61] SummarizedExperiment_1.20.0 spam_2.5-1                 
#>  [63] beeswarm_0.2.3              modeltools_0.2-23          
#>  [65] xfun_0.19                   zoo_1.8-8                  
#>  [67] tidyselect_1.1.0            ZIM_1.1.0                  
#>  [69] reshape2_1.4.4              purrr_0.3.4                
#>  [71] kernlab_0.9-29              EDASeq_2.24.0              
#>  [73] viridisLite_0.3.0           snow_0.4-3                 
#>  [75] rtracklayer_1.50.0          rlang_0.4.9                
#>  [77] hexbin_1.28.1               glue_1.4.2                 
#>  [79] RColorBrewer_1.1-2          fpc_2.2-8                  
#>  [81] matrixStats_0.57.0          MatrixGenerics_1.2.0       
#>  [83] stringr_1.4.0               fields_11.6                
#>  [85] ggsignif_0.6.0              DESeq2_1.30.0              
#>  [87] SparseM_1.78                httpuv_1.5.4               
#>  [89] class_7.3-17                BPSC_0.99.2                
#>  [91] BiocNeighbors_1.8.1         annotate_1.68.0            
#>  [93] jsonlite_1.7.1              XVector_0.30.0             
#>  [95] tmvnsim_1.0-2               bit_4.0.4                  
#>  [97] mime_0.9                    gridExtra_2.3              
#>  [99] gplots_3.1.1                Rsamtools_2.6.0            
#> [101] zingeR_0.1.0                stringi_1.5.3              
#> [103] gmodels_2.18.1              rhdf5filters_1.2.0         
#> [105] bitops_1.0-6                maps_3.3.0                 
#> [107] RSQLite_2.2.1               tidyr_1.1.2                
#> [109] pheatmap_1.0.12             data.table_1.13.2          
#> [111] rstudioapi_0.13             GenomicAlignments_1.26.0   
#> [113] nlme_3.1-150                qvalue_2.22.0              
#> [115] scran_1.18.1                fastcluster_1.1.25         
#> [117] locfit_1.5-9.4              scone_1.14.0               
#> [119] listenv_0.8.0               cobs_1.3-4                 
#> [121] R.oo_1.24.0                 prabclus_2.3-2             
#> [123] segmented_1.3-0             dbplyr_2.0.0               
#> [125] BiocGenerics_0.36.0         readxl_1.3.1               
#> [127] lifecycle_0.2.0             ROTS_1.18.0                
#> [129] munsell_0.5.0               cellranger_1.1.0           
#> [131] R.methodsS3_1.8.1           moments_0.14               
#> [133] hwriter_1.3.2               caTools_1.18.0             
#> [135] codetools_0.2-18            coda_0.19-4                
#> [137] Biobase_2.50.0              GenomeInfoDb_1.26.1        
#> [139] vipor_0.4.5                 htmlTable_2.1.0            
#> [141] bayNorm_1.8.0               rARPACK_0.11-0             
#> [143] xtable_1.8-4                SAVER_1.1.2                
#> [145] ROCR_1.0-11                 diptest_0.75-7             
#> [147] formatR_1.7                 lpsymphony_1.18.0          
#> [149] abind_1.4-5                 FNN_1.1.3                  
#> [151] parallelly_1.21.0           RANN_2.6.1                 
#> [153] askpass_1.1                 sparsesvd_0.2              
#> [155] CompQuadForm_1.4.3          GenomicRanges_1.42.0       
#> [157] tibble_3.0.4                ggdendro_0.1.22            
#> [159] cluster_2.1.0               future.apply_1.6.0         
#> [161] Matrix_1.2-18               ellipsis_0.3.1             
#> [163] prettyunits_1.1.1           shinyBS_0.61               
#> [165] NOISeq_2.34.0               shinydashboard_0.7.1       
#> [167] mclust_5.4.7                igraph_1.2.6               
#> [169] ggstance_0.3.4              slam_0.1-47                
#> [171] testthat_3.0.0              doSNOW_1.0.19              
#> [173] htmltools_0.5.0             BiocFileCache_1.14.0       
#> [175] yaml_2.2.1                  GenomicFeatures_1.42.1     
#> [177] XML_3.99-0.5                ggpubr_0.4.0               
#> [179] DrImpute_1.0                foreign_0.8-80             
#> [181] scuttle_1.0.3               fitdistrplus_1.1-1         
#> [183] BiocParallel_1.24.1         aroma.light_3.20.0         
#> [185] bit64_4.0.5                 foreach_1.5.1              
#> [187] robustbase_0.93-6           outliers_0.14              
#> [189] Biostrings_2.58.0           combinat_0.0-8             
#> [191] rsvd_1.0.3                  iCOBRA_1.18.0              
#> [193] memoise_1.1.0               evaluate_0.14              
#> [195] VGAM_1.1-4                  nonnest2_0.5-5             
#> [197] forcats_0.5.0               rio_0.5.16                 
#> [199] geneplotter_1.68.0          permute_0.9-5              
#> [201] curl_4.3                    fdrtool_1.2.15             
#> [203] conquer_1.0.2               edgeR_3.32.0               
#> [205] checkmate_2.0.0             truncnorm_1.0-8            
#> [207] DECENT_1.1.0                tensorA_0.36.2             
#> [209] ellipse_0.4.2               ggplot2_3.3.2              
#> [211] openxlsx_4.2.3              rstatix_0.6.0              
#> [213] ggrepel_0.8.2               scDD_1.14.0                
#> [215] tools_4.0.3                 sandwich_3.0-0             
#> [217] magrittr_2.0.1              RCurl_1.98-1.2             
#> [219] car_3.0-10                  pbivnorm_0.6.0             
#> [221] bayesm_3.1-4                xml2_1.3.2                 
#> [223] EBSeq_1.30.0                httr_1.4.2                 
#> [225] assertthat_0.2.1            rmarkdown_2.5              
#> [227] Rhdf5lib_1.12.0             boot_1.3-25                
#> [229] globals_0.14.0              R6_2.5.0                   
#> [231] nnet_7.3-14                 progress_1.2.2             
#> [233] genefilter_1.72.0           gtools_3.8.2               
#> [235] statmod_1.4.35              beachmat_2.6.1             
#> [237] BiocSingular_1.6.0          rhdf5_2.34.0               
#> [239] splines_4.0.3               carData_3.0-4              
#> [241] colorspace_2.0-0            amap_0.8-18                
#> [243] generics_0.1.0              stats4_4.0.3               
#> [245] NBPSeq_0.3.0                base64enc_0.1-3            
#> [247] compositions_2.0-0          baySeq_2.24.0              
#> [249] pillar_1.4.7                HSMMSingleCell_1.10.0      
#> [251] GenomeInfoDbData_1.2.4      plyr_1.8.6                 
#> [253] dotCall64_1.0-0             gtable_0.3.0               
#> [255] zip_2.1.1                   SCnorm_1.12.0              
#> [257] monocle_2.18.0              knitr_1.30                 
#> [259] RcppArmadillo_0.10.1.2.0    latticeExtra_0.6-29        
#> [261] biomaRt_2.46.0              IRanges_2.24.0             
#> [263] fastmap_1.0.1               doParallel_1.0.16          
#> [265] pscl_1.5.5                  flexmix_2.3-17             
#> [267] quantreg_5.75               AnnotationDbi_1.52.0       
#> [269] broom_0.7.2                 openssl_1.4.3              
#> [271] scales_1.1.1                arm_1.11-2                 
#> [273] backports_1.2.0             plotrix_3.7-8              
#> [275] IHW_1.18.0                  S4Vectors_0.28.0           
#> [277] densityClust_0.3            lme4_1.1-25                
#> [279] hms_0.5.3                   Rtsne_0.15                 
#> [281] dplyr_1.0.2                 shiny_1.5.0                
#> [283] grid_4.0.3                  Formula_1.2-4              
#> [285] blockmodeling_1.0.0         crayon_1.3.4               
#> [287] MAST_1.16.0                 RUVSeq_1.24.0              
#> [289] sparseMatrixStats_1.2.0     viridis_0.5.1              
#> [291] rpart_4.1-15                zinbwave_1.12.0            
#> [293] compiler_4.0.3
```

# Version 1.1.2 (2018-06-20)

* `estimateParam()` error fixed concerning expression cleanup.
* precompiled vignette in inst/doc/.

# Version 1.1.1 (2018-04-19)

* `simulateDE()` now with the option to perform DE testing on filtered/imputed counts (option `DEFilter`)

# Version 1.1.0 (2018-03-29)

* simulation of batch effects (see options `p.B`, `bLFC` and `bPattern` in `DESetup()` and `simulateCounts()`)
* simulation of spike-in expression (see `estimateSpike` , `plotSpike` and option `spikeIns` in `simulateDE` and `simulateCounts()`)
* simulation of multiple sample groups (e.g. single cell populations) with `simulateCounts()`
* imputation and prefiltering options prior to normalisation in DE power simulations added (scImpute, scone, Seurat, DrImpute, SAVER)
* additional normalisation options and DE tools (esp. for single cells) included in `simulateDE()`
* evaluation of simulation setup using estimated versus true value comparisons of library size factors and log2 fold changes in `evaluateSim()` and `plotEvalSim()`
* increased flexibility in preprocessing for distribution evaluation in `evaluateDist()`


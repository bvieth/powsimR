
# estimateUnwantedVariation -----------------------------------------------


# estimateCellType --------------------------------------------------------

estimateCellType <- function(countData,
                             readData = NULL,
                             annotData = NULL,
                             spikeData = NULL,
                             spikeInfo = NULL,
                             Lengths = NULL,
                             MeanFragLengths = NULL,
                             RNAseq = c('bulk', 'singlecell'),
                             Protocol = c('UMI', 'Read'),
                             Distribution = c('NB', 'ZINB'),
                             Normalisation = c("TMM", "MR", "PosCounts", "UQ",
                                               "scran", "Linnorm", "sctransform",
                                               "SCnorm", "Census", "depth", "none"),
                             GeneFilter = 0.05,
                             SampleFilter = 5,
                             sigma = 1.96,
                             NCores = NULL,
                             verbose = TRUE){

}

# simulateCounts ----------------------------------------------------------





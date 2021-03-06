% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimate.R
\name{estimateParam}
\alias{estimateParam}
\title{Estimate simulation parameters}
\usage{
estimateParam(countData,
readData = NULL,
batchData = NULL,
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
verbose=TRUE)
}
\arguments{
\item{countData}{is a (UMI) count \code{matrix} of gene expression.
Rows correspond to genes, columns to samples.
The gene names should be given as \code{rownames} without "_" in the names.
The samples names should be given as \code{colnames}.
The count matrix should only contain the expression of one group, e.g. wildtype / untreated control / one cell type population.}

\item{readData}{is a the matching read count \code{matrix} of gene expression if \code{countData} is UMI
and the same formatting should be applied. Default is \code{NULL} and
users only need to supply the read count matrix if they plan to apply downsampling of UMI counts for simulations, see \code{\link{Setup}}.}

\item{batchData}{is a \code{data.frame} for batch annotation.
Rows correspond to samples. The first column should contain the batches, e.g. 'a', 'b', 'c', etc.}

\item{spikeData}{is a count \code{matrix}.
Rows correspond to spike-ins, columns to samples.
The order of columns should be the same as in the \code{countData}.
This is only needed for spike-in aware normalisation methods ('MR', 'Linnorm', 'scran', 'SCnorm', 'bayNorm', 'Census'), see Details.}

\item{spikeInfo}{is a molecule count \code{matrix} of spike-ins.
Rows correspond to spike-ins. The order of rows should be the same as in the \code{spikeData}.
The column names should be 'SpikeID' and 'SpikeInput' for molecule counts of spike-ins.
This is only needed for spike-in aware normalisation methods (), see Details.}

\item{Lengths}{is a numeric vector of transcript lengths with the same length and order as the rows in countData.
This variable is only needed for internal gene length corrections (TPM), see Details.}

\item{MeanFragLengths}{is a numeric vector of mean fragment lengths with the same length as columns in countData.
This variable is only needed for internal gene length corrections (TPM), see Details.}

\item{RNAseq}{is a character value: "bulk" or "singlecell".}

\item{Protocol}{is a character value defining the type of counts given in \code{countData}.
Options are "UMI" (e.g. 10X Genomics, CEL-seq2) or "Read" (e.g. Smart-seq2).}

\item{Distribution}{is a character value: "NB" for negative binomial or "ZINB" for zero-inflated negative binomial distribution fitting.}

\item{Normalisation}{is a character value: 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm',
'SCnorm', 'Census', 'depth', 'none'.
For more information, please consult the Details section.}

\item{GeneFilter}{is a numeric vector indicating the minimal proportion of nonzero expression values
for a gene across all samples to be considered expressed and used for normalisation and parameter estimation.
The default is \code{0.05}, i.e. at least 5\% of the expression values per gene need to be nonzero.}

\item{SampleFilter}{is a numeric vector indicating the minimal number of MADs (median absolute deviation)
away from the median number of features detected as well as sequencing depth across all samples
so that outlying samples are removed prior to normalisation and parameter estimation.
The default is \code{5}, i.e. at least 5 MADs away from the median.
Choose higher values if you want to filter out less samples.
This parameter is particularly important for single cells to ensure reliable parameter estimation.
For more information, please consult \code{\link[scater]{isOutlier}}.}

\item{sigma}{The variability band width for mean-dispersion loess fit defining the prediction interval for read count simulation. Default is 1.96, i.e. 95\% interval. For more information see \code{\link[msir]{loess.sd}}.}

\item{NCores}{The number of cores for normalisation method SCnorm and Census.
The default \code{NULL} means 1 core.}

\item{verbose}{Logical value to indicate whether to print function information.
Default is \code{TRUE}.}
}
\value{
List object with the following entries:
\item{Parameters}{A list object containing the estimated moments for the full, dropped out genes, dropped out samples and filtered normalized count matrix. For more information please consult the details section and the plot made with \code{\link{plotParam}}.}
\item{Fit}{A list object containing the fitting results of the mean-dispersion and mean-dropout relation as well as the estimated parameter data used for the fits. For more information please consult the details section and the plot made with \code{\link{plotParam}}.}
\item{totalS,totalG}{Number of samples and genes provided with at least one read count.}
\item{DropOuts}{List object containing logical vectors for gene and sample dropouts after applying gene frequency and sample outlier filtering.}
\item{sf}{The estimated library size factor per sample.}
\item{Lengths -..- SampleFilter}{The chosen parameters settings.}
}
\description{
This function estimates and returns parameters needed for power simulations.\cr
The user needs to choose the following options at least: specify a gene expression matrix; the type of RNA-seq experiment, i.e. bulk or single cell; the recommended distribution is negative binomial (NB) except for single-cell full-length Smart-seq2 read data where we recommend zero-inflated NB (ZINB); the preferred normalisation method, we recommend scran for single cell and TMM or MR for bulk.\cr
The other parameters are optional (additional data) or have preset values (gene and sample filtering). Please consult the detailed arguments description.
}
\details{
\strong{Normalisation Methods}
\describe{
\item{TMM, UQ}{employ the edgeR style normalization of weighted trimmed mean of M-values and upperquartile
as implemented in \code{\link[edgeR]{calcNormFactors}}, respectively.}
\item{MR, PosCounts}{employ the DESeq2 style normalization of median ratio method and a modified geometric mean method
as implemented in \code{\link[DESeq2]{estimateSizeFactors}}, respectively.}
\item{scran, SCnorm}{apply the deconvolution and quantile regression normalization methods developed for sparse RNA-seq data
as implemented in \code{\link[scran]{calculateSumFactors}} and \code{\link[SCnorm]{SCnorm}}, respectively. Spike-ins can also be supplied for both methods via \code{spikeData}. Note, however that this means for scran that the normalisation as implemented in \code{\link[scuttle]{computeSpikeFactors}} is also applied to genes (\code{general.use=TRUE}).}
\item{Linnorm}{apply the normalization method for sparse RNA-seq data
as implemented in \code{\link[Linnorm]{Linnorm.Norm}}.
For \code{Linnorm}, the user can also supply \code{spikeData}.}
\item{sctransform}{apply the normalization method developed for single-cell
UMI RNA-seq data as implemented in \code{\link[sctransform]{vst}}. }
\item{Census}{converts relative measures of TPM/FPKM values into mRNAs per cell (RPC) without the need of spike-in standards.
Census at least needs \code{Lengths} for single-end data and preferably \code{MeanFragLengths} for paired-end data.
The authors state that Census should not be used for UMI data.}
\item{depth}{Sequencing depth normalisation.}
\item{none}{No normalisation is applied. This approach can be used for prenormalized expression estimates, e.g. cufflinks, RSEM or salmon.}
}
}
\examples{
\dontrun{
# Single Cells
data("SmartSeq2_Gene_Read_Counts")
Batches <- data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
                      stringsAsFactors = FALSE, row.names = colnames(SmartSeq2_Gene_Read_Counts))
data("GeneLengths_mm10")
estparam <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
                          readData = NULL,
                          batchData = Batches,
                          spikeData = SmartSeq2_SpikeIns_Read_Counts,
                          spikeInfo = SmartSeq2_SpikeInfo,
                          Lengths = GeneLengths, MeanFragLengths = NULL,
                          RNAseq = 'singlecell', Protocol = 'Read',
                          Distribution = 'ZINB', Normalisation = "scran",
                          GeneFilter = 0.1, SampleFilter = 3,
                          sigma = 1.96, NCores = NULL, verbose = TRUE)

# Bulk
data("Bulk_Read_Counts")
data("GeneLengths_hg19")
estparam <- estimateParam(countData = Bulk_Read_Counts,
                          readData = NULL,
                          batchData = NULL,
                          spikeData = NULL,
                          spikeInfo = NULL,
                         Lengths = GeneLengths_hg19,
                          MeanFragLengths = NULL,
                          RNAseq = 'bulk', Protocol = 'Read',
                          Distribution = 'NB', Normalisation = "MR",
                          GeneFilter = 0.1, SampleFilter = 3,
                          sigma = 1.96, NCores = NULL, verbose = TRUE)

# plot the results of estimation
plotParam(estparam, Annot = FALSE)
}
}
\author{
Beate Vieth
}

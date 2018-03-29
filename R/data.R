#' Kolodziejczk et al. 2015: Gene Counts
#'
#' 869 mouse embryonic stem cells under three different conditions (cnts data object):
#' serum + LIF (242 cells), 2i + LIF (433 cells) and alternative 2i + LIF (194 cells).
#' Single cell RNA-Seq was performed using Fluidigm C1 system and
#' libraries were generated using Nextera XT (Illumina) kit.
#' The read count matrix was downsampled to 10000 genes.
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/}.
"kolodziejczk_cnts"


#' Kolodziejczk et al. 2015: Parameter Estimation
#'
#' 869 mouse embryonic stem cells under three different conditions (cnts data object):
#' serum + LIF (242 cells), 2i + LIF (433 cells) and alternative 2i + LIF (194 cells).
#' Single cell RNA-Seq was performed using Fluidigm C1 system and
#' libraries were generated using Nextera XT (Illumina) kit.
#' We reduced the data set to standard serum + LIF cultured cells for
#' parameter estimation.
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/}.
"kolodziejczk_param"

#' Kolodziejczk et al. 2015: Simulation Result
#'
#' 869 mouse embryonic stem cells under three different conditions (cnts data object):
#' serum + LIF (242 cells), 2i + LIF (433 cells) and alternative 2i + LIF (194 cells).
#' Single cell RNA-Seq was performed using Fluidigm C1 system and
#' libraries were generated using Nextera XT (Illumina) kit.
#' We reduced the data set to standard serum + LIF cultured cells for
#' parameter estimation. For the simulation 20% of 10000 simulated genes were differentially expressed with a broad symmetrical gamma distribution of log2 fold changes. For DE testing, limma trend was used.
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/}.
"kolodziejczk_simDE"


#' Ziegenhain et al. 2017: Smartseq2 Gene Counts
#'
#' J1 mouse embryonic stem cells (cultured in medium) were facs sorted and processed in two batches.
#' Single cell RNA-Seq was performed using Smartseq2 protocol and libraries were generated using Nextera XT (Illumina) kit.
#' Furthermore, spike-ins (ERCC) were added.
#' The read count matrix was downsampled to 10000 genes.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"smartseq2_gene_cnts"

#' Ziegenhain et al. 2017: SCRBseq Gene Counts
#'
#' J1 mouse embryonic stem cells (cultured in medium) were facs sorted and processed in two batches.
#' Single cell RNA-Seq was performed using Single Cell RNA Barcoding and Sequencing (SCRB-Seq) with unique molecular identifiers (UMI) and libraries were generated using Nextera XT (Illumina) kit.
#' Furthermore, spike-ins (ERCC) were added.
#' The read count matrix was downsampled to 10000 genes.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"scrbseq_gene_cnts"

#' Ziegenhain et al. 2017: Gene Lengths
#'
#' ENSEMBL gene length annotation for SCRBseq and Smartseq2 library preparation.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"gene_lengths"

#' Ziegenhain et al. 2017: SCRBseq Spike-Ins Counts
#'
#' Spike-ins read counts for SCRBseq library preparation.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"scrbseq_spike_cnts"

#' Ziegenhain et al. 2017: Smartseq2 Spike-Ins Counts
#'
#' Spike-ins read counts for Smartseq2 library preparation.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"smartseq2_spike_cnts"

#' Ziegenhain et al. 2017: SCRBseq Spike-Ins Information
#'
#' Spike-ins information table with molecules and spike-ins transcript lengths for SCRBseq library preparation.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"scrbseq_spike_info"

#' Ziegenhain et al. 2017: Smartseq2 Spike-Ins Information
#'
#' Spike-ins information table with molecules and spike-ins transcript lengths for Smartseq2 library preparation.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"smartseq2_spike_info"


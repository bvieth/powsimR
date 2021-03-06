% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SpikeCounts}
\alias{SpikeCounts}
\title{Ziegenhain et al. 2017: Spike-Ins Expression}
\format{
An object of class \code{"data.frame"}, with spike-ins per row and cells per column.
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
}
\usage{
data(CELseq2_SpikeIns_Read_Counts)
data(CELseq2_SpikeIns_UMI_Counts)
data(DropSeq_SpikeIns_Read_Counts)
data(DropSeq_SpikeIns_UMI_Counts)
data(MARSseq_SpikeIns_Read_Counts)
data(MARSseq_SpikeIns_UMI_Counts)
data(SCRBseq_SpikeIns_Read_Counts)
data(SCRBseq_SpikeIns_UMI_Counts)
data(SmartSeq_SpikeIns_Read_Counts)
data(SmartSeq2_SpikeIns_Read_Counts)
}
\description{
Spike-ins counts for Ziegenhain et al. 2017 library preparation protocols.\cr
Either Reads (Smartseq, Smartseq2) or UMI Counts (CELseq2, Dropseq, MARSseq, SCRBseq) per ERCC Spike-In and Cell.
}
\keyword{datasets}

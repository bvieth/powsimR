
#' Parekh et al. 2016: Gene Expression Matrices
#'
#' 250 ng of Universal Human Reference RNA (UHRR; Agilent Technologies; catalog #740000) and ERCC spike-in control mix I (Life Technologies) were used and cDNA was synthesized as described in the Smart-Seq2 protocol from Picelli et al. 2013.\cr
#' Base-calls were performed with Illumina bcl2fastq.\cr
#' All cDNA reads were aligned with NextGenMap 0.4.12. Smart-seq2 reads were mapped to GRCh37 (hg19) and ERCC synthetic spike-in reference.
#' Reads were assigned to genes using FeatureCount method in Rsubread package.\cr
#'
#' @docType data
#'
#' @usage
#' data(Bulk_Read_Counts)
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75823}.
"Bulk_Read_Counts"

#' Parekh et al. 2016: Gene Lengths
#'
#' ENSEMBL gene length annotation for gene expression quantification in Parekh et al. 2016.
#'
#' @docType data
#'
#' @usage
#' data(GeneLengths_hg19)
#'
#' @format A named vector of gene lengths.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75823}.
"GeneLengths_hg19"

#' Ziegenhain et al. 2017: Gene Lengths
#'
#' ENSEMBL gene length annotation for gene expression quantification in Ziegenhain et al. 2017.
#'
#' @docType data
#'
#' @usage
#' data(GeneLengths_mm10)
#'
#' @format A named vector of gene lengths.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"GeneLengths_mm10"

#' Ziegenhain et al. 2017: Gene Expression Matrices
#'
#' J1 mESC were cultured in ground pluripotency states by 2i/LIF culture as described previously.\cr
#' The cells were sorted (FACS) and processed in two batches.\cr
#' Library preparations were performed with different protocols. Furthermore, spike-ins (ERCC) were added.\cr
#' "CELseq2" : Library preparations were performed with CEL-seq2 on the Fluidigm C1, as described in Hashimshony et al. 2016.\cr
#' "DropSeq": Drop-seq protocol by emulsion of single-cells with primer-microbreads, as described in Macosko et al. 2015.\cr
#' "MARSseq" : MARS-seq protocol from FACS-sorted single cells as described in Jaitin et al. 2014.\cr
#' "SCRBseq" : SCRB-seq protocol from FACS-sorted single cells as described in Soumillon et al. 2014.\cr
#' "SmartSeq" : Library preparation was performed with the Fluidigm C1 platform as per manufacturers' protocol.\cr
#' "SmartSeq2" : Library preparations were performed using the Smart-seq2 protocol, as described in Picelli et al. 2013.\cr
#' \cr
#' Data Processing:\cr
#' base-calls were performed with Illumina bcl2fastq.\cr
#' All cDNA reads were aligned with Star v2.4.0.\cr
#' SCRB-seq and Drop-seq protocol reads were tagged with its mate-pair UMI and cell barcode; unique reads were obtained.\cr
#' Smart-seq and Smart-seq2 reads were assigned to genes using Rsubread.\cr
#' Only cells were considered with at least 1,000,000 reads.\cr
#' The R Data files are once for Read Counts and once for UMI Counts per Library Preparation Protocol.
#'
#' @docType data
#'
#' @usage
#' data(CELseq2_Gene_Read_Counts)
#' data(CELseq2_Gene_UMI_Counts)
#' data(DropSeq_Gene_Read_Counts)
#' data(DropSeq_Gene_UMI_Counts)
#' data(MARSseq_Gene_Read_Counts)
#' data(MARSseq_Gene_UMI_Counts)
#' data(SCRBseq_Gene_Read_Counts)
#' data(SCRBseq_Gene_UMI_Counts)
#' data(SmartSeq_Gene_Read_Counts)
#' data(SmartSeq2_Gene_Read_Counts)
#'
#' @format An object of class \code{"data.frame"}, with genes per row and cells per column.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"GeneCounts"

#' Ziegenhain et al. 2017: Spike-Ins Expression
#'
#' Spike-ins counts for Ziegenhain et al. 2017 library preparation protocols.\cr
#' Either Reads (Smartseq, Smartseq2) or UMI Counts (CELseq2, Dropseq, MARSseq, SCRBseq) per ERCC Spike-In and Cell.
#'
#' @docType data
#'
#' @usage
#' data(CELseq2_SpikeIns_Read_Counts)
#' data(CELseq2_SpikeIns_UMI_Counts)
#' data(DropSeq_SpikeIns_Read_Counts)
#' data(DropSeq_SpikeIns_UMI_Counts)
#' data(MARSseq_SpikeIns_Read_Counts)
#' data(MARSseq_SpikeIns_UMI_Counts)
#' data(SCRBseq_SpikeIns_Read_Counts)
#' data(SCRBseq_SpikeIns_UMI_Counts)
#' data(SmartSeq_SpikeIns_Read_Counts)
#' data(SmartSeq2_SpikeIns_Read_Counts)
#'
#' @format An object of class \code{"data.frame"}, with spike-ins per row and cells per column.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
"SpikeCounts"

#' Ziegenhain et al. 2017: Spike-Ins Information
#'
#' Spike-ins information table with molecules and spike-ins transcript lengths per library preparation protocol.
#'
#' @docType data
#'
#' @usage
#' data(CELseq2_SpikeInfo)
#' data(MARSseq_SpikeInfo)
#' data(SCRBseq_SpikeInfo)
#' data(SmartSeq_SpikeInfo)
#' data(SmartSeq2_SpikeInfo)
#'
#' @format An object of class \code{"data.frame"}, with spike-in IDs, number of spike-in molecules and lengths of spike-ins.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790}.
'SpikeInfo'

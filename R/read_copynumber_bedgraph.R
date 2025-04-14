#' Import HiFiCNV output copynumber.bedgraph and convert to Granges
#'
#' A wrapper of plyranges::read_bed_graph() to read HiFiCNV output copynumber.bedgraph file, convert to GRanges, and annotation the copy number score
#'
#' @details
#' Copy number variant type (SVTYPE) = "DUP" if score > 2 and = "DEL" if score = 0. NA otherwise.
#' 
#' @param file A path to the VCF (GZ) file 
#' @param keep_seq_levels A character vector indicating which seqlevels to keep
#'
#' @return A \code{GRanges} object with copy number variant type annotation
#'
#' @examples
#' 
#' @export
read_copynumber_bedgraph <- function(file, keep_seq_levels) {
  stopifnot(file.exists(file))
  
  copy_number <- plyranges::read_bed_graph(file) %>%
    plyranges::mutate(SVTYPE = case_when(
      score < 2 ~ 'DEL',
      score > 2 ~ 'DUP',
      .default = NA
    ))
  
  if (!is.null(keep_seq_levels))
    copy_number <- keepSeqlevels(copy_number,
                                 value=keep_seq_levels,
                                 pruning.mode='coarse')
  
  return(copy_number)
}
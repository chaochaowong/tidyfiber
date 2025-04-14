#' Import VCF (gz) to Granges object
#'
#' A wrapper VCF import function of VariantAnnotation to tidy up the VCF and convert to GRanges
#'
#' @param file A path to the VCF (GZ) file 
#' @param keep_seq_levels A character vector indicating which seqlevels to keep
#'
#' @return A \code{GRanges} object embedding variant information in the metadata columns
#'
#' @examples
#' 
#' @export
read_vcf <- function(file, keep_seq_levels=NULL) {
  
  stopifnot(file.exists(file))
  
  vcf <- VariantAnnotation::readVcf(file, genome="hg38")
  vcf_gr <- rowRanges(vcf)
  mcols(vcf_gr) <- append(mcols(vcf_gr), VariantAnnotation::info(vcf))
  end(vcf_gr) <- cf_gr$END
  
  if (!is.null(keep_seql_levels))
    vcf_gr <- GenomeInfoDb::keepSeqlevels(vcf_gr,
                                          value=keep_seq_levels,
                                          pruning.mode='coarse')
  
  return(vcf_gr)
}
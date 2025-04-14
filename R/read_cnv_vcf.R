#' Import VCF (gz) to Granges object
#'
#' A wrapper \code{VariantAnnotation::readVcf} to tidy up the VCF and convert to GRanges
#'
#' @param file A path to the VCF (GZ) file 
#' @param genome A character string of Genome identifier corresponding to chromosome names in the file. This identifier replaces the genome information in the VCF Seqinfo (i.e., \code{seqinfo(vcf)}). When not provided, genome is taken from the VCF file header
#' @param keep_seq_levels A character vector indicating which \code{seqlevels} to keep
#'
#' @return A \code{GRanges} object embedding variant information in the metadata columns
#'
#' @examples
#' @seealso [VariantAnnotation::readVcf()]
#' 
#' @export
read_cnv_vcf <- function(file, genome=null, 
                         keep_seq_levels=NULL) {
  
  stopifnot(file.exists(file))
  
  vcf <- VariantAnnotation::readVcf(file, genome=genome)
  vcf_gr <- rowRanges(vcf)
  mcols(vcf_gr) <- append(mcols(vcf_gr), VariantAnnotation::info(vcf))
  end(vcf_gr) <- cf_gr$END
  
  if (!is.null(keep_seql_levels))
    vcf_gr <- GenomeInfoDb::keepSeqlevels(vcf_gr,
                                          value=keep_seq_levels,
                                          pruning.mode='coarse')
  
  return(vcf_gr)
}
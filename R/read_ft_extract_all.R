#' Import ft extract --all BED GZ
#'
#' Tool to import the output of ft extract <input.bam> --all. It is wrapper to the import family of function in \code{rtracklayer}
#'
#' @param file A path to the BED GZ file yielded by ft extract --all
#' @param col_names An optional characger vector for including additional columns in file that is not part of the BED format
#' @param genome_info An optional character or Ranges object that contains information about the genome (e.g., 'hg38')
#' @param overlap_ragnes An optional \code{GRanges} in the file that overlap the ranges will be return.
#'
#' @return A \code{GRanges} object with 26 additional columns specified in the fibertools-rs `extract` documentation.
#'
#' @examples
#' 
#' @references 
#' The computational guide to Fiber-seq \url{https://fiberseq.github.io/fibertools/extracting/extract.html}
#' 
#' @export
read_ft_extract_all <- function(file, col_names = NULL,
                                genome_info = NULL,
                                overlap_ranges = NULL) {
  args <- plyranges:::norm_args_reader(genome_info)
  extra_cols <- c(
    sam_flag = "character",
    HP = "character",
    RG = "character",
    fiber_length = "numeric",
    fiber_sequence = "DNAStringSet",
    ec = "character",
    rq = "character",
    total_AT_bp = "numeric",
    total_m6a_bp = "numeric",
    total_nuc_bp = "numeric",
    total_msp_bp = "numeric",
    total_5mC_bp = "numeric",
    nuc_starts = "character",
    nuc_lengths = "character",
    ref_nuc_starts = "character",
    ref_nuc_lengths = "character",
    msp_starts = "character",
    msp_lengths = "character",
    fire = "character",
    ref_msp_starts = "character",
    ref_msp_lengths = "character",
    m6a = "character",
    ref_m6a = "character",
    m6a_qual = "character",
    `5mC` = "character",
    ref_5mC = "character",
    `5mC_qual` = "character"
  )

  rtracklayer::import.bed(file,
    colnames = col_names,
    extraCols = extra_cols,
    genome = args$genome_info,
    seqinfo = args$seq_info,
    which = overlap_ranges
  )
}

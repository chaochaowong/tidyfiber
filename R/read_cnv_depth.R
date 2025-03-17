#' import copy number depth bedgraph and process to 5k-binned
#'
#' This function takes a path of the depth's bedgraph
#'
#' @param file A numeric value to be doubled.
#' @return A data.frame of the copy number depth information from hificnv
#' @examples
#' my_function(5)  # Returns 10
#' my_function(10) # Returns 20
#' @export
binned_cnv <- function(file, keep_seq_level = paste0('chr', 1:22),
                       bin_size = 5e3L) {
  require(plyranges)

  # assert file exist
  assert <- all(file.exists(file))
  if (!assert)
    stop(paste0(file, ' does not exist.'))

  # import bedgraph
  # Read the bedGraph file using readr
  depth_gr <- plyranges::read_bed_graph(path)

  # keep standard chromosome levels
  # check if keep_seq_levels exists
  depth_gr <- GenomeInfoDb::keepSeqlevels(depth_gr,
                                          value = keep_seq_levels,
                                          pruning.mode = 'coarse')
  .binned_depth(depth_gr, bin_size = bin_size)
}

.binned_depth <- function(depth_gr, bin_size = 5e3L) {

  # Process data using dplyr pipeline
  binned_gr <- as.data.frame(depth_gr) %>%
    dplyr::mutate(bin_start = (start %/% bin_size) * bin_size) %>%
    dplyr::group_by(seqnames, bin_start) %>%
    dplyr::summarize(score = mean(score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(bin_end = bin_start + bin_size) %>%
    dplyr::select(seqnames, bin_start, bin_end, score) %>%
    dplyr::rename(start = bin_start, end = bin_end) %>%
    dplyr::mutate(delta = score - median(score)) %>%
    plyranges::as_granges()
}

define_pallette <- function(name) {
  crayon_palette <- c(
    "chr1" = "#F37735", "chr2" = "#7BC8A4", "chr3" = "#9B59B6", "chr4" = "#E74C3C",
    "chr5" = "#3498DB", "chr6" = "#2ECC71", "chr7" = "#F1C40F", "chr8" = "#E67E22",
    "chr9" = "#1ABC9C", "chr10" = "#9B59B6", "chr11" = "#34495E", "chr12" = "#D35400",
    "chr13" = "#E84393", "chr14" = "#16A085", "chr15" = "#27AE60", "chr16" = "#2980B9",
    "chr17" = "#8E44AD", "chr18" = "#C0392B", "chr19" = "#F39C12", "chr20" = "#BDC3C7",
    "chr21" = "#95A5A6", "chr22" = "#2C3E50", "chrX" = "#E91E63", "chrY" = "#607D8B"
  )
}


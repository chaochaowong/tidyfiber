#' import copy number depth bedgraph and process to 5k-binned
#'
#' This function takes a path of the depth's bedgraph
#'
#' @param file A numeric value to be doubled.
#' @return A data.frame of the copy number depth information from hificnv
#' @examples
#' my_function(5) # Returns 10
#' my_function(10) # Returns 20
#' @export
binned_depth <- function(file,
                         keep_seq_levels = paste0("chr", 1:22),
                         bin_size = 5e3L) {
  require(plyranges)
  # file -> bigwig file
  # assert file exist
  assert <- all(file.exists(file))
  if (!assert) {
    stop(paste0(file, " does not exist."))
  }

  # import bedgraph
  # Read the bedGraph file using readr
  depth_gr <- plyranges::read_bigwig(file)

  # keep standard chromosome levels
  # check if keep_seq_levels exists
  depth_gr <- GenomeInfoDb::keepSeqlevels(depth_gr,
    value = keep_seq_levels,
    pruning.mode = "coarse"
  )
  .binned_depth(depth_gr, bin_size = bin_size)
}

.binned_depth <- function(depth_gr, bin_size = 5e3L) {
  # Process data using dplyr pipeline
  binned_gr <- as.data.frame(depth_gr) %>%
    dplyr::mutate(bin_start = (start %/% bin_size) * bin_size) %>%
    dplyr::group_by(seqnames, bin_start) %>%
    dplyr::summarize(
      score = mean(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(bin_end = bin_start + bin_size) %>%
    dplyr::select(seqnames, bin_start, bin_end, score) %>%
    dplyr::rename(start = bin_start, end = bin_end) %>%
    dplyr::mutate(delta = score - median(score)) %>%
    plyranges::as_granges()
}

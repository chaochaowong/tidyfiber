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
  depth_gr <- plyranges::read_bed_graph(file)

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

plot_cna <- function(binned_gr, chrom_sizes_gr,
                     where=NULL) {
  require(RColorBrewer)
  
  # filter
  if (!is.null(where)) {
    binned_gr <- binned_gr %>%
      plyranges::filter_by_overlaps(where)
  }
  
  # sanity check
  if (length(binned_gr) < 1) 
    stop(paste0(as.character(where), ' is not within the binned_gr.'))
  
  
  keep_seq_levels <- seqlevels(binned_gr)
  
  # define color palette
  chrom_palette <-
    colorRampPalette(brewer.pal(12,"Set3"))(length(keep_seq_levels))
  
  # Assign colors to chromosome names
  chrom_colors <- setNames(chrom_palette, keep_seq_levels)
  
  # must have a sanity check on the seq-levels
  chrom_sizes_gr <- keepSeqlevels(chrom_sizes_gr,
                                  value=keep_seq_levels,
                                  pruning.mode='coarse')
  
  df <- as.data.frame(binned_gr) %>%
    dplyr::inner_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_bin', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_bin + cum_start) 
  
  df %>%
    ggplot(aes(y=score, x=global_start, color=seqnames)) +
    geom_point(size=0.2, alpha=0.7) +
    #geom_vline(xintercep = chrom_sizes_gr$cum_start,
    #           linetype = 'dashed', color='grey75',
    #           alpha=0.5) +
    scale_x_continuous(
      breaks = chrom_sizes_gr$cum_midpoint,  
      labels = seqlevels(chrom_sizes_gr) 
    ) +
    labs(title = "5k Binned HiFiCNV Depth", 
         x = "Chromosome", y = "Depth") +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.position = "none")  +
    scale_y_continuous(limits = c(0, 25)) + 
    scale_color_manual(values = chrom_colors)
}




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
binned_depth <- function(file, 
                         keep_seq_levels = paste0('chr', 1:22),
                         bin_size = 5e3L) {
  require(plyranges)
  # file -> bigwig file
  # assert file exist
  assert <- all(file.exists(file))
  if (!assert)
    stop(paste0(file, ' does not exist.'))

  # import bedgraph
  # Read the bedGraph file using readr
  depth_gr <- plyranges::read_bigwig(file)

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
                     where=NULL,
                     ylim = c(0, 75),
                     annotate_cnv = FALSE,
                     cnv_gr = NULL,
                     cnv_ypos = ylim[2]-5) {
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
    colorRampPalette(brewer.pal(12, 'Set3'))(length(keep_seq_levels))
  
  # Assign colors to chromosome names
  chrom_colors <- setNames(chrom_palette, keep_seq_levels)
  
  # must have a sanity check on the seq-levels
  chrom_sizes_gr <- keepSeqlevels(chrom_sizes_gr,
                                  value=keep_seq_levels,
                                  pruning.mode='coarse')
  
  # prepare data.frame for binned depth: get global starts
  df <- as.data.frame(binned_gr) %>%
    dplyr::inner_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_bin', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_bin + cum_start) 
  
  # plot CNA
  gg <- df %>%
    ggplot(aes(y=score, x=global_start, color=seqnames)) +
    geom_point(size=0.2, alpha=0.8, show.legend = FALSE) +
    scale_x_continuous(
      breaks = chrom_sizes_gr$cum_midpoint,  
      labels = seqlevels(chrom_sizes_gr) 
    ) +
    labs(title = "5k Binned HiFiCNV Depth", 
         x = "Chromosome", y = "Depth") +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.position = "bottom")  +
    scale_y_continuous(limits = ylim) + 
    scale_color_manual(values = chrom_colors)
  
  if (annotate_cnv) 
    gg <- .annotate_cnv(gg, keep_seq_levels, 
                        chrom_sizes_gr, 
                        cnv_gr, 
                        cnv_ypos=ylim[2]-5)
  return(gg)
}

.annotate_cnv <- function(gg, 
                          keep_seq_levels,
                          chrom_sizes_gr, 
                          cnv_gr, 
                          cnv_ypos=ylim[2]-5) {
  
  # sanity check
  # tidy up cnv_gr
  cnv_gr <- keepSeqlevels(cnv_gr,
                          value=keep_seq_levels,
                          pruning.mode='coarse')
  
  cnv_gr <- cnv_gr %>% 
    plyranges::filter(!is.na(SVTYPE))
  
  # convert to data.frame and get global_start and global_end
  cnv_df <- as.data.frame(cnv_gr) %>%
    dplyr::inner_join(as.data.frame(chrom_sizes_gr), by='seqnames',
                      suffix=c('_cnv', '_chrom_sizes')) %>%
    dplyr::mutate(global_start = start_cnv + cum_start,
                  global_end = end_cnv + cum_start) 
  gg <- gg +
    geom_rect(data = cnv_df,
              aes(xmin = global_start, xmax = global_end,
                  ymin = cnv_ypos, ymax = cnv_ypos+5, fill = SVTYPE),
              inherit.aes = FALSE,
              alpha = 0.5) +
    #guides(fill = guide_legend(title = "SVTYPE")) +
    theme(
      legend.key.size = unit(0.4, "cm"),       # smaller boxes
      legend.text = element_text(size = 8),    # smaller text
      legend.title = element_text(size = 9),   # smaller title
      legend.spacing.y = unit(0.1, "cm")       # less vertical spacing
    )
  
  return(gg)
}


